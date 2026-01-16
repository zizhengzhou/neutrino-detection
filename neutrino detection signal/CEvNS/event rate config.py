import numpy as np
import pandas as pd
import matplotlib

# [设置] Matplotlib 后端为 Agg (非交互模式)
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.integrate import simpson, quad
from scipy.special import spherical_jn
import os
import json
import glob
from multiprocessing import Pool

# --- PHYSICAL CONSTANTS ---
TARGET_MATERIAL = 'CaWO4'  # Material Name
ROI_MIN_KEV = 0.1  # [keV]
ROI_MAX_KEV = 1.0  # [keV]

SIN2_THETA_W = 0.23867
G_F = 1.1663787e-11  # [MeV^-2] Fermi coupling constant
H_BAR_C_CM = 197.327e-13  # [MeV*cm]
H_BAR_C_FM = 197.327  # [MeV*fm] for form factor
SECONDS_PER_DAY = 86400.0
DAYS_PER_YEAR = 365.25
MEV_PER_KEV = 1000.0
AMU_MEV = 931.4941

# --- TARGET PROPERTIES ---
TARGET_PROPERTIES = {
    'Si': {
        'molar_mass_kg': 28.0855 / 1000.0,
        'elements': [{'Z': 14, 'A': 28, 'n': 1}]
    },
    'Xe': {
        'molar_mass_kg': 131.293 / 1000.0,
        'elements': [{'Z': 54, 'A': 131, 'n': 1}]
    },
    'Ge': {
        'molar_mass_kg': 72.630 / 1000.0,
        'elements': [{'Z': 32, 'A': 73, 'n': 1}]
    },
    'GaAs': {
        'molar_mass_kg': 144.64 / 1000.0,
        'elements': [
            {'Z': 31, 'A': 70, 'n': 1},
            {'Z': 33, 'A': 75, 'n': 1}
        ]
    },
    'CaWO4': {
        'molar_mass_kg': 287.92 / 1000.0,
        'elements': [
            {'Z': 20, 'A': 40, 'n': 1},  # Ca
            {'Z': 74, 'A': 184, 'n': 1},  # W
            {'Z': 8, 'A': 16, 'n': 4}  # O
        ]
    },
    'Li2MoO4': {
        'molar_mass_kg': 173.83 / 1000.0,
        'elements': [
            {'Z': 3, 'A': 7, 'n': 2},
            {'Z': 42, 'A': 96, 'n': 1},
            {'Z': 8, 'A': 16, 'n': 4}
        ]
    },
}


# --- PHYSICS FUNCTIONS ---

def get_target_nuclei_per_kg(material):
    """Return molecules per kg."""
    props = TARGET_PROPERTIES.get(material)
    if not props:
        raise ValueError(f"Material {material} not found")
    return 6.02214076e23 / props['molar_mass_kg']


def get_atomic_masses(material):
    """Return list of masses [MeV] for all elements in the molecule."""
    props = TARGET_PROPERTIES.get(material)
    return [elem['A'] * AMU_MEV for elem in props['elements']]


def T_max_nucleus(E_nu, M_nucleus):
    """Max recoil energy for a given neutrino energy and nucleus mass."""
    # Prevent division by zero if E_nu is 0
    return np.divide(2 * E_nu ** 2, (M_nucleus + 2 * E_nu), where=(E_nu > 0))


def helm_form_factor(Q_mev, A):
    """
    Calculates Helm Form Factor F(Q^2). Vectorized for numpy arrays.
    Q_mev: Momentum transfer in MeV (can be array)
    A: Mass number
    """
    Q_fm_inv = Q_mev / H_BAR_C_FM
    s = 0.9
    R_A = (1.23 * A ** (1 / 3) - 0.6)
    qr = Q_fm_inv * R_A

    # Handle qr ~ 0 to avoid singularity
    # mask_small = np.abs(qr) < 1e-4
    # mask_large = ~mask_small

    # Using scipy's spherical_jn which handles small arguments well usually,
    # but explicit handling is safer for 0.

    with np.errstate(divide='ignore', invalid='ignore'):
        bessel_term = 3 * spherical_jn(n=1, z=qr) / qr

    # Fix division by zero (lim x->0 of 3*j1(x)/x = 1)
    bessel_term = np.where(np.abs(qr) < 1e-4, 1.0 - (qr ** 2) / 10.0, bessel_term)

    gaussian_term = np.exp(-(Q_fm_inv * s) ** 2 / 2)
    return bessel_term * gaussian_term


def calculate_dsigmadT_single_element_vectorized(E_nu_grid, T_recoil_vals, Z, A):
    """
    Calculates dSigma/dT for a SINGLE element type on a grid.
    Optimized for broadcasting:
    - E_nu_grid: Shape (N_E,)
    - T_recoil_vals: Shape (N_T, 1) or (N_T,)

    Returns: Matrix of shape (N_T, N_E) representing cross-section values.
    """
    M_nucleus = A * AMU_MEV
    N = A - Z
    Q_w = N - (1 - 4 * SIN2_THETA_W) * Z

    # Ensure T is column vector for broadcasting: (N_T, 1)
    T = T_recoil_vals.reshape(-1, 1)
    # Ensure E is row vector: (1, N_E)
    E = E_nu_grid.reshape(1, -1)

    # 1. Kinematic Cutoff: E_nu must be > E_min(T)
    # Inverting T_max = 2E^2 / (M + 2E) gives:
    # E_min = T/2 + 0.5 * sqrt(T^2 + 2*M*T)
    E_min = 0.5 * T + 0.5 * np.sqrt(T ** 2 + 2 * M_nucleus * T)

    # Mask where E_nu is physically capable of producing recoil T
    # Shape: (N_T, N_E)
    mask_kinematic = E > E_min

    # 2. Calculate Cross Section Terms
    # Q only depends on T
    Q = np.sqrt(2 * M_nucleus * T)
    F_Q2 = helm_form_factor(Q, A)  # Shape (N_T, 1)

    prefactor = (G_F ** 2 * M_nucleus) / (4 * np.pi)

    # Kinematic term depends on both T and E
    # (1 - MT / 2E^2)
    # Use 'where' to avoid division by zero or invalid values outside mask
    kinematic_term = np.zeros_like(mask_kinematic, dtype=float)

    valid_indices = mask_kinematic
    if np.any(valid_indices):
        # Calculation only needed where E > E_min
        # We broadcast T (column) and E (row)
        term = 1 - (M_nucleus * T) / (2 * E ** 2)
        # Apply strict kinematic cutoff (term must be >= 0, though E > E_min ensures this)
        np.maximum(term, 0, out=term)
        kinematic_term = np.where(valid_indices, term, 0.0)

    # Combine
    # shape: scalar * (N_T, 1) * (N_T, 1) * (N_T, N_E) -> (N_T, N_E)
    dsigma = prefactor * (Q_w ** 2) * (F_Q2 ** 2) * kinematic_term

    return dsigma


# --- [优化] 高速率计算函数 ---
def calculate_event_rate_optimized(T_recoil_mev, flux_function, E_nu_max, material):
    """
    High-performance calculation of dR/dT using vectorized numpy operations.
    Implements the formula:
      dR/dT = Sum_i [ w_i * N_i * Integral(Flux(E) * dSigma_i/dT(E) dE) ]
    """
    # 1. Prepare Integration Grid for Neutrino Energy
    # Create a fine grid for integration (Simpsons rule)
    # 1000 points is usually sufficient for smooth fluxes, 2000 for safety
    N_E_POINTS = 20000
    E_grid = np.linspace(1e-5, E_nu_max, N_E_POINTS)

    # Evaluate Flux on this grid once
    Flux_grid = flux_function(E_grid)  # Shape (N_E,)
    # Reshape for broadcasting: (1, N_E)
    Flux_grid_row = Flux_grid.reshape(1, -1)

    # Conversion factor: Natural units to Rate
    # [MeV^-3] * [Flux] * [dE] -> Counts...
    # Need to multiply by (hbar*c)^2 * Nuclei_per_kg * Time_units
    conversion = (H_BAR_C_CM ** 2) * get_target_nuclei_per_kg(material) * SECONDS_PER_DAY * DAYS_PER_YEAR / MEV_PER_KEV

    total_rate = np.zeros_like(T_recoil_mev)

    props = TARGET_PROPERTIES.get(material)

    # Loop over elements (Sum_i)
    for elem in props['elements']:
        Z = elem['Z']
        A = elem['A']
        n_atoms = elem['n']

        # Calculate dSigma matrix: Shape (N_T, N_E)
        # This handles the E_cutoff(M_i) internally via masking (returning 0)
        dSigma_matrix = calculate_dsigmadT_single_element_vectorized(E_grid, T_recoil_mev, Z, A)

        # Integrand = Flux(E) * dSigma(E, T)
        # Broadcasting: (1, N_E) * (N_T, N_E) -> (N_T, N_E)
        integrand = Flux_grid_row * dSigma_matrix

        # Integrate over E (axis 1)
        # using Simpson's rule is better than trapz for smooth functions
        integral_results = simpson(y=integrand, x=E_grid, axis=1)

        # Add to total rate: w_i * N_i * Integral
        # (Note: get_target_nuclei_per_kg gives molecules/kg, so we multiply by n_atoms)
        total_rate += n_atoms * integral_results

    return total_rate * conversion


# --- Helper for Total Cross Section (Fig 3) ---
def calculate_total_sigma_at_E_compound(E_nu, material):
    """
    Calculates total cross section at a fixed Energy E_nu.
    Sum of integrals: Sigma_tot = Sum_i ( n_i * Integral(dSigma_i/dT dT) )
    """
    props = TARGET_PROPERTIES.get(material)
    total_sigma = 0.0

    for elem in props['elements']:
        Z = elem['Z']
        A = elem['A']
        n_atoms = elem['n']
        M_nucleus = A * AMU_MEV

        T_max = T_max_nucleus(E_nu, M_nucleus)

        # Use simple quad for 1D integral as this is not the bottleneck
        # dSigma/dT function for quad (scalar T)
        def integrand(T):
            # Re-use the vector logic but for scalar inputs or simplify
            # Inline simplified version for speed
            Q = np.sqrt(2 * M_nucleus * T)
            if Q == 0: return 0
            F = helm_form_factor(Q, A)
            kin = 1 - (M_nucleus * T) / (2 * E_nu ** 2)
            N_n = A - Z
            Qw = N_n - (1 - 4 * SIN2_THETA_W) * Z
            pre = (G_F ** 2 * M_nucleus) / (4 * np.pi)
            return pre * (Qw ** 2) * (F ** 2) * kin

        val, _ = quad(integrand, 0, T_max, limit=100)
        total_sigma += n_atoms * val

    return total_sigma * (H_BAR_C_CM ** 2)


# --- CORE LOGIC WRAPPER ---
def run_analysis(config):
    file_path = "../" + config['flux_file_path']
    save_name = config['figure_save_name'].replace("NMM", "CEvNS")
    raw_title = config['title']
    title_str = f"{TARGET_MATERIAL} : CEvNS : {raw_title}"
    flux_x_limit = config.get('plot_flux_x_limit', 15.0)

    print(f"--> [PID {os.getpid()}] Processing: {file_path}")

    # 1. Load Data
    try:
        flux_data = pd.read_csv(file_path, comment='#')
        flux_data = flux_data.sort_values(by=flux_data.columns[0]).reset_index(drop=True)
        E_nu_data, phi_data = flux_data.iloc[:, 0].values, flux_data.iloc[:, 1].values
    except Exception as e:
        print(f"    [Error] Reading CSV file '{file_path}' failed: {e}")
        return

    # 2. Interpolate Flux
    flux_function = interp1d(E_nu_data, phi_data, kind='cubic', bounds_error=False, fill_value=0.0)

    E_nu_max_data = E_nu_data[-1]
    E_nu_peak = E_nu_data[np.argmax(phi_data)]

    # 3. Define Recoil Energy Range (T)
    # Determine plotting range based on Lightest nucleus (extends to highest T)
    masses = get_atomic_masses(TARGET_MATERIAL)
    M_LIGHTEST = min(masses)
    T_max_abs_kev = T_max_nucleus(E_nu_max_data, M_LIGHTEST) * MEV_PER_KEV

    # Generate T points (log scale)
    T_recoil_kev = np.logspace(np.log10(1e-3), np.log10(T_max_abs_kev * 1.05), 300)
    T_recoil_mev = T_recoil_kev / MEV_PER_KEV

    # 4. Calculate Rates (Vectorized & Physics Corrected)
    rate_cevns_yearly = calculate_event_rate_optimized(T_recoil_mev, flux_function, E_nu_max_data, TARGET_MATERIAL)

    # 5. ROI Integration
    mask_roi = (T_recoil_kev >= ROI_MIN_KEV) & (T_recoil_kev <= ROI_MAX_KEV)
    if np.any(mask_roi):
        total_counts_roi = simpson(y=rate_cevns_yearly[mask_roi], x=T_recoil_kev[mask_roi])
    else:
        total_counts_roi = 0.0

    # 6. Calculate Other Plots Data (Total Sigma & Diff Sigma at Max E)
    E_sigma_plot = np.linspace(0.1, E_nu_max_data, 100)
    sigma_total_vals = np.array([calculate_total_sigma_at_E_compound(E, TARGET_MATERIAL) for E in E_sigma_plot])
    sigma_at_peak = calculate_total_sigma_at_E_compound(E_nu_peak, TARGET_MATERIAL)

    # Diff Sigma at Max Energy (Weighted sum of elements)
    # Using the same vectorization logic but with single E point
    diff_sigma_max = np.zeros_like(T_recoil_mev)
    props = TARGET_PROPERTIES.get(TARGET_MATERIAL)
    sigma_conversion = (H_BAR_C_CM ** 2) / MEV_PER_KEV  # to cm^2/keV

    E_max_grid = np.array([E_nu_max_data])
    for elem in props['elements']:
        dsig = calculate_dsigmadT_single_element_vectorized(E_max_grid, T_recoil_mev, elem['Z'], elem['A'])
        diff_sigma_max += elem['n'] * dsig.flatten()

    diff_sigma_max *= sigma_conversion

    # --- PLOTTING ---
    plt.style.use('seaborn-v0_8-talk')
    plt.rcParams.update({
        'font.family': 'serif', 'mathtext.fontset': 'cm', 'font.size': 16,
        'axes.labelsize': 18, 'axes.titlesize': 16, 'xtick.labelsize': 16,
        'ytick.labelsize': 16, 'legend.fontsize': 15, 'figure.titlesize': 24
    })

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(title_str)

    # Flux
    ax_flux = axes[0, 0]
    ax_flux.plot(E_nu_data, phi_data, 'o', ms=4, alpha=0.5)
    E_plot = np.linspace(0, E_nu_max_data, 500)
    ax_flux.plot(E_plot, flux_function(E_plot), '-', lw=2)
    ax_flux.set(xlabel='Neutrino Energy [MeV]', ylabel='Flux', title='Flux Input')
    ax_flux.set_xlim(0, flux_x_limit if flux_x_limit else None)
    ax_flux.grid(True)

    # Diff Sigma
    ax_diff = axes[0, 1]
    ax_diff.loglog(T_recoil_kev, diff_sigma_max, 'r-', lw=2)
    ax_diff.set(xlabel='Recoil T [keV]', ylabel=r'd$\sigma$/dT [cm$^2$/keV]',
                title=f'Diff. Cross Section @ {E_nu_max_data:.1f} MeV')
    ax_diff.grid(True, which='both')
    ax_diff.set_ylim(bottom=1e-50)  # simple catch

    # Total Sigma
    ax_tot = axes[1, 0]
    ax_tot.semilogy(E_sigma_plot, sigma_total_vals, 'purple', lw=2)
    ax_tot.plot(E_nu_peak, sigma_at_peak, 'ko')
    ax_tot.text(E_nu_peak, sigma_at_peak, f" {sigma_at_peak:.1e}", va='bottom')
    ax_tot.set(xlabel='Neutrino Energy [MeV]', ylabel=r'$\sigma_{tot}$ [cm$^2$]', title='Total Cross Section')
    ax_tot.grid(True, which='both')

    # Rate
    ax_rate = axes[1, 1]
    ax_rate.loglog(T_recoil_kev, rate_cevns_yearly, 'g-', lw=2)
    ax_rate.fill_between(T_recoil_kev, rate_cevns_yearly, 0, where=mask_roi, color='green', alpha=0.3)
    ax_rate.text(0.05, 0.1, f"ROI Sum: {total_counts_roi:.2f} cts/kg/yr", transform=ax_rate.transAxes,
                 bbox=dict(facecolor='wheat', alpha=0.8))
    ax_rate.set(xlabel='Recoil T [keV]', ylabel='Rate [cts/kg/yr/keV]', title='Differential Rate')
    ax_rate.grid(True, which='both')

    # Fix limits for Rate if empty
    if np.max(rate_cevns_yearly) > 0:
        ax_rate.set_ylim(bottom=np.max(rate_cevns_yearly) * 1e-6)

    plt.tight_layout()
    os.makedirs(os.path.dirname(save_name), exist_ok=True)
    plt.savefig(save_name, dpi=300)
    print(f"    [Success] Saved: {save_name}")
    plt.close(fig)


# --- MAIN EXECUTION ---
if __name__ == "__main__":
    setting_files = glob.glob("../data/**/setting.json", recursive=True)

    all_tasks = []
    for s_file in setting_files:
        try:
            with open(s_file, 'r', encoding='utf-8') as f:
                content = json.load(f)
                if isinstance(content, list): all_tasks.extend(content)
        except:
            pass

    if not all_tasks:
        print("No tasks found.")
        exit()

    # 根据你的机器调整进程数，向量化计算后CPU占用较高，建议设置为核心数-1
    num_processes = 8
    print(f"Starting {len(all_tasks)} jobs for {TARGET_MATERIAL}...")

    with Pool(processes=num_processes) as pool:
        pool.map(run_analysis, all_tasks)

    print("Done.")