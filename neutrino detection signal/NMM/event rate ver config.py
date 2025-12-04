import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.integrate import quad, simpson
import os
import json
import glob

# --- PHYSICAL CONSTANTS ---
NMM_IN_BOHR_MAGNETONS = 2.8e-11
TARGET_MATERIAL = 'Si'
ROI_MIN_KEV = 0.010
ROI_MAX_KEV = 0.10
SIN2_THETA_W = 0.23867
G_F = 1.1663787e-11
M_e = 0.51099895
ALPHA = 1 / 137.036
H_BAR_C = 197.327e-13
SECONDS_PER_DAY = 86400.0
DAYS_PER_YEAR = 365.25
MEV_PER_KEV = 1000.0


# --- HELPER FUNCTIONS ---
def get_target_electrons_per_kg(material):
    if material == 'Si':
        z = 14
        molar_mass_kg = 28.0855 / 1000.0
    else:
        raise ValueError(f"Unsupported material: {material}")
    avogadro_number = 6.02214076e23
    return z * (avogadro_number / molar_mass_kg)


def T_max(E_nu):
    return E_nu / (1 + M_e / (2 * E_nu))


def calculate_dsigmadT_sm(E_nu, T, neutrino_type):
    if T <= 0 or T > T_max(E_nu): return 0.0
    flavor = neutrino_type.replace('_bar', '')
    if flavor == 'e':
        gV_base, gA_base = 0.5 + 2 * SIN2_THETA_W, 0.5
    elif flavor in ['mu', 'tau']:
        gV_base, gA_base = -0.5 + 2 * SIN2_THETA_W, -0.5
    else:
        raise ValueError(f"Invalid neutrino flavor in '{neutrino_type}'")
    gV = gV_base
    gA = -gA_base if 'bar' in neutrino_type else gA_base
    term1 = (gV + gA) ** 2
    term2 = (gV - gA) ** 2 * (1 - T / E_nu) ** 2
    term3 = (gV ** 2 - gA ** 2) * M_e * T / E_nu ** 2
    prefactor = G_F ** 2 * M_e / (2 * np.pi)
    return prefactor * (term1 + term2 - term3)


def dsigmadT_nmm(E_nu, T):
    if T <= 0 or T > T_max(E_nu): return 0.0
    mu_nu_sq = NMM_IN_BOHR_MAGNETONS ** 2
    prefactor = np.pi * ALPHA ** 2 / M_e ** 2
    kinematic_term = (1 / T - 1 / E_nu - M_e / (2 * E_nu ** 2))
    if kinematic_term < 0: return 0.0
    return prefactor * mu_nu_sq * kinematic_term


def calculate_event_rate(T_recoil_mev, flux_function, E_nu_max, neutrino_type):
    rate_sm, rate_nmm = np.zeros_like(T_recoil_mev), np.zeros_like(T_recoil_mev)
    conversion = (H_BAR_C ** 2) * get_target_electrons_per_kg(
        TARGET_MATERIAL) * SECONDS_PER_DAY * DAYS_PER_YEAR / MEV_PER_KEV

    for i, T in enumerate(T_recoil_mev):
        E_nu_min = 0.5 * (T + np.sqrt(T ** 2 + 2 * T * M_e))
        if E_nu_min >= E_nu_max: continue
        integrand_sm = lambda E: flux_function(E) * calculate_dsigmadT_sm(E, T, neutrino_type)
        integrand_nmm = lambda E: flux_function(E) * dsigmadT_nmm(E, T)
        try:
            integral_sm, _ = quad(integrand_sm, E_nu_min, E_nu_max, limit=100)
            integral_nmm, _ = quad(integrand_nmm, E_nu_min, E_nu_max, limit=100)
        except:
            integral_sm, integral_nmm = 0.0, 0.0
        rate_sm[i] = integral_sm * conversion
        rate_nmm[i] = integral_nmm * conversion
    return rate_sm, rate_nmm


# --- CORE LOGIC WRAPPER ---
def run_analysis(config):
    file_path = "../" + config['flux_file_path']
    save_name = config['figure_save_name']

    raw_title = config['title']
    title_str = f"NMM : {raw_title}"

    nu_type = config['neutrino_type']
    flux_x_limit = config.get('plot_flux_x_limit', 15.0)

    print(f"--> Processing: {file_path}")
    print(f"    Target: {save_name}")

    try:
        flux_data = pd.read_csv(file_path, comment='#')
        flux_data = flux_data.sort_values(by=flux_data.columns[0]).reset_index(drop=True)
        E_nu_data, phi_data = flux_data.iloc[:, 0].values, flux_data.iloc[:, 1].values
    except Exception as e:
        print(f"    [Error] Reading CSV file '{file_path}' failed: {e}")
        return

    flux_function = interp1d(E_nu_data, phi_data, kind='cubic', bounds_error=False, fill_value=0.0)
    E_nu_max = E_nu_data[-1]
    E_nu_peak = E_nu_data[np.argmax(phi_data)]

    T_max_plot_mev = T_max(E_nu_max)
    T_max_plot_kev = T_max_plot_mev * MEV_PER_KEV
    T_peak_limit_kev = T_max(E_nu_peak) * MEV_PER_KEV

    T_recoil_kev = np.logspace(np.log10(1e-3), np.log10(T_max_plot_kev * 1.1), 300)
    T_recoil_mev = T_recoil_kev / MEV_PER_KEV

    rate_sm, rate_nmm = calculate_event_rate(T_recoil_mev, flux_function, E_nu_max, nu_type)

    mask_roi = (T_recoil_kev >= ROI_MIN_KEV) & (T_recoil_kev <= ROI_MAX_KEV)
    if np.any(mask_roi):
        counts_sm = simpson(y=rate_sm[mask_roi], x=T_recoil_kev[mask_roi])
        counts_nmm = simpson(y=rate_nmm[mask_roi], x=T_recoil_kev[mask_roi])
    else:
        counts_sm, counts_nmm = 0.0, 0.0

    # --- Plotting Configuration ---
    plt.style.use('seaborn-v0_8-talk')

    plt.rcParams.update({
        'font.family': 'serif',
        'mathtext.fontset': 'cm',
        'font.size': 16,
        'axes.labelsize': 18,
        'axes.titlesize': 16,  # 子图标题小一号
        'xtick.labelsize': 16,
        'ytick.labelsize': 16,
        'legend.fontsize': 15,
        'figure.titlesize': 24
    })

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(f'{title_str}')

    ax_flux, ax_sigma_sm, ax_sigma_nmm, ax_rate = axes.flatten()

    # [Plot 1: Flux]
    E_plot = np.linspace(0, E_nu_max, 500)
    ax_flux.plot(E_nu_data, phi_data, 'o', label='Data Points', markersize=6)
    ax_flux.plot(E_plot, flux_function(E_plot), '-', label='Interpolation', lw=3)

    # [修复1] 恢复 Peak 数值显示
    ax_flux.axvline(E_nu_peak, color='k', ls=':', label=f'Peak $E_\\nu$ ({E_nu_peak:.2f} MeV)')

    x_lim_actual = flux_x_limit if flux_x_limit else E_nu_max * 1.1
    # [修复2] 恢复原始 Axis Labels
    ax_flux.set(xlim=(0, x_lim_actual),
                xlabel='Neutrino Energy $E_\\nu$ [MeV]',
                ylabel=r'Flux [$\nu$ / cm$^2$ / s / MeV]',
                title='Input Neutrino Flux')

    # [修复3] 强制使用科学计数法
    ax_flux.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

    ax_flux.legend(loc='upper right')
    ax_flux.grid(True)

    # [Sigma Setup]
    sigma_conversion = H_BAR_C ** 2 / MEV_PER_KEV
    sigma_sm_max = np.array([calculate_dsigmadT_sm(E_nu_max, T, nu_type) for T in T_recoil_mev]) * sigma_conversion
    sigma_nmm_max = np.array([dsigmadT_nmm(E_nu_max, T) for T in T_recoil_mev]) * sigma_conversion

    # [Plot 2: SM]
    ax_sigma_sm.loglog(T_recoil_kev, sigma_sm_max, color='C0', lw=3, label=f'Limit ($E_\\nu$={E_nu_max:.1f} MeV)')
    # [修复1] 恢复 Cutoff 数值显示
    ax_sigma_sm.axvline(T_peak_limit_kev, color='k', ls=':', label=f'Peak Flux ({T_peak_limit_kev:.2f} keV)')
    # [修复2] 恢复原始 Labels
    ax_sigma_sm.set(xlabel='Recoil Energy T [keV]', ylabel=r'd$\sigma$/dT [cm$^2$/keV]',
                    title='SM Differential Cross Section')
    if np.any(sigma_sm_max > 0): ax_sigma_sm.set_ylim(bottom=sigma_sm_max[sigma_sm_max > 0].min() / 2)
    ax_sigma_sm.legend()
    ax_sigma_sm.grid(True, which='both')

    # [Plot 3: NMM]
    ax_sigma_nmm.loglog(T_recoil_kev, sigma_nmm_max, color='C1', lw=3, label=f'Limit ($E_\\nu$={E_nu_max:.1f} MeV)')
    # [修复1] 恢复 Cutoff 数值显示
    ax_sigma_nmm.axvline(T_peak_limit_kev, color='k', ls=':', label=f'Peak Flux ({T_peak_limit_kev:.2f} keV)')
    # [修复2] 恢复原始 Labels
    ax_sigma_nmm.set(xlabel='Recoil Energy T [keV]', ylabel=r'd$\sigma$/dT [cm$^2$/keV]',
                     title='NMM Differential Cross Section')
    if np.any(sigma_nmm_max > 0): ax_sigma_nmm.set_ylim(bottom=sigma_nmm_max[sigma_nmm_max > 0].min() / 2)
    ax_sigma_nmm.legend()
    ax_sigma_nmm.grid(True, which='both')

    # [Plot 4: Rate]
    ax_rate.loglog(T_recoil_kev, rate_sm, label='Standard Model (SM)', color='blue', lw=3)
    # [修复4] 恢复原始图例格式 (包含 mu_nu 数值)
    ax_rate.loglog(T_recoil_kev, rate_nmm, label=f'NMM ($\\mu_\\nu={NMM_IN_BOHR_MAGNETONS:.1e} \\mu_B$)', color='red',
                   lw=3)

    ax_rate.fill_between(T_recoil_kev, rate_sm, 0, where=mask_roi, color='blue', alpha=0.1)
    ax_rate.fill_between(T_recoil_kev, rate_nmm, 0, where=mask_roi, color='red', alpha=0.1)

    stats_text = (f"ROI: {ROI_MIN_KEV}-{ROI_MAX_KEV} keV\n"
                  f"SM Counts: {counts_sm:.1f} /yr\n"
                  f"NMM Counts: {counts_nmm:.1f} /yr")
    ax_rate.text(0.05, 0.1, stats_text, transform=ax_rate.transAxes, fontsize=14,
                 verticalalignment='bottom', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # [修复2] 恢复原始 Labels
    ax_rate.set(xlabel='Recoil Energy T [keV]', ylabel='Event Rate [counts / kg / year / keV]',
                title='Final Differential Event Rate')

    max_val = max(rate_sm.max(), rate_nmm.max()) if len(rate_sm) > 0 else 1.0
    if max_val > 0: ax_rate.set_ylim(bottom=max_val * 1e-6, top=max_val * 5)
    ax_rate.legend()
    ax_rate.grid(True, which='both')

    # Save
    os.makedirs(os.path.dirname(save_name), exist_ok=True)
    plt.tight_layout()
    try:
        plt.savefig(save_name, dpi=300)
        print(f"    [Success] Saved to: {save_name}")
    except Exception as e:
        print(f"    [Error] Saving figure failed: {e}")
    plt.close(fig)


# --- MAIN EXECUTION ---
if __name__ == "__main__":
    print("Searching for setting.json files...")

    setting_files = glob.glob("../data/**/setting.json", recursive=True)

    if not setting_files:
        print("No setting.json files found in '../data/'.")
        exit()

    print(f"Found {len(setting_files)} setting file(s).")
    all_tasks = []
    for s_file in setting_files:
        try:
            with open(s_file, 'r', encoding='utf-8') as f:
                configs = json.load(f)
                if isinstance(configs, list):
                    all_tasks.extend(configs)
        except Exception as e:
            print(f"Error reading {s_file}: {e}")

    print(f"Total plots to generate: {len(all_tasks)}")
    print("-" * 40)

    for i, config in enumerate(all_tasks):
        print(f"Job {i + 1}/{len(all_tasks)}")
        run_analysis(config)
        print("-" * 20)

    print("All jobs completed.")