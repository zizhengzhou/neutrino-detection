import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d
from scipy.integrate import cumulative_trapezoid
from matplotlib.lines import Line2D

# ==========================================
# 0. Configuration & Style (Paper Quality)
# ==========================================
plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 10,
    'figure.figsize': (7, 8),  # A4 Width fit
    'mathtext.fontset': 'stix',
    'font.family': 'STIXGeneral',
    'lines.linewidth': 2
})

FILENAME = '../data/Daya Bay/calc_DayaBay_10m_with_uncertainty.csv'

# Define cases
cases = [
    {'label': 'Si', 'T_eV': 10, 'mass': 28.0855, 'color': '#1f77b4', 'ls': '-'},
    {'label': 'Ge', 'T_eV': 10, 'mass': 72.64, 'color': '#2ca02c', 'ls': '--'},
    {'label': 'Xe', 'T_eV': 630, 'mass': 131.293, 'color': '#d62728', 'ls': '-.'},
    # {'label': 'Xe', 'T_eV': 1360, 'mass': 131.293, 'color': '#9467bd', 'ls': ':'}
]


# ==========================================
# 1. Physics Functions
# ==========================================
def calculate_Enu_min(T_eV, mass_amu):
    T_MeV = T_eV * 1.0e-6
    M_MeV = mass_amu * 931.494
    E_nu = (T_MeV + np.sqrt(T_MeV ** 2 + 2 * T_MeV * M_MeV)) / 2
    return E_nu


# ==========================================
# 2. Data Loading & Processing
# ==========================================
try:
    flux_data = pd.read_csv(FILENAME, comment='#')
    flux_data = flux_data.sort_values(by=flux_data.columns[0]).reset_index(drop=True)
    E_nu_data = flux_data.iloc[:, 0].values
    phi_data = flux_data.iloc[:, 1].values
    print(f"Successfully loaded data from {FILENAME}")
except Exception as e:
    print(f"[Warning] Could not read '{FILENAME}': {e}")
    print("[Info] Generating DUMMY data (0-13 MeV)...")
    E_nu_data = np.linspace(0.01, 13, 3000)
    phi_data = (E_nu_data ** 2) * np.exp(-0.8 * E_nu_data) * 1e9
    phi_data[E_nu_data > 12.5] = 0

# --- Calculate Integrated Flux ---
integral_reversed = -cumulative_trapezoid(phi_data[::-1], E_nu_data[::-1], initial=0)
cumulative_flux = integral_reversed[::-1]
total_flux = cumulative_flux[0]
print("total_flux:", total_flux)

integral_func = interp1d(E_nu_data, cumulative_flux, kind='linear', fill_value="extrapolate")

# ==========================================
# 3. Plotting
# ==========================================
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [1, 1]})
plt.subplots_adjust(hspace=0.08)

# --- Top Panel: Differential Flux (Linear) ---
ax1.plot(E_nu_data, phi_data, color='k', lw=2.5, label='Total Flux Spectrum')
ax1.set_ylabel(r'Diff. Flux $\phi(E_\nu)$' + '\n' + r'$[\mathrm{cm}^{-2}\,\mathrm{s}^{-1}\,\mathrm{MeV}^{-1}]$')
ax1.set_title('Usage of Neutrino Flux Comparison in CENvS', pad=15, fontweight='bold')
ax1.grid(True, linestyle=':', alpha=0.5)

# --- Bottom Panel: Integrated Flux (Log Scale) ---
ax2.plot(E_nu_data, cumulative_flux, color='k', lw=2.5)
ax2.set_yscale('log')  # Log scale
ax2.set_ylabel(r'Integrated Flux $\Phi(>E_\nu)$' + '\n' + r'$[\mathrm{cm}^{-2}\,\mathrm{s}^{-1}]$')
ax2.set_xlabel(r'Neutrino Energy $E_\nu$ [MeV]')
ax2.grid(True, linestyle=':', alpha=0.5, which='both')

# --- Plot Cases and Legends ---
top_handles = []
bottom_handles = []

for case in cases:
    # Calc
    E_min = calculate_Enu_min(case['T_eV'], case['mass'])
    flux_avail = integral_func(E_min)
    pct_avail = (flux_avail / total_flux) * 100

    # Top Panel Plot
    ax1.axvline(E_min, color=case['color'], linestyle=case['ls'], lw=1.5)
    ax1.fill_between(E_nu_data, 0, phi_data, where=(E_nu_data >= E_min),
                     color=case['color'], alpha=0.08)

    top_label = f"{case['label']} ($T^\\mathrm{{cutoff}}={case['T_eV']}$ eV)"
    top_handles.append(Line2D([0], [0], color=case['color'], lw=2, linestyle=case['ls'],
                              label=top_label))

    # Bottom Panel Plot
    ax2.axvline(E_min, color=case['color'], linestyle=case['ls'], lw=1.5)
    ax2.plot(E_min, flux_avail, marker='o', color=case['color'], markersize=6, zorder=5)

    # Bottom Legend Label: E^cutoff and Flux usage
    # Using \mathrm for non-italic 'cutoff'
    label_stats = (f"{case['label']}: "
                   f"$E_\\nu^\\mathrm{{cutoff}}$: {E_min:.2f} MeV\n"
                   f"Flux usage: {pct_avail:.2f}%")
    bottom_handles.append(Line2D([0], [0], color=case['color'], lw=2, linestyle=case['ls'],
                                 label=label_stats))

# --- Finalizing Top Panel ---
top_flux_handle = Line2D([0], [0], color='k', lw=2, label='Neutrino Flux')
ax1.legend(handles=[top_flux_handle] + top_handles, loc='upper right', frameon=True, fancybox=True)
ax1.set_ylim(bottom=0)

# --- Finalizing Bottom Panel ---
bottom_flux_handle = Line2D([0], [0], color='k', lw=2, label='Cumulated Neutrino Flux')
ax2.legend(handles=[bottom_flux_handle] + bottom_handles, loc='upper right', frameon=True, fancybox=True,
           bbox_to_anchor=(1.0, 1.0), borderpad=0.8, labelspacing=0.8)

# --- Set Axis Limits ---
ax2.set_xlim(0, 10)

# Set Y-limit for Log Scale
# Lower limit: 0.01% of total flux or similar small number to avoid log(0)
y_min_limit = total_flux * 1e-5
ax2.set_ylim(y_min_limit, total_flux * 5)

plt.tight_layout()
plt.savefig('neutrino_flux_usage_log.png', dpi=300, bbox_inches='tight')
print("Plot saved as 'neutrino_flux_usage_log.png'")
plt.show()
