import numpy as np
import pandas as pd
import matplotlib

# [设置] Matplotlib 后端为 Agg (非交互模式)
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.integrate import simpson
from scipy.special import spherical_jn
import os
import json
import glob
from multiprocessing import Pool
import math

# --- PHYSICAL CONSTANTS ---
SIN2_THETA_W = 0.23867
G_F = 1.1663787e-11  # [MeV^-2] Fermi coupling constant
H_BAR_C_CM = 197.327e-13  # [MeV*cm]
H_BAR_C_FM = 197.327  # [MeV*fm] for form factor
SECONDS_PER_DAY = 86400.0
DAYS_PER_YEAR = 365.25
MEV_PER_KEV = 1000.0
AMU_MEV = 931.4941

# --- [配置] 靶材注册表 ---
TARGET_REGISTRY = {
    'Si': {
        'label': r'Si',
        'roi_kev': (0.1, 1.0),
        'color': 'tab:blue',
        'molar_mass_kg': 28.0855 / 1000.0,
        'elements': [{'Z': 14, 'A': 28, 'n': 1}]
    },
    'Ge': {
        'label': r'Ge',
        'roi_kev': (0.1, 1.0),
        'color': 'tab:green',
        'molar_mass_kg': 72.630 / 1000.0,
        'elements': [{'Z': 32, 'A': 73, 'n': 1}]
    },
    # 'Xe': {
    #     'label': r'Xe',
    #     'roi_kev': (0.63, 1.36),
    #     'color': 'tab:cyan',
    #     'molar_mass_kg': 131.293 / 1000.0,
    #     'elements': [{'Z': 54, 'A': 131, 'n': 1}]
    # },
    'CaWO4': {
        'label': r'CaWO$_4$',
        'roi_kev': (0.1, 1.0),
        'color': 'tab:purple',
        'molar_mass_kg': 287.92 / 1000.0,
        'elements': [{'Z': 20, 'A': 40, 'n': 1}, {'Z': 74, 'A': 184, 'n': 1}, {'Z': 8, 'A': 16, 'n': 4}]
    },
    'PbWO4': {
        'label': r'PbWO$_4$',
        'roi_kev': (0.1, 1.0),
        'color': 'tab:brown',
        'molar_mass_kg': 455.0 / 1000.0,
        'elements': [{'Z': 82, 'A': 208, 'n': 1}, {'Z': 74, 'A': 184, 'n': 1}, {'Z': 8, 'A': 16, 'n': 4}]
    },
    'Li2MoO4': {
        'label': r'Li$_2$MoO$_4$',
        'roi_kev': (0.1, 1.0),
        'color': 'tab:red',
        'molar_mass_kg': 173.83 / 1000.0,
        'elements': [{'Z': 3, 'A': 7, 'n': 2}, {'Z': 42, 'A': 96, 'n': 1}, {'Z': 8, 'A': 16, 'n': 4}]
    },
    # 'GaAs': {
    #     'label': r'GaAs',
    #     'roi_kev': (0.1, 1.0),
    #     'color': 'tab:orange',
    #     'molar_mass_kg': 144.64 / 1000.0,
    #     'elements': [{'Z': 31, 'A': 70, 'n': 1}, {'Z': 33, 'A': 75, 'n': 1}]
    # }
}


# --- PHYSICS HELPER FUNCTIONS ---

def get_target_nuclei_per_kg(target_key):
    props = TARGET_REGISTRY[target_key]
    return 6.02214076e23 / props['molar_mass_kg']


def get_atomic_masses(target_key):
    props = TARGET_REGISTRY[target_key]
    return [elem['A'] * AMU_MEV for elem in props['elements']]


def T_max_nucleus(E_nu, M_nucleus):
    return np.divide(2 * E_nu ** 2, (M_nucleus + 2 * E_nu), where=(E_nu > 0))


def helm_form_factor(Q_mev, A):
    Q_fm_inv = Q_mev / H_BAR_C_FM
    s = 0.9
    R_A = (1.23 * A ** (1 / 3) - 0.6)
    qr = Q_fm_inv * R_A
    with np.errstate(divide='ignore', invalid='ignore'):
        bessel_term = 3 * spherical_jn(n=1, z=qr) / qr
    bessel_term = np.where(np.abs(qr) < 1e-4, 1.0 - (qr ** 2) / 10.0, bessel_term)
    gaussian_term = np.exp(-(Q_fm_inv * s) ** 2 / 2)
    return bessel_term * gaussian_term


def calculate_dsigmadT_matrix(E_nu_grid, T_recoil_vals, Z, A):
    M_nucleus = A * AMU_MEV
    N = A - Z
    Q_w = N - (1 - 4 * SIN2_THETA_W) * Z

    T = T_recoil_vals.reshape(-1, 1)
    E = E_nu_grid.reshape(1, -1)

    E_min = 0.5 * T + 0.5 * np.sqrt(T ** 2 + 2 * M_nucleus * T)
    mask_kinematic = E > E_min

    Q = np.sqrt(2 * M_nucleus * T)
    F_Q2 = helm_form_factor(Q, A)

    prefactor = (G_F ** 2 * M_nucleus) / (4 * np.pi)

    kinematic_term = np.zeros_like(mask_kinematic, dtype=float)
    if np.any(mask_kinematic):
        term = 1 - (M_nucleus * T) / (2 * E ** 2)
        np.maximum(term, 0, out=term)
        kinematic_term = np.where(mask_kinematic, term, 0.0)

    return prefactor * (Q_w ** 2) * (F_Q2 ** 2) * kinematic_term


def calculate_event_rate_vectorized(T_recoil_mev, flux_function, E_nu_max, target_key):
    N_E_POINTS = 5000
    E_grid = np.linspace(1e-5, E_nu_max, N_E_POINTS)

    Flux_grid_row = flux_function(E_grid).reshape(1, -1)

    conversion = (H_BAR_C_CM ** 2) * get_target_nuclei_per_kg(
        target_key) * SECONDS_PER_DAY * DAYS_PER_YEAR / MEV_PER_KEV
    total_rate = np.zeros_like(T_recoil_mev)
    props = TARGET_REGISTRY[target_key]

    for elem in props['elements']:
        Z, A, n_atoms = elem['Z'], elem['A'], elem['n']
        dSigma_matrix = calculate_dsigmadT_matrix(E_grid, T_recoil_mev, Z, A)
        integrand = Flux_grid_row * dSigma_matrix
        integral_results = simpson(y=integrand, x=E_grid, axis=1)
        total_rate += n_atoms * integral_results

    return total_rate * conversion


# --- PLOTTING FUNCTIONS ---

def setup_plot_style():
    """配置全局绘图样式"""
    plt.style.use('seaborn-v0_8-paper')
    plt.rcParams.update({
        'font.family': 'serif',
        'mathtext.fontset': 'cm',
        'font.size': 14,
        'axes.labelsize': 16,
        'axes.titlesize': 16,
        'xtick.labelsize': 14,
        'ytick.labelsize': 14,
        'legend.fontsize': 13,
        'axes.grid': True,
        'grid.alpha': 0.5,
        'grid.linestyle': '--'
    })


def plot_flux_standalone(E_nu_vals, phi_vals, flux_function, task_info):
    """图1：单独的中微子通量"""
    setup_plot_style()
    source_title = task_info['title']

    fname, fext = os.path.splitext(os.path.basename(task_info['save_name_base']))
    save_dir = os.path.dirname(task_info['save_name_base'])
    full_save_path = os.path.join(save_dir, f"{fname}_Flux{fext}")

    E_nu_max = E_nu_vals[-1]
    E_plot = np.linspace(0, E_nu_max, 500)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(E_nu_vals, phi_vals, 'o', color='gray', alpha=0.6, label='Data', markersize=5)
    ax.plot(E_plot, flux_function(E_plot), 'k-', lw=2.5, label='Interpolation')

    # idx_max = np.argmax(phi_vals)
    # E_peak = E_nu_vals[idx_max]
    # phi_peak = phi_vals[idx_max]
    # ax.text(E_peak, phi_peak * 1.05, f' Peak E: {E_peak:.2f} MeV',
    #         color='tab:red', fontsize=12, fontweight='bold', ha='center')

    ax.set_title(f"Flux: {source_title}", fontsize=18)
    ax.set_xlabel(r"$E_\nu$ [MeV]")
    ax.set_ylabel(r"Flux [$\nu$ / cm$^2$ / s / MeV]")
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

    xlim = task_info.get('flux_x_limit', E_nu_max * 1.1)
    ax.set_xlim(0, xlim)
    ax.set_ylim(bottom=0)

    plt.tight_layout()
    try:
        plt.savefig(full_save_path, dpi=300)
        print(f"    [Flux] Saved: {full_save_path}")
    except Exception as e:
        print(f"    [Error] Saving Flux plot: {e}")
    plt.close(fig)


def plot_rates_multipanel(E_nu_max, flux_function, task_info):
    """图2：多面板分别展示详细信息（含 ROI 区域）"""
    setup_plot_style()
    source_title = task_info['title']

    fname, fext = os.path.splitext(os.path.basename(task_info['save_name_base']))
    save_dir = os.path.dirname(task_info['save_name_base'])
    full_save_path = os.path.join(save_dir, f"{fname}_Rates_Panels{fext}")

    targets = list(TARGET_REGISTRY.keys())
    n_targets = len(targets)
    cols = 3
    rows = math.ceil(n_targets / cols)

    fig, axes = plt.subplots(rows, cols, figsize=(cols * 5.5, rows * 4.5))
    fig.suptitle(f"Event Rates Breakdown: {source_title}", fontsize=20, y=0.98)

    axes_flat = axes.flatten()

    for i, ax in enumerate(axes_flat):
        if i >= n_targets:
            ax.axis('off')
            continue

        target_name = targets[i]
        props = TARGET_REGISTRY[target_name]

        masses = get_atomic_masses(target_name)
        T_max_abs_kev = T_max_nucleus(E_nu_max, min(masses)) * MEV_PER_KEV
        T_kev = np.logspace(np.log10(1e-3), np.log10(T_max_abs_kev), 300)
        T_mev = T_kev / MEV_PER_KEV
        rate = calculate_event_rate_vectorized(T_mev, flux_function, E_nu_max, target_name)

        roi_min, roi_max = props['roi_kev']
        mask_roi = (T_kev >= roi_min) & (T_kev <= roi_max)
        total_counts = simpson(y=rate[mask_roi], x=T_kev[mask_roi]) if np.any(mask_roi) else 0.0

        color = props.get('color', 'black')
        ax.loglog(T_kev, rate, color=color, lw=2.5)
        ax.fill_between(T_kev, rate, 0, where=mask_roi, color=color, alpha=0.25)

        info_str = (f"ROI: {roi_min}-{roi_max} keV\n"
                    f"Rate: {total_counts:.2f}")
        ax.text(0.96, 0.96, info_str, transform=ax.transAxes,
                ha='right', va='top', fontsize=12,
                bbox=dict(facecolor='white', alpha=0.9, boxstyle='round,pad=0.3'))

        ax.set_title(props['label'], fontsize=16, fontweight='bold')
        ax.set_xlabel("Recoil Energy [keV]")
        if i % cols == 0:
            ax.set_ylabel("Rate [cts/kg/yr/keV]")

        ax.set_xlim(1e-2, T_max_abs_kev * 1.5)
        if np.max(rate) > 0:
            ax.set_ylim(bottom=np.max(rate) * 1e-5)

    plt.tight_layout(rect=[0, 0.02, 1, 0.96])
    try:
        plt.savefig(full_save_path, dpi=300)
        print(f"    [Panel] Saved: {full_save_path}")
    except Exception as e:
        print(f"    [Error] Saving Panel plot: {e}")
    plt.close(fig)


def plot_rates_overlay(E_nu_max, flux_function, task_info):
    """图3：[叠加] 所有材料的 Rate 叠加在同一张图，包含阴影 Highlight"""
    setup_plot_style()
    source_title = task_info['title']

    fname, fext = os.path.splitext(os.path.basename(task_info['save_name_base']))
    save_dir = os.path.dirname(task_info['save_name_base'])
    full_save_path = os.path.join(save_dir, f"{fname}_Rates_Overlay{fext}")

    fig, ax = plt.subplots(figsize=(10, 8))

    targets = list(TARGET_REGISTRY.keys())
    max_T_endpoint = 0

    for target_name in targets:
        props = TARGET_REGISTRY[target_name]
        masses = get_atomic_masses(target_name)
        T_max_abs_kev = T_max_nucleus(E_nu_max, min(masses)) * MEV_PER_KEV
        max_T_endpoint = max(max_T_endpoint, T_max_abs_kev)

        T_kev = np.logspace(np.log10(1e-3), np.log10(T_max_abs_kev), 400)
        T_mev = T_kev / MEV_PER_KEV

        rate = calculate_event_rate_vectorized(T_mev, flux_function, E_nu_max, target_name)

        roi_min, roi_max = props['roi_kev']
        mask_roi = (T_kev >= roi_min) & (T_kev <= roi_max)

        if np.any(mask_roi):
            total_counts = simpson(y=rate[mask_roi], x=T_kev[mask_roi])
        else:
            total_counts = 0.0

        # 1. 绘制曲线
        label_str = f"{props['label']} : {total_counts:.2f} cts/kg/yr"
        ax.loglog(T_kev, rate, color=props['color'], lw=3, label=label_str)

        # 2. [新增] 绘制阴影 Highlight
        # alpha 设置得比 Panel 图低一些，防止叠加时颜色太深
        ax.fill_between(T_kev, rate, 0, where=mask_roi,
                        color=props['color'], alpha=0.15)

    ax.set_title(f"CEvNS Event Rates Comparison\nSource: {source_title}", fontsize=20, pad=15)
    ax.set_xlabel("Recoil Energy T [keV]", fontsize=18)
    ax.set_ylabel(r"Differential Event Rate [counts $\cdot$ kg$^{-1} \cdot$ yr$^{-1} \cdot$ keV$^{-1}$]", fontsize=18)

    ax.set_xlim(1e-2, max_T_endpoint * 1.5)
    ax.set_ylim(bottom=1e-3)

    ax.legend(loc='best', fontsize=14, frameon=True, fancybox=True, framealpha=0.9, shadow=True)

    ax.tick_params(which='major', length=8, width=1.5)
    ax.tick_params(which='minor', length=4, width=1)

    plt.tight_layout()
    try:
        plt.savefig(full_save_path, dpi=300)
        print(f"    [Overlay] Saved: {full_save_path}")
    except Exception as e:
        print(f"    [Error] Saving Overlay plot: {e}")
    plt.close(fig)


def process_source_group(task_info):
    file_path = task_info['flux_file_path']
    print(f"--> [PID {os.getpid()}] Processing: {task_info['title']}")

    try:
        real_path = "../" + file_path if not os.path.exists(file_path) else file_path
        flux_data = pd.read_csv(real_path, comment='#')
        flux_data = flux_data.sort_values(by=flux_data.columns[0]).reset_index(drop=True)
        E_nu_vals, phi_vals = flux_data.iloc[:, 0].values, flux_data.iloc[:, 1].values

        flux_function = interp1d(E_nu_vals, phi_vals, kind='cubic', bounds_error=False, fill_value=0.0)
        E_nu_max = E_nu_vals[-1]
    except Exception as e:
        print(f"    [Error] Failed to load flux {file_path}: {e}")
        return

    # 1. Flux Plot
    plot_flux_standalone(E_nu_vals, phi_vals, flux_function, task_info)

    # 2. Panel Plot
    plot_rates_multipanel(E_nu_max, flux_function, task_info)

    # 3. Overlay Plot (With Highlights)
    plot_rates_overlay(E_nu_max, flux_function, task_info)


# --- MAIN EXECUTION ---
if __name__ == "__main__":
    print("Searching for setting.json files...")
    setting_files = glob.glob("../data/**/setting.json", recursive=True)

    if not setting_files:
        print("No setting.json files found.")
        exit()

    unique_tasks = {}
    for s_file in setting_files:
        try:
            with open(s_file, 'r', encoding='utf-8') as f:
                configs = json.load(f)
                if not isinstance(configs, list): configs = [configs]
                for cfg in configs:
                    key = (cfg['flux_file_path'], cfg['title'])
                    if key not in unique_tasks:
                        unique_tasks[key] = {
                            'flux_file_path': cfg['flux_file_path'],
                            'title': cfg['title'],
                            'save_name_base': cfg['figure_save_name'],
                            'plot_flux_x_limit': cfg.get('plot_flux_x_limit', None)
                        }
        except Exception:
            pass

    task_list = list(unique_tasks.values())
    print(f"Found {len(task_list)} unique Source configurations.")
    print(f"Targets: {list(TARGET_REGISTRY.keys())}")
    print("-" * 40)

    with Pool(processes=min(len(task_list), 8)) as pool:
        pool.map(process_source_group, task_list)

    print("-" * 40)
    print("All jobs completed.")