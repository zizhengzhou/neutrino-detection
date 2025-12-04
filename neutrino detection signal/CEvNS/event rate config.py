import numpy as np
import pandas as pd
import matplotlib

# [修改 1] 设置 Matplotlib 后端为 Agg (非交互模式)
# 必须在 import pyplot 之前设置，防止多进程绘图时 GUI 冲突
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.integrate import quad, simpson
from scipy.special import spherical_jn
import os
import json
import glob

# [修改 2] 引入多进程库
from multiprocessing import Pool, cpu_count

# --- PHYSICAL CONSTANTS ---
TARGET_MATERIAL = 'Si'
ROI_MIN_KEV = 0.010  # [keV]
ROI_MAX_KEV = 0.100  # [keV]

SIN2_THETA_W = 0.23867
G_F = 1.1663787e-11  # [MeV^-2] Fermi coupling constant
H_BAR_C_CM = 197.327e-13  # [MeV*cm]
H_BAR_C_FM = 197.327  # [MeV*fm] for form factor
SECONDS_PER_DAY = 86400.0
DAYS_PER_YEAR = 365.25
MEV_PER_KEV = 1000.0
AMU_MEV = 931.4941

TARGET_PROPERTIES = {
    'Si': {'Z': 14, 'A': 28, 'molar_mass_kg': 28.0855 / 1000.0},
    'Xe': {'Z': 54, 'A': 131, 'molar_mass_kg': 131.293 / 1000.0},
    'Ge': {'Z': 32, 'A': 73, 'molar_mass_kg': 72.630 / 1000.0},
    # 砷化镓 (Gallium Arsenide) - 等效处理
    # GaAs 是化合物，Ga(Z=31, A~69.7) 和 As(Z=33, A~74.9)
    # 平均 Z = 32, 平均质量 = (69.723 + 74.922)/2 = 72.3225
    'GaAs': {'Z': 32, 'A': 72, 'molar_mass_kg': 72.3225 / 1000.0}
}


# --- PHYSICS FUNCTIONS (CEvNS Specific) ---

def get_target_nuclei_per_kg(material):
    """Calculates number of nuclei per kg for the target."""
    props = TARGET_PROPERTIES.get(material)
    return 6.02214076e23 / props['molar_mass_kg']


def get_nucleus_mass_mev(material):
    """Returns nucleus mass in MeV."""
    props = TARGET_PROPERTIES.get(material)
    return props['A'] * AMU_MEV


def T_max_nucleus(E_nu, M_nucleus):
    """Calculates maximum kinematic recoil energy."""
    return (2 * E_nu ** 2) / (M_nucleus + 2 * E_nu)


def helm_form_factor(Q, A):
    """Calculates Helm Form Factor F(Q^2)."""
    Q_fm_inv = Q / H_BAR_C_FM
    s = 0.9
    R_A = (1.23 * A ** (1 / 3) - 0.6)
    qr = Q_fm_inv * R_A
    if qr == 0: return 1.0
    if abs(qr) < 1e-4:
        # 使用泰勒展开近似: 1 - x^2/10
        bessel_term = 1.0 - (qr ** 2) / 10.0
    else:
        bessel_term = 3 * spherical_jn(n=1, z=qr, derivative=False) / qr
    gaussian_term = np.exp(-(Q_fm_inv * s) ** 2 / 2)
    return bessel_term * gaussian_term


def calculate_dsigmadT_cevns(E_nu, T, material):
    """Calculates differential cross-section dSigma/dT [MeV^-3]."""
    M_nucleus = get_nucleus_mass_mev(material)
    if T <= 0 or T > T_max_nucleus(E_nu, M_nucleus):
        return 0.0
    props = TARGET_PROPERTIES.get(material)
    Z, N = props['Z'], props['A'] - props['Z']
    Q_w = N - (1 - 4 * SIN2_THETA_W) * Z
    Q = np.sqrt(2 * M_nucleus * T)
    F_Q2 = helm_form_factor(Q, props['A'])
    prefactor = (G_F ** 2 * M_nucleus) / (4 * np.pi)
    kinematic_term = (1 - (M_nucleus * T) / (2 * E_nu ** 2))
    return prefactor * (Q_w ** 2) * (F_Q2 ** 2) * kinematic_term


def calculate_total_sigma_at_E(E_nu, material):
    """Calculates total cross-section sigma_tot [cm^2] at a specific E_nu."""
    M_nucleus = get_nucleus_mass_mev(material)
    T_max = T_max_nucleus(E_nu, M_nucleus)
    integrand = lambda T: calculate_dsigmadT_cevns(E_nu, T, material)
    try:
        integral_val, _ = quad(integrand, 0, T_max, limit=100)
    except:
        integral_val = 0.0
    sigma_cm2 = integral_val * (H_BAR_C_CM ** 2)
    return sigma_cm2


def calculate_event_rate(T_recoil_mev, flux_function, E_nu_max, material):
    """Calculates differential event rate [counts/kg/year/keV]."""
    rate = np.zeros_like(T_recoil_mev)
    conversion = (H_BAR_C_CM ** 2) * get_target_nuclei_per_kg(material) * SECONDS_PER_DAY * DAYS_PER_YEAR / MEV_PER_KEV
    M_nucleus = get_nucleus_mass_mev(material)

    for i, T in enumerate(T_recoil_mev):
        # Kinematics: E_nu must be large enough to produce recoil T
        E_nu_min = np.sqrt(M_nucleus * T / 2)
        # Refined E_nu_min calculation to match T_max inverse more accurately if needed,
        # but sqrt(M*T/2) is the standard approximation for E << M.
        # Exact inverse: E = T/2 + sqrt(T^2/4 + M*T/2)
        E_nu_min_exact = 0.5 * T + np.sqrt(0.25 * T ** 2 + M_nucleus * T / 2)

        if E_nu_min_exact >= E_nu_max:
            continue

        integrand = lambda E: flux_function(E) * calculate_dsigmadT_cevns(E, T, material)
        try:
            integral, _ = quad(integrand, E_nu_min_exact, E_nu_max, limit=100)
        except:
            integral = 0.0
        rate[i] = integral * conversion
    return rate


# --- CORE LOGIC WRAPPER ---
def run_analysis(config):
    # 1. 路径修正：添加 "../"
    file_path = "../" + config['flux_file_path']

    # 2. 保存文件名修正：把 "NMM" 替换为 "CEvNS" 以免覆盖
    save_name = config['figure_save_name'].replace("NMM", "CEvNS")

    # 3. 标题修正
    raw_title = config['title']
    title_str = f"CEvNS : {raw_title}"

    flux_x_limit = config.get('plot_flux_x_limit', 15.0)

    # 打印进程信息
    print(f"--> [PID {os.getpid()}] Processing: {file_path}")

    # Load Data
    try:
        flux_data = pd.read_csv(file_path, comment='#')
        flux_data = flux_data.sort_values(by=flux_data.columns[0]).reset_index(drop=True)
        E_nu_data, phi_data = flux_data.iloc[:, 0].values, flux_data.iloc[:, 1].values
    except Exception as e:
        print(f"    [Error] Reading CSV file '{file_path}' failed: {e}")
        return

    # Interpolate
    flux_function = interp1d(E_nu_data, phi_data, kind='cubic', bounds_error=False, fill_value=0.0)

    # Physical parameters
    M_NUCLEUS = get_nucleus_mass_mev(TARGET_MATERIAL)
    E_nu_max_data = E_nu_data[-1]
    E_nu_peak = E_nu_data[np.argmax(phi_data)]

    # Kinematic limits
    T_max_abs_kev = T_max_nucleus(E_nu_max_data, M_NUCLEUS) * MEV_PER_KEV
    T_peak_limit_kev = T_max_nucleus(E_nu_peak, M_NUCLEUS) * MEV_PER_KEV

    # Calculate Event Rates (Fig 4 Data)
    T_recoil_kev = np.logspace(np.log10(1e-3), np.log10(T_max_abs_kev * 1.1), 300)
    T_recoil_mev = T_recoil_kev / MEV_PER_KEV

    rate_cevns_yearly = calculate_event_rate(T_recoil_mev, flux_function, E_nu_max_data, TARGET_MATERIAL)

    # ROI Integration
    mask_roi = (T_recoil_kev >= ROI_MIN_KEV) & (T_recoil_kev <= ROI_MAX_KEV)
    if np.any(mask_roi):
        T_roi = T_recoil_kev[mask_roi]
        Rate_roi = rate_cevns_yearly[mask_roi]
        total_counts_roi = simpson(y=Rate_roi, x=T_roi)
    else:
        total_counts_roi = 0.0

    # Calculate Total Cross Section (Fig 3 Data)
    # This is specific to CEvNS script, distinct from NMM script
    E_sigma_plot = np.linspace(0.1, E_nu_max_data, 100)
    sigma_total_vals = np.array([calculate_total_sigma_at_E(E, TARGET_MATERIAL) for E in E_sigma_plot])
    sigma_at_peak = calculate_total_sigma_at_E(E_nu_peak, TARGET_MATERIAL)

    # Calculate Diff Cross Section at Max Energy (Fig 2 Data)
    sigma_conversion = (H_BAR_C_CM ** 2) / MEV_PER_KEV
    diff_sigma_max = np.array(
        [calculate_dsigmadT_cevns(E_nu_max_data, T, TARGET_MATERIAL) for T in T_recoil_mev]) * sigma_conversion

    # --- PLOTTING ---
    plt.style.use('seaborn-v0_8-talk')

    plt.rcParams.update({
        'font.family': 'serif',
        'mathtext.fontset': 'cm',
        'font.size': 16,
        'axes.labelsize': 18,
        'axes.titlesize': 16,  # 子图标题
        'xtick.labelsize': 16,
        'ytick.labelsize': 16,
        'legend.fontsize': 15,
        'figure.titlesize': 24
    })

    # A4 适配尺寸
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(title_str)

    ax_flux = axes[0, 0]
    ax_diff_sigma = axes[0, 1]
    ax_tot_sigma = axes[1, 0]
    ax_rate = axes[1, 1]

    # === FIG 1: Flux ===
    E_plot = np.linspace(0, E_nu_max_data, 500)
    ax_flux.plot(E_nu_data, phi_data, 'o', label='Data Points', markersize=6)
    ax_flux.plot(E_plot, flux_function(E_plot), '-', label='Interpolation', lw=3)

    # 恢复 Peak 标记
    ax_flux.axvline(E_nu_peak, color='k', ls=':', label=f'Peak $E_\\nu$ ({E_nu_peak:.2f} MeV)')

    x_max_limit = flux_x_limit if flux_x_limit else E_nu_max_data * 1.1
    ax_flux.set(xlim=(0, x_max_limit),
                xlabel='Neutrino Energy $E_\\nu$ [MeV]',
                ylabel=r'Flux [$\nu$ / cm$^2$ / s / MeV]',
                title='Input Neutrino Flux')

    # 强制科学计数法
    ax_flux.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    ax_flux.legend(loc='upper right')
    ax_flux.grid(True)

    # === FIG 2: Diff Sigma ===
    ax_diff_sigma.loglog(T_recoil_kev, diff_sigma_max, lw=3, color='tab:red')
    # 恢复 Cutoff 标记
    ax_diff_sigma.axvline(T_peak_limit_kev, color='k', ls=':', label=f'Peak Flux ({T_peak_limit_kev:.1e} keV)')

    ax_diff_sigma.set(xlabel='Recoil Energy T [keV]',
                      ylabel=r'd$\sigma$/dT [cm$^2$/keV]',
                      title='Differential Cross Section')

    if np.any(diff_sigma_max > 0):
        ax_diff_sigma.set_ylim(bottom=diff_sigma_max[diff_sigma_max > 0].min() / 2)

    ax_diff_sigma.grid(True, which='both')
    ax_diff_sigma.legend()

    # === FIG 3: Total Sigma (CEvNS Specific) ===
    ax_tot_sigma.plot(E_sigma_plot, sigma_total_vals, lw=3, color='tab:purple')

    # Annotate Peak

    # === 原来的代码 ===
    # ax_tot_sigma.plot(E_nu_peak, sigma_at_peak, 'ko')
    # ax_tot_sigma.text(E_nu_peak, sigma_at_peak * 1.3,
    #                   f'$\\sigma_{{tot}}(E_{{peak}}) \\approx {sigma_at_peak:.1e}$ cm$^2$',
    #                   fontsize=14, ha='right', bbox=dict(facecolor='white', alpha=0.7))
    #
    # === 修改后的自适应代码 ===
    ax_tot_sigma.plot(E_nu_peak, sigma_at_peak, 'ko')  # 画黑点

    # 1. 获取当前 X 轴的范围
    x_min, x_max = ax_tot_sigma.get_xlim()

    # 2. 计算 Peak 位置在图中的相对比例 (0.0 - 1.0)
    relative_pos = (E_nu_peak - x_min) / (x_max - x_min)

    # 3. 智能判断对齐方式
    if relative_pos > 0.6:
        # 如果点在右边 (超过60%的位置)，文字就要“右对齐”，往左延伸，避免超出右边界
        my_ha = 'right'
        my_offset = (-10, 10)  # (x偏移, y偏移) 单位是点(points)，向左10，向上10
    else:
        # 如果点在左边，文字“左对齐”，往右延伸
        my_ha = 'left'
        my_offset = (10, 10)  # 向右10，向上10

    # 4. 使用 annotate 进行标注
    annotation_text = f'$\\sigma_{{tot}}(E_{{peak}}) \\approx {sigma_at_peak:.1e}$ cm$^2$'

    ax_tot_sigma.annotate(annotation_text,
                          xy=(E_nu_peak, sigma_at_peak),  # 数据点坐标
                          xytext=my_offset,  # 文字偏移量
                          textcoords='offset points',  # 偏移量单位：点 (不受数据坐标系影响)
                          ha=my_ha,  # 智能水平对齐
                          va='bottom',  # 垂直对齐：文字在点上方
                          fontsize=14,
                          bbox=dict(facecolor='white', alpha=0.7, boxstyle='round,pad=0.3'))
    ax_tot_sigma.set(xlabel='Neutrino Energy $E_\\nu$ [MeV]',
                     ylabel=r'Total Cross Section $\sigma$ [cm$^2$]',
                     title=r'Total Cross Section $\sigma_{tot}(E_\nu)$',
                     xlim=(0, x_max_limit))
    ax_tot_sigma.set_yscale('log')
    ax_tot_sigma.grid(True, which='both')

    # === FIG 4: Rate ===
    ax_rate.loglog(T_recoil_kev, rate_cevns_yearly, lw=3, color='tab:green')

    # ROI Shading
    ax_rate.fill_between(T_recoil_kev, rate_cevns_yearly, 0,
                         where=mask_roi, color='tab:green', alpha=0.3)

    text_str = (f"ROI: {ROI_MIN_KEV}-{ROI_MAX_KEV} keV\n"
                f"Total: {total_counts_roi:.1f} cts/kg/yr")

    ax_rate.text(0.05, 0.1, text_str,
                 transform=ax_rate.transAxes, fontsize=14, verticalalignment='bottom',
                 bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    ax_rate.set(xlabel='Recoil Energy T [keV]',
                ylabel='Event Rate [counts / kg / year / keV]',
                title='Expected Event Rate')

    valid_rates = rate_cevns_yearly[rate_cevns_yearly > 0]
    if len(valid_rates) > 0:
        ax_rate.set_ylim(bottom=valid_rates.min(), top=valid_rates.max() * 5)

    # ax_rate.legend(loc='upper right')
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

    total_jobs = len(all_tasks)
    num_processes = cpu_count() - 1

    print(f"Total CEvNS jobs to process: {total_jobs}")
    print(f"Starting parallel execution with {num_processes} processes...")
    print("-" * 40)

    # [修改 3] 使用 Pool 进行并行计算
    # 使用上下文管理器自动关闭 Pool
    with Pool(processes=num_processes) as pool:
        # map 将列表中的任务分配给不同的进程
        pool.map(run_analysis, all_tasks)

    print("-" * 40)
    print("All jobs completed.")