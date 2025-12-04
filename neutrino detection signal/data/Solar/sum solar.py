import pandas as pd
import numpy as np
from scipy.interpolate import PchipInterpolator
import matplotlib.pyplot as plt
import glob
import os

# --- Configuration ---
DATA_DIR = './'  # 数据所在的文件夹
OUTPUT_FILENAME = 'equal_flavor_solar_flux.csv'  # 输出文件名
OUTPUT_PATH = os.path.join(DATA_DIR, OUTPUT_FILENAME)

FIXED_DIVISOR = 3
NUM_POINTS_IN_GRID = 2000  # 增加网格点数以提高精度


# --- Helper: Safe Interpolation ---
def safe_interpolate(x_new, x_data, y_data):
    """
    使用 PCHIP 插值 (单调三次插值)，避免出现负值振荡。
    同时处理边界情况：超出原始数据范围的部分填 0。
    """
    # 1. 创建插值器
    interpolator = PchipInterpolator(x_data, y_data)

    # 2. 计算插值
    y_new = interpolator(x_new)

    # 3. 处理边界 (Extrapolation handling)
    # Pchip 默认会外推，我们需要手动把超出 x_data 范围的部分设为 0
    mask_out_of_bounds = (x_new < x_data.min()) | (x_new > x_data.max())
    y_new[mask_out_of_bounds] = 0.0

    # 4. 物理约束：通量不能为负
    y_new = np.clip(y_new, 0.0, None)

    return y_new


# --- Main Script ---

# 1. 查找文件 (Search Pattern)
search_pattern = os.path.join(DATA_DIR, '*.csv')
all_files = glob.glob(search_pattern)

# 2. 过滤掉输出文件本身 (Exclusion Logic)
input_files = [f for f in all_files if os.path.basename(f) != OUTPUT_FILENAME]

if not input_files:
    print(f"Error: No input files found in '{DATA_DIR}'.")
    print(f"Make sure you are running this script from the correct directory.")
else:
    print(f"Target Directory: {DATA_DIR}")
    print(f"Output File will be: {OUTPUT_PATH}")
    print(f"Found {len(input_files)} input files to process.")

    # 存储读取的数据，避免重复读取
    loaded_data = []
    global_min_energy = float('inf')
    global_max_energy = float('-inf')

    plt.figure(figsize=(10, 7))

    # 3. 读取并分析范围
    for filepath in input_files:
        try:
            # 读取 CSV (假设第一列是能量，第二列是通量)
            df = pd.read_csv(filepath, comment='#')
            # 排序并去空值
            df = df.sort_values(by=df.columns[0]).dropna()

            energy = df.iloc[:, 0].values
            flux = df.iloc[:, 1].values

            # 简单的单位清洗（防止数据里全是0导致报错）
            if len(energy) < 2:
                print(f"Skipping {os.path.basename(filepath)}: Not enough data points.")
                continue

            # 更新全局能量范围
            global_min_energy = min(global_min_energy, energy[0])
            global_max_energy = max(global_max_energy, energy[-1])

            loaded_data.append({
                'name': os.path.basename(filepath),
                'energy': energy,
                'flux': flux
            })

            # 画个图预览一下
            plt.plot(energy, flux, '.', markersize=2, alpha=0.5, label=f'Raw: {os.path.basename(filepath)}')

        except Exception as e:
            print(f"Error reading {os.path.basename(filepath)}: {e}")

    if loaded_data:
        print(f"\nGlobal Energy Range: {global_min_energy:.4f} - {global_max_energy:.4f} MeV")

        # 4. 创建统一网格
        common_energy_grid = np.linspace(global_min_energy, global_max_energy, NUM_POINTS_IN_GRID)
        summed_flux = np.zeros_like(common_energy_grid)

        # 5. 插值并求和
        print("Interpolating and Summing...")
        for item in loaded_data:
            interpolated_flux = safe_interpolate(common_energy_grid, item['energy'], item['flux'])
            summed_flux += interpolated_flux

        # 6. 除以 3 (Equal Flavor Assumption)
        final_flux = summed_flux / FIXED_DIVISOR

        # 7. 保存结果
        output_df = pd.DataFrame({
            'Energy_MeV': common_energy_grid,
            'Flux': final_flux
        })

        # 去掉通量为 0 的行（通常是超出所有谱范围的高能或低能端）
        # output_df = output_df[output_df['Flux'] > 0]
        # 注：为了保持网格完整性，有时候保留0也可以，这里选择保留以便画图连续

        output_df.to_csv(OUTPUT_PATH, index=False, float_format='%.6e')
        print(f"-" * 30)
        print(f"SUCCESS! File saved to: {OUTPUT_PATH}")
        print(f"-" * 30)

        # 8. 最终画图确认
        plt.plot(common_energy_grid, final_flux, 'k-', linewidth=2, label=f'Total / {FIXED_DIVISOR} (Result)')

        plt.xlabel('Energy [MeV]')
        plt.ylabel('Flux [cm$^{-2}$ s$^{-1}$ MeV$^{-1}$]')
        plt.title('Solar Neutrino Flux Summation')
        plt.yscale('log')
        plt.legend(loc='best', fontsize='small', ncol=2)
        plt.grid(True, which='both', alpha=0.3)

        # 限制一下Y轴下限，防止log图因为0报错
        valid_flux = final_flux[final_flux > 0]
        if len(valid_flux) > 0:
            plt.ylim(bottom=valid_flux.min() * 0.5)

        plt.tight_layout()
        plt.show()

    else:
        print("No valid data loaded.")