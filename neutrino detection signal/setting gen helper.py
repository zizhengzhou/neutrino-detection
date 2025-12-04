import os
import json

# --- 配置 ---
DATA_ROOT = "data"
OUTPUT_DIR = "output"


# 默认参数猜测逻辑 (保持不变)
def guess_params(folder_name, file_name):
    lower_name = file_name.lower()
    lower_folder = folder_name.lower()

    # 1. 猜测中微子类型
    nu_type = "e"
    if "anti" in lower_name or "bar" in lower_name:
        if "mu" in lower_name:
            nu_type = "mu_bar"
        elif "tau" in lower_name:
            nu_type = "tau_bar"
        else:
            nu_type = "e_bar"
    else:
        if "mu" in lower_name:
            nu_type = "mu"
        elif "tau" in lower_name:
            nu_type = "tau"
        else:
            nu_type = "e"

    # 2. 猜测画图 X 轴上限 (MeV)
    if lower_folder in ["daya bay"]:
        x_limit = 10.0
    elif lower_folder in ["csns"]:
        x_limit = 50.0
    elif lower_folder in ["t2k"]:
        x_limit = 15000.0
    elif lower_folder in ["test"]:
        x_limit = 20.0
    else:
        x_limit = 15.0

    return nu_type, x_limit


def generate_json():
    # 确保 output 目录存在
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    if not os.path.exists(DATA_ROOT):
        print(f"Error: 目录 '{DATA_ROOT}' 不存在。")
        return

    for folder_name in os.listdir(DATA_ROOT):
        folder_path = os.path.join(DATA_ROOT, folder_name)

        if not os.path.isdir(folder_path):
            continue

        configs = []
        files = sorted(os.listdir(folder_path))

        for file_name in files:
            if not file_name.endswith(".csv"):
                continue

            # --- 特殊过滤: Solar 只取一个 ---
            if folder_name == "solar" and file_name != "equal_flavor_solar_flux.csv":
                continue

            # 获取猜测参数
            nu_type, x_limit = guess_params(folder_name, file_name)

            # --- 核心修改部分 ---

            # 1. 原始文件名（去掉后缀）
            raw_name = os.path.splitext(file_name)[0]

            # 2. 处理 TITLE：替换下划线为空格，移除前缀，防止 LaTeX 报错
            # 例如: "equal_flavor_solar_flux" -> "equal flavor solar flux"
            title_folder = folder_name.replace("_", " ")
            title_file = raw_name.replace("_", " ")

            # 最终 Title 格式: "solar - equal flavor solar flux (e)"
            title = f"{title_folder} - {title_file} ({nu_type})"

            # 3. 处理文件保存路径：文件名里最好用下划线，避免空格
            save_name_clean = raw_name.replace(" ", "_")
            save_folder_clean = folder_name.replace(" ", "_")

            file_path = f"{DATA_ROOT}/{folder_name}/{file_name}"
            save_path = f"{OUTPUT_DIR}/NMM_{save_folder_clean}_{save_name_clean}.png"

            config_entry = {
                "flux_file_path": file_path,
                "figure_save_name": save_path,
                "title": title,
                "neutrino_type": nu_type,
                "plot_flux_x_limit": x_limit
            }
            configs.append(config_entry)

        if configs:
            json_path = os.path.join(folder_path, "setting_disable.json")
            with open(json_path, 'w', encoding='utf-8') as f:
                json.dump(configs, f, indent=4, ensure_ascii=False)
            print(f"[OK] Generated {json_path}")
        else:
            print(f"[SKIP] No matching CSVs in {folder_name}")


if __name__ == "__main__":
    generate_json()