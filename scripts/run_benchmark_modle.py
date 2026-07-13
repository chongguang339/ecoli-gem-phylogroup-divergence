import os
import glob
import time
import subprocess
import re
import csv
import zipfile

# ======================= 1. 严格对照：对齐您的集群路径与版本 =======================
BASE_DIR = "/home/u2024110181382/Project/04.pipeline/modle_compare"
INPUT_DIR = os.path.join(BASE_DIR, "00.data")

# 🌟【版本 1：原版默认】直接调用环境变量/Conda环境中的原始 gapseq 
GAPSEQ_DEFAULT = "gapseq"

# 🌟【版本 2：高效优化】统一使用你通过 Git 安装、并且修改了源码的绝对路径版
GAPSEQ_OPTIMIZED = "/home/u2024110181382/Project/02.E.coli/biosoft/gapseq/gapseq"

# 优化管线特有的依赖库与培养基
M9_CSV = "/home/u2024110181382/Project/02.E.coli/biosoft/gapseq/server_deploy/M9_glucose_medium.csv"
FASTDB = "/home/u2024110181382/Project/02.E.coli/biosoft/gapseq/enterobacterales_gapseq_fastdb.faa"

# ======================= 2. 初始化/读取追踪表格 =======================
csv_path = os.path.join(BASE_DIR, "benchmark_metrics.csv")
fieldnames = ["Accession", "Actual_N50", "FNA_Path", "Default_Time_Sec", "Optimized_Time_Sec", "Default_Reactions", "Optimized_Reactions"]

# --- 纯 Python 的 FASTA 读取和 N50 计算函数（彻底移除 Bio 依赖） ---
def calculate_n50_from_fasta(fna_path):
    lengths = []
    current_len = 0
    with open(fna_path, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('>'):
                if current_len > 0:
                    lengths.append(current_len)
                    current_len = 0
            else:
                current_len += len(line.strip())
        if current_len > 0:
            lengths.append(current_len)
            
    if not lengths:
        return 0
    lengths.sort(reverse=True)
    total_len = sum(lengths)
    cumulative_len = 0
    for l in lengths:
        cumulative_len += l
        if cumulative_len >= total_len / 2:
            return l
    return 0

rows = []
if os.path.exists(csv_path):
    print(f"📖 检测到已存在的基准表格，执行断点续传模式。")
    with open(csv_path, 'r', newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
else:
    print(f"🔍 未检测到表格，正在自动扫描 {INPUT_DIR} 下的子文件夹与压缩包...")
    fna_files_info = []
    
    for item in os.listdir(INPUT_DIR):
        item_path = os.path.join(INPUT_DIR, item)
        if os.path.isdir(item_path):
            accession = item
            zip_path = os.path.join(item_path, "ncbi_dataset.zip")
            extract_path = os.path.join(item_path, "ncbi_dataset_extracted")
            
            if os.path.exists(zip_path):
                if not os.path.exists(extract_path):
                    print(f"📦 正在为您解压 {accession} 的基因组数据包...")
                    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                        zip_ref.extractall(extract_path)
            
            search_path = os.path.join(extract_path, "**", "*_genomic.fna")
            found_fnas = glob.glob(search_path, recursive=True)
            
            if found_fnas:
                fna_files_info.append((accession, found_fnas[0]))
    
    if not fna_files_info:
        print("❌ 错误：未能在指定路径下找到解压的 *_genomic.fna 文件！请检查数据目录。")
        exit(1)
        
    for accession, fna in sorted(fna_files_info, key=lambda x: x[0]):
        n50 = calculate_n50_from_fasta(fna)
        rows.append({
            "Accession": accession,
            "Actual_N50": str(n50),
            "FNA_Path": os.path.abspath(fna),
            "Default_Time_Sec": "",
            "Optimized_Time_Sec": "",
            "Default_Reactions": "",
            "Optimized_Reactions": ""
        })
    
    with open(csv_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

def save_csv(path, data_rows):
    with open(path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(data_rows)

def count_sbml_reactions(folder_path):
    xml_files = glob.glob(os.path.join(folder_path, "*.xml"))
    if not xml_files:
        return 0
    try:
        with open(xml_files[0], 'r', encoding='utf-8') as f:
            content = f.read()
        return len(re.findall(r'<reaction\s', content))
    except Exception as e:
        print(f"⚠️ 解析模型反应数失败: {e}")
        return 0

def check_and_print_error_log(log_path):
    if os.path.exists(log_path):
        print(f"\n❌ [ERROR DETECTED] 检测到程序异常终止，正在为您提取错误日志 ({os.path.basename(log_path)})：")
        print("----------------------------------------------------------------------")
        try:
            with open(log_path, 'r', encoding='utf-8') as f:
                lines = f.readlines()
                for line in lines[-25:]:
                    print(f"  > {line.strip()}")
        except Exception as e:
            print(f"  读取日志内容失败: {e}")
        print("----------------------------------------------------------------------\n")

# ======================= 3. 循环迭代目标基因组 =======================
for index, row in enumerate(rows):
    accession = row['Accession']
    fna_path = row.get('FNA_Path', '').strip()
    actual_n50 = int(row['Actual_N50']) if row['Actual_N50'] else 0
    
    if not fna_path or not os.path.exists(fna_path):
        fallback_search = os.path.join(INPUT_DIR, accession, "ncbi_dataset_extracted", "**", "*_genomic.fna")
        fallback_files = glob.glob(fallback_search, recursive=True)
        if fallback_files:
            fna_path = fallback_files[0]
        else:
            print(f"⚠️ 警告：无法为样本 {accession} 恢复基因组文件路径，自动跳过。")
            continue

    fna_path = os.path.abspath(fna_path)
    row['FNA_Path'] = fna_path
    save_csv(csv_path, rows)

    sample_out_root = os.path.join(BASE_DIR, accession)
    sample_default_dir = os.path.join(sample_out_root, "results_default")
    sample_opt_dir = os.path.join(sample_out_root, "results_optimized")
    
    os.makedirs(sample_default_dir, exist_ok=True)
    os.makedirs(sample_opt_dir, exist_ok=True)
    
    print("\n" + "="*70)
    print(f"🔥 正在测试样本 [{index+1}/{len(rows)}]: {accession} (N50: {actual_n50:,} bp)")
    print(f"📍 基因组绝对路径: {fna_path}")
    print(f"📂 输出目录准备完毕: {sample_out_root}")
    print("="*70)
    
    # ------------------ 核心对照组 1：原版默认 gapseq 流程 ------------------
    has_default_xml = len(glob.glob(os.path.join(sample_default_dir, "*.xml"))) > 0
    is_valid_default = row.get('Default_Time_Sec') and row['Default_Time_Sec'].strip() != "" and int(row['Default_Time_Sec']) > 5
    
    if is_valid_default and has_default_xml:
        print(f"⏭️ [SKIP] 原版默认流程已成功跑完，自动跳过。")
    else:
        print(f"👉 [1/2] 启动【原版默认】管线中... (调用 Conda 环境原版，不加额外参数)")
        log_file = os.path.join(sample_default_dir, f"{accession}_default_run.log")
        cmd_default = f"{GAPSEQ_DEFAULT} doall {fna_path}"
        
        t0 = time.time()
        with open(log_file, "w") as f_log:
            subprocess.run(cmd_default, shell=True, cwd=sample_default_dir, executable='/bin/bash', stdout=f_log, stderr=f_log)
        t1 = time.time()
        
        elapsed = int(t1 - t0)
        rxn_count = count_sbml_reactions(sample_default_dir)
        
        row['Default_Time_Sec'] = str(elapsed)
        row['Default_Reactions'] = str(rxn_count)
        save_csv(csv_path, rows)
        print(f"✅ 默认流程完成：耗时 {elapsed} 秒 | 反应数: {rxn_count}")
        
        if elapsed <= 5:
            check_and_print_error_log(log_file)
        
    # ------------------ 核心实验组 2：高效优化管线 ------------------
    has_opt_xml = len(glob.glob(os.path.join(sample_opt_dir, "*.xml"))) > 0
    is_valid_opt = row.get('Optimized_Time_Sec') and row['Optimized_Time_Sec'].strip() != "" and int(row['Optimized_Time_Sec']) > 5
    
    if is_valid_opt and has_opt_xml:
        print(f"⏭️ [SKIP] 高效优化流程已成功跑完，自动跳过。")
    else:
        print(f"👉 [2/2] 启动【高效优化】管线中... (调用 Git 修改版绝对路径 + 组合加速参数)")
        log_file = os.path.join(sample_opt_dir, f"{accession}_opt_pipeline.log")
        cmd_opt = f"{GAPSEQ_OPTIMIZED} doall -A diamond -S identical -E {FASTDB} -m {M9_CSV} -f {sample_opt_dir} -K 4 {fna_path}"
        
        t0 = time.time()
        with open(log_file, "w") as f_log:
            subprocess.run(cmd_opt, shell=True, cwd=sample_opt_dir, executable='/bin/bash', stdout=f_log, stderr=f_log)
        t1 = time.time()
        
        elapsed = int(t1 - t0)
        rxn_count = count_sbml_reactions(sample_opt_dir)
        
        row['Optimized_Time_Sec'] = str(elapsed)
        row['Optimized_Reactions'] = str(rxn_count)
        save_csv(csv_path, rows)
        print(f"✅ 优化流程完成：耗时 {elapsed} 秒 | 反应数: {rxn_count}")
        
        if elapsed <= 5:
            check_and_print_error_log(log_file)

print("\n🎉 基准压力测试监控程序执行结束。")
