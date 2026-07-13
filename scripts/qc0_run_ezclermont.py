import os
import subprocess
import csv
import shutil
import re
from concurrent.futures import ProcessPoolExecutor, as_completed

# ==================== 【绝对路径配置区】 ====================
# 只专注这一个基因组文件夹
GENOME_DIR = "/home/u2024110181382/Project/02.E.coli/drep_result/high_quality_genomes" 
OUTPUT_FILE = "/home/u2024110181382/Project/04.pipeline/B2/Phylogroup_Results.csv"
MAX_WORKERS = 10
# ============================================================

def process_single_genome(filename, ez_executable):
    fasta_path = os.path.join(GENOME_DIR, filename)
    
    try:
        # 运行 ezclermont
        res = subprocess.run(
            [ez_executable, fasta_path], 
            capture_output=True, 
            text=True, 
            check=False, 
            env=os.environ
        )
        
        # 融合 stdout 和 stderr 捕获所有输出
        full_output = (res.stdout or "") + "\n" + (res.stderr or "")
        
        # 利用你原本成功的正则逻辑强行搜寻 Clermont type
        match = re.search(r"Clermont type(?:\s+is)?:\s*([A-Za-z0-9/_\-]+)", full_output, re.IGNORECASE)
        if match:
            return {"Genome": filename, "Phylogroup": match.group(1).strip().upper()}
            
        if res.returncode == 0:
            stdout_clean = res.stdout.strip()
            if stdout_clean:
                phylogroup = stdout_clean.split()[-1]
                return {"Genome": filename, "Phylogroup": phylogroup.upper()}
                
        err_info = res.stderr.strip().replace('\n', ' ')[:50]
        return {"Genome": filename, "Phylogroup": f"Failed_Exit_{res.returncode}: {err_info}"}
        
    except Exception as e:
        return {"Genome": filename, "Phylogroup": f"SysErr: {type(e).__name__}"}

def main():
    print(f"【纯净版】EzClermont 52核心全量基因组跑批启动...")
    
    ez_executable = shutil.which("ezclermont")
    if not ez_executable:
        ez_executable = "ezclermont"
        
    if not os.path.exists(GENOME_DIR):
        print(f"❌ 错误：基因组目录不存在：{GENOME_DIR}")
        return
        
    # 直接获取该目录下所有的 .fna 文件
    all_fna_files = [f for f in os.listdir(GENOME_DIR) if f.endswith('.fna')]
    total_genomes = len(all_fna_files)
    print(f"📂 锁定待测 .fna 基因组共 {total_genomes} 个。开始多进程并行计算...")
    
    results = []
    completed = 0
    
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(process_single_genome, fname, ez_executable): fname for fname in all_fna_files}
        for future in as_completed(futures):
            results.append(future.result())
            completed += 1
            if completed % 100 == 0 or completed == total_genomes:
                print(f" ⏳ 分型进度: {completed} / {total_genomes}...")
                
    with open(OUTPUT_FILE, mode='w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['Genome', 'Phylogroup'])
        for row in results:
            writer.writerow([row['Genome'], row['Phylogroup']])
            
    print(f"\n🎉 运行完成！结果已写入：{OUTPUT_FILE}")

if __name__ == "__main__":
    main()