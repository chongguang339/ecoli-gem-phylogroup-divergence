#!/usr/bin/env bash
#SBATCH --job-name=ecoli_broad_focus_180
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=52
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --output=/home/u2024110181382/Project/04.pipeline/revision4_finalcohort_broad_focused_20261030/logs/slurm_%j.out
set -euo pipefail
root=/home/u2024110181382/Project/04.pipeline/revision4_finalcohort_broad_focused_20261030
source /home/u2024110181382/miniconda3/etc/profile.d/conda.sh
conda activate gapseq
python "$root/scripts/prepare_sample.py"
exe=/home/u2024110181382/Project/02.E.coli/biosoft/gapseq/gapseq
medium=/home/u2024110181382/Project/02.E.coli/biosoft/gapseq/server_deploy/M9_glucose_medium.csv
fastdb=/home/u2024110181382/Project/02.E.coli/biosoft/gapseq/enterobacterales_gapseq_fastdb.faa
for p in "$exe" "$medium" "$fastdb"; do test -s "$p" || { echo "missing $p"; exit 2; }; done
sha256sum "$exe" "$medium" "$fastdb" "$root/inputs/sample_180.csv" > "$root/inputs/run_inputs.sha256"
awk -F, 'NR>1 && !seen[$2]++ {print $1}' "$root/inputs/sample_180.csv" > "$root/inputs/pilot_ids.txt"
parallel --will-cite -j 12 "$root/scripts/run_pair.sh" {} :::: "$root/inputs/pilot_ids.txt"
pilot_ok=$(grep -l $'\tfocused\tSUCCESS\|\tfocused\tSKIPPED_EXISTING' "$root"/logs/*.status.tsv 2>/dev/null | wc -l)
echo "pilot_focused_success=$pilot_ok"
if [[ "$pilot_ok" -lt 8 ]]; then echo 'Pilot failed; stopping expansion'; exit 3; fi
grep -vxFf "$root/inputs/pilot_ids.txt" "$root/inputs/sample_180_ids.txt" > "$root/inputs/remainder_ids.txt"
parallel --will-cite -j 12 "$root/scripts/run_pair.sh" {} :::: "$root/inputs/remainder_ids.txt"
python - "$root" <<'PY'
import pathlib,sys,csv,collections
r=pathlib.Path(sys.argv[1]); rows=[]
for p in sorted((r/'logs').glob('*.status.tsv')):
 v=p.read_text().strip().split('\t');rows.append(v)
with (r/'summary/run_status.tsv').open('w') as h:
 h.write('Genome_ID\tConfig\tStatus\tReturn_Code\tElapsed_Seconds\tFinal_XML\n')
 for v in rows:h.write('\t'.join(v)+'\n')
print('run statuses',dict(collections.Counter((v[1],v[2]) for v in rows)),flush=True)
PY
