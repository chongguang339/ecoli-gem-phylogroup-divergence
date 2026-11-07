#!/usr/bin/env bash
set -u
id="$1"
root=/home/u2024110181382/Project/04.pipeline/revision4_finalcohort_broad_focused_20261030
pipeline=/home/u2024110181382/Project/04.pipeline
exe=/home/u2024110181382/Project/02.E.coli/biosoft/gapseq/gapseq
medium=/home/u2024110181382/Project/02.E.coli/biosoft/gapseq/server_deploy/M9_glucose_medium.csv
fastdb=/home/u2024110181382/Project/02.E.coli/biosoft/gapseq/enterobacterales_gapseq_fastdb.faa
fasta=$(awk -F, -v id="$id" '$1==id {print $3; exit}' "$root/inputs/sample_180.csv")
if [[ ! -s "$fasta" ]]; then echo "missing FASTA $id" >&2; exit 2; fi
for config in broad focused; do
 out="$root/models/$id/$config"; mkdir -p "$out"
 final="$out/$id.xml"; log="$root/logs/${id}_${config}.log"; status="$root/logs/${id}_${config}.status.tsv"
 if [[ -s "$final" ]]; then printf '%s\t%s\tSKIPPED_EXISTING\t0\t0\t%s\n' "$id" "$config" "$final" > "$status"; continue; fi
 args=(-A diamond -S identical -m "$medium" -K 4 -f "$out")
 if [[ "$config" == focused ]]; then args+=(-E "$fastdb"); fi
 start=$(date +%s)
 (cd "$out" && /usr/bin/time -v "$exe" doall "${args[@]}" "$fasta") > "$log" 2>&1
 rc=$?; elapsed=$(($(date +%s)-start))
 if [[ $rc -eq 0 && -s "$final" ]]; then state=SUCCESS; else state=FAILED; fi
 printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$id" "$config" "$state" "$rc" "$elapsed" "$final" > "$status"
done
