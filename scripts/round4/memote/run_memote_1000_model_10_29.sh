#!/usr/bin/env bash
set -euo pipefail

# Run: bash run_memote_1000_model_10_29.sh
# Optional: WORKERS=8 bash run_memote_1000_model_10_29.sh

BASE=/home/u2024110181382/Project/04.pipeline/model_10.29
ARCHIVE="$BASE/Supplementary_Models_1000_gapfill_annotated.tar.gz"
OUT="$BASE/memote_1000_output"
MODELS="$OUT/models"
REPORTS="$OUT/reports"
LOGS="$OUT/logs"
STATUS="$OUT/status"
WORKERS="${WORKERS:-4}"
EXPECTED_SHA256=e8f33b10f9bb5dec4f476cd46180c0cdc930b34de11d5b56c188e43777a93f62

[[ -s "$ARCHIVE" ]] || { echo "Missing archive: $ARCHIVE" >&2; exit 1; }
[[ "$WORKERS" =~ ^[1-9][0-9]*$ ]] || { echo "WORKERS must be a positive integer" >&2; exit 1; }
actual_sha256=$(sha256sum "$ARCHIVE" | awk '{print $1}')
[[ "$actual_sha256" == "$EXPECTED_SHA256" ]] || {
    echo "Archive SHA256 differs from the frozen submission archive." >&2
    echo "Expected: $EXPECTED_SHA256" >&2
    echo "Actual:   $actual_sha256" >&2
    exit 1
}

mkdir -p "$MODELS" "$REPORTS" "$LOGS" "$STATUS"
if [[ ! -e "$OUT/extraction_complete" ]]; then
    echo "Extracting frozen annotated models to $OUT"
    tar -xzf "$ARCHIVE" -C "$OUT"
    touch "$OUT/extraction_complete"
fi

mapfile -d '' model_files < <(find "$MODELS" -maxdepth 1 -type f -name '*.xml' -print0 | sort -z)
[[ ${#model_files[@]} -eq 1000 ]] || {
    echo "Expected exactly 1000 XML files in $MODELS; found ${#model_files[@]}" >&2
    exit 1
}
printf '%s\n' "${model_files[@]}" > "$OUT/model_paths.txt"

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate memote-conda
command -v memote >/dev/null || { echo "memote not found in memote-conda" >&2; exit 1; }
memote --version > "$OUT/memote_version.txt" 2>&1
date -u +'%Y-%m-%dT%H:%M:%SZ' > "$OUT/run_started_utc.txt"
printf 'Archive SHA256: %s\nModels: %s\nWorkers: %s\nOutput: %s\n' \
    "$actual_sha256" "${#model_files[@]}" "$WORKERS" "$OUT"

run_one() {
    local model="$1" name report log status_file tmp start end code
    name=$(basename "$model" .xml)
    report="$REPORTS/$name.html"
    log="$LOGS/$name.log"
    status_file="$STATUS/$name.tsv"
    if [[ -s "$report" && -s "$status_file" ]] &&
       [[ $(cut -f2 "$status_file") == 0 ]] &&
       grep -q '"score":{"sections":' "$report"; then
        return 0
    fi
    start=$(date +%s)
    if timeout 3h memote report snapshot --filename "$report" "$model" > "$log" 2>&1; then
        code=0
    else
        code=$?
    fi
    end=$(date +%s)
    tmp="$status_file.tmp.$$"
    printf '%s\t%s\t%s\t%s\t%s\n' "$name" "$code" "$((end-start))" "$model" "$report" > "$tmp"
    mv "$tmp" "$status_file"
    return 0
}
export REPORTS LOGS STATUS
export -f run_one

printf '%s\0' "${model_files[@]}" | xargs -0 -P "$WORKERS" -n 1 bash -c 'run_one "$1"' _

python - "$STATUS" "$OUT/memote_all_models.csv" <<'PY'
import csv, json, pathlib, re, sys
status_dir, output = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
rows = []
for item in sorted(status_dir.glob('*.tsv')):
    name, code, seconds, model, report = item.read_text().rstrip('\n').split('\t')
    row = dict(model_id=name, exit_code=int(code), wall_seconds=int(seconds),
               model_xml=model, report_html=report)
    if row['exit_code'] == 0:
        html = pathlib.Path(report).read_bytes()
        match = re.search(rb'"score":(\{"sections":.*?\}),"cards"', html)
        if match is None:
            row['exit_code'] = 125
            row['error'] = 'score object missing in HTML'
        else:
            score = json.loads(match.group(1))
            row['total_score'] = score.get('total_score')
            for section in score.get('sections', []):
                row['section_' + str(section.get('section'))] = section.get('score')
    rows.append(row)
fields = sorted({key for row in rows for key in row})
with output.open('w', newline='', encoding='utf-8-sig') as handle:
    writer = csv.DictWriter(handle, fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)
passed = sum(row['exit_code'] == 0 for row in rows)
print(f'MEMOTE completed successfully: {passed}/{len(rows)}; summary: {output}')
if passed != 1000:
    print('Inspect failed status/*.tsv and logs/*.log, then rerun the same shell script.', file=sys.stderr)
    sys.exit(1)
PY

date -u +'%Y-%m-%dT%H:%M:%SZ' > "$OUT/run_finished_utc.txt"
