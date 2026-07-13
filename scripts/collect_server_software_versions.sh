#!/usr/bin/env bash
set -euo pipefail

OUTDIR="${1:-server_software_versions_$(date +%Y%m%d_%H%M%S)}"
GAPSEQ_DIR="${GAPSEQ_DIR:-/home/u2024110181382/Project/02.E.coli/biosoft/gapseq}"
GAPSEQ_BIN="${GAPSEQ_BIN:-${GAPSEQ_DIR}/gapseq}"
FASTDB="${FASTDB:-${GAPSEQ_DIR}/enterobacterales_gapseq_fastdb.faa}"
M9_MEDIUM="${M9_MEDIUM:-${GAPSEQ_DIR}/server_deploy/M9_glucose_medium.csv}"

mkdir -p "$OUTDIR"
TXT="${OUTDIR}/server_software_versions.txt"
CSV="${OUTDIR}/server_software_versions.csv"
CHECKS="${OUTDIR}/workflow_file_checksums.tsv"

run_capture() {
  local name="$1"
  shift
  {
    echo "### ${name}"
    if command -v "$1" >/dev/null 2>&1 || [[ -x "$1" ]]; then
      "$@" 2>&1 || true
    else
      echo "NOT_FOUND: $1"
    fi
    echo
  } >> "$TXT"
}

csv_row() {
  local tool="$1"
  local command="$2"
  local version="$3"
  printf '"%s","%s","%s"\n' "$tool" "$command" "$version" >> "$CSV"
}

capture_version_line() {
  local tool="$1"
  local command_label="$2"
  shift 2
  local output
  if command -v "$1" >/dev/null 2>&1 || [[ -x "$1" ]]; then
    output=$("$@" 2>&1 | head -n 5 | tr '\n' '; ' | sed 's/"/'\''/g') || true
  else
    output="NOT_FOUND: $1"
  fi
  csv_row "$tool" "$command_label" "$output"
}

{
  echo "Server software/version audit"
  echo "Generated: $(date -Iseconds)"
  echo "Hostname: $(hostname)"
  echo "User: $(whoami)"
  echo "Working directory: $(pwd)"
  echo "GAPSEQ_DIR: ${GAPSEQ_DIR}"
  echo "GAPSEQ_BIN: ${GAPSEQ_BIN}"
  echo "FASTDB: ${FASTDB}"
  echo "M9_MEDIUM: ${M9_MEDIUM}"
  echo
} > "$TXT"

printf '"Tool","Command","Version_Output"\n' > "$CSV"

capture_version_line "bash" "bash --version" bash --version
capture_version_line "python" "python --version" python --version
capture_version_line "python3" "python3 --version" python3 --version
capture_version_line "Rscript" "Rscript --version" Rscript --version
capture_version_line "conda" "conda --version" conda --version
capture_version_line "mamba" "mamba --version" mamba --version
capture_version_line "gapseq_optimized" "${GAPSEQ_BIN} --version" "${GAPSEQ_BIN}" --version
capture_version_line "gapseq_path" "gapseq --version" gapseq --version
capture_version_line "diamond" "diamond version" diamond version
capture_version_line "blastp" "blastp -version" blastp -version
capture_version_line "fastANI" "fastANI --version" fastANI --version
capture_version_line "prokka" "prokka --version" prokka --version
capture_version_line "checkm2" "checkm2 --version" checkm2 --version
capture_version_line "eggnog-mapper" "emapper.py --version" emapper.py --version
capture_version_line "signalp6" "signalp6 --version" signalp6 --version
capture_version_line "signalp" "signalp --version" signalp --version
capture_version_line "tmhmm" "tmhmm -h" tmhmm -h
capture_version_line "CRISPRCasTyper" "CRISPRCasTyper --version" CRISPRCasTyper --version
capture_version_line "minced" "minced -version" minced -version
capture_version_line "seqkit" "seqkit version" seqkit version

{
  echo "### PATH"
  echo "$PATH"
  echo
  echo "### Conda environments"
  command -v conda >/dev/null 2>&1 && conda info --envs 2>&1 || echo "conda not found"
  echo
  echo "### Python package versions"
  python - <<'PY' 2>&1 || true
import importlib.metadata as md
packages = [
    "cobra", "memote", "pandas", "numpy", "scipy", "sklearn",
    "matplotlib", "seaborn", "Bio", "openpyxl", "docx"
]
for pkg in packages:
    query = {"sklearn": "scikit-learn", "Bio": "biopython", "docx": "python-docx"}.get(pkg, pkg)
    try:
        print(f"{query}=={md.version(query)}")
    except Exception as exc:
        print(f"{query}: NOT_FOUND ({exc})")
PY
  echo
  echo "### R session info"
  Rscript -e 'sessionInfo()' 2>&1 || echo "Rscript not available"
  echo
  echo "### gapseq executable and source provenance"
  ls -l "${GAPSEQ_BIN}" 2>&1 || true
  if [[ -d "${GAPSEQ_DIR}/.git" ]]; then
    git -C "${GAPSEQ_DIR}" rev-parse HEAD 2>&1 || true
    git -C "${GAPSEQ_DIR}" status --short 2>&1 || true
    git -C "${GAPSEQ_DIR}" log -1 --oneline 2>&1 || true
  else
    echo "No .git directory detected at ${GAPSEQ_DIR}"
  fi
} >> "$TXT"

{
  printf "file\tbytes\tsha256\n"
  for f in \
    "${GAPSEQ_BIN}" \
    "${GAPSEQ_DIR}/src/doall.sh" \
    "${GAPSEQ_DIR}/src/gapseq_find.sh" \
    "${GAPSEQ_DIR}/src/transporter.sh" \
    "${GAPSEQ_DIR}/src/gf.suite.R" \
    "${FASTDB}" \
    "${M9_MEDIUM}"; do
    if [[ -f "$f" ]]; then
      bytes=$(wc -c < "$f" | tr -d ' ')
      hash=$(sha256sum "$f" | awk '{print $1}')
      printf "%s\t%s\t%s\n" "$f" "$bytes" "$hash"
    else
      printf "%s\tNOT_FOUND\tNOT_FOUND\n" "$f"
    fi
  done
} > "$CHECKS"

echo "Wrote:"
echo "  ${TXT}"
echo "  ${CSV}"
echo "  ${CHECKS}"
