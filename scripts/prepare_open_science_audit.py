from __future__ import annotations

import csv
import shutil
import zipfile
from pathlib import Path


ROOT = Path(r"D:\pipeline\analysis_v5")
OUT = ROOT / "06.open_science_package"
CODE_OUT = OUT / "code_repository_candidate"


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="", encoding="utf-8-sig") as fh:
        writer = csv.DictWriter(fh, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def collect_files(base: Path, patterns: tuple[str, ...]) -> list[Path]:
    files: list[Path] = []
    if not base.exists():
        return files
    for pattern in patterns:
        files.extend(base.rglob(pattern))
    return sorted({p for p in files if p.is_file()})


def main() -> None:
    OUT.mkdir(exist_ok=True)

    code_files = sorted((ROOT / "scripts").glob("*.py"))
    ez_script = ROOT / "脚本" / "qc0_run_ezclermont.py"
    if ez_script.exists():
        code_files.append(ez_script)

    code_rows = []
    for f in code_files:
        code_rows.append(
            {
                "file": rel(f),
                "bytes": f.stat().st_size,
                "public_status": "include",
                "note": "Analysis, figure-generation, manuscript-update, or wrapper script used during this study.",
            }
        )
    write_csv(OUT / "code_manifest.csv", code_rows)

    data_rows = []
    supp_root = ROOT / "05.supplementary_data"
    for f in sorted(supp_root.rglob("*")):
        if f.is_file():
            data_rows.append(
                {
                    "file": rel(f),
                    "bytes": f.stat().st_size,
                    "upload_target": "Zenodo/Figshare supplementary data",
                    "status": "available",
                }
            )
    write_csv(OUT / "supplementary_data_upload_manifest.csv", data_rows)

    model_patterns = ("*.xml", "*.sbml", "*.json")
    candidate_model_files = [
        p
        for p in collect_files(ROOT, model_patterns)
        if ".python_packages" not in str(p)
        and "02.论文" not in str(p)
        and ("results_optimized" in str(p) or "mg1655_ascii_models" in str(p) or "tmp_repaired_sbml" in str(p))
    ]
    model_rows = []
    optimized_draft = [p for p in candidate_model_files if "results_optimized" in str(p) and p.name.endswith("-draft.xml")]
    for f in candidate_model_files:
        model_rows.append(
            {
                "file": rel(f),
                "bytes": f.stat().st_size,
                "detected_category": (
                    "20-genome/reference optimized model"
                    if "results_optimized" in str(f)
                    else "MG1655 or repaired benchmark model"
                ),
                "status": "available locally but not the full 1,000-GEM cohort",
            }
        )
    write_csv(OUT / "detected_model_files_not_full_1000.csv", model_rows)

    missing_rows = [
        {
            "item": "Repository URL",
            "required_for_submission": "yes",
            "current_status": "missing placeholder in manuscript",
            "action": "Create GitHub repository or Zenodo software record, upload code, then replace manuscript placeholder with DOI/URL.",
        },
        {
            "item": "1,000 optimized GEM model files",
            "required_for_submission": "yes",
            "current_status": f"not found as a complete model package in analysis_v5; detected {len(optimized_draft)} optimized draft XML files, mainly benchmark/reference models",
            "action": "Collect all 1,000 optimized SBML/XML or JSON files from the production gapseq output directory and package as Data S10 or external repository archive.",
        },
        {
            "item": "Software environment file",
            "required_for_submission": "strongly recommended",
            "current_status": "not found as conda/environment.yml or requirements.txt in the project root",
            "action": "Add requirements.txt or environment.yml listing Python, cobra, pandas, scipy, scikit-learn, matplotlib/seaborn, python-docx, and external tools with versions.",
        },
        {
            "item": "License",
            "required_for_submission": "strongly recommended",
            "current_status": "not found in open-science package",
            "action": "Add LICENSE, preferably MIT/BSD-3-Clause for code; use CC BY 4.0 for data where allowed.",
        },
        {
            "item": "Model-to-genome manifest",
            "required_for_submission": "yes",
            "current_status": "not found for 1,000 model files",
            "action": "Create CSV with Genome_ID, phylogroup, source assembly accession, optimized model filename, file size, MD5/SHA256, and COBRApy load status.",
        },
    ]
    write_csv(OUT / "missing_items_for_submission.csv", missing_rows)

    CODE_OUT.mkdir(parents=True, exist_ok=True)
    (CODE_OUT / "scripts").mkdir(exist_ok=True)
    for f in sorted((ROOT / "scripts").glob("*.py")):
        shutil.copy2(f, CODE_OUT / "scripts" / f.name)
    if ez_script.exists():
        (CODE_OUT / "scripts" / "qc0_run_ezclermont.py").write_bytes(ez_script.read_bytes())

    readme = """# E. coli phylogroup-resolved GEM workflow

This repository candidate contains the analysis and figure-generation scripts used for the manuscript:

High-Throughput GEM Workflow Optimization for Phylogroup-Resolved Functional Divergence Analysis in Escherichia coli

## Contents

- `scripts/`: Python wrappers, statistical analyses, FBA benchmarking, figure generation, and manuscript update utilities.
- `../supplementary_data_upload_manifest.csv`: files already available for upload as Supplementary Data.
- `../missing_items_for_submission.csv`: items that must be completed before journal submission.

## Reproducibility note

Production scripts were run on local high-performance-computing storage during analysis. For public reuse, users should replace local input/output paths with project-relative paths or command-line arguments, while preserving the documented input table names and output file structure.

## Data and models

Processed analysis tables are available in the Supplementary Data package. The complete set of 1,000 optimized GEM SBML/JSON files should be deposited separately as an archival data record because model files are large and are not currently included in this code package.

## Suggested citation

Please cite the manuscript and the archived Zenodo/Figshare record DOI after deposition.
"""
    (CODE_OUT / "README.md").write_text(readme, encoding="utf-8")

    upload_plan = """# Open-science upload plan

## Recommended deposition structure

1. GitHub repository: code only
   - Include `scripts/`, `README.md`, `requirements.txt` or `environment.yml`, and `LICENSE`.
   - Do not upload large model archives to GitHub.

2. Zenodo or Figshare: data and model archives
   - Upload `05.supplementary_data/supplementary_data_full.zip`.
   - Upload a separate archive containing the 1,000 optimized GEM SBML/XML or JSON model files.
   - Upload a model manifest with Genome_ID, phylogroup, filename, checksum, model format, and COBRApy load status.

3. Manuscript replacement text
   - Replace `Repository URL: pending insertion before submission` with the GitHub URL and/or Zenodo DOI.
   - Add a sentence: `The 1,000 optimized GEMs are available as Supplementary Data S10 at [DOI/URL].`

## Suggested archive names

- `ecoli_gem_phylogroup_code_v1.0.0.zip`
- `ecoli_gem_phylogroup_supplementary_data_v1.0.0.zip`
- `ecoli_1000_optimized_gem_models_sbml_v1.0.0.zip`
- `ecoli_1000_optimized_gem_models_manifest.csv`
"""
    (OUT / "OPEN_SCIENCE_UPLOAD_PLAN.md").write_text(upload_plan, encoding="utf-8")

    zip_path = OUT / "ecoli_gem_phylogroup_code_candidate.zip"
    if zip_path.exists():
        zip_path.unlink()
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for f in CODE_OUT.rglob("*"):
            if f.is_file():
                zf.write(f, f.relative_to(CODE_OUT.parent))

    print(f"Wrote {OUT}")
    print(f"Code files: {len(code_rows)}")
    print(f"Supplementary files: {len(data_rows)}")
    print(f"Detected candidate model files: {len(model_rows)}")
    print(f"Detected optimized draft XML files: {len(optimized_draft)}")
    print(f"Code zip: {zip_path}")


if __name__ == "__main__":
    main()
