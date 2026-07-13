from __future__ import annotations

import csv
import hashlib
import shutil
from pathlib import Path

import pandas as pd


ROOT = Path(r"D:\pipeline\analysis_v5")
MODEL_DIR = ROOT / "exported_xml_models"
MODEL_ARCHIVE = ROOT / "exported_models.tar.gz"
SUPP = ROOT / "05.supplementary_data"
DATA_S10 = SUPP / "Data_S10_optimized_gem_models"
OPEN = ROOT / "06.open_science_package"


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


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


def main() -> None:
    DATA_S10.mkdir(parents=True, exist_ok=True)
    OPEN.mkdir(exist_ok=True)

    cohort = pd.read_csv(ROOT / "00.data" / "genome_assembly_and_cog_metrics.csv")
    qc = pd.read_csv(ROOT / "memote_reports" / "cobrapy_qc_results.csv")
    def normalize_model_id(value: object) -> str:
        text = str(value)
        if text.endswith(".xml"):
            text = text[:-4]
        if text.endswith(".fna"):
            text = text[:-4]
        return text

    cohort["_norm_id"] = cohort["Genome_ID"].map(normalize_model_id)
    qc["_norm_id"] = qc["Genome_ID"].map(normalize_model_id)
    cohort_map = cohort.set_index("_norm_id")["Phylogroup"].to_dict()
    cohort_original_id = cohort.set_index("_norm_id")["Genome_ID"].to_dict()
    qc_map = qc.set_index("_norm_id").to_dict("index")

    model_files = sorted(MODEL_DIR.glob("*.xml"))
    rows = []
    unmatched = []
    for f in model_files:
        genome_id = f.stem
        if genome_id not in cohort_map:
            unmatched.append(f.name)
        qc_row = qc_map.get(genome_id, {})
        rows.append(
            {
                "Genome_ID": genome_id,
                "Cohort_Genome_ID": cohort_original_id.get(genome_id, ""),
                "Phylogroup": cohort_map.get(genome_id, ""),
                "Model_File": f.name,
                "Format": "SBML/XML",
                "Relative_Path_In_Archive": f"exported_xml_models/{f.name}",
                "Bytes": f.stat().st_size,
                "SHA256": sha256_file(f),
                "COBRApy_QC_Status": "available in cobrapy_qc_results.csv" if genome_id in qc_map else "not found in QC table",
                "Total_Reactions": qc_row.get("Total_Reactions", ""),
                "Total_Metabolites": qc_row.get("Total_Metabolites", ""),
                "Mass_Balanced_Ratio(%)": qc_row.get("Mass_Balanced_Ratio(%)", ""),
                "Blocked_Reactions_Ratio(%)": qc_row.get("Blocked_Reactions_Ratio(%)", ""),
                "Deadend_Metabolites_Ratio(%)": qc_row.get("Deadend_Metabolites_Ratio(%)", ""),
            }
        )

    manifest = DATA_S10 / "ecoli_1000_optimized_gem_models_manifest.csv"
    write_csv(manifest, rows)

    archive_sha = sha256_file(MODEL_ARCHIVE)
    archive_rows = [
        {
            "Archive_File": MODEL_ARCHIVE.name,
            "Bytes": MODEL_ARCHIVE.stat().st_size,
            "SHA256": archive_sha,
            "Model_Count_In_Manifest": len(rows),
            "Model_Format": "SBML/XML",
            "Directory_When_Unpacked": "exported_xml_models/",
        }
    ]
    write_csv(DATA_S10 / "ecoli_1000_optimized_gem_models_archive_checksum.csv", archive_rows)

    readme = f"""Data S10. Optimized GEM model archive

This supplementary dataset contains the optimized gapseq genome-scale metabolic models generated for the 1,000-genome Escherichia coli analytical cohort.

Files:
- exported_models.tar.gz: compressed archive containing 1,000 SBML/XML GEM files.
- ecoli_1000_optimized_gem_models_manifest.csv: per-model manifest with Genome_ID, phylogroup, filename, byte size, SHA256 checksum, and corresponding COBRApy QC metrics.
- ecoli_1000_optimized_gem_models_archive_checksum.csv: archive-level byte size and SHA256 checksum.

Archive SHA256:
{archive_sha}

Model count:
{len(rows)} XML files.

Notes:
- Model filenames are keyed by Genome_ID and match the 1,000-genome cohort table in Data S1.
- COBRApy-derived quality-control metrics are also summarized in Data S9.
- For journal submission, this archive should be deposited in Zenodo, Figshare, BioModels, or the journal's data repository. The resulting DOI/URL should be inserted into the Data availability statement.
"""
    (DATA_S10 / "README.txt").write_text(readme, encoding="utf-8")

    dest_archive = DATA_S10 / MODEL_ARCHIVE.name
    if not dest_archive.exists() or dest_archive.stat().st_size != MODEL_ARCHIVE.stat().st_size:
        shutil.copy2(MODEL_ARCHIVE, dest_archive)

    shutil.copy2(manifest, OPEN / manifest.name)
    shutil.copy2(DATA_S10 / "ecoli_1000_optimized_gem_models_archive_checksum.csv", OPEN / "ecoli_1000_optimized_gem_models_archive_checksum.csv")

    missing_rows = [
        {
            "check": "model_file_count",
            "status": "pass" if len(rows) == 1000 else "fail",
            "value": len(rows),
            "expected": 1000,
        },
        {
            "check": "all_models_match_cohort_table",
            "status": "pass" if not unmatched else "fail",
            "value": len(unmatched),
            "expected": 0,
        },
        {
            "check": "archive_present",
            "status": "pass" if MODEL_ARCHIVE.exists() else "fail",
            "value": str(MODEL_ARCHIVE),
            "expected": "exported_models.tar.gz",
        },
        {
            "check": "archive_size_bytes",
            "status": "pass" if MODEL_ARCHIVE.stat().st_size > 0 else "fail",
            "value": MODEL_ARCHIVE.stat().st_size,
            "expected": ">0",
        },
    ]
    write_csv(DATA_S10 / "data_s10_integrity_checks.csv", missing_rows)
    shutil.copy2(DATA_S10 / "data_s10_integrity_checks.csv", OPEN / "data_s10_integrity_checks.csv")

    print(f"Models: {len(rows)}")
    print(f"Unmatched to cohort: {len(unmatched)}")
    print(f"Archive bytes: {MODEL_ARCHIVE.stat().st_size}")
    print(f"Archive SHA256: {archive_sha}")
    print(f"Wrote {DATA_S10}")


if __name__ == "__main__":
    main()
