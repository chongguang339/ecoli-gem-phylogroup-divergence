from pathlib import Path
import zipfile

import pandas as pd


ROOT = Path(r"D:\pipeline\analysis_v5")
SUPP = ROOT / "05.supplementary_data"
MANIFEST = SUPP / "supplementary_data_manifest.csv"

NEW_ROWS = [
    {
        "Supplementary_Data": "Data S1",
        "Description": "Phylogroup proportional allocation plan from 10,951 EzClermont-typed source genomes.",
        "Source_Path": str(ROOT / "03.analysis" / "phylogroup_proportional_sampling_plan.csv"),
        "Packaged_Path": r"Data_S1_cohort_and_phylogroup\phylogroup_proportional_sampling_plan.csv",
        "Status": "copied",
    },
    {
        "Supplementary_Data": "Data S1",
        "Description": "Available cohort-level GEM QC and topology metrics; uncomputed advanced QC fields explicitly marked.",
        "Source_Path": str(ROOT / "03.analysis" / "gem_model_qc_available_metrics.csv"),
        "Packaged_Path": r"Data_S1_cohort_and_phylogroup\gem_model_qc_available_metrics.csv",
        "Status": "copied",
    },
]


def main():
    df = pd.read_csv(MANIFEST)
    existing = set(df["Packaged_Path"].astype(str))
    rows = []
    for row in NEW_ROWS:
        packaged = SUPP / row["Packaged_Path"]
        row = row.copy()
        row["Bytes"] = packaged.stat().st_size if packaged.exists() else 0
        if row["Packaged_Path"] not in existing:
            rows.append(row)
    if rows:
        df = pd.concat([df, pd.DataFrame(rows)], ignore_index=True)
        df.to_csv(MANIFEST, index=False)

    zip_path = SUPP / "supplementary_data_full.zip"
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in SUPP.rglob("*"):
            if path.is_file() and path.name not in {"supplementary_data_full.zip", "supplementary_data_package_index_only.zip"}:
                zf.write(path, path.relative_to(SUPP))
    print(MANIFEST)
    print(zip_path)


if __name__ == "__main__":
    main()
