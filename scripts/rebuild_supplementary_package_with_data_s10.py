from __future__ import annotations

import csv
import zipfile
from pathlib import Path


ROOT = Path(r"D:\pipeline\analysis_v5")
SUPP = ROOT / "05.supplementary_data"
MANIFEST = SUPP / "supplementary_data_manifest.csv"
ZIP = SUPP / "supplementary_data_full_with_models.zip"


def main() -> None:
    rows = []
    data_dirs = sorted([p for p in SUPP.iterdir() if p.is_dir() and p.name.startswith("Data_S")])
    for d in data_dirs:
        data_id = d.name.split("_", 2)[0] + " " + d.name.split("_", 2)[1] if "_" in d.name else d.name
        for f in sorted(d.rglob("*")):
            if f.is_file():
                rows.append(
                    {
                        "Supplementary_Data": data_id,
                        "Packaged_Path": str(f.relative_to(SUPP)),
                        "Bytes": f.stat().st_size,
                        "Status": "available",
                    }
                )

    with MANIFEST.open("w", newline="", encoding="utf-8-sig") as fh:
        writer = csv.DictWriter(fh, fieldnames=["Supplementary_Data", "Packaged_Path", "Bytes", "Status"])
        writer.writeheader()
        writer.writerows(rows)

    if ZIP.exists():
        ZIP.unlink()
    with zipfile.ZipFile(ZIP, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
        zf.write(MANIFEST, MANIFEST.relative_to(SUPP))
        for d in data_dirs:
            for f in sorted(d.rglob("*")):
                if f.is_file():
                    zf.write(f, f.relative_to(SUPP))

    print(f"Manifest rows: {len(rows)}")
    print(f"Wrote {MANIFEST}")
    print(f"Wrote {ZIP}")
    print(f"Zip bytes: {ZIP.stat().st_size}")


if __name__ == "__main__":
    main()
