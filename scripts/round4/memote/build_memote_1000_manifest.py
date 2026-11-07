import csv
import hashlib
import tarfile
from pathlib import Path

root = Path(r"D:\pipeline\analysis_v5")
archive_dir = root / "18.mSystems_revision3_final_20260911" / "05_model_archive"
model_manifest = archive_dir / "Supplementary_Models_1000_gapfill_annotated_manifest.csv"
crosswalk = root / "09.mSystems_revision_data_freeze" / "03_phylogroup_crosswalk" / "ezclermont_to_final_phylogroup_crosswalk.csv"
source_dir = root / "exported_xml_models_remote_verified"
output = root / "19.mSystems_revision4_20261028" / "06_reproducibility" / "MEMOTE_1000_annotated_model_manifest.tsv"

with model_manifest.open(encoding="utf-8-sig", newline="") as fh:
    models = list(csv.DictReader(fh))
with crosswalk.open(encoding="utf-8-sig", newline="") as fh:
    labels = {row["Genome_ID"]: row for row in csv.DictReader(fh) if row["Included_In_Final_Cohort"] == "True"}
assert len(models) == len(labels) == 1000
assert {r["Genome_ID"] for r in models} == set(labels)

tar_path = archive_dir / "Supplementary_Models_1000_gapfill_annotated.tar.gz"
with tarfile.open(tar_path, "r:gz") as tar:
    members = {m.name: m for m in tar if m.isfile()}
assert {r["Archive_Member"] for r in models} == set(members)

rows = []
for row in models:
    name = row["Genome_ID"]
    source = source_dir / (name + ".xml")
    if not source.is_file():
        raise RuntimeError(f"Missing source model: {source}")
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    if digest != row["Source_SHA256"]:
        raise RuntimeError(f"Source hash mismatch: {name}")
    label = labels[name]
    rows.append({
        "Genome_UID": label["Genome_UID"], "Genome_ID": name,
        "Accession": "_".join(name.split("_")[:2]),
        "Final_Phylogroup": label["Final_Phylogroup"],
        "Archive_Member": row["Archive_Member"],
        "Annotated_SHA256": row["Annotated_SHA256"],
        "Source_SHA256": row["Source_SHA256"],
        "Source_Reaction_Count": row["Reaction_Count"],
        "Gapfilled_Reactions_Marked": row["Gapfilled_Reactions_Marked"],
        "Local_Source_XML": str(source),
    })
with output.open("w", newline="", encoding="utf-8") as fh:
    writer = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
print(f"Validated {len(rows)} cohort labels, archive members and source model hashes; wrote {output}")
