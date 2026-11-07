"""Independently audit the two downloaded raw-input archives and their indices."""

import csv
import hashlib
import json
import tarfile
from pathlib import Path

root = Path(__file__).resolve().parent
raw = root / "data_download"
upload = root / "01_zenodo_upload"
qa = root / "03_qa"
server = json.loads((qa / "raw_input_package_results.json").read_text(encoding="utf-8"))
expected = {Path(row["archive"]).name: row for row in server}

with (raw / "raw_cog_audit" / "frozen_1000_file_checksum_audit.csv").open(encoding="utf-8-sig", newline="") as f:
    frozen = list(csv.DictReader(f))
annotations = {row["Genome_ID"]: row for row in frozen if row["File_Type"] == "EggNOG_Annotations"}
fastas = {row["Genome_ID"]: row for row in frozen if row["File_Type"] == "FASTA"}
with (raw / "revision4_finalcohort_broad_focused_20261030" / "inputs" / "sample_180.csv").open(encoding="utf-8-sig", newline="") as f:
    sample = {row["Genome_ID"]: row for row in csv.DictReader(f)}


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


result = []
for name, kind, count, source in (
    ("round4_eggnog_annotations_1000.tar.gz", "raw_annotations", 1000, annotations),
    ("round4_sample180_fastas.tar.gz", "sample180", 180, fastas),
):
    path = upload / name
    assert path.is_file() and path.stat().st_size == expected[name]["bytes"]
    digest = sha256(path)
    assert digest == expected[name]["sha256"]
    with tarfile.open(path, "r:gz") as tar:
        members = tar.getmembers()
        file_members = [m for m in members if m.name.startswith(kind + "/")]
        assert len(file_members) == count and len({m.name for m in members}) == len(members)
        index = json.load(tar.extractfile("SHA256_INPUT_INDEX.json"))
        assert len(index) == count
        assert {r["archive_member"] for r in index} == {m.name for m in file_members}
        for row in index:
            genome = Path(row["archive_member"]).name.split(".emapper.annotations")[0] if kind == "raw_annotations" else Path(row["archive_member"]).name.removesuffix(".fna")
            assert source[genome]["SHA256"].lower() == row["sha256"].lower()
            assert int(source[genome]["Bytes"]) == row["bytes"]
            if kind == "sample180":
                assert genome in sample
    result.append({"archive": name, "bytes": path.stat().st_size, "sha256": digest,
                   "indexed_input_files": count, "metadata_matches_frozen_manifest": True})

(qa / "DOWNLOADED_RAW_ARCHIVE_VALIDATION.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
with (upload / "RAW_INPUT_SHA256SUMS.txt").open("w", encoding="utf-8") as f:
    for row in result:
        f.write(f"{row['sha256']}  {row['archive']}\n")
print(json.dumps(result, indent=2))
