"""Hash-verify and package the 1,000 raw eggNOG annotations and 180 FASTAs.

This script reads existing files only. It performs no biological re-analysis.
"""

import csv
import hashlib
import io
import json
import tarfile
from pathlib import Path

ROOT = Path("/home/u2024110181382/Project/04.pipeline")
OUT = ROOT / "round4_public_missing_inputs_20261107"
OUT.mkdir(exist_ok=True)
FREEZE = ROOT / "revision_20260905_data_freeze" / "frozen_1000_file_checksum_audit.csv"
SAMPLE = ROOT / "revision4_finalcohort_broad_focused_20261030" / "inputs" / "sample_180.csv"
COG_DIR = ROOT / "revision4_source_inventory_20261029"


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


with FREEZE.open(encoding="utf-8-sig", newline="") as stream:
    freeze = list(csv.DictReader(stream))
ann = [row for row in freeze if row["File_Type"] == "EggNOG_Annotations"]
fasta = {row["Genome_ID"]: row for row in freeze if row["File_Type"] == "FASTA"}
with SAMPLE.open(encoding="utf-8-sig", newline="") as stream:
    sample = list(csv.DictReader(stream))
assert len(ann) == 1000 and len({r["Genome_ID"] for r in ann}) == 1000
assert len(sample) == 180 and len({r["Genome_ID"] for r in sample}) == 180


def verify(source, expected_sha, expected_size):
    path = Path(source)
    assert path.is_file(), f"MISSING: {path}"
    assert path.stat().st_size == int(expected_size), f"SIZE MISMATCH: {path}"
    actual = sha256(path)
    assert actual.lower() == expected_sha.lower(), f"SHA256 MISMATCH: {path}: {actual} != {expected_sha}"
    return path, actual


annotation_files = []
for row in ann:
    path, actual = verify(row["Path"], row["SHA256"], row["Bytes"])
    annotation_files.append((path, f"raw_annotations/{row['Genome_ID']}.emapper.annotations", actual))
print("VERIFIED 1000 eggNOG annotation SHA256", flush=True)

fasta_files = []
for row in sample:
    frozen = fasta[row["Genome_ID"]]
    assert frozen["Path"] == row["FASTA"], row["Genome_ID"]
    assert frozen["SHA256"].lower() == row["Frozen_FASTA_SHA256"].lower(), row["Genome_ID"]
    path, actual = verify(row["FASTA"], row["Frozen_FASTA_SHA256"], frozen["Bytes"])
    fasta_files.append((path, f"sample180/{row['Genome_ID']}.fna", actual))
print("VERIFIED 180 FASTA SHA256", flush=True)


def package(name, files, supplements):
    target = OUT / name
    if target.exists():
        raise FileExistsError(f"Refusing to overwrite existing archive: {target}")
    index = [{"archive_member": member, "source": str(path), "bytes": path.stat().st_size,
              "sha256": digest} for path, member, digest in files]
    with tarfile.open(target, "w:gz", compresslevel=3) as tar:
        for path, member, _ in files:
            tar.add(path, arcname=member, recursive=False)
        for path, member in supplements:
            tar.add(path, arcname=member, recursive=False)
        content = json.dumps(index, ensure_ascii=False, indent=2).encode("utf-8")
        info = tarfile.TarInfo("SHA256_INPUT_INDEX.json")
        info.size = len(content)
        tar.addfile(info, io.BytesIO(content))
    with tarfile.open(target, "r:gz") as tar:
        assert len(tar.getmembers()) == len(files) + len(supplements) + 1
        for member in tar:
            if member.isfile():
                tar.extractfile(member).read(1)
    record = {"archive": str(target), "file_count": len(files), "bytes": target.stat().st_size,
              "sha256": sha256(target)}
    print(json.dumps(record), flush=True)
    return record


results = [
    package("round4_eggnog_annotations_1000.tar.gz", annotation_files, [
        (FREEZE, "provenance/frozen_1000_file_checksum_audit.csv"),
        (COG_DIR / "eggnog_1000_raw_cog_audit.tsv", "provenance/eggnog_1000_raw_cog_audit.tsv"),
        (COG_DIR / "eggnog_1000_raw_cog_summary.json", "provenance/eggnog_1000_raw_cog_summary.json"),
    ]),
    package("round4_sample180_fastas.tar.gz", fasta_files, [
        (SAMPLE, "provenance/sample_180.csv"),
    ]),
]
(OUT / "raw_input_package_results.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
print("COMPLETE", flush=True)
