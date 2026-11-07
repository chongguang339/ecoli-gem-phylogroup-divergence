"""Independently rescan raw COG_category fields inside the downloaded TAR.GZ."""

import csv
import hashlib
import json
import tarfile
from pathlib import Path

root = Path(__file__).resolve().parent
archive = root / "data_download" / "raw_cog_audit" / "round4_eggnog_annotations_1000.tar.gz"
freeze_file = root / "data_download" / "raw_cog_audit" / "frozen_1000_file_checksum_audit.csv"
audit_file = root / "data_download" / "raw_cog_audit" / "eggnog_1000_raw_cog_audit.tsv"
freeze = {r["Genome_ID"]: r for r in csv.DictReader(freeze_file.open(encoding="utf-8-sig"))
          if r["File_Type"] == "EggNOG_Annotations"}
audit = {r["Genome_ID"]: r for r in csv.DictReader(audit_file.open(encoding="utf-8"), delimiter="\t")}
assert len(freeze) == len(audit) == 1000

counts = {"files": 0, "protein_rows": 0, "R": 0, "X": 0, "Y": 0, "sha256_matched": 0,
          "audit_rows_matched": 0}
with tarfile.open(archive, "r|gz") as tar:
    for member in tar:
        if not member.name.startswith("raw_annotations/"):
            continue
        genome = Path(member.name).name.removesuffix(".emapper.annotations")
        assert genome in freeze and genome in audit, genome
        source = tar.extractfile(member)
        assert source is not None
        digest = hashlib.sha256()
        column = None
        rows = 0
        r_count = x_count = y_count = 0
        for raw in source:
            digest.update(raw)
            line = raw.decode("utf-8", errors="replace")
            if column is None:
                if line.startswith("#query"):
                    headers = line.lstrip("#").rstrip("\r\n").split("\t")
                    column = headers.index("COG_category")
                continue
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            rows += 1
            category = fields[column] if len(fields) > column else "-"
            r_count += "R" in category
            x_count += "X" in category
            y_count += "Y" in category
        assert column is not None, genome
        assert member.size == int(freeze[genome]["Bytes"]), genome
        assert digest.hexdigest() == freeze[genome]["SHA256"].lower(), genome
        assert (rows, r_count, x_count, y_count) == tuple(int(audit[genome][key]) for key in
                                                         ("Annotated_Rows", "COG_R", "COG_X", "COG_Y")), genome
        counts["files"] += 1
        counts["protein_rows"] += rows
        counts["R"] += r_count
        counts["X"] += x_count
        counts["Y"] += y_count
        counts["sha256_matched"] += 1
        counts["audit_rows_matched"] += 1
assert counts["files"] == 1000
(root / "03_qa" / "INDEPENDENT_RAW_COG_RESCAN.json").write_text(json.dumps(counts, indent=2), encoding="utf-8")
print(json.dumps(counts, indent=2))
