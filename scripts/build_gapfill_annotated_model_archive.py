from __future__ import annotations

import csv
import hashlib
import io
import os
import tarfile
from collections import defaultdict
from pathlib import Path

from lxml import etree


ROOT = Path(os.environ.get("ECOLI_GEM_ROOT", Path(__file__).resolve().parents[1]))
SOURCE_ARCHIVE = ROOT / "exported_models.tar.gz"
EVIDENCE = ROOT / "10.mSystems_revision2_model_evidence" / "02_model_evidence" / "gapfilled_reactions_long.csv"
OUT_DIR = ROOT / "18.mSystems_revision3_final_20260911" / "05_model_archive"
ARCHIVE = OUT_DIR / "Supplementary_Models_1000_gapfill_annotated.tar.gz"
MANIFEST = OUT_DIR / "Supplementary_Models_1000_gapfill_annotated_manifest.csv"
ARCHIVE_SHA = OUT_DIR / "Supplementary_Models_1000_gapfill_annotated.tar.gz.sha256"
NS = "https://github.com/chongguang339/ecoli-gem-phylogroup-divergence/ns/gapfill"


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def load_gapfilled() -> dict[str, set[str]]:
    result: dict[str, set[str]] = defaultdict(set)
    with EVIDENCE.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            if row["Gap_Filled"].lower() == "true":
                result[row["Genome_ID"]].add(row["Reaction_ID"])
    return result


def annotate_model(source_name: str, source_bytes: bytes, expected: set[str]) -> tuple[bytes, int, int]:
    parser = etree.XMLParser(remove_blank_text=False, huge_tree=True)
    tree = etree.parse(io.BytesIO(source_bytes), parser)
    reactions = tree.xpath("//*[local-name()='listOfReactions']/*[local-name()='reaction']")
    found: set[str] = set()
    for reaction in reactions:
        rid = reaction.get("id")
        if rid not in expected:
            continue
        found.add(rid)
        annotation = next((child for child in reaction if etree.QName(child).localname == "annotation"), None)
        if annotation is None:
            annotation = etree.Element("annotation")
            insert_at = 1 if len(reaction) and etree.QName(reaction[0]).localname == "notes" else 0
            reaction.insert(insert_at, annotation)
        marker = etree.SubElement(annotation, etree.QName(NS, "gapFill"), nsmap={"msystems": NS})
        marker.set("status", "gapfilled")
        marker.set("definition", "present_in_final_model_and_absent_from_matched_draft_model")
    missing = expected - found
    if missing:
        raise RuntimeError(f"{source_name}: {len(missing)} expected gap-filled reactions were not found")
    data = etree.tostring(tree, xml_declaration=True, encoding="UTF-8", pretty_print=False)
    return data, len(reactions), len(found)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    gapfilled = load_gapfilled()
    rows: list[dict[str, object]] = []
    with tarfile.open(SOURCE_ARCHIVE, "r:gz") as source_archive:
        source_members = sorted(
            (member for member in source_archive.getmembers() if member.isfile() and member.name.endswith(".xml")),
            key=lambda member: member.name,
        )
        if len(source_members) != 1000:
            raise RuntimeError(f"Expected 1000 source models, found {len(source_members)}")
        with tarfile.open(ARCHIVE, "w:gz", compresslevel=6) as archive:
            for index, source_member in enumerate(source_members, start=1):
                source_name = Path(source_member.name).name
                genome_id = Path(source_name).stem
                expected = gapfilled.get(genome_id)
                if expected is None:
                    raise RuntimeError(f"No gap-fill evidence for {genome_id}")
                source_handle = source_archive.extractfile(source_member)
                if source_handle is None:
                    raise RuntimeError(f"Cannot read {source_member.name}")
                original = source_handle.read()
                annotated, reaction_count, marked_count = annotate_model(source_name, original, expected)
                info = tarfile.TarInfo(name=f"models/{source_name}")
                info.size = len(annotated)
                info.mtime = 0
                info.mode = 0o644
                archive.addfile(info, io.BytesIO(annotated))
                rows.append(
                    {
                        "Genome_ID": genome_id,
                        "Archive_Member": info.name,
                        "Reaction_Count": reaction_count,
                        "Gapfilled_Reactions_Marked": marked_count,
                        "Source_SHA256": sha256(original),
                        "Annotated_SHA256": sha256(annotated),
                        "Annotation_Namespace": NS,
                        "Annotation_Definition": "present in final model and absent from matched draft model",
                    }
                )
                if index % 100 == 0:
                    print(f"processed={index}", flush=True)
    with MANIFEST.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    digest = hashlib.sha256(ARCHIVE.read_bytes()).hexdigest()
    ARCHIVE_SHA.write_text(f"{digest}  {ARCHIVE.name}\n", encoding="ascii")
    print(f"archive={ARCHIVE}")
    print(f"sha256={digest}")


if __name__ == "__main__":
    main()
