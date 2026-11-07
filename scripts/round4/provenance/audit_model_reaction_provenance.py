from __future__ import annotations

import json
import tarfile
from collections import Counter
from pathlib import Path
from xml.etree import ElementTree as ET

import pandas as pd


ROOT = Path(r"D:\pipeline\analysis_v5")
OUT = ROOT / "19.mSystems_revision4_20261028/01_analysis"
ARCHIVE = ROOT / "18.mSystems_revision3_final_20260911/05_model_archive/Supplementary_Models_1000_gapfill_annotated.tar.gz"
MANIFEST = ROOT / "10.mSystems_revision2_model_evidence/02_model_evidence/reaction_evidence_counts_by_genome.csv"


def local(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def reaction_kind(rid: str, name: str) -> str:
    s = f"{rid} {name}".lower()
    if rid.startswith("R_EX_") or " exchange" in s:
        return "exchange"
    if "biomass" in s or rid.lower().startswith("r_bio"):
        return "biomass"
    if "transport" in s:
        return "transport_named"
    if "spontaneous" in s:
        return "spontaneous_named"
    return "other_or_unresolved"


def main() -> None:
    rows = []
    files = 0
    with tarfile.open(ARCHIVE, "r|gz") as tar:
        for member in tar:
            if not member.isfile() or not member.name.endswith(".xml"):
                continue
            file = tar.extractfile(member)
            assert file
            genome = Path(member.name).stem
            counter = Counter()
            for _, element in ET.iterparse(file, events=("end",)):
                if local(element.tag) != "reaction":
                    continue
                rid = element.attrib.get("id", "")
                name = element.attrib.get("name", "")
                has_gene = any(local(c.tag) == "geneProductAssociation" for c in element.iter())
                gapfilled = any(local(c.tag) == "gapFill" and c.attrib.get("status") == "gapfilled" for c in element.iter())
                kind = reaction_kind(rid, name)
                counter[(has_gene, gapfilled, kind)] += 1
                element.clear()
            for (gene, gap, kind), n in counter.items():
                rows.append({"Genome_ID": genome, "gene_association": gene, "gapfilled": gap, "reaction_class_heuristic": kind, "n": n})
            files += 1
            if files % 100 == 0:
                print(f"Parsed {files} models", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "reaction_provenance_cross_classification_by_genome.csv", index=False)
    totals = df.groupby(["gene_association", "gapfilled", "reaction_class_heuristic"], dropna=False).n.sum().reset_index()
    totals.to_csv(OUT / "reaction_provenance_cross_classification_total.csv", index=False)
    per = df.groupby(["Genome_ID", "gene_association", "gapfilled"]).n.sum().unstack([1, 2], fill_value=0)
    per.columns = [f"gene_{str(a).lower()}_gapfilled_{str(b).lower()}" for a, b in per.columns]
    per = per.reset_index()
    per.to_csv(OUT / "reaction_provenance_four_way_by_genome.csv", index=False)
    original = pd.read_csv(MANIFEST)
    merged = original.merge(per, on="Genome_ID", validate="one_to_one")
    assert len(merged) == 1000
    assert (merged[[c for c in per if c != "Genome_ID"]].sum(axis=1) == merged.Final_Reactions).all()
    assert (merged[[c for c in per if c.startswith("gene_true")]].sum(axis=1) == merged.Gene_Evidenced_Reactions).all()
    assert (merged[[c for c in per if c.endswith("gapfilled_true")]].sum(axis=1) == merged.Gap_Filled_Reactions).all()
    summary = {"models": files, "reaction_occurrences": int(df.n.sum()), "four_way_totals": {c: int(per[c].sum()) for c in per if c != "Genome_ID"}}
    (OUT / "reaction_provenance_audit_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
