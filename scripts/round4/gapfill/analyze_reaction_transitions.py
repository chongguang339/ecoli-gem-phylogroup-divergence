from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from xml.etree import ElementTree as ET

import pandas as pd


ROOT = Path(r"D:\pipeline\analysis_v5")
OUT = ROOT / "19.mSystems_revision4_20261028" / "01_analysis"
OUT.mkdir(parents=True, exist_ok=True)
REPORTS = ROOT / "11.benchmark_20x6_provenance/03_reports"
RUNS = ROOT / "11.benchmark_20x6_provenance/remote_downloaded_20x6_runs"
DEFAULT = ROOT / "04.基准测试/test_20/results_default"


def reaction_ids(path: Path) -> set[str]:
    ids: set[str] = set()
    for _, elem in ET.iterparse(path, events=("end",)):
        if elem.tag.rsplit("}", 1)[-1] == "reaction":
            rid = elem.attrib.get("id")
            if rid:
                ids.add(rid)
        elem.clear()
    return ids


def locate(acc: str, version: str) -> tuple[Path, Path]:
    if version == "V1":
        folder = DEFAULT / acc
    else:
        corrected = RUNS / "runs_corrected_sanitized" / acc / (version + "s")
        folder = corrected if corrected.exists() else RUNS / "runs" / acc / version
    finals = [p for p in folder.glob("*.xml") if not p.name.endswith("-draft.xml")]
    drafts = list(folder.glob("*-draft.xml"))
    if len(finals) != 1 or len(drafts) != 1:
        raise FileNotFoundError(f"Expected one final/draft for {acc} {version}: {folder} {finals} {drafts}")
    return finals[0], drafts[0]


def stage(v_a_draft: bool, v_b_draft: bool) -> str:
    if v_a_draft and not v_b_draft:
        return "draft_prediction_or_evidence_stage_loss"
    if not v_a_draft and not v_b_draft:
        return "A_gapfill_or_finalization_only"
    if v_b_draft:
        return "present_in_B_draft_but_absent_B_final"
    return "other"


def main() -> None:
    delta = pd.read_csv(REPORTS / "benchmark_20x6_pairwise_reaction_delta_evidence_long.csv", low_memory=False)
    delta = delta[((delta.Comparison == "V1_vs_V2") | (delta.Comparison == "V4_vs_V5")) & (delta.Direction == "A_only")].copy()
    accessions = sorted(delta.Accession.unique())
    assert len(accessions) == 20
    caches: dict[tuple[str, str], tuple[set[str], set[str]]] = {}
    paths: list[dict] = []
    for acc in accessions:
        for v in ["V1", "V2", "V4", "V5"]:
            final, draft = locate(acc, v)
            caches[(acc, v)] = reaction_ids(final), reaction_ids(draft)
            paths.append({"Accession": acc, "Version": v, "final_xml": str(final), "draft_xml": str(draft), "final_n": len(caches[(acc, v)][0]), "draft_n": len(caches[(acc, v)][1])})
    pd.DataFrame(paths).to_csv(OUT / "paired_SBML_source_manifest.csv", index=False)
    records = []
    for row in delta.itertuples(index=False):
        va, vb = ("V1", "V2") if row.Comparison == "V1_vs_V2" else ("V4", "V5")
        rid = row.Reaction_ID
        a_fin, a_dra = caches[(row.Accession, va)]
        b_fin, b_dra = caches[(row.Accession, vb)]
        if rid not in a_fin or rid in b_fin:
            raise AssertionError(f"Delta source mismatch {row.Accession} {row.Comparison} {rid}")
        is_exchange = str(row.Base_Reaction_ID).startswith("EX_")
        records.append({
            "Accession": row.Accession,
            "Comparison": row.Comparison,
            "Reaction_ID": rid,
            "Base_Reaction_ID": row.Base_Reaction_ID,
            "Reaction_Name": row.Reaction_Name,
            "A_draft": rid in a_dra,
            "A_final": True,
            "B_draft": rid in b_dra,
            "B_final": False,
            "stage_class": stage(rid in a_dra, rid in b_dra),
            "is_exchange": is_exchange,
            "gene_product_count": row.XML_GeneProduct_Count,
            "Evidence_Class": row.Evidence_Class,
            "EC_List": row.EC_List,
            "MetaCyc_Pathway_List": row.MetaCyc_Pathway_List,
        })
    detail = pd.DataFrame(records)
    detail.to_csv(OUT / "V1_V2_and_V4_V5_reaction_stage_audit.csv", index=False)
    summary = detail.groupby(["Comparison", "stage_class", "is_exchange"], dropna=False).size().rename("n").reset_index()
    summary.to_csv(OUT / "reaction_stage_summary.csv", index=False)
    evidence = detail.groupby(["Comparison", "Evidence_Class"], dropna=False).size().rename("n").reset_index()
    evidence.to_csv(OUT / "reaction_evidence_class_summary.csv", index=False)

    # Independent protein-level alignments are not reaction-level causes. Aggregate them
    # separately and retain mixed classes rather than forcing one explanation per reaction.
    hit = pd.read_csv(ROOT / "18.mSystems_revision3_final_20260911/02_analysis_evidence/02_database_comparison/v4_only_reaction_linked_protein_hits.csv", low_memory=False)
    class_counts = hit.groupby(["Accession", "V4_Only_Reaction_ID", "Paired_Hit_Class"]).size().unstack(fill_value=0).reset_index()
    class_counts.to_csv(OUT / "V4_only_reaction_protein_hit_classes.csv", index=False)
    hit_summary = hit.groupby("Paired_Hit_Class").size().rename("n_protein_reaction_records").reset_index()
    hit_summary.to_csv(OUT / "V4_only_protein_hit_class_summary.csv", index=False)
    print(json.dumps({"accessions": len(accessions), "reaction_rows": len(detail), "stage_summary": summary.to_dict("records"), "hit_summary": hit_summary.to_dict("records")}, indent=2), flush=True)


if __name__ == "__main__":
    main()
