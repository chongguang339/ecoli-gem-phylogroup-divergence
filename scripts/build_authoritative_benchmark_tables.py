from __future__ import annotations

import csv
import os
import re
import statistics
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path


ROOT = Path(os.environ.get("ECOLI_GEM_ROOT", Path(__file__).resolve().parents[1]))
OUT = ROOT / "18.mSystems_revision3_final_20260911" / "02_analysis_evidence"
RUN_ROOT = ROOT / "11.benchmark_20x6_provenance" / "remote_downloaded_20x6_runs"
SOURCE_MANIFEST = (
    ROOT
    / "11.benchmark_20x6_provenance"
    / "03_reports"
    / "benchmark_20x6_runtime_reaction_provenance_corrected_integrated.csv"
)
REACTION_EVIDENCE = (
    ROOT
    / "11.benchmark_20x6_provenance"
    / "03_reports"
    / "benchmark_20x6_reaction_presence_evidence_long.csv"
)
EC_RECURRENCE = (
    ROOT / "11.benchmark_20x6_provenance" / "03_reports" / "benchmark_20x6_ec_recurrence_by_run.csv"
)
METACYC_RECURRENCE = (
    ROOT / "11.benchmark_20x6_provenance" / "03_reports" / "benchmark_20x6_metacyc_recurrence_by_run.csv"
)
HIT_TABLE = (
    ROOT
    / "17.broad_focused_protein_evidence"
    / "remote_downloaded"
    / "outputs_20260910"
    / "tables"
    / "broad_vs_focused_best_hit_per_query.csv"
)


def q1(values: list[float]) -> float:
    return statistics.quantiles(values, n=4, method="inclusive")[0]


def q3(values: list[float]) -> float:
    return statistics.quantiles(values, n=4, method="inclusive")[2]


def find_final_xml(accession: str, version: str, corrected: bool) -> Path:
    candidates: list[Path] = []
    if version == "V1":
        candidates = list((ROOT / "04.基准测试" / "test_20" / "results_default" / accession).glob("*.xml"))
    elif version == "V6":
        candidates = list((ROOT / "04.基准测试" / "test_20" / "results_optimized" / accession).glob("*.xml"))
    elif corrected:
        candidates = list((RUN_ROOT / "runs_corrected_sanitized" / accession / f"{version}s").glob("*.xml"))
    else:
        candidates = list((RUN_ROOT / "runs" / accession / version).glob("*.xml"))
    finals = [p for p in candidates if not p.name.endswith("-draft.xml")]
    if len(finals) != 1:
        raise RuntimeError(f"Expected one final XML for {accession} {version}; found {finals}")
    return finals[0]


def find_draft_xml(final_xml: Path) -> Path:
    draft = final_xml.with_name(final_xml.stem + "-draft.xml")
    if not draft.exists():
        raise FileNotFoundError(draft)
    return draft


def reaction_ids(path: Path) -> set[str]:
    root = ET.parse(path).getroot()
    reaction_list = next(node for node in root.iter() if node.tag.endswith("listOfReactions"))
    return {
        node.attrib["id"]
        for node in list(reaction_list)
        if node.tag.endswith("reaction") and "id" in node.attrib
    }


def base_reaction_id(reaction_id: str) -> str:
    value = reaction_id.removeprefix("R_")
    return re.sub(r"_[a-z]\d+$", "", value)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_annotation_sets() -> dict[tuple[str, str], tuple[set[str], set[str]]]:
    ec_pattern = re.compile(r"(?<![0-9])(?:[0-9]+|-)\.(?:[0-9]+|-)\.(?:[0-9]+|-)\.(?:[0-9]+|-)(?![0-9])")
    values: dict[tuple[str, str], tuple[set[str], set[str]]] = defaultdict(lambda: (set(), set()))
    with EC_RECURRENCE.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            key = (row["Accession"], row["Version_ID"])
            values[key][0].update(ec_pattern.findall(row.get("EC", "") or ""))
    with METACYC_RECURRENCE.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            key = (row["Accession"], row["Version_ID"])
            raw = (row.get("MetaCyc_Pathway") or "").strip().strip("|")
            if raw and raw not in {"/", "//", "NA"}:
                values[key][1].add(raw)
    return values


def build_ablation_tables() -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    annotations = parse_annotation_sets()
    runs: list[dict[str, object]] = []
    for row in read_csv(SOURCE_MANIFEST):
        accession = row["Accession"]
        version = row["Version_ID"]
        corrected = row.get("Corrected_Rerun", "").lower() == "true"
        final_xml = find_final_xml(accession, version, corrected)
        draft_xml = find_draft_xml(final_xml)
        final_ids = reaction_ids(final_xml)
        draft_ids = reaction_ids(draft_xml)
        ecs, pathways = annotations.get((accession, version), (set(), set()))
        runs.append(
            {
                "Accession": accession,
                "Version_ID": version,
                "Version_Label": row["Version_Label"],
                "Wallclock_Seconds": float(row["Wallclock_Seconds"]),
                "Final_Reactions_SBML": len(final_ids),
                "Draft_Reactions_SBML": len(draft_ids),
                "Gapfilled_Reactions_DraftFinal": len(final_ids - draft_ids),
                "Unique_EC_Numbers": len(ecs),
                "Unique_MetaCyc_Pathways": len(pathways),
                "Corrected_Rerun": corrected,
                "Final_XML_Relative_Path": str(final_xml.relative_to(ROOT)),
                "Reaction_Count_Method": "direct child reaction elements in SBML listOfReactions",
            }
        )

    summaries: list[dict[str, object]] = []
    for version in ["V1", "V2", "V3", "V4", "V5", "V6"]:
        subset = [r for r in runs if r["Version_ID"] == version]
        summary: dict[str, object] = {
            "Version_ID": version,
            "Version_Label": subset[0]["Version_Label"],
            "N": len(subset),
        }
        for column in [
            "Wallclock_Seconds",
            "Final_Reactions_SBML",
            "Gapfilled_Reactions_DraftFinal",
            "Unique_EC_Numbers",
            "Unique_MetaCyc_Pathways",
        ]:
            vals = [float(r[column]) for r in subset]
            summary[f"{column}_Median"] = statistics.median(vals)
            summary[f"{column}_Q1"] = q1(vals)
            summary[f"{column}_Q3"] = q3(vals)
        summaries.append(summary)

    write_csv(OUT / "01_ablation" / "benchmark_20x6_authoritative_run_table.csv", runs)
    write_csv(OUT / "01_ablation" / "benchmark_20x6_authoritative_version_summary.csv", summaries)
    return runs, summaries


def load_hits() -> dict[str, dict[str, str]]:
    return {row["Query_ID"]: row for row in read_csv(HIT_TABLE)}


def locate_v4_table(accession: str, corrected: bool) -> Path:
    if corrected:
        folder = RUN_ROOT / "runs_corrected_sanitized" / accession / "V4s"
    else:
        folder = RUN_ROOT / "runs" / accession / "V4"
    tables = list(folder.glob("*-all-Reactions.tbl"))
    if len(tables) != 1:
        raise RuntimeError(f"Expected one V4 reaction table for {accession}; found {tables}")
    return tables[0]


def load_reaction_annotations(path: Path, accession: str) -> dict[str, dict[str, set[str]]]:
    mapping: dict[str, dict[str, set[str]]] = defaultdict(
        lambda: {"queries": set(), "ecs": set(), "pathways": set()}
    )
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        header = None
        for line in handle:
            if line.startswith("pathway\t"):
                header = line.rstrip("\n\r").split("\t")
                break
        if header is None:
            raise RuntimeError(f"No table header in {path}")
        reader = csv.DictReader(handle, fieldnames=header, delimiter="\t")
        for row in reader:
            protein = (row.get("stitle") or "").split(" # ", 1)[0].strip()
            if not protein:
                continue
            query_id = f"{accession}|{protein}"
            for rid in (row.get("dbhit") or "").split():
                mapping[rid]["queries"].add(query_id)
                ec = (row.get("ec") or "").strip().strip("/")
                if ec:
                    mapping[rid]["ecs"].update(part for part in re.split(r"[/;, ]+", ec) if part)
                pathway = (row.get("pathway") or "").strip().strip("|")
                if pathway:
                    mapping[rid]["pathways"].add(pathway)
    return mapping


def build_database_evidence(runs: list[dict[str, object]]) -> None:
    hit_by_query = load_hits()
    run_lookup = {(str(r["Accession"]), str(r["Version_ID"])): r for r in runs}
    corrected_accessions = {
        str(r["Accession"])
        for r in runs
        if r["Version_ID"] == "V4" and bool(r["Corrected_Rerun"])
    }
    detailed: list[dict[str, object]] = []
    reaction_rows: list[dict[str, object]] = []
    all_accessions = sorted({str(r["Accession"]) for r in runs})
    for accession in all_accessions:
        v4 = Path(ROOT / str(run_lookup[(accession, "V4")]["Final_XML_Relative_Path"]))
        v5 = Path(ROOT / str(run_lookup[(accession, "V5")]["Final_XML_Relative_Path"]))
        v4_only = sorted(reaction_ids(v4) - reaction_ids(v5))
        reaction_map = load_reaction_annotations(locate_v4_table(accession, accession in corrected_accessions), accession)
        for reaction_id in v4_only:
            base_id = base_reaction_id(reaction_id)
            annotations = reaction_map.get(base_id, {"queries": set(), "ecs": set(), "pathways": set()})
            queries = sorted(annotations["queries"])
            reaction_rows.append(
                {
                    "Accession": accession,
                    "V4_Only_Reaction_ID": reaction_id,
                    "Base_Reaction_ID": base_id,
                    "Direct_Protein_Query_Count": len(queries),
                    "EC_Numbers": ";".join(sorted(annotations["ecs"])),
                    "MetaCyc_Pathways": ";".join(sorted(annotations["pathways"])),
                }
            )
            for query_id in queries:
                hit = hit_by_query.get(query_id)
                detailed.append(
                    {
                        "Accession": accession,
                        "V4_Only_Reaction_ID": reaction_id,
                        "Base_Reaction_ID": base_id,
                        "Query_ID": query_id,
                        "Protein_Token": query_id.split("|", 1)[1],
                        "Paired_Hit_Class": hit.get("Paired_Hit_Class", "not_searched") if hit else "not_searched",
                        "Broad_HighConfidence": hit.get("Broad_HighConfidence", "") if hit else "",
                        "Broad_pident": hit.get("Broad_pident", "") if hit else "",
                        "Broad_qcovhsp": hit.get("Broad_qcovhsp", "") if hit else "",
                        "Broad_scovhsp": hit.get("Broad_scovhsp", "") if hit else "",
                        "Broad_evalue": hit.get("Broad_evalue", "") if hit else "",
                        "Focused_HighConfidence": hit.get("Focused_HighConfidence", "") if hit else "",
                        "Focused_pident": hit.get("Focused_pident", "") if hit else "",
                        "Focused_qcovhsp": hit.get("Focused_qcovhsp", "") if hit else "",
                        "Focused_scovhsp": hit.get("Focused_scovhsp", "") if hit else "",
                        "Focused_evalue": hit.get("Focused_evalue", "") if hit else "",
                    }
                )

    class_counts: dict[str, set[str]] = defaultdict(set)
    for row in detailed:
        class_counts[str(row["Paired_Hit_Class"])].add(str(row["Query_ID"]))
    unique_queries = {str(row["Query_ID"]) for row in detailed}
    summary = [
        {"Metric": "V4_only_reaction_occurrences", "Value": len(reaction_rows)},
        {"Metric": "V4_only_unique_reactions", "Value": len({r["Base_Reaction_ID"] for r in reaction_rows})},
        {"Metric": "V4_only_reaction_occurrences_with_direct_protein", "Value": sum(int(r["Direct_Protein_Query_Count"]) > 0 for r in reaction_rows)},
        {"Metric": "V4_only_unique_direct_proteins", "Value": len(unique_queries)},
        {
            "Metric": "V4_only_unique_EC_numbers",
            "Value": len({ec for row in reaction_rows for ec in str(row["EC_Numbers"]).split(";") if ec}),
        },
        {
            "Metric": "V4_only_unique_MetaCyc_pathways",
            "Value": len({p for row in reaction_rows for p in str(row["MetaCyc_Pathways"]).split(";") if p}),
        },
    ]
    for cls in ["both_high_confidence", "broad_only_high_confidence", "focused_only_high_confidence", "neither_high_confidence", "not_searched"]:
        n = len(class_counts.get(cls, set()))
        summary.append({"Metric": f"reaction_linked_proteins_{cls}", "Value": n})
        summary.append({"Metric": f"reaction_linked_proteins_{cls}_percent", "Value": 100 * n / len(unique_queries) if unique_queries else 0})
    summary.append({"Metric": "high_confidence_threshold", "Value": "identity>=40%; query coverage>=70%; subject coverage>=50%; E-value<=1e-10"})

    write_csv(OUT / "02_database_comparison" / "v4_only_reactions_by_genome.csv", reaction_rows)
    write_csv(OUT / "02_database_comparison" / "v4_only_reaction_linked_protein_hits.csv", detailed)
    write_csv(OUT / "02_database_comparison" / "v4_v5_reaction_linked_database_summary.csv", summary)


def main() -> None:
    runs, _ = build_ablation_tables()
    build_database_evidence(runs)
    print(OUT)


if __name__ == "__main__":
    main()
