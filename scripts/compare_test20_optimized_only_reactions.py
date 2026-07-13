from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd


ROOT = Path(r"D:\pipeline\analysis_v5")
TEST20 = ROOT / "04.基准测试" / "test_20"
OUT = TEST20 / "reaction_set_delta"
OUT.mkdir(exist_ok=True)


def find_final_xml(folder: Path) -> Path:
    xmls = [p for p in folder.glob("*.xml") if "-draft" not in p.name and "_draft" not in p.name]
    if len(xmls) != 1:
        raise FileNotFoundError(f"Expected exactly one final XML in {folder}, found {len(xmls)}: {xmls}")
    return xmls[0]


def reaction_base_id(rid: str) -> str:
    rid = re.sub(r"^R_", "", rid)
    rid = re.sub(r"_[a-z]\d+$", "", rid)
    return rid


def read_reactions(xml_path: Path) -> dict[str, str]:
    root = ET.parse(xml_path).getroot()
    reactions: dict[str, str] = {}
    for elem in root.iter():
        if elem.tag.split("}")[-1] == "reaction":
            rid = elem.attrib.get("id", "").strip()
            if not rid:
                continue
            reactions[rid] = elem.attrib.get("name", "").strip()
    return reactions


def matrix_from_long(long_df: pd.DataFrame, id_col: str, name_col: str, genomes: list[str], out_path: Path):
    if long_df.empty:
        pd.DataFrame(columns=[id_col, name_col, "n_genomes", *genomes]).to_csv(out_path, index=False)
        return
    mat = (
        long_df.assign(value=1)
        .pivot_table(index=[id_col, name_col], columns="Genome_ID", values="value", aggfunc="max", fill_value=0)
        .reindex(columns=genomes, fill_value=0)
        .astype(int)
    )
    mat.insert(0, "n_genomes", mat.sum(axis=1))
    mat = mat.reset_index().sort_values(["n_genomes", id_col], ascending=[False, True])
    mat.to_csv(out_path, index=False)


def main():
    metrics = pd.read_csv(TEST20 / "benchmark_metrics.csv")
    genomes = metrics["Accession"].astype(str).tolist()
    summary_rows = []
    opt_only_rows = []
    default_only_rows = []
    opt_only_base_rows = []

    for genome in genomes:
        default_xml = find_final_xml(TEST20 / "results_default" / genome)
        optimized_xml = find_final_xml(TEST20 / "results_optimized" / genome)
        default = read_reactions(default_xml)
        optimized = read_reactions(optimized_xml)
        default_ids = set(default)
        optimized_ids = set(optimized)
        opt_only = sorted(optimized_ids - default_ids)
        default_only = sorted(default_ids - optimized_ids)
        shared = sorted(optimized_ids & default_ids)

        default_base = {reaction_base_id(x) for x in default_ids}
        optimized_base = {reaction_base_id(x) for x in optimized_ids}
        opt_only_base = sorted(optimized_base - default_base)
        base_name = {}
        for rid, name in optimized.items():
            base_name.setdefault(reaction_base_id(rid), name)

        summary_rows.append(
            {
                "Genome_ID": genome,
                "Default_XML": str(default_xml),
                "Optimized_XML": str(optimized_xml),
                "Default_Reactions_XML": len(default_ids),
                "Optimized_Reactions_XML": len(optimized_ids),
                "Shared_Reactions": len(shared),
                "Optimized_only_Reactions": len(opt_only),
                "Default_only_Reactions": len(default_only),
                "Optimized_only_Base_Reactions": len(opt_only_base),
                "Jaccard_similarity_exact": len(shared) / len(default_ids | optimized_ids),
            }
        )
        for rid in opt_only:
            opt_only_rows.append(
                {
                    "Genome_ID": genome,
                    "Reaction_ID": rid,
                    "Base_Reaction_ID": reaction_base_id(rid),
                    "Reaction_Name": optimized.get(rid, ""),
                }
            )
        for rid in default_only:
            default_only_rows.append(
                {
                    "Genome_ID": genome,
                    "Reaction_ID": rid,
                    "Base_Reaction_ID": reaction_base_id(rid),
                    "Reaction_Name": default.get(rid, ""),
                }
            )
        for bid in opt_only_base:
            opt_only_base_rows.append(
                {
                    "Genome_ID": genome,
                    "Base_Reaction_ID": bid,
                    "Representative_Reaction_Name": base_name.get(bid, ""),
                }
            )

    summary = pd.DataFrame(summary_rows).merge(metrics, left_on="Genome_ID", right_on="Accession", how="left")
    summary["Runtime_speedup"] = summary["Default_Time_Sec"] / summary["Optimized_Time_Sec"]
    summary["Optimized_only_fraction_of_optimized"] = (
        summary["Optimized_only_Reactions"] / summary["Optimized_Reactions_XML"]
    )
    summary["Default_only_fraction_of_default"] = (
        summary["Default_only_Reactions"] / summary["Default_Reactions_XML"]
    )
    summary.to_csv(OUT / "reaction_set_delta_summary.csv", index=False)

    opt_long = pd.DataFrame(opt_only_rows)
    def_long = pd.DataFrame(default_only_rows)
    opt_base_long = pd.DataFrame(opt_only_base_rows)
    opt_long.to_csv(OUT / "optimized_only_reactions_long.csv", index=False)
    def_long.to_csv(OUT / "default_only_reactions_long.csv", index=False)
    opt_base_long.to_csv(OUT / "optimized_only_base_reactions_long.csv", index=False)

    matrix_from_long(
        opt_long,
        "Reaction_ID",
        "Reaction_Name",
        genomes,
        OUT / "optimized_only_reaction_matrix.csv",
    )
    matrix_from_long(
        opt_base_long,
        "Base_Reaction_ID",
        "Representative_Reaction_Name",
        genomes,
        OUT / "optimized_only_base_reaction_matrix.csv",
    )

    # A compact recurrence table is convenient for supplement text.
    recurrence = (
        opt_long.groupby(["Reaction_ID", "Reaction_Name"], dropna=False)["Genome_ID"]
        .agg(n_genomes="nunique", genomes=lambda s: ";".join(sorted(s)))
        .reset_index()
        .sort_values(["n_genomes", "Reaction_ID"], ascending=[False, True])
    )
    recurrence.to_csv(OUT / "optimized_only_reaction_recurrence.csv", index=False)

    print(f"outputs: {OUT}")
    print("summary")
    cols = [
        "Genome_ID",
        "Default_Reactions_XML",
        "Optimized_Reactions_XML",
        "Shared_Reactions",
        "Optimized_only_Reactions",
        "Default_only_Reactions",
        "Optimized_only_Base_Reactions",
        "Jaccard_similarity_exact",
        "Runtime_speedup",
    ]
    print(summary[cols].to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    print("\noptimized-only recurrence top 20")
    if recurrence.empty:
        print("No optimized-only reactions.")
    else:
        print(recurrence.head(20).to_string(index=False))


if __name__ == "__main__":
    main()
