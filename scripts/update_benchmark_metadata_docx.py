from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd
from docx import Document


ROOT = Path(r"D:\pipeline\analysis_v5")
DOCX = ROOT / "02.论文" / "Evolutionary Divergence of Genomic Defense.docx"
BACKUP = ROOT / "02.论文" / "Evolutionary Divergence of Genomic Defense.before_benchmark_metadata_update.docx"
META = ROOT / "00.data" / "all_20_genomes_metadata.csv"
DELTA = ROOT / "04.基准测试" / "test_20" / "reaction_set_delta" / "reaction_set_delta_summary.csv"
OUT = ROOT / "03.analysis" / "benchmark_20_with_phylogroup_cds_reaction_delta.csv"
CHECKLIST = ROOT / "03.analysis" / "reviewer_issue_resolution_checklist.csv"
OUT_TEXT_INDEX = ROOT / "02.论文" / "docx_text_index_after_benchmark_metadata_update.txt"


def replace_exact(doc: Document, old: str, new: str) -> None:
    for p in doc.paragraphs:
        if p.text == old:
            p.text = new
            return
    raise ValueError(f"paragraph not found: {old[:120]}")


def main() -> None:
    meta = pd.read_csv(META)
    delta = pd.read_csv(DELTA)
    merged = delta.merge(meta, on="Accession", how="left")
    cols = [
        "Accession",
        "Genome_ID",
        "Phylogroup",
        "CDS_Count",
        "Actual_N50",
        "Default_Time_Sec",
        "Optimized_Time_Sec",
        "Runtime_speedup",
        "Default_Reactions_XML",
        "Optimized_Reactions_XML",
        "Shared_Reactions",
        "Optimized_only_Reactions",
        "Default_only_Reactions",
        "Jaccard_similarity_exact",
        "Optimized_only_fraction_of_optimized",
        "Default_only_fraction_of_default",
    ]
    merged[cols].to_csv(OUT, index=False, encoding="utf-8-sig")

    phy_counts = meta["Phylogroup"].value_counts().to_dict()
    explicit_order = ["A", "B1", "B2", "C", "E"]
    count_text_parts = [f"{g} (n = {phy_counts[g]})" for g in explicit_order if g in phy_counts]
    cryptic_n = int(sum(v for k, v in phy_counts.items() if str(k).lower() in {"cryptic", "u/cryptic", "u"}))
    if cryptic_n:
        count_text_parts.append(f"cryptic/U (n = {cryptic_n})")
    count_text = ", ".join(count_text_parts)
    cds_min = int(meta["CDS_Count"].min())
    cds_max = int(meta["CDS_Count"].max())
    cds_med = float(meta["CDS_Count"].median())

    if not BACKUP.exists():
        shutil.copy2(DOCX, BACKUP)

    doc = Document(DOCX)
    replace_exact(
        doc,
        "Computational benchmarking was expanded from the original single MG1655 run to an N50-stratified stress-test subset of 20 E. coli assemblies selected to represent a broad range of assembly contiguity (Actual_N50 = 22,535-675,909 bp). Each genome was reconstructed with the gapseq default workflow and the optimized workflow, and wall-clock time, reaction count, speedup (default runtime / optimized runtime), reaction delta, and reaction retention were summarized as distributions. To test whether model compaction reflected simple truncation or a systematic change in reaction content, final SBML/XML `listOfReactions` entries were parsed for each paired reconstruction, and optimized-only, default-only, shared-reaction, and exact Jaccard matrices were generated. The benchmark table currently contains N50 and reaction metrics but does not encode phylogroup or CDS metadata; therefore, the analysis is reported as assembly-quality stratified rather than as a definitive phylogroup/CDS-balanced benchmark.",
        f"Computational benchmarking was expanded from the original single MG1655 run to an N50-stratified stress-test subset of 20 E. coli assemblies selected to represent a broad range of assembly contiguity (Actual_N50 = 22,535-675,909 bp). Post hoc metadata annotation showed that this subset covered {count_text}, with CDS counts ranging from {cds_min:,} to {cds_max:,} (median {cds_med:,.1f}). Each genome was reconstructed with the gapseq default workflow and the optimized workflow, and wall-clock time, reaction count, speedup (default runtime / optimized runtime), reaction delta, and reaction retention were summarized as distributions. To test whether model compaction reflected simple truncation or a systematic change in reaction content, final SBML/XML `listOfReactions` entries were parsed for each paired reconstruction, and optimized-only, default-only, shared-reaction, and exact Jaccard matrices were generated. Because D, F, and G were not represented in the 20-genome stress-test subset, the benchmark is reported as an N50-stratified, phylogroup/CDS-annotated performance test rather than a fully balanced all-phylogroup benchmark.",
    )
    replace_exact(
        doc,
        "Together, the stratified assembly benchmark, paired reaction-set matrices, and reference-model comparison support the optimized workflow as a scalable reconstruction strategy for large cohorts. The evidence now rests on a distribution of runtime effects, explicit optimized-only/default-only reaction matrices, and multiple reference genomes, not solely on MG1655. At the same time, the current benchmark metadata do not yet allow a strict claim that the 20-genome subset is balanced across all phylogroups and CDS ranges; this should be completed by linking benchmark accessions to phylogroup and CDS annotations in a future revision or supplementary table.",
        "Together, the stratified assembly benchmark, paired reaction-set matrices, and reference-model comparison support the optimized workflow as a scalable reconstruction strategy for large cohorts. The evidence now rests on a distribution of runtime effects, explicit optimized-only/default-only reaction matrices, accession-level phylogroup/CDS metadata, and multiple reference genomes, not solely on MG1655. The 20-genome stress-test subset should still be interpreted as an N50-stratified benchmark with phylogroup/CDS annotation, not as a fully balanced all-phylogroup benchmark, because D, F, and G are absent from this small performance subset.",
    )
    replace_exact(
        doc,
        "Figure 1. Stratified runtime and reference-model benchmarking. (A) Speedup distribution across 20 assemblies spanning Actual_N50 = 22,535-675,909 bp. (B) Default versus optimized wall-clock runtime distributions. (C) Reaction-count compaction and reaction retention in the 20-genome stress-test subset; paired XML reaction matrices further identify optimized-only and default-only reactions for every genome. (D) Four available curated-reference genomes show that most reactions removed by the optimized workflow were not present in the matched BiGG reference model. The current benchmark table lacks phylogroup/CDS metadata, and the available SD01 growth matrix represents in silico screens rather than experimental Biolog truth data.",
        f"Figure 1. Stratified runtime and reference-model benchmarking. (A) Speedup distribution across 20 assemblies spanning Actual_N50 = 22,535-675,909 bp and CDS_Count = {cds_min:,}-{cds_max:,}; the subset includes {count_text}. (B) Default versus optimized wall-clock runtime distributions. (C) Reaction-count compaction and reaction retention in the 20-genome stress-test subset; paired XML reaction matrices further identify optimized-only and default-only reactions for every genome. (D) Four available curated-reference genomes show that most reactions removed by the optimized workflow were not present in the matched BiGG reference model. The available SD01 growth matrix represents in silico screens rather than experimental Biolog truth data.",
    )
    replace_exact(
        doc,
        "Our updated benchmark replaces the original single-genome acceleration claim with a distribution-based assessment. Across 20 N50-stratified low-contiguity assemblies, optimized gapseq reduced median runtime from 19,821.5 s to 1,260.5 s, corresponding to a median 15.81-fold speedup and a narrow IQR of 15.40-16.35. Four curated-reference genomes showed a median 12.76-fold speedup. This stronger runtime result is more useful for cohort-scale planning than the former 4.64-fold MG1655 point estimate, but it also changes the interpretation: the benchmark currently supports assembly-quality-stratified acceleration and reference-aware structural compaction, not a fully phylogroup/CDS-balanced validation unless accession-level metadata are linked.",
        f"Our updated benchmark replaces the original single-genome acceleration claim with a distribution-based assessment. Across 20 N50-stratified low-contiguity assemblies annotated with phylogroup and CDS metadata ({count_text}; CDS_Count range {cds_min:,}-{cds_max:,}), optimized gapseq reduced median runtime from 19,821.5 s to 1,260.5 s, corresponding to a median 15.81-fold speedup and a narrow IQR of 15.40-16.35. Four curated-reference genomes showed a median 12.76-fold speedup. This stronger runtime result is more useful for cohort-scale planning than the former 4.64-fold MG1655 point estimate, but it should be interpreted as an N50-stratified performance benchmark with post hoc lineage/CDS annotation rather than as a deliberately balanced all-phylogroup benchmark.",
    )
    doc.save(DOCX)

    if CHECKLIST.exists():
        df = pd.read_csv(CHECKLIST)
        mask = df["Reviewer_issue"].eq("2.2 Benchmark distribution")
        df.loc[mask, "Status"] = "resolved with scope qualification"
        df.loc[mask, "Current_resolution"] = (
            "20 N50-stratified genomes now report speedup distribution, reaction-set deltas, "
            "and accession-level phylogroup/CDS metadata."
        )
        df.loc[mask, "Data_or_action_still_needed"] = (
            "No further data needed for an N50-stratified annotated benchmark; do not claim a fully "
            "balanced all-phylogroup benchmark because D, F, and G are absent from the 20-genome subset."
        )
        df.to_csv(CHECKLIST, index=False, encoding="utf-8-sig")

    with OUT_TEXT_INDEX.open("w", encoding="utf-8") as f:
        for i, p in enumerate(Document(DOCX).paragraphs):
            f.write(f"{i:04d}: {p.text}\n")

    print(f"Wrote {OUT}")
    print(f"Updated {DOCX}")
    print(f"Backup {BACKUP}")
    print(f"Phylogroup counts: {phy_counts}")
    print(f"CDS range: {cds_min}-{cds_max}, median {cds_med}")


if __name__ == "__main__":
    main()
