from __future__ import annotations

import csv
import shutil
from pathlib import Path

import pandas as pd
from docx import Document


ROOT = Path(r"D:\pipeline\analysis_v5")
DOCX = ROOT / "02.论文" / "Evolutionary Divergence of Genomic Defense.docx"
BACKUP = ROOT / "02.论文" / "Evolutionary Divergence of Genomic Defense.before_reviewer_issue_cleanup.docx"
DELTA = ROOT / "04.基准测试" / "test_20" / "reaction_set_delta" / "reaction_set_delta_summary.csv"
RECURRENCE = ROOT / "04.基准测试" / "test_20" / "reaction_set_delta" / "optimized_only_reaction_recurrence.csv"
OUT_CHECKLIST = ROOT / "03.analysis" / "reviewer_issue_resolution_checklist.csv"
OUT_TEXT_INDEX = ROOT / "02.论文" / "docx_text_index_after_reviewer_issue_cleanup.txt"


def fmt(x: float, digits: int = 2) -> str:
    return f"{x:.{digits}f}"


def replace_exact(doc: Document, old: str, new: str) -> None:
    for p in doc.paragraphs:
        if p.text == old:
            p.text = new
            return
    raise ValueError(f"paragraph not found: {old[:80]}")


def main() -> None:
    if not BACKUP.exists():
        shutil.copy2(DOCX, BACKUP)

    delta = pd.read_csv(DELTA)
    recurrence = pd.read_csv(RECURRENCE)
    recurrence_count_col = "Genome_Count" if "Genome_Count" in recurrence.columns else "n_genomes"
    universal = int((recurrence[recurrence_count_col] == len(delta)).sum())
    unique_opt = int(recurrence.shape[0])

    opt_only_min = int(delta["Optimized_only_Reactions"].min())
    opt_only_med = float(delta["Optimized_only_Reactions"].median())
    opt_only_max = int(delta["Optimized_only_Reactions"].max())
    def_only_min = int(delta["Default_only_Reactions"].min())
    def_only_med = float(delta["Default_only_Reactions"].median())
    def_only_max = int(delta["Default_only_Reactions"].max())
    jacc_med = float(delta["Jaccard_similarity_exact"].median())
    opt_frac_med = 100 * float(delta["Optimized_only_fraction_of_optimized"].median())
    def_frac_med = 100 * float(delta["Default_only_fraction_of_default"].median())

    doc = Document(DOCX)

    replace_exact(
        doc,
        "A total of 10,951 Escherichia coli genome assemblies were retrieved from the NCBI Assembly database. Genome quality was assessed using CheckM2 v1.0.2 with the neural network-specific completeness model [12]. Only genomes passing strict quality thresholds (completeness ≥ 95%, contamination ≤ 5%, coding density ≥ 0.85) were retained. This resulted in a high-quality cohort of 10,951 genomes with a median genome size of 5.13 Mb, median GC content of 50.7%, and a median of 4,863 coding sequences (Table 2). Representative genome quality metrics are available in `quality_report.tsv`.",
        "The quality-filtered source dataset used for sampling comprised 10,951 Escherichia coli genome assemblies from the NCBI Assembly database. Genome quality was assessed using CheckM2 v1.0.2 with the neural network-specific completeness model [12], and only genomes passing strict quality thresholds (completeness >= 95%, contamination <= 5%, coding density >= 0.85) were retained for downstream sampling. This source set had a median genome size of 5.13 Mb, median GC content of 50.7%, and a median of 4,863 coding sequences (Table 2). Representative genome quality metrics are available in `quality_report.tsv`; the pre-filter retrieval count should be reported if it differs from this quality-filtered source set.",
    )

    replace_exact(
        doc,
        "Computational benchmarking was expanded from the original single MG1655 run to an N50-stratified stress-test subset of 20 E. coli assemblies selected to represent a broad range of assembly contiguity (Actual_N50 = 22,535-675,909 bp). Each genome was reconstructed with the gapseq default workflow and the optimized workflow, and wall-clock time, reaction count, speedup (default runtime / optimized runtime), reaction delta, and reaction retention were summarized as distributions. The benchmark table currently contains N50 and reaction metrics but does not encode phylogroup or CDS metadata; therefore, the analysis is reported as assembly-quality stratified rather than as a definitive phylogroup/CDS-balanced benchmark.",
        "Computational benchmarking was expanded from the original single MG1655 run to an N50-stratified stress-test subset of 20 E. coli assemblies selected to represent a broad range of assembly contiguity (Actual_N50 = 22,535-675,909 bp). Each genome was reconstructed with the gapseq default workflow and the optimized workflow, and wall-clock time, reaction count, speedup (default runtime / optimized runtime), reaction delta, and reaction retention were summarized as distributions. To test whether model compaction reflected simple truncation or a systematic change in reaction content, final SBML/XML `listOfReactions` entries were parsed for each paired reconstruction, and optimized-only, default-only, shared-reaction, and exact Jaccard matrices were generated. The benchmark table currently contains N50 and reaction metrics but does not encode phylogroup or CDS metadata; therefore, the analysis is reported as assembly-quality stratified rather than as a definitive phylogroup/CDS-balanced benchmark.",
    )

    replace_exact(
        doc,
        "The optimized workflow also produced substantially more compact reconstructions in the stress-test subset (Figure 1C). Median reaction count decreased from 3,020.0 reactions under default gapseq to 2,373.5 reactions under the optimized workflow, corresponding to a median reaction retention of 79.1% (range 76.3-80.1%). These values indicate a reproducible speed-and-compactness shift, but they should not be interpreted alone as proof of functional equivalence; reaction-set and phenotype-level validation are required for that claim.",
        f"The optimized workflow also produced substantially more compact reconstructions in the stress-test subset (Figure 1C). Median reaction count decreased from 3,020.0 reactions under default gapseq to 2,373.5 reactions under the optimized workflow, corresponding to a median reaction retention of 79.1% (range 76.3-80.1%). Reaction-set comparison showed that this compaction was not a simple deletion-only process: optimized models contained {opt_only_min}-{opt_only_max} optimized-only reactions per genome (median {fmt(opt_only_med, 1)}; median {fmt(opt_frac_med, 2)}% of optimized reactions), including {unique_opt} unique optimized-only reactions across the subset and {universal} reactions recurrent in all 20 genomes. Conversely, default models contained {def_only_min}-{def_only_max} default-only reactions per genome (median {fmt(def_only_med, 1)}; median {fmt(def_frac_med, 2)}% of default reactions), and the median exact reaction-set Jaccard similarity between paired models was {fmt(jacc_med, 3)}. These values support a reproducible speed-and-compactness shift with systematic reaction redistribution, but they should not be interpreted alone as proof of functional equivalence; phenotype-level validation remains required for that claim.",
    )

    replace_exact(
        doc,
        "Together, the stratified assembly benchmark and the reference-model comparison support the optimized workflow as a scalable reconstruction strategy for large cohorts. The evidence now rests on a distribution of runtime effects and on multiple reference genomes, not solely on MG1655. At the same time, the current benchmark metadata do not yet allow a strict claim that the 20-genome subset is balanced across all phylogroups and CDS ranges; this should be completed by linking benchmark accessions to phylogroup and CDS annotations in a future revision or supplementary table.",
        "Together, the stratified assembly benchmark, paired reaction-set matrices, and reference-model comparison support the optimized workflow as a scalable reconstruction strategy for large cohorts. The evidence now rests on a distribution of runtime effects, explicit optimized-only/default-only reaction matrices, and multiple reference genomes, not solely on MG1655. At the same time, the current benchmark metadata do not yet allow a strict claim that the 20-genome subset is balanced across all phylogroups and CDS ranges; this should be completed by linking benchmark accessions to phylogroup and CDS annotations in a future revision or supplementary table.",
    )

    replace_exact(
        doc,
        "Figure 1. Stratified runtime and reference-model benchmarking. (A) Speedup distribution across 20 assemblies spanning Actual_N50 = 22,535-675,909 bp. (B) Default versus optimized wall-clock runtime distributions. (C) Reaction-count compaction and reaction retention in the 20-genome stress-test subset. (D) Four available curated-reference genomes show that most reactions removed by the optimized workflow were not present in the matched BiGG reference model. The current benchmark table lacks phylogroup/CDS metadata, and the available SD01 growth matrix represents in silico screens rather than experimental Biolog truth data.",
        "Figure 1. Stratified runtime and reference-model benchmarking. (A) Speedup distribution across 20 assemblies spanning Actual_N50 = 22,535-675,909 bp. (B) Default versus optimized wall-clock runtime distributions. (C) Reaction-count compaction and reaction retention in the 20-genome stress-test subset; paired XML reaction matrices further identify optimized-only and default-only reactions for every genome. (D) Four available curated-reference genomes show that most reactions removed by the optimized workflow were not present in the matched BiGG reference model. The current benchmark table lacks phylogroup/CDS metadata, and the available SD01 growth matrix represents in silico screens rather than experimental Biolog truth data.",
    )

    replace_exact(
        doc,
        "Several limitations remain. First, pathway names are represented by KEGG identifiers only; pathway-name annotation would improve biological interpretability of the 15 corrected signals. Second, the residual phylogroup analysis uses ANI PCoA covariates and label-permutation sensitivity tests rather than a full core-gene tree with formal phylogenetic comparative modeling; it therefore distinguishes continuous genomic relatedness from categorical labels but does not replace tree-based evolutionary inference. Third, although TMHMM coverage was improved from 846 to 999 genomes after reprocessing, one B1 genome still lacks detailed TMHMM output and was treated as missing in membrane-protein analyses. Fourth, CRISPR spacer sharing is sequence-level but not yet matched to phage/protospacer databases, so spacer targets and ecological exposures remain unresolved. Finally, all GEM predictions remain in silico and should be experimentally validated in representative strains.",
        "Several limitations remain. First, pathway names are represented by KEGG identifiers only; pathway-name annotation would improve biological interpretability of the 15 corrected signals. Second, phylogroup assignment currently relies on ClermonTyping output, and an independent concordance analysis using EzClermont, MLST, Mash/ANI placement against typed references, or a core-gene tree should be added before making strong taxonomic claims. Third, the residual phylogroup analysis uses ANI PCoA covariates and label-permutation sensitivity tests rather than a full core-gene tree with formal phylogenetic comparative modeling; it therefore distinguishes continuous genomic relatedness from categorical labels but does not replace tree-based evolutionary inference. Fourth, although TMHMM coverage was improved from 846 to 999 genomes after reprocessing, one B1 genome still lacks detailed TMHMM output and was treated as missing in membrane-protein analyses; a DeepTMHMM or equivalent modern topology reanalysis would further strengthen membrane-protein inference. Fifth, CRISPR spacer sharing is sequence-level but not yet matched to phage/protospacer databases, so spacer targets and ecological exposures remain unresolved. Finally, all GEM predictions remain in silico and should be experimentally validated in representative strains.",
    )

    replace_exact(
        doc,
        "In conclusion, this revised analysis presents an optimized and more carefully validated E. coli GEM reconstruction framework applied to a corrected 1,000-genome phylogroup cohort. Stratified benchmarking across 20 assemblies showed a median 15.81-fold runtime acceleration, and four available curated-reference genomes showed a median 12.76-fold acceleration with mostly non-reference-supported reactions removed by optimization. The rebalanced cohort eliminated the earlier B2 sampling defect and substantially weakened the claim of D-specific virulence enrichment. FastANI residual analysis further showed that KEGG pathway composition is strongly tied to continuous core-genome distance, whereas discrete phylogroup labels retain smaller residual effects for KEGG, COG, GEM topology, CRISPR spacer composition, and SignalP/TMHMM surface architecture. Remaining publication-critical work is concentrated in adding true multi-strain experimental Biolog or carbon-source data for functional accuracy validation.",
        f"In conclusion, this revised analysis presents an optimized and more carefully validated E. coli GEM reconstruction framework applied to a corrected 1,000-genome phylogroup cohort. Stratified benchmarking across 20 assemblies showed a median 15.81-fold runtime acceleration, paired XML reaction matrices identified {fmt(opt_only_med, 1)} optimized-only reactions per genome at the median, and four available curated-reference genomes showed a median 12.76-fold acceleration with mostly non-reference-supported reactions removed by optimization. The rebalanced cohort eliminated the earlier B2 sampling defect and substantially weakened the claim of D-specific virulence enrichment. FastANI residual analysis further showed that KEGG pathway composition is strongly tied to continuous core-genome distance, whereas discrete phylogroup labels retain smaller residual effects for KEGG, COG, GEM topology, CRISPR spacer composition, and SignalP/TMHMM surface architecture. The main remaining validation priority is true multi-strain experimental Biolog or carbon-source data for functional accuracy benchmarking.",
    )

    doc.save(DOCX)

    issues = [
        ("2.1 Cohort design and B2 underrepresentation", "resolved in analysis/manuscript", "Re-run cohort now includes B2 n=244 and D-specific virulence claims are weakened.", "Independent phylogroup concordance using EzClermont/MLST/Mash/core-gene placement is still needed for the strongest taxonomic claims."),
        ("2.2 Benchmark distribution", "partially resolved", "20 N50-stratified genomes now report speedup distribution and four reference-model genomes.", "Benchmark accession-to-phylogroup/CDS metadata are needed to claim formal phylogroup/CDS-balanced benchmarking."),
        ("2.3 No-content-loss overclaim", "resolved with caution", "Paired XML reaction-set matrices quantify optimized-only, default-only, shared reactions, and Jaccard similarity.", "Functional equivalence still requires phenotype truth data."),
        ("2.4 Biolog/carbon-source truth data", "not resolved", "MG1655 M9 panel retained only as proof-of-concept; in silico SD01 table not used as experimental truth.", "Experimental Biolog or curated carbon-source growth matrices for at least 3-5 strains such as W3110, BW25113, CFT073, O157:H7 Sakai."),
        ("2.5 PERMANOVA with PERMDISP failure", "resolved", "Balanced PERMANOVA, PERMDISP warning, and cautious interpretation added; R2 not treated as pure explained variance.", "None for current scope."),
        ("2.6 Circularity of phylogroup/pathway structure", "largely resolved", "FastANI distance matrix used for Mantel and residual phylogroup models.", "A formal core-gene tree model would further strengthen evolutionary inference."),
        ("2.7 CRISPR count underpowering", "partially resolved", "Spacer-level sequence analysis replaces simple array/Cas counts for multivariate structure.", "Protospacer/phage target matching is still needed for ecological interpretation."),
        ("2.8 Coevolution/tradeoff overinterpretation", "resolved by text weakening", "Discussion now uses modular, conservative interpretation rather than a single host-pathogen axis.", "None for current scope."),
        ("3.1 Internal consistency", "partially resolved", "Old B2 counts and old 56.1% claim removed; 10,951 wording changed to quality-filtered source set.", "Need original pre-QC retrieval count if it differs from 10,951."),
        ("3.2 ClermonTyping validation", "not resolved", "Limitation now explicitly states independent concordance is required.", "EzClermont/MLST/Mash/core-gene typed-reference output."),
        ("3.3 SignalP/TMHMM normalization", "partially resolved", "Per-1000-CDS normalization and effect sizes are reported.", "Optional per-Mb and genome-size residual models can be added if requested."),
        ("3.4 TMHMM v2.0 depreciation", "not resolved", "TMHMM missingness corrected to 999/1000 and limitation now recommends DeepTMHMM/equivalent reanalysis.", "DeepTMHMM or equivalent modern topology outputs."),
        ("3.5 KEGG pathway count/source", "resolved", "Manuscript now states 264 non-redundant KO columns after duplicate map-column removal.", "None for current scope."),
        ("3.6 COG multiple testing", "resolved", "COG tests now report BH and Bonferroni correction across 26 categories.", "None for current scope."),
        ("3.7 BiGG structural comparison", "partially resolved", "Four available curated-reference genomes are analyzed and reaction-set deltas are provided.", "For stronger claims, add cross-namespace reaction mapping and manual audit of reference-supported removed reactions."),
        ("4 Reporting placeholders", "not resolved", "Manuscript still requires final repository URL/DOI, author details, funding, and final supplement packaging.", "Repository URL/DOI, funding statement, author contributions, conflict statement, and final supplementary file list."),
    ]
    OUT_CHECKLIST.parent.mkdir(parents=True, exist_ok=True)
    with OUT_CHECKLIST.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow(["Reviewer_issue", "Status", "Current_resolution", "Data_or_action_still_needed"])
        writer.writerows(issues)

    with OUT_TEXT_INDEX.open("w", encoding="utf-8") as f:
        for i, p in enumerate(Document(DOCX).paragraphs):
            f.write(f"{i:04d}: {p.text}\n")

    print(f"Updated {DOCX}")
    print(f"Backup {BACKUP}")
    print(f"Wrote {OUT_CHECKLIST}")
    print(f"Optimized-only reactions per genome: {opt_only_min}-{opt_only_max}, median {opt_only_med}")
    print(f"Default-only reactions per genome: {def_only_min}-{def_only_max}, median {def_only_med}")
    print(f"Unique optimized-only reactions: {unique_opt}; universal: {universal}")


if __name__ == "__main__":
    main()
