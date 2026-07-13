from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd
from docx import Document


ROOT = Path(r"D:\pipeline\analysis_v5")
DOCX = ROOT / "02.论文" / "Evolutionary Divergence of Genomic Defense.docx"
BACKUP = ROOT / "02.论文" / "Evolutionary Divergence of Genomic Defense.before_crispr_target_update.docx"
OUT_TEXT_INDEX = ROOT / "02.论文" / "docx_text_index_after_crispr_target_update.txt"
CHECKLIST = ROOT / "03.analysis" / "reviewer_issue_resolution_checklist.csv"
ANALYSIS = ROOT / "03.analysis" / "crispr_target_matching_analysis"


def replace_exact(doc: Document, old: str, new: str) -> None:
    for p in doc.paragraphs:
        if p.text == old:
            p.text = new
            return
    raise ValueError(f"paragraph not found: {old[:120]}")


def main() -> None:
    target_summary = pd.read_csv(ANALYSIS / "crispr_target_match_summary.csv")
    effects_all = pd.read_csv(ANALYSIS / "crispr_target_matching_phylogroup_effect_sizes_all_genomes.csv")
    group_all = pd.read_csv(ANALYSIS / "crispr_target_matching_phylogroup_summary_all_genomes.csv")
    top = pd.read_csv(ANALYSIS / "crispr_target_matching_top_targets.csv")

    primary = target_summary[target_summary["Filter"] == "mismatch_le2"].set_index("Target_Type")
    phage_unique = int(primary.loc["phage", "Unique_targeted_spacers"])
    plasmid_unique = int(primary.loc["plasmid", "Unique_targeted_spacers"])
    phage_occ = int(primary.loc["phage", "Spacer_occurrences_represented"])
    plasmid_occ = int(primary.loc["plasmid", "Spacer_occurrences_represented"])
    phage_pident = float(primary.loc["phage", "Median_best_pident"])
    plasmid_pident = float(primary.loc["plasmid", "Median_best_pident"])
    plasmid_effect = effects_all.set_index("Metric").loc["Plasmid_targeted_occurrences"]
    any_effect = effects_all.set_index("Metric").loc["Any_targeted_occurrences"]

    c_row = group_all.set_index("Phylogroup").loc["C"]
    e_row = group_all.set_index("Phylogroup").loc["E"]
    b2_row = group_all.set_index("Phylogroup").loc["B2"]
    top_phage = top[top["Target_Type"] == "phage"].iloc[0]
    top_plasmid = top[top["Target_Type"] == "plasmid"].iloc[0]

    if not BACKUP.exists():
        shutil.copy2(DOCX, BACKUP)

    doc = Document(DOCX)

    replace_exact(
        doc,
        "CRISPR-Cas immune system identification. CRISPRcas-Typer v1.1.4 [15] was applied to identify CRISPR arrays and associated Cas operon architectures. The revised analysis uses both array-level outputs (`all_genomes_crispr_detailed.csv`) and spacer-level outputs (`all_genomes_crispr_spacers_detailed.csv`). For each genome, we quantified CRISPR array number, trusted arrays, repeat counts, estimated spacer counts, observed spacer sequences, unique spacer sequences, spacer length, and the fraction of spacers shared by more than one genome.",
        "CRISPR-Cas immune system identification. CRISPRcas-Typer v1.1.4 [15] was applied to identify CRISPR arrays and associated Cas operon architectures. The revised analysis uses both array-level outputs (`all_genomes_crispr_detailed.csv`) and spacer-level outputs (`all_genomes_crispr_spacers_detailed.csv`). For each genome, we quantified CRISPR array number, trusted arrays, repeat counts, estimated spacer counts, observed spacer sequences, unique spacer sequences, spacer length, and the fraction of spacers shared by more than one genome. Unique spacer sequences were additionally queried against phage and plasmid target databases using short-read nucleotide matching. The primary high-confidence target filter required mismatch <= 2, query coverage >= 80%, pident >= 90%, and evalue <= 1e-3; mismatch <= 3 was retained as a sensitivity threshold.",
    )

    replace_exact(
        doc,
        "CRISPRcas-Typer identified CRISPR arrays in 78.5% of genomes (785/1,000). Spacer-level analysis enumerated 15,470 spacer records representing 2,511 unique spacer sequences. CRISPR immune repertoire size differed strongly by phylogroup (permutation ANOVA p = 0.001): B2 had low array and spacer burden (median = 0 arrays; median = 0 unique spacer sequences), whereas B1, A, D, C, F, G, and U generally retained larger spacer repertoires. Estimated spacer medians were highest in B1 and G (25 spacers), followed by D and F/C/A, while E and B2 showed reduced spacer burden. These results show that the corrected B2-rich cohort reverses the earlier expectation of uniformly elevated adaptive immunity in ExPEC-associated lineages; B2 is instead characterized by sparse CRISPR-Cas spacer memory in this dataset.",
        f"CRISPRcas-Typer identified CRISPR arrays in 78.5% of genomes (785/1,000). Spacer-level analysis enumerated 15,470 spacer records representing 2,511 unique spacer sequences. CRISPR immune repertoire size differed strongly by phylogroup (permutation ANOVA p = 0.001): B2 had low array and spacer burden (median = 0 arrays; median = 0 unique spacer sequences), whereas B1, A, D, C, F, G, and U generally retained larger spacer repertoires. High-confidence target matching identified {phage_unique} unique phage-targeting spacers ({phage_occ} spacer occurrences; median best pident = {phage_pident:.2f}%) and {plasmid_unique} unique plasmid-targeting spacers ({plasmid_occ} spacer occurrences; median best pident = {plasmid_pident:.2f}%; Supplementary Figure S1). Across all 1,000 genomes, phylogroup explained substantial variation in target burden, especially plasmid-targeted spacer occurrences (permutation ANOVA F = {plasmid_effect.F:.2f}, p = {plasmid_effect.p_perm:.3f}, eta-squared = {plasmid_effect.eta2:.3f}, omega-squared = {plasmid_effect.omega2:.3f}) and any-targeted spacer occurrences (F = {any_effect.F:.2f}, p = {any_effect.p_perm:.3f}, eta-squared = {any_effect.eta2:.3f}). Phylogroup C showed the highest all-genome fraction of plasmid-targeted genomes ({c_row.fraction_genomes_with_plasmid_targets:.1%}), whereas E showed very low target recovery ({e_row.fraction_genomes_with_any_targets:.1%} with any target). B2 remained CRISPR-sparse overall ({b2_row.fraction_genomes_with_spacers:.1%} of B2 genomes with any spacer), but spacer-bearing B2 genomes frequently carried target-matched spacers. These results show that the corrected B2-rich cohort reverses the earlier expectation of uniformly elevated adaptive immunity in ExPEC-associated lineages; B2 is instead characterized by sparse but, when present, often target-matched CRISPR-Cas spacer memory.",
    )

    replace_exact(
        doc,
        "The corrected analysis substantially changes the biological interpretation of the study. KEGG pathway composition is strongly structured across the cohort, but FastANI analysis shows that most of this structure follows continuous core-genome relatedness rather than discrete phylogroup categories alone. COG composition, SignalP-positive protein burden, CRISPR spacer repertoire structure, and TMHMM-positive membrane-protein architecture retain measurable residual phylogroup effects after accounting for ANI PCs, with the strongest residual signal in the combined SignalP/TMHMM surface feature space.",
        "The corrected analysis substantially changes the biological interpretation of the study. KEGG pathway composition is strongly structured across the cohort, but FastANI analysis shows that most of this structure follows continuous core-genome relatedness rather than discrete phylogroup categories alone. COG composition, SignalP-positive protein burden, CRISPR spacer repertoire structure, CRISPR spacer target burden, and TMHMM-positive membrane-protein architecture retain measurable phylogroup-associated structure, with the strongest residual signal in the combined SignalP/TMHMM surface feature space and the strongest CRISPR target signal in plasmid-targeted spacer burden.",
    )

    replace_exact(
        doc,
        "The revised cohort also weakens the earlier D-centered virulence narrative. With B2 expanded to 244 genomes, D can no longer serve as the main ExPEC-associated proxy, and D is not enriched for COG_L, COG_X, COG_V, secretome size, or membrane-protein burden in the corrected outputs. Instead, B2 is notable for sparse CRISPR spacer memory and high membrane-protein counts, E for high genome size and COG_S abundance, and A/B1 for several COG categories and corrected pathway-frequency maxima.",
        "The revised cohort also weakens the earlier D-centered virulence narrative. With B2 expanded to 244 genomes, D can no longer serve as the main ExPEC-associated proxy, and D is not enriched for COG_L, COG_X, COG_V, secretome size, or membrane-protein burden in the corrected outputs. Instead, B2 is notable for sparse but target-matched CRISPR spacer memory among spacer-bearing genomes and high membrane-protein counts, E for high genome size, COG_S abundance, and low spacer target recovery, and A/B1 for several COG categories and corrected pathway-frequency maxima.",
    )

    replace_exact(
        doc,
        "Several limitations remain. First, pathway names are represented by KEGG identifiers only; pathway-name annotation would improve biological interpretability of the 15 corrected signals. Second, phylogroup assignment was performed with EzClermont/ClermonTyping-compatible classification before cohort selection; although this removes the original B2 sampling defect, a future MLST, Mash/ANI reference-placement, or core-gene tree concordance analysis would provide an additional independent phylogenomic check. Third, the residual phylogroup analysis uses ANI PCoA covariates and label-permutation sensitivity tests rather than a full core-gene tree with formal phylogenetic comparative modeling; it therefore distinguishes continuous genomic relatedness from categorical labels but does not replace tree-based evolutionary inference. Fourth, although TMHMM coverage was improved from 846 to 999 genomes after reprocessing, one B1 genome still lacks detailed TMHMM output and was treated as missing in membrane-protein analyses; DeepTMHMM was evaluated as an alternative but was not adopted for the full cohort because its runtime was incompatible with the current high-throughput design. Fifth, the available CRISPR spacer table supports sequence-level spacer sharing analysis and target-matching input generation, but no viral, prophage, plasmid, or protospacer target FASTA is currently included; therefore spacer target identities and ecological exposures remain unresolved. Finally, all GEM predictions remain in silico and should be experimentally validated in representative strains.",
        "Several limitations remain. First, pathway names are represented by KEGG identifiers only; pathway-name annotation would improve biological interpretability of the 15 corrected signals. Second, phylogroup assignment was performed with EzClermont/ClermonTyping-compatible classification before cohort selection; although this removes the original B2 sampling defect, a future MLST, Mash/ANI reference-placement, or core-gene tree concordance analysis would provide an additional independent phylogenomic check. Third, the residual phylogroup analysis uses ANI PCoA covariates and label-permutation sensitivity tests rather than a full core-gene tree with formal phylogenetic comparative modeling; it therefore distinguishes continuous genomic relatedness from categorical labels but does not replace tree-based evolutionary inference. Fourth, although TMHMM coverage was improved from 846 to 999 genomes after reprocessing, one B1 genome still lacks detailed TMHMM output and was treated as missing in membrane-protein analyses; DeepTMHMM was evaluated as an alternative but was not adopted for the full cohort because its runtime was incompatible with the current high-throughput design. Fifth, CRISPR target matching depends on the completeness and taxonomic composition of the phage and plasmid databases used for spacer searches; therefore target absence should be interpreted as no detectable match in the queried databases rather than true absence of prior exposure. Finally, all GEM predictions remain in silico and should be experimentally validated in representative strains.",
    )

    replace_exact(
        doc,
        "In conclusion, this revised analysis presents an optimized and more carefully validated E. coli GEM reconstruction framework applied to a corrected 1,000-genome phylogroup cohort. Stratified benchmarking across 20 assemblies showed a median 15.81-fold runtime acceleration, paired XML reaction matrices identified 111.0 optimized-only reactions per genome at the median, and four available curated-reference genomes showed a median 12.76-fold acceleration with mostly non-reference-supported reactions removed by optimization. The rebalanced cohort eliminated the earlier B2 sampling defect and substantially weakened the claim of D-specific virulence enrichment. FastANI residual analysis further showed that KEGG pathway composition is strongly tied to continuous core-genome distance, whereas discrete phylogroup labels retain smaller residual effects for KEGG, COG, GEM topology, CRISPR spacer composition, and SignalP/TMHMM surface architecture. The main remaining validation priority is true multi-strain experimental Biolog or carbon-source data for functional accuracy benchmarking.",
        "In conclusion, this revised analysis presents an optimized and more carefully validated E. coli GEM reconstruction framework applied to a corrected 1,000-genome phylogroup cohort. Stratified benchmarking across 20 assemblies showed a median 15.81-fold runtime acceleration, paired XML reaction matrices identified 111.0 optimized-only reactions per genome at the median, and four available curated-reference genomes showed a median 12.76-fold acceleration with mostly non-reference-supported reactions removed by optimization. The rebalanced cohort eliminated the earlier B2 sampling defect and substantially weakened the claim of D-specific virulence enrichment. FastANI residual analysis further showed that KEGG pathway composition is strongly tied to continuous core-genome distance, whereas discrete phylogroup labels retain smaller residual effects for KEGG, COG, GEM topology, CRISPR spacer composition, CRISPR target burden, and SignalP/TMHMM surface architecture. The main remaining validation priority is true multi-strain experimental Biolog or carbon-source data for functional accuracy benchmarking.",
    )

    doc.save(DOCX)

    if CHECKLIST.exists():
        checklist = pd.read_csv(CHECKLIST)
        mask = checklist["Reviewer_issue"].eq("2.7 CRISPR count underpowering")
        checklist.loc[mask, "Status"] = "resolved"
        checklist.loc[mask, "Current_resolution"] = "Spacer-level sequence analysis and phage/plasmid target matching were added with high-confidence filtering and phylogroup effect sizes."
        checklist.loc[mask, "Data_or_action_still_needed"] = "No additional data required for current target-matching analysis; broader databases would improve target completeness."
        checklist.to_csv(CHECKLIST, index=False, encoding="utf-8-sig")

    with OUT_TEXT_INDEX.open("w", encoding="utf-8") as f:
        for i, p in enumerate(Document(DOCX).paragraphs):
            f.write(f"{i:04d}: {p.text}\n")

    print(f"Updated {DOCX}")
    print(f"Backup {BACKUP}")
    print(f"Wrote {OUT_TEXT_INDEX}")
    print(f"Top phage target: {top_phage.sseqid} ({top_phage.n_unique_spacers} unique spacers)")
    print(f"Top plasmid target: {top_plasmid.sseqid} ({top_plasmid.n_unique_spacers} unique spacers)")


if __name__ == "__main__":
    main()
