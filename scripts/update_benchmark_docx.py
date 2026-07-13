from pathlib import Path
from shutil import copy2

from docx import Document
from docx.shared import Inches


ROOT = Path(r"D:\pipeline\analysis_v5")
DOCX = ROOT / "02.论文" / "Evolutionary Divergence of Genomic Defense.docx"
BACKUP = ROOT / "02.论文" / "Evolutionary Divergence of Genomic Defense.before_benchmark_update.docx"
FIG1 = ROOT / "01.figures" / "Figure1_benchmark_distribution.png"
INDEX = ROOT / "02.论文" / "docx_text_index_after_benchmark_update.txt"


REPLACEMENTS = {
    17: "Results: We developed and deployed an optimized high-throughput annotation-and-reconstruction pipeline on 10,951 high-quality E. coli genomes passing CheckM2 quality thresholds, from which 1,000 representative genomes were selected for multi-omic functional profiling and GEM reconstruction. The corrected phylogroup design substantially increased B2 representation (A = 298, B2 = 244, B1 = 200, D = 103, E = 73, C = 36, F = 31, G = 11, U = 4). In the updated benchmark, an N50-stratified 20-genome stress-test subset showed a median 15.81-fold acceleration of optimized gapseq relative to the default workflow (IQR 15.40-16.35; range 14.78-17.77), reducing median wall-clock time from 19,821.5 s to 1,260.5 s. Four available curated-reference genomes showed a median 12.76-fold acceleration (range 11.30-15.52). Because the currently available benchmark table records N50 and reaction metrics but not phylogroup/CDS metadata, we treat this as an assembly-quality-stratified benchmark and avoid claiming complete phylogroup/CDS coverage until metadata linkage is available. Functional analysis of the revised cohort showed that phylogroup membership was associated with genome architecture, secretion, membrane localization, and CRISPR spacer structure, whereas KEGG pathway composition showed weak global phylogroup separation after correction (Jaccard PERMANOVA R2 = 0.0078, p = 0.447; Hellinger PERMANOVA R2 = 0.0067, p = 0.793; PERMDISP p = 0.956).",
    32: "Here we address these challenges in a unified pipeline deployed on 1,000 representative E. coli genomes. We report: (i) a multi-layer computational optimization of gapseq evaluated as a distribution across an N50-stratified 20-genome stress-test subset, rather than as a single-genome point estimate; (ii) a rebalanced phylogroup cohort with adequate B2 representation; (iii) GEM- and pathway-level analyses that separate centroid structure from dispersion effects; and (iv) corrected COG, CRISPR spacer, SignalP, and TMHMM analyses designed to reduce over-interpretation under multiple testing and unbalanced sampling.",
    78: "Computational benchmarking was expanded from the original single MG1655 run to an N50-stratified stress-test subset of 20 E. coli assemblies selected to represent a broad range of assembly contiguity (Actual_N50 = 22,535-675,909 bp). Each genome was reconstructed with the gapseq default workflow and the optimized workflow, and wall-clock time, reaction count, speedup (default runtime / optimized runtime), reaction delta, and reaction retention were summarized as distributions. The benchmark table currently contains N50 and reaction metrics but does not encode phylogroup or CDS metadata; therefore, the analysis is reported as assembly-quality stratified rather than as a definitive phylogroup/CDS-balanced benchmark.",
    80: "Structural reference-model benchmarking was performed for the four curated-reference genomes available in the current `modle_compare` directory: iB21_1397, iEC1364_W, iECDH10B_1368, and iJN1463. For each genome, default and optimized gapseq reconstructions were compared with the corresponding BiGG reference model by reaction counts and by the fraction of reactions removed by optimization that were nevertheless present in the curated reference. The Biolog benchmarking directory also contains an extracted SD01 Table 5 binary growth matrix for BL21(DE3), W3110, and DH10B, with Nissle 1917 not recovered from the table; however, the accompanying README identifies this table as in silico growth screens rather than experimental Biolog truth data. Consequently, these files were used as a template for future multi-strain phenotype benchmarking and were not interpreted as independent experimental validation. The MG1655 M9 carbon-source analysis is retained only as a single-strain proof-of-concept.",
    97: " 3.1 Stratified Benchmarking Shows 11-18-Fold Runtime Acceleration with Compact, Reference-Aware Models",
    101: "Across the 20-genome N50-stratified stress-test subset, the optimized gapseq workflow consistently outperformed the default workflow (Figure 1A-B). Median default runtime was 19,821.5 s per genome, whereas median optimized runtime was 1,260.5 s. The resulting speedup distribution had a median of 15.81-fold (IQR 15.40-16.35), with all genomes falling between 14.78-fold and 17.77-fold. This distribution-based benchmark directly addresses the concern that pipeline performance should not be represented by a single point estimate.",
    103: "The optimized workflow also produced substantially more compact reconstructions in the stress-test subset (Figure 1C). Median reaction count decreased from 3,020.0 reactions under default gapseq to 2,373.5 reactions under the optimized workflow, corresponding to a median reaction retention of 79.1% (range 76.3-80.1%). These values indicate a reproducible speed-and-compactness shift, but they should not be interpreted alone as proof of functional equivalence; reaction-set and phenotype-level validation are required for that claim.",
    105: " 3.1.2 Reference-model structural comparison",
    107: "The four available curated-reference genomes provided an additional structural benchmark outside MG1655 (Figure 1D). Runtime acceleration remained substantial but slightly lower than in the low-contiguity stress-test subset, with a median speedup of 12.76-fold and a range of 11.30-15.52-fold. Across iB21_1397, iEC1364_W, iECDH10B_1368, and iJN1463, optimized reaction retention ranged from 77.0% to 87.5%. When reactions present in the default reconstruction but absent after optimization were compared against the matched BiGG reference models, only 10.55-17.53% of these removed reactions were present in the curated references. Thus, most removed reactions were not supported by the corresponding reference reconstructions, although the reference-supported subset of removed reactions should be audited before making strain-specific functional claims.",
    109: "Together, the stratified assembly benchmark and the reference-model comparison support the optimized workflow as a scalable reconstruction strategy for large cohorts. The evidence now rests on a distribution of runtime effects and on multiple reference genomes, not solely on MG1655. At the same time, the current benchmark metadata do not yet allow a strict claim that the 20-genome subset is balanced across all phylogroups and CDS ranges; this should be completed by linking benchmark accessions to phylogroup and CDS annotations in a future revision or supplementary table.",
    111: " 3.1.3 Phenotype benchmarking and current limits",
    113: "Phenotype benchmarking remains the most important unresolved validation layer. The previous MG1655 single-carbon-source panel under M9 constraints is retained as a proof-of-concept because it tests whether medium-aware gap-filling can recover expected substrate utilization in a well-characterized strain. However, it does not satisfy the stronger requirement for multi-strain experimental Biolog validation across W3110, BW25113, CFT073, O157:H7 Sakai, or comparable strains.",
    115: "For MG1655, the M9-constrained optimized model achieved accuracy = 82.86% and sensitivity = 96.55% (TP = 28, TN = 1, FP = 5, FN = 1; 35-carbon-source panel). In the newly supplied multi-model folder, SD01 Table 5 could be extracted for BL21(DE3), W3110, and DH10B, whereas Nissle 1917 was not found; importantly, the README identifies these entries as in silico growth screens rather than experimental Biolog measurements. We therefore do not report multi-strain experimental accuracy from these files. A publication-ready phenotype benchmark will require a true experimental truth table for at least three to five strains, after which sensitivity, specificity, balanced accuracy, MCC, and per-substrate error classes should be reported.",
    117: "Figure 1. Stratified runtime and reference-model benchmarking. (A) Speedup distribution across 20 assemblies spanning Actual_N50 = 22,535-675,909 bp. (B) Default versus optimized wall-clock runtime distributions. (C) Reaction-count compaction and reaction retention in the 20-genome stress-test subset. (D) Four available curated-reference genomes show that most reactions removed by the optimized workflow were not present in the matched BiGG reference model. The current benchmark table lacks phylogroup/CDS metadata, and the available SD01 growth matrix represents in silico screens rather than experimental Biolog truth data.",
    183: "Our updated benchmark replaces the original single-genome acceleration claim with a distribution-based assessment. Across 20 N50-stratified low-contiguity assemblies, optimized gapseq reduced median runtime from 19,821.5 s to 1,260.5 s, corresponding to a median 15.81-fold speedup and a narrow IQR of 15.40-16.35. Four curated-reference genomes showed a median 12.76-fold speedup. This stronger runtime result is more useful for cohort-scale planning than the former 4.64-fold MG1655 point estimate, but it also changes the interpretation: the benchmark currently supports assembly-quality-stratified acceleration and reference-aware structural compaction, not a fully phylogroup/CDS-balanced validation unless accession-level metadata are linked.",
    194: "In conclusion, this revised analysis presents an optimized and more carefully validated E. coli GEM reconstruction framework applied to a corrected 1,000-genome phylogroup cohort. Stratified benchmarking across 20 assemblies showed a median 15.81-fold runtime acceleration, and four available curated-reference genomes showed a median 12.76-fold acceleration with mostly non-reference-supported reactions removed by optimization. The rebalanced cohort eliminated the earlier B2 sampling defect and substantially weakened the claim of D-specific virulence enrichment, while corrected KEGG, COG, secretion, membrane, and CRISPR spacer analyses provide a more conservative picture of phylogroup-associated functional divergence. Remaining publication-critical work is concentrated in two areas: linking benchmark genomes to phylogroup/CDS metadata for a formally balanced benchmark, and adding true multi-strain experimental Biolog or carbon-source data for functional accuracy validation.",
}


def main():
    copy2(DOCX, BACKUP)
    doc = Document(DOCX)
    for idx, text in REPLACEMENTS.items():
        doc.paragraphs[idx].text = text

    if not doc.inline_shapes:
        raise RuntimeError("No inline image found to replace")
    shape = doc.inline_shapes[0]
    blip = shape._inline.graphic.graphicData.pic.blipFill.blip
    image_part = doc.part.related_parts[blip.embed]
    image_part._blob = FIG1.read_bytes()
    shape.width = Inches(6.7)
    shape.height = Inches(4.88)

    doc.save(DOCX)
    reloaded = Document(DOCX)
    with INDEX.open("w", encoding="utf-8") as f:
        for i, paragraph in enumerate(reloaded.paragraphs):
            f.write(f"{i:04d}: {paragraph.text}\n")

    print(f"saved={DOCX}")
    print(f"backup={BACKUP}")
    print(f"index={INDEX}")
    print(f"inline_shapes={len(reloaded.inline_shapes)}")


if __name__ == "__main__":
    main()
