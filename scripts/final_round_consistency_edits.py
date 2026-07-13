from pathlib import Path

from docx import Document


DOCX = Path(r"D:\pipeline\analysis_v5\02.论文\Evolutionary Divergence of Genomic Defense.first_version_clean.docx")


def set_text(paragraph, text):
    paragraph.clear()
    paragraph.add_run(text)


def main():
    doc = Document(DOCX)
    for para in doc.paragraphs:
        text = para.text
        if text.startswith("For each of the 1,000 representative genomes"):
            set_text(
                para,
                "For each of the 1,000 phylogroup-stratified genomes, both the optimized (M9-constrained) model and its pre-filling draft model were assessed for three topological metrics: number of genes (Genes), reactions (Reactions), and metabolites (Metabolites). Biomass flux (growth rate under LB-equivalent rich medium) was computed via FBA to classify models as Viable (biomass flux > 0) or Lethal (zero growth). Phylogroup-stratified topology distributions were visualized as violin plots.",
            )
        elif text.startswith("Figure 1. Stratified runtime"):
            set_text(
                para,
                "Figure 1. Stratified runtime and reference-model and SD01 benchmarking. (A) Speedup distribution across 20 assemblies spanning Actual_N50 = 22,535-675,909 bp and CDS_Count = 4,252-5,584; the subset includes A (n = 4), B1 (n = 7), B2 (n = 4), C (n = 1), E (n = 2), cryptic/U (n = 2). (B) Default versus optimized wall-clock runtime distributions. (C) Reaction-count compaction and reaction retention in the 20-genome stress-test subset; paired XML reaction matrices further identify optimized-only and default-only reactions for every genome. (D) Five reference model systems, including MG1655/iML1515 with the benchmark-average runtime basis, show that most reactions removed by the optimized workflow were not present in the matched BiGG reference model. SD01 Table 5 carbon-source mapping and 770-condition FBA outputs are provided in Supplementary Data S5-S6 as an in silico model-consistency benchmark, not an experimental Biolog validation dataset.",
            )
        elif text.startswith("Across the 1,000 viable models in the current re-run tables"):
            set_text(
                para,
                "Across the 1,000 viable models in the cohort-level topology tables, topology statistics were: median GEM_Genes = 8,007, median GEM_Reactions = 2,403, and median GEM_Metabolites = 12,354. Reaction counts were compact across the cohort (IQR: 2,380-2,421; range: 2,285-2,480), while metabolite counts ranged from 11,756 to 12,734. These values describe the aggregated GEM topology fields in `recalculated_annotations.csv` and are interpreted as cohort-level model topology metrics rather than one-to-one equivalents of manually curated iML1515 identifier counts.",
            )
        elif text.startswith("The computational benchmark uses a distribution-based assessment"):
            set_text(
                para,
                "The computational benchmark uses a distribution-based assessment across 20 N50-stratified low-contiguity assemblies annotated with phylogroup and CDS metadata (A (n = 4), B1 (n = 7), B2 (n = 4), C (n = 1), E (n = 2), cryptic/U (n = 2); CDS_Count range 4,252-5,584). Optimized gapseq reduced median runtime from 19,821.5 s to 1,260.5 s, corresponding to a median 15.81-fold speedup and a narrow IQR of 15.40-16.35. Five reference model systems were organized for structure and phenotype benchmarking; MG1655/iML1515 used the same benchmark-average runtime basis, giving a five-model reference interpretation centered on the 12.76-fold speedup estimate. This distributional benchmark is useful for cohort-scale planning.",
            )
        elif text.startswith("In conclusion, this study presents"):
            set_text(
                para,
                "In conclusion, this study presents an optimized and systematically benchmarked E. coli GEM reconstruction framework applied to a 1,000-genome phylogroup-matched cohort. Stratified benchmarking across 20 assemblies showed a median 15.81-fold runtime acceleration, paired XML reaction matrices identified 111.0 optimized-only reactions per genome at the median, and five reference model systems supported reaction-set and phenotype-benchmark organization, with MG1655/iML1515 using the same benchmark-average runtime basis. Automated SD01 carbon-source mapping translated 738/770 aerobic/anaerobic environments into model-readable exchanges and produced reproducible FBA prediction matrices for curated BiGG and gapseq models, providing a transparent in silico consistency benchmark while preserving the need for independent experimental validation. FastANI residual analysis showed that KEGG pathway composition is strongly tied to continuous core-genome distance, whereas discrete phylogroup labels retain smaller residual effects for KEGG, COG, GEM topology, CRISPR spacer composition, CRISPR target burden, and SignalP/TMHMM surface architecture. The main validation priority is true multi-strain experimental Biolog or carbon-source data for functional accuracy benchmarking.",
            )
    doc.save(DOCX)
    print(f"saved {DOCX}")


if __name__ == "__main__":
    main()
