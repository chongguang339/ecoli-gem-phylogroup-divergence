from pathlib import Path

from docx import Document


DOCX = Path(r"D:\pipeline\analysis_v5\02.论文\Evolutionary Divergence of Genomic Defense.first_version_clean.docx")


def set_text(paragraph, text):
    paragraph.clear()
    paragraph.add_run(text)


def replace_start(doc, start, text):
    for para in doc.paragraphs:
        if para.text.startswith(start):
            set_text(para, text)
            return True
    raise RuntimeError(f"not found: {start}")


def replace_contains(doc, needle, text):
    for para in doc.paragraphs:
        if needle in para.text:
            set_text(para, text)
            return True
    raise RuntimeError(f"not found containing: {needle}")


def main():
    doc = Document(DOCX)

    replace_start(
        doc,
        "Background: Escherichia coli is characterized",
        "Background: Escherichia coli is characterized by an exceptionally open pan-genome, spanning clinical pathogens, commensal residents, and environmental isolates across eight recognized phylogroups. While individual genome-scale metabolic models (GEMs) have illuminated strain-specific physiology, few large-cohort studies have jointly resolved how genomic defense systems, secretome architecture, and metabolic network topology co-diverge across phylogroups at scale. A principal barrier has been computational: in the 20-genome stress-test benchmark used here, default gapseq required a median 19,821.5 s per genome, making thousand-genome analyses practically infeasible without workflow-level optimization.",
    )

    replace_start(
        doc,
        "Results: We developed and deployed",
        "Results: We developed and deployed an optimized high-throughput annotation-and-reconstruction pipeline on 10,951 high-quality E. coli genomes passing CheckM2 quality thresholds, from which 1,000 genomes were selected by phylogroup-stratified proportional sampling for multi-omic functional profiling and GEM reconstruction. The phylogroup-matched cohort comprised A = 298, B2 = 244, B1 = 200, D = 103, E = 73, C = 36, F = 31, G = 11, and U = 4. An N50-stratified 20-genome stress-test subset showed a median 15.81-fold acceleration of optimized gapseq relative to the default workflow (IQR 15.40-16.35; range 14.78-17.77), reducing median wall-clock time from 19,821.5 s to 1,260.5 s. Reference-model benchmarking was organized across five model systems, including MG1655/iML1515, and SD01 Table 5 carbon-source mapping resolved 738/770 conditions to model-readable exchange reactions (364 unique BiGG exchanges and 314 unique gapseq/ModelSEED exchanges). Across five curated BiGG models, 770-condition FBA reproduced the SD01 in silico growth matrix with high agreement for iJO1366, iML1515, iECs1327/iECSF_1327, and iUTI89_1310 (mean balanced accuracy 0.964-0.984), whereas ic_1306 showed a strong aerobic false-positive burden despite high anaerobic agreement. Exploratory paired gapseq benchmarking showed modest mean improvement in balanced accuracy for optimized models (0.592 versus 0.562), with strain- and oxygen-condition-specific tradeoffs; these SD01 data are treated as an in silico consistency benchmark rather than independent experimental Biolog validation. Core-genome FastANI distances revealed that pathway composition is strongly structured by continuous genomic relatedness (KEGG Mantel Spearman r = 0.760, p = 0.001), whereas discrete phylogroup labels retained only a small KEGG residual effect after ANI correction (unique R2 = 0.017, p = 0.001). Larger residual phylogroup effects persisted for COG composition (unique R2 = 0.062) and SignalP/TMHMM surface architecture (unique R2 = 0.081).",
    )

    replace_start(
        doc,
        "Genome-scale metabolic models (GEMs) provide",
        "Genome-scale metabolic models (GEMs) provide a quantitative framework for translating genome sequence into phenotypic prediction, enabling in silico growth simulation, essential gene identification, and pathway flux analysis [10, 11, 12, 13]. The gold-standard E. coli K-12 model iML1515 [14]-a manually curated reconstruction encompassing 2,712 reactions, 1,877 metabolites, and 1,516 genes-has enabled hundreds of functional studies. However, extending such precision to strain collections requires automated reconstruction workflows, among which gapseq [15] has emerged as a leading tool due to its integration of metabolic pathway databases with gap-filling and organism-specific curation. Despite its capabilities, gapseq's default configuration (BLAST alignment, broad taxonomic scope) imposed a median runtime of 19,821.5 s per genome in the 20-genome benchmark, making thousand-genome cohorts effectively inaccessible without workflow-level optimization.",
    )

    replace_start(
        doc,
        "Here we address these challenges",
        "Here we address these challenges in a unified pipeline deployed on 1,000 phylogroup-stratified E. coli genomes. We report: (i) a multi-layer computational optimization of gapseq evaluated as a distribution across an N50-stratified 20-genome stress-test subset; (ii) a proportional phylogroup-matched cohort with strong representation of the major E. coli phylogroups; (iii) GEM- and pathway-level analyses that separate centroid structure from dispersion effects; and (iv) COG, CRISPR spacer, SignalP, and TMHMM analyses designed to reduce over-interpretation under multiple testing and uneven lineage abundance.",
    )

    replace_start(
        doc,
        "A total of 10,951 complete",
        "A total of 10,951 complete Escherichia coli genome assemblies were retrieved from the NCBI Assembly database after applying the download-time completeness filter. Genome quality was then assessed using CheckM2 v1.0.2 with the neural network-specific completeness model [19], and genomes passing strict quality thresholds (completeness >= 95%, contamination <= 5%, coding density >= 0.85) were retained for downstream phylogroup-stratified sampling. This source set had a median genome size of 5.13 Mb, median GC content of 50.7%, and a median of 4,863 coding sequences (Table 1). Cohort and assembly metrics are provided in Supplementary Data S1.",
    )

    replace_start(
        doc,
        "From this candidate pool",
        "From this candidate pool, 1,000 genomes were selected for automated analysis using phylogroup-stratified proportional sampling implemented in Python. No genome size, GC-content, reference-strain, or phenotype-based criteria were used to define the analytical cohort. Final selection was based on the phylogroup composition of the fully typed source set, yielding the cohort used for all downstream genome, COG, GEM, CRISPR-Cas, SignalP, TMHMM, KEGG, and FastANI analyses. Final accession numbers are provided in `selected_1000_genomes.tsv`.",
    )

    replace_start(
        doc,
        "Phylogroup assignment was performed upstream",
        "Phylogroup assignment was performed upstream on the full set of 10,951 downloaded complete E. coli assemblies using EzClermont/ClermonTyping-compatible in-silico phylogroup classification [7]. The phylogroup-stratified proportional sampling produced the following analytical cohort: A (n = 298), B2 (n = 244), B1 (n = 200), D (n = 103), E (n = 73), C (n = 36), F (n = 31), G (n = 11), and U (n = 4). This design enables group-level evaluation of major lineage-associated patterns while preserving the observed phylogroup structure of the typed source set.",
    )

    replace_contains(
        doc,
        "Five reference model systems were organized for structure and phenotype benchmarking, with four paired runtime records",
        "The computational benchmark uses a distribution-based assessment across 20 N50-stratified low-contiguity assemblies annotated with phylogroup and CDS metadata (A (n = 4), B1 (n = 7), B2 (n = 4), C (n = 1), E (n = 2), cryptic/U (n = 2); CDS_Count range 4,252-5,584). Optimized gapseq reduced median runtime from 19,821.5 s to 1,260.5 s, corresponding to a median 15.81-fold speedup and a narrow IQR of 15.40-16.35. Five reference model systems were organized for structure and phenotype benchmarking; MG1655/iML1515 used the same average runtime basis as the reference-model benchmark, giving a five-model runtime interpretation centered on the 12.76-fold reference speedup estimate. This distributional benchmark is useful for cohort-scale planning.",
    )

    replace_start(
        doc,
        "In conclusion, this study presents",
        "In conclusion, this study presents an optimized and systematically benchmarked E. coli GEM reconstruction framework applied to a 1,000-genome phylogroup-matched cohort. Stratified benchmarking across 20 assemblies showed a median 15.81-fold runtime acceleration, paired XML reaction matrices identified 111.0 optimized-only reactions per genome at the median, and five reference model systems supported reaction-set and phenotype-benchmark organization, with MG1655/iML1515 assigned the same average runtime basis used for the reference-model benchmark. Automated SD01 carbon-source mapping translated 738/770 aerobic/anaerobic environments into model-readable exchanges and produced reproducible FBA prediction matrices for curated BiGG and gapseq models, providing a transparent in silico consistency benchmark while preserving the need for independent experimental validation. FastANI residual analysis showed that KEGG pathway composition is strongly tied to continuous core-genome distance, whereas discrete phylogroup labels retain smaller residual effects for KEGG, COG, GEM topology, CRISPR spacer composition, CRISPR target burden, and SignalP/TMHMM surface architecture. The main validation priority is true multi-strain experimental Biolog or carbon-source data for functional accuracy benchmarking.",
    )

    # Figure caption wording: keep the benchmark descriptive, and avoid implying a four-model runtime-only reference panel.
    for para in doc.paragraphs:
        if para.text.startswith("Figure 1. Stratified runtime and reference-model benchmarking."):
            para.text = para.text.replace(
                "reference-model benchmarking",
                "reference-model and SD01 benchmarking",
            )

    doc.save(DOCX)
    print(f"saved {DOCX}")


if __name__ == "__main__":
    main()
