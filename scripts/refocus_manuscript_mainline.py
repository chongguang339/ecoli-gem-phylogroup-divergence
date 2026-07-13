from pathlib import Path
import shutil

from docx import Document


DOCX = Path(r"D:\pipeline\analysis_v5\02.论文\Evolutionary Divergence of Genomic Defense.docx")
BACKUP = DOCX.with_name("Evolutionary Divergence of Genomic Defense.before_refocus_mainline.docx")


def set_text(paragraph, text):
    paragraph.clear()
    paragraph.add_run(text)


def replace_start(doc, start, text):
    for para in doc.paragraphs:
        if para.text.strip().startswith(start):
            set_text(para, text)
            return True
    raise RuntimeError(f"not found: {start}")


def main():
    shutil.copy2(DOCX, BACKUP)
    doc = Document(DOCX)

    # Title: single main line.
    set_text(
        doc.paragraphs[0],
        "High-Throughput GEM Workflow Optimization for Phylogroup-Resolved Functional Divergence Analysis in Escherichia coli",
    )

    replace_start(
        doc,
        "Background:",
        "Background: Escherichia coli is characterized by an exceptionally open pan-genome spanning commensal, environmental, and pathogenic lineages. Genome-scale metabolic models (GEMs) provide a framework for translating this genomic diversity into functional hypotheses, but default automated reconstruction remains computationally expensive at cohort scale. The central question of this study is whether workflow-level GEM optimization can enable reproducible, phylogroup-resolved functional analysis across a 1,000-genome E. coli cohort, while keeping auxiliary defense and surface-feature annotations subordinate to the GEM-centered analysis.",
    )

    replace_start(
        doc,
        "Results:",
        "Results: We developed and deployed an optimized high-throughput annotation-and-reconstruction pipeline on 10,951 high-quality E. coli genomes passing CheckM2 quality thresholds, from which 1,000 genomes were selected by phylogroup-stratified proportional sampling for GEM reconstruction and functional profiling. The final analytical cohort comprised A = 298, B2 = 244, B1 = 200, D = 103, E = 73, C = 36, F = 31, G = 11, and U = 4. An N50-stratified 20-genome stress-test subset showed a median 15.81-fold acceleration of optimized gapseq relative to the default workflow (IQR 15.40-16.35; range 14.78-17.77), reducing median wall-clock time from 19,821.5 s to 1,260.5 s. Reference-model benchmarking was organized across five model systems, including MG1655/iML1515, with measured 20-genome runtime distributions reported separately from the benchmark-average runtime basis used for iML1515 harmonization. SD01 Table 5 carbon-source mapping resolved 738/770 conditions to model-readable exchange reactions. Because SD01 Table 5 is an in silico literature matrix rather than independent experimental Biolog data, these simulations are interpreted as model-consistency evidence rather than phenotype-accuracy validation. Core-genome FastANI distances revealed that pathway composition is strongly structured by continuous genomic relatedness (KEGG Mantel Spearman r = 0.760, p = 0.001), whereas discrete phylogroup labels retained only a small KEGG residual effect after ANI correction (unique R2 = 0.017, p = 0.001). Auxiliary COG, CRISPR spacer, SignalP, and TMHMM layers were used to contextualize GEM and pathway divergence rather than as independent mechanistic endpoints.",
    )

    replace_start(
        doc,
        "Conclusions:",
        "Conclusions: This work establishes a GEM-centered, phylogroup-resolved workflow for large-scale E. coli functional divergence analysis. The optimized reconstruction pipeline enables cohort-scale GEM generation, while FastANI residual analysis shows that pathway composition is driven largely by continuous core-genome relatedness with smaller residual contributions from discrete phylogroup labels. CRISPR, secretion-signal, membrane-protein, and COG annotations provide supporting functional context, but substrate-level predictive performance still requires independent experimental carbon-source validation.",
    )

    replace_start(
        doc,
        "Here we address these challenges",
        "Here we address these challenges in a unified pipeline deployed on 1,000 phylogroup-stratified E. coli genomes. The central aim is to test whether workflow-level optimization of automated GEM reconstruction can support cohort-scale, phylogroup-aware metabolic and pathway analysis without overstating substrate-level predictive performance. CRISPR-Cas, SignalP, TMHMM, and COG analyses are included as auxiliary annotation layers to contextualize the GEM and KEGG results, not as separate mechanistic claims. We report: (i) a multi-layer computational optimization of gapseq evaluated as a distribution across an N50-stratified 20-genome stress-test subset; (ii) a proportional phylogroup-matched cohort generated from the full EzClermont-typed source set; (iii) GEM- and pathway-level analyses that separate centroid structure from dispersion effects and continuous ANI relatedness; and (iv) focused auxiliary annotation results used to interpret lineage-associated functional context.",
    )

    replace_start(
        doc,
        "The proportional allocation was calculated",
        "The proportional allocation was calculated from the complete typed source set before analytical-label harmonization using n_g = round(N_g / 10,951 x 1,000), followed by final adjustment of rounding residuals to produce exactly 1,000 sampled genomes. The allocation plan was: A = 3,228 -> 296 planned samples (29.48%), B2 = 2,661 -> 243 (24.30%), B1 = 2,183 -> 199 (19.93%), D = 1,123 -> 103 (10.25%), E = 791 -> 72 (7.22%), C = 396 -> 36 (3.62%), F = 335 -> 31 (3.06%), G = 125 -> 11 (1.14%), U = 48 -> 4 (0.44%), CRYPTIC = 29 -> 3 (0.26%), U/CRYPTIC = 16 -> 1 (0.15%), and EC_CONTROL_FAIL = 16 -> 1 (0.15%). The sampled list was written to `Sampled_Phylogroup_Results.csv`, and the proportional allocation table is provided as `phylogroup_proportional_sampling_plan.csv`. For the final analytical cohort, non-standard or ambiguous labels were harmonized into the nine reported analytical categories, producing A (n = 298), B2 (n = 244), B1 (n = 200), D (n = 103), E (n = 73), C (n = 36), F (n = 31), G (n = 11), and U (n = 4).",
    )

    replace_start(
        doc,
        "CRISPR-Cas immune system identification.",
        "Auxiliary CRISPR-Cas immune-system annotation. CRISPRcas-Typer v1.1.4 [22] was applied to identify CRISPR arrays and associated Cas operon architectures. These outputs were treated as supporting functional context for lineage-associated genome dynamics rather than as a primary mechanistic endpoint. The analysis uses both array-level outputs (`all_genomes_crispr_detailed.csv`) and spacer-level outputs (`all_genomes_crispr_spacers_detailed.csv`). For each genome, we quantified CRISPR array number, trusted arrays, repeat counts, estimated spacer counts, observed spacer sequences, unique spacer sequences, spacer length, and the fraction of spacers shared by more than one genome. Unique spacer sequences were queried against phage and plasmid target databases using short-read nucleotide matching. The primary high-confidence target filter required mismatch <= 2, query coverage >= 80%, pident >= 90%, and evalue <= 1e-3; mismatch <= 3 was retained as a sensitivity threshold. The exported best-hit tables collapse multiple hits per spacer to the best available target under each mismatch threshold. Phage and plasmid target categories were treated as database-match categories rather than mutually exclusive biological classes because homologous mobile-element sequences can occur in both databases. No randomized shuffled-spacer negative-control analysis is included in the current outputs; therefore, spacer-target results are interpreted conservatively as database-supported target matches.",
    )

    replace_start(doc, "2.6 Benchmarking: Runtime, Structure, and Model Consistency", "2.6 Benchmarking: Runtime, Model Structure, and Consistency")

    replace_start(
        doc,
        "For the SD01 770-condition",
        "For the SD01 770-condition BiGG-model benchmark, five curated models were evaluated: iJO1366, iML1515, ic_1306, iUTI89_1310, and the 1327-gene proxy iECSF_1327 for the requested iECs1327 comparison. Public BiGG SBML files were used for these curated models [29, 30]. Paired gapseq reconstructions were evaluated for the available matching strain systems. All model inclusion, mapped-exchange coverage, simulation status, and confusion-matrix outputs are reported in the supplementary benchmark tables. The SD01 benchmark is therefore described throughout as model-consistency testing rather than independent phenotype validation.",
    )

    replace_start(
        doc,
        "For each of the 1,000 phylogroup-stratified genomes",
        "For each of the 1,000 phylogroup-stratified genomes, both the optimized (M9-constrained) model and its pre-filling draft model were assessed for three available topological metrics: number of genes (Genes), reactions (Reactions), and metabolites (Metabolites). Biomass flux (growth rate under LB-equivalent rich medium) was computed via FBA to classify models as Viable (biomass flux > 0) or Lethal (zero growth). These outputs provide cohort-level topology and solvability checks, not a full MEMOTE-style model-quality audit. Advanced QC metrics requested by GEM reviewers, including blocked-reaction fraction, dead-end metabolites, mass-balance status, exchange-reaction count, biomass-reaction identity, and gapfilled-reaction count, are tabulated as not computed from the current cohort outputs in `gem_model_qc_available_metrics.csv` to avoid treating uncomputed metrics as validated results.",
    )

    replace_start(doc, "3.1 Stratified Benchmarking Shows", "3.1 GEM Workflow Benchmarking Shows")
    replace_start(doc, "3.2 Phylogroup-Stratified Genomic Architecture Reveals Divergent Immune Defense and Secretory Strategies", "3.2 Auxiliary Annotation Layers Contextualize Phylogroup-Resolved GEM Variation")
    replace_start(doc, "3.2.2 CRISPR-Cas immune defense divergence", "3.2.2 CRISPR-Cas spacer context")
    replace_start(doc, "3.2.3 Secretome and transmembrane proteome divergence", "3.2.3 SignalP and TMHMM surface-feature context")
    replace_start(doc, "3.2.4 COG functional category divergence", "3.2.4 COG functional-category context")
    replace_start(doc, "3.4 Cohort-Level GEM Topology Shows", "3.4 Cohort-Level GEM Topology Provides the Primary Functional Comparison")

    replace_start(
        doc,
        "The KEGG matrix contained",
        "The KEGG matrix contained 264 non-redundant ko pathways after removal of duplicated map columns. Hellinger-transformed pathway abundance PCA explained 19.1% and 9.0% of variance on PC1 and PC2, respectively, while binary pathway PCA explained 15.0% and 8.8%. PERMANOVA and PERMDISP together indicated phylogroup-associated pathway structure jointly driven by centroid separation and heterogeneous within-group dispersion: binary Jaccard PERMANOVA gave F = 288.60, R2 = 0.6997, p = 0.001, and Hellinger abundance PERMANOVA gave F = 142.58, R2 = 0.5351, p = 0.001, while PERMDISP was also significant for binary (F = 26.57, p = 0.001) and Hellinger analyses (F = 48.65, p = 0.001). These R2 values are therefore not interpreted as the fraction of pathway variation explained by discrete phylogroup identity alone.",
    )

    replace_start(
        doc,
        "Figure 3. KEGG pathway landscape",
        "Figure 3. KEGG pathway landscape in the analyzed cohort. Hellinger-transformed pathway PCA and complementary PERMANOVA/PERMDISP testing indicate phylogroup-associated pathway structure jointly driven by centroid separation and dispersion heterogeneity; PERMANOVA R2 values are not interpreted as pure phylogroup-explained variance. The heatmap shows the 15 pathways significant after BH and Bonferroni correction.",
    )

    replace_start(
        doc,
        "Figure 2. Phylogroup-stratified genome",
        "Figure 2. Auxiliary phylogroup-stratified annotation layers supporting GEM interpretation. Boxplots show group distributions across genome size, GEM reactions/metabolites, SignalP-positive proteins per 1,000 CDS, TMHMM-positive proteins per 1,000 CDS (n = 999; one unresolved B1 genome treated as missing), and unique CRISPR spacer sequences. Supplementary Figure S1 summarizes phage- and plasmid-targeted CRISPR spacer matches.",
    )

    replace_start(
        doc,
        "The computational benchmark uses",
        "The computational benchmark anchors the manuscript's main contribution: a distribution-based assessment of optimized gapseq performance across 20 N50-stratified low-contiguity assemblies annotated with phylogroup and CDS metadata (A (n = 4), B1 (n = 7), B2 (n = 4), C (n = 1), E (n = 2), cryptic/U (n = 2); CDS_Count range 4,252-5,584). Optimized gapseq reduced median runtime from 19,821.5 s to 1,260.5 s, corresponding to a median 15.81-fold speedup and a narrow IQR of 15.40-16.35. Five reference model systems were organized for structure and model-consistency benchmarking; MG1655/iML1515 used the same benchmark-average runtime basis, giving a five-model reference interpretation centered on the 12.76-fold speedup estimate while keeping measured 20-genome runtime distributions separate. This benchmark supports cohort-scale GEM construction, but it should not be treated as direct evidence of universal substrate-level predictive performance.",
    )

    replace_start(
        doc,
        "The integrated analysis supports",
        "The integrated analysis is therefore best read as a GEM-centered workflow and application study. The optimized reconstruction pipeline makes 1,000-genome GEM analysis practical; the phylogroup-stratified cohort then provides a controlled use case for asking how GEM topology and KEGG pathway composition vary with both categorical phylogroup labels and continuous core-genome relatedness. COG, CRISPR spacer, SignalP, and TMHMM results are retained as contextual annotation layers that help interpret the GEM/pathway findings, but they are not presented as separate proof of lineage-specific virulence mechanisms.",
    )

    replace_start(
        doc,
        "Several caveats should guide interpretation.",
        "Several caveats should guide interpretation. First, pathway names are represented primarily by KEGG identifiers; pathway-name annotation, representative KO/gene members, Cramer's V, and effect-size tables are needed for detailed biological interpretation of the 15 adjusted signals. Second, phylogroup assignment was performed with EzClermont/ClermonTyping-compatible classification before cohort selection; the exact run script and proportional allocation are reported, but the sampling log does not currently expose a fixed random seed or full within-group ordering rule beyond phylogroup-proportional allocation, so the final repository should include the sampling script and `Sampled_Phylogroup_Results.csv`. Third, the residual phylogroup analysis uses ANI PCoA covariates and label-permutation sensitivity tests rather than a full core-gene tree with formal phylogenetic comparative modeling; it therefore distinguishes continuous genomic relatedness from categorical labels but does not replace tree-based evolutionary inference. Fourth, the current cohort-level GEM QC includes topology and solvability checks but does not compute blocked-reaction fraction, dead-end metabolites, mass-balance status, exchange-reaction count, biomass-reaction identity, or gapfilled-reaction count across all 1,000 models. Fifth, CRISPR target matching depends on the completeness and taxonomic composition of the phage and plasmid databases used for spacer searches, and no shuffled-spacer negative-control table is included in the current outputs; therefore target absence should be interpreted as no detectable match in the queried databases rather than true absence of prior exposure. Sixth, the SD01 carbon-source benchmark uses the Monk et al. Table 5 in silico growth-screen matrix, not independent experimental Biolog data; 32/770 conditions remained unmapped automatically, iECs1327 required the iECSF_1327/SE15 proxy, and several gapseq reconstructions showed oxygen-condition-specific sensitivity losses. Reaction compaction should not be interpreted by itself as proof of functional equivalence; phenotype-level and substrate-level validation remain necessary in representative strains.",
    )

    replace_start(
        doc,
        "In conclusion, this study presents",
        "In conclusion, this study presents an optimized and systematically benchmarked E. coli GEM reconstruction framework applied to a 1,000-genome phylogroup-matched cohort. Stratified benchmarking across 20 assemblies showed a median 15.81-fold runtime acceleration, paired XML reaction matrices identified 111.0 optimized-only reactions per genome at the median, and five reference model systems supported reaction-set and model-consistency benchmarking, with MG1655/iML1515 using the same benchmark-average runtime basis. Automated SD01 carbon-source mapping translated 738/770 aerobic/anaerobic environments into model-readable exchanges and produced reproducible FBA prediction matrices for curated BiGG and gapseq models, providing a transparent in silico consistency benchmark rather than independent experimental phenotype validation. FastANI residual analysis showed that KEGG pathway composition is strongly tied to continuous core-genome distance, whereas discrete phylogroup labels retain smaller residual effects for KEGG, COG, GEM topology, CRISPR spacer composition, CRISPR target burden, and SignalP/TMHMM surface architecture. Overall, the manuscript's primary contribution is a scalable GEM workflow and its conservative application to phylogroup-resolved functional divergence; the main validation priority remains true multi-strain experimental Biolog or carbon-source data for substrate-level predictive benchmarking.",
    )

    replace_start(
        doc,
        "Availability of data and materials:",
        "Availability of data and materials: All genome accession numbers, phylogroup assignments, annotation outputs, GEM topology metrics, available GEM QC fields, benchmark metadata, reaction-set matrices, SD01 carbon-source exchange mappings, and FBA prediction/model-consistency matrices are organized in `05.supplementary_data` and will be deposited with the manuscript as Supplementary Data Files. Genome sequences were obtained from NCBI Assembly (https://www.ncbi.nlm.nih.gov/assembly/). The EzClermont run script (`qc0_run_ezclermont.py`), proportional sampling plan (`phylogroup_proportional_sampling_plan.csv`), available GEM QC table (`gem_model_qc_available_metrics.csv`), modified gapseq scripts, sampling outputs, and pipeline configuration files are available in the project script directory and will be deposited at [repository URL to be supplied before submission].",
    )

    # Table 1 updates.
    table = doc.tables[0]
    table.cell(1, 3).text = "quality_report.tsv"
    table.cell(2, 3).text = "selected_1000_genomes.tsv; phylogroup_proportional_sampling_plan.csv; recalculated_annotations.csv; gem_model_qc_available_metrics.csv"
    table.cell(6, 0).text = "SD01 curated-model consistency FBA"

    doc.save(DOCX)
    print(f"saved {DOCX}")
    print(f"backup {BACKUP}")


if __name__ == "__main__":
    main()
