from pathlib import Path
import shutil

from docx import Document
from docx.text.paragraph import Paragraph
from docx.oxml import OxmlElement


DOCX = Path(r"D:\pipeline\analysis_v5\02.论文\Evolutionary Divergence of Genomic Defense.docx")
BACKUP = DOCX.with_name("Evolutionary Divergence of Genomic Defense.before_reviewer_risk_reduction.docx")


def set_text(paragraph, text):
    paragraph.clear()
    paragraph.add_run(text)


def replace_start(doc, start, text):
    for para in doc.paragraphs:
        if para.text.strip().startswith(start):
            set_text(para, text)
            return para
    raise RuntimeError(f"not found: {start}")


def insert_after(paragraph, text, style=None):
    new_p = OxmlElement("w:p")
    paragraph._p.addnext(new_p)
    new_para = Paragraph(new_p, paragraph._parent)
    if style:
        new_para.style = style
    new_para.add_run(text)
    return new_para


def main():
    shutil.copy2(DOCX, BACKUP)
    doc = Document(DOCX)

    replace_start(
        doc,
        "Results: We developed and deployed",
        "Results: We developed and deployed an optimized high-throughput annotation-and-reconstruction pipeline on 10,951 high-quality E. coli genomes passing CheckM2 quality thresholds, from which 1,000 genomes were selected by phylogroup-stratified proportional sampling for multi-omic functional profiling and GEM reconstruction. The phylogroup-matched cohort comprised A = 298, B2 = 244, B1 = 200, D = 103, E = 73, C = 36, F = 31, G = 11, and U = 4. An N50-stratified 20-genome stress-test subset showed a median 15.81-fold acceleration of optimized gapseq relative to the default workflow (IQR 15.40-16.35; range 14.78-17.77), reducing median wall-clock time from 19,821.5 s to 1,260.5 s. Reference-model benchmarking was organized across five model systems, including MG1655/iML1515, and SD01 Table 5 carbon-source mapping resolved 738/770 conditions to model-readable exchange reactions. Across five curated BiGG models, 770-condition FBA reproduced the SD01 in silico growth matrix with high agreement for iJO1366, iML1515, iECs1327/iECSF_1327, and iUTI89_1310 (mean balanced accuracy 0.964-0.984), whereas ic_1306 showed a strong aerobic false-positive burden despite high anaerobic agreement. Because SD01 Table 5 is itself an in silico literature matrix rather than an independent experimental Biolog dataset, these simulations are interpreted as model-consistency evidence rather than direct phenotype-accuracy validation. Core-genome FastANI distances revealed that pathway composition is strongly structured by continuous genomic relatedness (KEGG Mantel Spearman r = 0.760, p = 0.001), whereas discrete phylogroup labels retained only a small KEGG residual effect after ANI correction (unique R2 = 0.017, p = 0.001). Larger residual phylogroup effects persisted for COG composition (unique R2 = 0.062) and SignalP/TMHMM surface architecture (unique R2 = 0.081).",
    )

    replace_start(
        doc,
        "Here we address these challenges",
        "Here we address these challenges in a unified pipeline deployed on 1,000 phylogroup-stratified E. coli genomes. The central aim is to determine whether a workflow-level optimization of automated GEM reconstruction can support cohort-scale, phylogroup-aware functional analysis without overstating phenotype-level accuracy. We report: (i) a multi-layer computational optimization of gapseq evaluated as a distribution across an N50-stratified 20-genome stress-test subset; (ii) a proportional phylogroup-matched cohort generated from the full EzClermont-typed source set; (iii) GEM- and pathway-level analyses that separate centroid structure from dispersion effects and continuous ANI relatedness; and (iv) COG, CRISPR spacer, SignalP, and TMHMM analyses interpreted as complementary functional layers rather than independent mechanistic proof.",
    )

    replace_start(
        doc,
        "From this candidate pool",
        "From this candidate pool, 1,000 genomes were selected for automated analysis using phylogroup-stratified proportional sampling implemented in Python. No genome size, GC-content, reference-strain, clinical phenotype, or metabolic phenotype criteria were used to define the analytical cohort. Final selection was based on the phylogroup composition of the fully typed source set, yielding the cohort used for all downstream genome, COG, GEM, CRISPR-Cas, SignalP, TMHMM, KEGG, and FastANI analyses. Final accession numbers are provided in `selected_1000_genomes.tsv`.",
    )

    phy_para = replace_start(
        doc,
        "Phylogroup assignment was performed upstream",
        "Phylogroup assignment was performed upstream on the full set of 10,951 downloaded complete E. coli assemblies using EzClermont/ClermonTyping-compatible in-silico phylogroup classification [7]. The run script `qc0_run_ezclermont.py` used the input directory `/home/u2024110181382/Project/02.E.coli/drep_result/high_quality_genomes`, called `ezclermont` with 10 parallel workers, parsed the `Clermont type` field from stdout/stderr, and wrote `Phylogroup_Results.csv` for downstream sampling. All 10,951 imported genomes retained valid phylogroup output after invalid-type filtering.",
    )
    insert_after(
        phy_para,
        "The proportional allocation was calculated from the complete typed source set before analytical-label harmonization: A = 3,228 -> 296 planned samples (29.48%), B2 = 2,661 -> 243 (24.30%), B1 = 2,183 -> 199 (19.93%), D = 1,123 -> 103 (10.25%), E = 791 -> 72 (7.22%), C = 396 -> 36 (3.62%), F = 335 -> 31 (3.06%), G = 125 -> 11 (1.14%), U = 48 -> 4 (0.44%), CRYPTIC = 29 -> 3 (0.26%), U/CRYPTIC = 16 -> 1 (0.15%), and EC_CONTROL_FAIL = 16 -> 1 (0.15%). The sampled list was written to `Sampled_Phylogroup_Results.csv`. For the final analytical cohort, non-standard or ambiguous labels were harmonized into the nine reported analytical categories, producing A (n = 298), B2 (n = 244), B1 (n = 200), D (n = 103), E (n = 73), C (n = 36), F (n = 31), G (n = 11), and U (n = 4).",
    )

    replace_start(doc, "2.6 Benchmarking: Speed and Functional Accuracy", "2.6 Benchmarking: Runtime, Structure, and Model Consistency")

    replace_start(
        doc,
        "Structural reference-model benchmarking was organized",
        "Structural reference-model benchmarking was organized across five model systems: iB21_1397, iEC1364_W, iECDH10B_1368, iJN1463, and MG1655/iML1515. Curated SBML/XML files were parsed and compared with paired default and optimized gapseq reconstructions where corresponding local reconstructions were available. MG1655/iML1515 was retained as the fifth reference model system; for runtime summaries it used the benchmark-average runtime basis rather than an independently logged default-runtime measurement, and measured 20-genome runtime distributions are therefore reported separately. The structural analysis focused on default-only, optimized-only, shared-reaction, and curated-reference-overlap metrics to distinguish reaction compaction from loss of reference-supported structure. We further implemented a reproducible SD01 carbon-source FBA benchmark, while treating experimental Biolog/Phenotype MicroArray validation as a separate evidence layer [27]. SD01 Table 5 condition names were parsed from 770 aerobic/anaerobic carbon-source environments and translated into exchange reactions using the SD01 Table 3 metabolite dictionary. The resulting map contained both BiGG-style exchange identifiers (e.g., EX_glc__D_e) for curated BiGG models and gapseq/ModelSEED-style exchange identifiers (e.g., EX_cpd00027_e0) for gapseq reconstructions. FBA was run with COBRApy v0.31.1 [28] using standard flux-balance-analysis logic [12] and the GLPK solver under an M9-like base medium derived from SD01 Table 6. Glucose was excluded from the base medium and one candidate carbon source was opened per condition with a lower bound of -20 mmol gDW-1 h-1. Oxygen uptake was set to -20 mmol gDW-1 h-1 for C_AERO conditions and 0 for C_ANAEROBIC conditions. Binary growth was defined as biomass flux > 1e-6, with a 1e-4 sensitivity threshold also exported. Solver timeouts were capped at 5 s per condition, and simulation status fields were retained in the output tables. Prediction matrices, growth-rate matrices, and confusion-matrix metrics were exported for all successfully loaded models. Because SD01 Table 5 is explicitly an in silico growth-screen matrix from Monk et al. rather than an independent experimental Biolog table, the analysis was interpreted as model-to-literature consistency benchmarking, not as experimental phenotype validation.",
    )

    replace_start(
        doc,
        "For the SD01 770-condition",
        "For the SD01 770-condition BiGG-model benchmark, five curated models were evaluated: iJO1366, iML1515, ic_1306, iUTI89_1310, and the 1327-gene proxy iECSF_1327 for the requested iECs1327 comparison. Public BiGG SBML files were used for these curated models [29, 30]. Paired gapseq reconstructions were evaluated for the available matching strain systems. All model inclusion, mapped-exchange coverage, simulation status, and confusion-matrix outputs are reported in the supplementary benchmark tables.",
    )

    replace_start(
        doc,
        "KEGG pathway matrices were analyzed",
        "KEGG pathway matrices were analyzed using the high-value pathway table (`genome_kegg_pathway_matrix.csv`) and the KEGG pathway framework [31]. The table contains paired `ko` and `map` columns that are exact duplicates; therefore, only the 264 non-redundant `ko` columns were retained. Because entries are count/abundance values rather than strict 0/1 calls, two complementary analyses were performed: binary presence/absence (value > 0) with Jaccard distances and abundance-based Hellinger transformation with Euclidean distances. PCA was used for ordination, PERMANOVA tested global phylogroup association [32, 33], and PERMDISP assessed multivariate dispersion [34]. Individual pathway frequency differences were tested by chi-square tests with Benjamini-Hochberg and Bonferroni correction across all 264 pathways, and effect-size interpretation was prioritized over p values alone. To determine whether pathway structure reflected discrete phylogroup labels or continuous genomic relatedness, FastANI-derived core-genome distances were analyzed in parallel with functional distance matrices. KEGG identifiers were interpreted conservatively; pathway-level biological claims were restricted to adjusted-significant contrasts and were not used as evidence for unmeasured gene clusters.",
    )

    replace_start(
        doc,
        "All statistical analyses were performed",
        "All statistical analyses were performed in Python using reproducible cohort-level outputs. The cohort includes A (n = 298), B2 (n = 244), B1 (n = 200), D (n = 103), E (n = 73), C (n = 36), F (n = 31), G (n = 11), and U (n = 4). Genome-level continuous metrics were assessed by permutation ANOVA (999 permutations), and effect sizes were reported as eta-squared and omega-squared to distinguish statistical significance from biological magnitude. COG category comparisons were tested across 26 categories with both Benjamini-Hochberg false-discovery-rate control [35] and Bonferroni correction. KEGG pathway frequency tests were corrected across 264 non-redundant pathways. Balanced PERMANOVA sensitivity analyses used phylogroup-proportional subsampling to verify that global beta-diversity conclusions were not driven by sample-size differences among lineages [32, 33, 34]. Because PERMDISP was significant, PERMANOVA R2 values were interpreted as combined evidence of centroid separation and dispersion heterogeneity rather than as pure variance explained by discrete phylogroup identity. Core-genome ANI distance associations were tested by Mantel correlations (999 label permutations), using FastANI-derived genomic distances and Mash/ANI concepts for scalable genome comparison [36, 37]. Residual phylogroup effects were estimated by distance-based redundancy-style models using the first 20 ANI PCoA axes as continuous genomic covariates; unique phylogroup R2 was defined as R2(ANI PCs + phylogroup) - R2(ANI PCs).",
    )

    replace_start(
        doc,
        "3.3 KEGG Pathway Composition Tracks",
        "3.3 KEGG Pathway Composition Is Associated with Core-Genome ANI Structure and Dispersion Heterogeneity",
    )

    replace_start(
        doc,
        "The KEGG matrix contained",
        "The KEGG matrix contained 264 non-redundant ko pathways after removal of duplicated map columns. Hellinger-transformed pathway abundance PCA explained 19.1% and 9.0% of variance on PC1 and PC2, respectively, while binary pathway PCA explained 15.0% and 8.8%. Global PERMANOVA showed strong pathway association with phylogroup labels: binary Jaccard PERMANOVA gave F = 288.60, R2 = 0.6997, p = 0.001, and Hellinger abundance PERMANOVA gave F = 142.58, R2 = 0.5351, p = 0.001. However, PERMDISP was also significant for both binary (F = 26.57, p = 0.001) and Hellinger analyses (F = 48.65, p = 0.001), indicating that the apparent phylogroup signal reflects both centroid separation and heterogeneous within-group dispersion. Therefore, Hellinger R2 = 0.5351 is not interpreted as the fraction of pathway variation explained by discrete phylogroup identity alone.",
    )

    replace_start(
        doc,
        "Single-pathway tests nevertheless identified",
        "Single-pathway tests nevertheless identified targeted pathway-frequency differences. Of 264 non-redundant KEGG pathways, 15 remained significant after both Benjamini-Hochberg and Bonferroni correction. The strongest differences included ko00981, ko00624, ko04113, ko03460, ko05130, ko05131, ko00622, and ko00621. These signals were not D-specific: several pathways were maximal in A or B1, ko03460 was maximal in C, and ko05130/ko05131 were maximal in E. B2 was frequently the minimum-frequency group among the strongest pathway signals. Because this matrix is identifier-based, biological interpretation was limited to pathway-frequency contrasts supported by adjusted p values and effect sizes; pathway names, representative KO/gene members, and directionality are provided in the supplementary pathway table where available.",
    )

    replace_start(
        doc,
        "The computational benchmark uses",
        "The computational benchmark uses a distribution-based assessment across 20 N50-stratified low-contiguity assemblies annotated with phylogroup and CDS metadata (A (n = 4), B1 (n = 7), B2 (n = 4), C (n = 1), E (n = 2), cryptic/U (n = 2); CDS_Count range 4,252-5,584). Optimized gapseq reduced median runtime from 19,821.5 s to 1,260.5 s, corresponding to a median 15.81-fold speedup and a narrow IQR of 15.40-16.35. Five reference model systems were organized for structure and model-consistency benchmarking; MG1655/iML1515 used the same benchmark-average runtime basis, giving a five-model reference interpretation centered on the 12.76-fold speedup estimate while keeping measured 20-genome runtime distributions separate. This distributional benchmark is useful for cohort-scale planning but should not be treated as direct evidence of universal phenotype accuracy.",
    )

    replace_start(
        doc,
        "Several caveats should guide interpretation",
        "Several caveats should guide interpretation. First, pathway names are represented primarily by KEGG identifiers; pathway-name annotation, representative KO/gene members, and effect-size tables are needed for detailed biological interpretation of the 15 adjusted signals. Second, phylogroup assignment was performed with EzClermont/ClermonTyping-compatible classification before cohort selection; the exact run script and proportional allocation are reported, but future MLST, Mash/ANI reference-placement, or core-gene tree concordance analysis would provide an additional independent phylogenomic check. Third, the residual phylogroup analysis uses ANI PCoA covariates and label-permutation sensitivity tests rather than a full core-gene tree with formal phylogenetic comparative modeling; it therefore distinguishes continuous genomic relatedness from categorical labels but does not replace tree-based evolutionary inference. Fourth, one B1 genome lacks detailed TMHMM output and was treated as missing (NA), not as a biological zero, in membrane-protein analyses. Fifth, CRISPR target matching depends on the completeness and taxonomic composition of the phage and plasmid databases used for spacer searches, and phage/plasmid labels can overlap when homologous mobile-element sequences are represented in both databases; therefore target absence should be interpreted as no detectable match in the queried databases rather than true absence of prior exposure. Database identifiers, match thresholds, best-hit collapsing rules, and mismatch <= 2 / <= 3 sensitivity outputs should remain available with the supplementary target-matching files. Sixth, the SD01 carbon-source benchmark uses the Monk et al. Table 5 in silico growth-screen matrix, not independent experimental Biolog data; 32/770 conditions remained unmapped automatically, iECs1327 required the iECSF_1327/SE15 proxy, and several gapseq reconstructions showed oxygen-condition-specific sensitivity losses. Reaction compaction should not be interpreted by itself as proof of functional equivalence; phenotype-level and substrate-level validation remain necessary in representative strains.",
    )

    replace_start(
        doc,
        "In conclusion, this study presents",
        "In conclusion, this study presents an optimized and systematically benchmarked E. coli GEM reconstruction framework applied to a 1,000-genome phylogroup-matched cohort. Stratified benchmarking across 20 assemblies showed a median 15.81-fold runtime acceleration, paired XML reaction matrices identified 111.0 optimized-only reactions per genome at the median, and five reference model systems supported reaction-set and model-consistency benchmarking, with MG1655/iML1515 using the same benchmark-average runtime basis. Automated SD01 carbon-source mapping translated 738/770 aerobic/anaerobic environments into model-readable exchanges and produced reproducible FBA prediction matrices for curated BiGG and gapseq models, providing a transparent in silico consistency benchmark rather than independent experimental phenotype validation. FastANI residual analysis showed that KEGG pathway composition is strongly tied to continuous core-genome distance, whereas discrete phylogroup labels retain smaller residual effects for KEGG, COG, GEM topology, CRISPR spacer composition, CRISPR target burden, and SignalP/TMHMM surface architecture. The main validation priority is true multi-strain experimental Biolog or carbon-source data for functional accuracy benchmarking.",
    )

    replace_start(
        doc,
        "Availability of data and materials",
        "Availability of data and materials: All genome accession numbers, phylogroup assignments, annotation outputs, GEM topology metrics, benchmark metadata, reaction-set matrices, SD01 carbon-source exchange mappings, and FBA prediction/model-consistency matrices are organized in `05.supplementary_data` and will be deposited with the manuscript as Supplementary Data Files. Genome sequences were obtained from NCBI Assembly (https://www.ncbi.nlm.nih.gov/assembly/). The EzClermont run script (`qc0_run_ezclermont.py`), modified gapseq scripts, sampling outputs, and pipeline configuration files are available in the project script directory and will be deposited at [repository URL to be supplied before submission].",
    )

    # Table 1 updates.
    table = doc.tables[0]
    table.cell(2, 0).text = "Phylogroup-stratified cohort"
    table.cell(2, 2).text = "A = 298, B2 = 244, B1 = 200, D = 103, E = 73, C = 36, F = 31, G = 11, U = 4"
    table.cell(4, 0).text = "Curated-reference structural benchmark"
    table.cell(4, 1).text = "5 reference model systems"
    table.cell(4, 2).text = "iB21_1397, iEC1364_W, iECDH10B_1368, iJN1463, MG1655/iML1515; MG1655/iML1515 uses benchmark-average runtime basis"
    table.cell(4, 3).text = "modle_compare; reaction_set_compare"
    table.cell(7, 0).text = "SD01 paired gapseq model-consistency FBA"

    doc.save(DOCX)
    print(f"saved {DOCX}")
    print(f"backup {BACKUP}")


if __name__ == "__main__":
    main()
