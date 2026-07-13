from pathlib import Path
import shutil

from docx import Document


DOCX = Path(r"D:\pipeline\analysis_v5\02.论文\Evolutionary Divergence of Genomic Defense.docx")
BACKUP = DOCX.with_name("Evolutionary Divergence of Genomic Defense.before_first_version_cleanup.docx")
OUT = DOCX.with_name("Evolutionary Divergence of Genomic Defense.first_version_clean.docx")


def set_text(paragraph, text):
    style = paragraph.style
    alignment = paragraph.alignment
    paragraph.clear()
    paragraph.style = style
    paragraph.alignment = alignment
    paragraph.add_run(text)


def replace_by_start(doc, start, new_text):
    for para in doc.paragraphs:
        if para.text.startswith(start):
            set_text(para, new_text)
            return True
    raise RuntimeError(f"paragraph not found: {start}")


def replace_exact(doc, old_text, new_text):
    for para in doc.paragraphs:
        if para.text == old_text:
            set_text(para, new_text)
            return True
    raise RuntimeError(f"paragraph not found: {old_text[:80]}")


def replace_contains(doc, needle, new_text):
    for para in doc.paragraphs:
        if needle in para.text:
            set_text(para, new_text)
            return True
    raise RuntimeError(f"paragraph not found containing: {needle}")


def main():
    shutil.copy2(DOCX, BACKUP)
    doc = Document(DOCX)

    replace_by_start(
        doc,
        "Results: We developed and deployed",
        "Results: We developed and deployed an optimized high-throughput annotation-and-reconstruction pipeline on 10,951 high-quality E. coli genomes passing CheckM2 quality thresholds, from which 1,000 representative genomes were selected for multi-omic functional profiling and GEM reconstruction. The phylogroup-matched cohort comprised A = 298, B2 = 244, B1 = 200, D = 103, E = 73, C = 36, F = 31, G = 11, and U = 4. An N50-stratified 20-genome stress-test subset showed a median 15.81-fold acceleration of optimized gapseq relative to the default workflow (IQR 15.40-16.35; range 14.78-17.77), reducing median wall-clock time from 19,821.5 s to 1,260.5 s. Reference-model benchmarking was organized across five model systems, including MG1655/iML1515, and SD01 Table 5 carbon-source mapping resolved 738/770 conditions to model-readable exchange reactions (364 unique BiGG exchanges and 314 unique gapseq/ModelSEED exchanges). Across five curated BiGG models, 770-condition FBA reproduced the SD01 in silico growth matrix with high agreement for iJO1366, iML1515, iECs1327/iECSF_1327, and iUTI89_1310 (mean balanced accuracy 0.964-0.984), whereas ic_1306 showed a strong aerobic false-positive burden despite high anaerobic agreement. Exploratory paired gapseq benchmarking showed modest mean improvement in balanced accuracy for optimized models (0.592 versus 0.562), with strain- and oxygen-condition-specific tradeoffs; these SD01 data are treated as an in silico consistency benchmark rather than independent experimental Biolog validation. Core-genome FastANI distances revealed that pathway composition is strongly structured by continuous genomic relatedness (KEGG Mantel Spearman r = 0.760, p = 0.001), whereas discrete phylogroup labels retained only a small KEGG residual effect after ANI correction (unique R2 = 0.017, p = 0.001). Larger residual phylogroup effects persisted for COG composition (unique R2 = 0.062) and SignalP/TMHMM surface architecture (unique R2 = 0.081).",
    )

    replace_by_start(
        doc,
        "Conclusions: This work establishes",
        "Conclusions: This work establishes a phylogroup-resolved metabolic and functional atlas of E. coli based on a proportional phylogroup-matched cohort. The data support lineage-associated differences in COG profiles, secretion signals, CRISPR spacer repertoires, membrane-protein architecture, and KEGG pathway composition; however, FastANI residual analysis indicates that much of the pathway signal tracks continuous core-genome relatedness rather than discrete phylogroup categories alone. These results provide a conservative and reproducible framework for interpreting E. coli functional divergence in large genome cohorts.",
    )

    replace_by_start(
        doc,
        "Here we address these challenges",
        "Here we address these challenges in a unified pipeline deployed on 1,000 representative E. coli genomes. We report: (i) a multi-layer computational optimization of gapseq evaluated as a distribution across an N50-stratified 20-genome stress-test subset; (ii) a proportional phylogroup-matched cohort with strong representation of the major E. coli phylogroups; (iii) GEM- and pathway-level analyses that separate centroid structure from dispersion effects; and (iv) COG, CRISPR spacer, SignalP, and TMHMM analyses designed to reduce over-interpretation under multiple testing and uneven lineage abundance.",
    )

    replace_by_start(
        doc,
        "Phylogroup assignment was performed upstream",
        "Phylogroup assignment was performed upstream on the full set of 10,951 downloaded complete E. coli assemblies using EzClermont/ClermonTyping-compatible in-silico phylogroup classification [7]. The 1,000-genome analytical cohort was selected from this fully typed source set by proportional stratified sampling according to phylogroup abundance, producing a phylogroup-matched distribution: A (n = 298), B2 (n = 244), B1 (n = 200), D (n = 103), E (n = 73), C (n = 36), F (n = 31), G (n = 11), and U (n = 4). This design enables group-level evaluation of major lineage-associated patterns while preserving the observed phylogroup structure of the typed source set.",
    )

    replace_by_start(
        doc,
        "CRISPR-Cas immune system identification.",
        "CRISPR-Cas immune system identification. CRISPRcas-Typer v1.1.4 [22] was applied to identify CRISPR arrays and associated Cas operon architectures. The analysis uses both array-level outputs (`all_genomes_crispr_detailed.csv`) and spacer-level outputs (`all_genomes_crispr_spacers_detailed.csv`). For each genome, we quantified CRISPR array number, trusted arrays, repeat counts, estimated spacer counts, observed spacer sequences, unique spacer sequences, spacer length, and the fraction of spacers shared by more than one genome. Unique spacer sequences were additionally queried against phage and plasmid target databases using short-read nucleotide matching. The primary high-confidence target filter required mismatch <= 2, query coverage >= 80%, pident >= 90%, and evalue <= 1e-3; mismatch <= 3 was retained as a sensitivity threshold. Phage and plasmid target categories were treated as database-match categories rather than mutually exclusive biological classes because homologous mobile-element sequences can occur in both databases.",
    )

    replace_by_start(
        doc,
        "Transmembrane topology prediction.",
        "Transmembrane topology prediction. TMHMM v2.0 [25] was applied to all predicted protein sequences to enumerate proteins containing at least one transmembrane helix. The TMHMM-positive protein table contains 1,125,213 protein-level records across 999 of 1,000 genomes. Per-genome analyses used TMHMM-positive protein count, total predicted helices, mean helices per membrane protein, maximum helices, and TMHMM proteins per 1,000 CDS. The single genome without TMHMM detail output was retained as missing (NA), not as a biological zero.",
    )

    replace_by_start(
        doc,
        "Computational benchmarking was expanded",
        "Computational benchmarking used an N50-stratified stress-test subset of 20 E. coli assemblies selected to represent a broad range of assembly contiguity (Actual_N50 = 22,535-675,909 bp). Metadata annotation showed that this subset covered A (n = 4), B1 (n = 7), B2 (n = 4), C (n = 1), E (n = 2), and cryptic/U (n = 2), with CDS counts ranging from 4,252 to 5,584 (median 4,842.0). Each genome was reconstructed with the gapseq default workflow and the optimized workflow, and wall-clock time, reaction count, speedup (default runtime / optimized runtime), reaction delta, and reaction retention were summarized as distributions. To test whether model compaction reflected simple truncation or a systematic change in reaction content, final SBML/XML `listOfReactions` entries were parsed for each paired reconstruction, and optimized-only, default-only, shared-reaction, and exact Jaccard matrices were generated.",
    )

    replace_by_start(
        doc,
        "Structural reference-model benchmarking",
        "Structural reference-model benchmarking was organized across five model systems: iB21_1397, iEC1364_W, iECDH10B_1368, iJN1463, and MG1655/iML1515. Curated SBML/XML files were parsed and compared with paired default and optimized gapseq reconstructions where corresponding local reconstructions were available. This analysis focused on default-only, optimized-only, shared-reaction, and curated-reference-overlap metrics to distinguish reaction compaction from loss of reference-supported structure. We further implemented a reproducible SD01 carbon-source FBA benchmark, while treating experimental Biolog/Phenotype MicroArray validation as a separate evidence layer [27]. SD01 Table 5 condition names were parsed from 770 aerobic/anaerobic carbon-source environments and translated into exchange reactions using the SD01 Table 3 metabolite dictionary. The resulting map contained both BiGG-style exchange identifiers (e.g., EX_glc__D_e) for curated BiGG models and gapseq/ModelSEED-style exchange identifiers (e.g., EX_cpd00027_e0) for gapseq reconstructions. FBA was run with COBRApy v0.31.1 [28] using standard flux-balance-analysis logic [12] and the GLPK solver under an M9-like base medium derived from SD01 Table 6. Glucose was excluded from the base medium and one candidate carbon source was opened per condition with a lower bound of -20 mmol gDW-1 h-1. Oxygen uptake was set to -20 mmol gDW-1 h-1 for C_AERO conditions and 0 for C_ANAEROBIC conditions. Binary growth was defined as biomass flux > 1e-6, with a 1e-4 sensitivity threshold also exported. Solver timeouts were capped at 5 s per condition, and simulation/load status fields were retained in the output tables. Prediction matrices, growth-rate matrices, and confusion-matrix metrics were exported for all successfully loaded models. Because SD01 Table 5 is explicitly an in silico growth-screen matrix from Monk et al. rather than an independent experimental Biolog table, the analysis was interpreted as model-to-literature consistency benchmarking, not as experimental phenotype validation.",
    )

    replace_by_start(
        doc,
        "For the SD01 770-condition",
        "For the SD01 770-condition BiGG-model benchmark, five curated models were evaluated: iJO1366, iML1515, ic_1306, iUTI89_1310, and the 1327-gene proxy iECSF_1327 for the requested iECs1327 comparison. Public BiGG SBML files were used for these curated models [29, 30]. Paired gapseq reconstructions were evaluated for the available matching strain systems. All model-load status, mapped-exchange coverage, simulation status, and confusion-matrix outputs are reported in the supplementary benchmark tables.",
    )

    replace_by_start(
        doc,
        "KEGG pathway matrices were analyzed",
        "KEGG pathway matrices were analyzed using the high-value pathway table (`genome_kegg_pathway_matrix.csv`) and the KEGG pathway framework [31]. The table contains paired `ko` and `map` columns that are exact duplicates; therefore, only the 264 non-redundant `ko` columns were retained. Because entries are count/abundance values rather than strict 0/1 calls, two complementary analyses were performed: binary presence/absence (value > 0) with Jaccard distances and abundance-based Hellinger transformation with Euclidean distances. PCA was used for ordination, PERMANOVA tested global phylogroup association [32, 33], and PERMDISP assessed multivariate dispersion [34]. Individual pathway frequency differences were tested by chi-square tests with Benjamini-Hochberg and Bonferroni correction across all 264 pathways. To determine whether pathway structure reflected discrete phylogroup labels or continuous genomic relatedness, FastANI-derived core-genome distances were analyzed in parallel with functional distance matrices.",
    )

    replace_by_start(
        doc,
        "All statistical analyses were performed",
        "All statistical analyses were performed in Python using reproducible cohort-level outputs. The cohort includes A (n = 298), B2 (n = 244), B1 (n = 200), D (n = 103), E (n = 73), C (n = 36), F (n = 31), G (n = 11), and U (n = 4). Genome-level continuous metrics were assessed by permutation ANOVA (999 permutations), and effect sizes were reported as eta-squared and omega-squared to distinguish statistical significance from biological magnitude. COG category comparisons were tested across 26 categories with both Benjamini-Hochberg false-discovery-rate control [35] and Bonferroni correction. KEGG pathway frequency tests were corrected across 264 non-redundant pathways. Balanced PERMANOVA sensitivity analyses used phylogroup-proportional subsampling to verify that global beta-diversity conclusions were not driven by sample-size differences among lineages [32, 33, 34]. Core-genome ANI distance associations were tested by Mantel correlations (999 label permutations), using FastANI-derived genomic distances and Mash/ANI concepts for scalable genome comparison [36, 37]. Residual phylogroup effects were estimated by distance-based redundancy-style models using the first 20 ANI PCoA axes as continuous genomic covariates; unique phylogroup R2 was defined as R2(ANI PCs + phylogroup) - R2(ANI PCs).",
    )

    replace_by_start(
        doc,
        "Together, the stratified assembly benchmark",
        "Together, the stratified assembly benchmark, paired reaction-set matrices, and reference-model comparison support the optimized workflow as a scalable reconstruction strategy for large cohorts. The evidence rests on a distribution of runtime effects, explicit optimized-only/default-only reaction matrices, accession-level phylogroup/CDS metadata, and multiple reference model systems.",
    )

    replace_exact(doc, " 3.1.3 Phenotype benchmarking and current limits", " 3.1.3 SD01 carbon-source consistency benchmarking")

    replace_by_start(
        doc,
        "Phenotype benchmarking remains",
        "The SD01 carbon-source benchmark separates three evidence levels. First, the MG1655 35-carbon-source M9 panel provides a single-strain proof-of-concept. Second, the SD01 770-condition benchmark tests whether curated BiGG models and paired gapseq reconstructions reproduce the Monk et al. in silico growth-screen matrix after automated exchange-ID mapping [38]. Third, independent experimental Biolog or carbon-source datasets remain the appropriate validation layer for phenotype-level predictive accuracy across strains.",
    )

    replace_by_start(
        doc,
        "For MG1655, the M9-constrained",
        "For MG1655, the M9-constrained optimized model achieved accuracy = 82.86% and sensitivity = 96.55% (TP = 28, TN = 1, FP = 5, FN = 1; 35-carbon-source panel). The SD01 Table 5 carbon-source analysis evaluated 770 aerobic/anaerobic environments. Automated mapping resolved 738/770 conditions to exchange reactions, leaving 32 unmapped conditions that were retained as missing-exchange predictions. Across five curated BiGG models, carbon-source FBA showed high agreement with the SD01 in silico matrix for iECs1327/iECSF_1327 (balanced accuracy 0.989 aerobic and 0.979 anaerobic), iJO1366 (0.979 and 0.963), iML1515 (0.979 and 0.948), and iUTI89_1310 (0.989 and 0.979). The ic_1306 CFT073 model retained high anaerobic agreement (balanced accuracy 0.969) but showed poor aerobic specificity (balanced accuracy 0.557; FP = 183), identifying a model- or medium-specific false-positive burden that requires targeted curation rather than being hidden in an aggregate score.",
    )

    replace_by_start(
        doc,
        "The revised 1,000-genome cohort",
        "The 1,000-genome cohort spans nine assigned phylogroup labels: A (29.8%), B2 (24.4%), B1 (20.0%), D (10.3%), E (7.3%), C (3.6%), F (3.1%), G (1.1%), and U (0.4%). Genome sizes range from 4.49 to 5.86 Mb (mean +/- SD: 5.11 +/- 0.28 Mb). Phylogroup E shows the largest mean genome size (5.52 Mb), followed by F (5.31 Mb), G (5.27 Mb), D (5.22 Mb), and B2 (5.19 Mb), whereas A and U occupy the smaller end of the distribution. Coding sequence counts range from 4,154 to 6,352 per genome and remain strongly correlated with genome size (Pearson r = 0.97), as expected for bacteria with broadly conserved gene density.",
    )

    replace_by_start(
        doc,
        "CRISPRcas-Typer identified CRISPR arrays",
        "CRISPRcas-Typer identified CRISPR arrays in 78.5% of genomes (785/1,000). Spacer-level analysis enumerated 15,470 spacer records representing 2,511 unique spacer sequences. CRISPR immune repertoire size differed strongly by phylogroup (permutation ANOVA p = 0.001): B2 had low array and spacer burden (median = 0 arrays; median = 0 unique spacer sequences), whereas B1, A, D, C, F, G, and U generally retained larger spacer repertoires. High-confidence target matching identified 65 unique phage-matched spacers (450 spacer occurrences; median best pident = 96.88%) and 262 unique plasmid-matched spacers (2,246 spacer occurrences; median best pident = 100.00%; Supplementary Figure S1). These categories overlapped for 39 unique spacers, yielding 288 unique spacers with at least one high-confidence mobile-element target match. Across all 1,000 genomes, phylogroup was associated with substantial variation in target burden, especially plasmid-matched spacer occurrences (permutation ANOVA F = 133.57, p = 0.001, eta-squared = 0.519, omega-squared = 0.515) and any-targeted spacer occurrences (F = 136.03, p = 0.001, eta-squared = 0.523). In this dataset, phylogroup C showed the highest all-genome fraction of plasmid-matched genomes (100.0%; n = 36), whereas E showed very low target recovery (4.1% with any target). B2 remained CRISPR-sparse overall (18.4% of B2 genomes with any spacer), but spacer-bearing B2 genomes frequently carried target-matched spacers. These results indicate phylogroup-associated heterogeneity in adaptive immune memory rather than a uniformly elevated CRISPR-Cas burden in ExPEC-associated lineages.",
    )

    replace_by_start(
        doc,
        "SignalP-positive proteins were recovered",
        "SignalP-positive proteins were recovered for all 1,000 genomes after quality control of protein-description fields containing unescaped commas. The cohort contained 560,585 SignalP-positive protein records. SignalP burden differed by phylogroup both as total count and after normalization by CDS (permutation ANOVA p = 0.001 for both). Median SignalP-positive proteins per 1,000 CDS ranged from 107.8 in U and 110.8 in E to 120.0 in G, with B2 (118.9), B1 (118.1), C (117.9), F (116.6), A (116.5), and D (115.3) occupying intermediate positions. Thus, secretion-signal abundance is phylogroup-associated, but the data do not support a uniquely D-enriched secretome.",
    )

    replace_by_start(
        doc,
        "TMHMM-positive proteins were identified",
        "TMHMM-positive proteins were identified in 999 of 1,000 genomes, yielding 1,125,213 protein-level records. Only one B1 genome (GCF_045037575.1_ASM4503757v1_genomic.fna) lacked TMHMM detail output and was treated as missing rather than as a true zero. TMHMM burden differed significantly across phylogroups (permutation ANOVA p = 0.001), with substantial effect sizes for absolute TMHMM-positive protein count (eta-squared = 0.415; omega-squared = 0.410) and total predicted helices (eta-squared = 0.525; omega-squared = 0.521), and a moderate effect after normalization by CDS (TMHMM proteins per 1,000 CDS: eta-squared = 0.289; omega-squared = 0.283). Median TMHMM-positive protein counts were highest in G (1,171), E (1,167), F (1,162), C (1,161), and B2 (1,153), whereas A and U were lower. After CDS normalization, the strongest medians were observed in C (238.7), A (237.6), B1 (236.4), B2 (236.3), and G (235.9) TMHMM-positive proteins per 1,000 CDS, indicating phylogroup-associated membrane-protein architecture.",
    )

    replace_by_start(
        doc,
        "COG functional profiling",
        "COG functional profiling across 26 categories revealed broad phylogroup-associated shifts after formal multiple-testing correction: 22 of 26 categories were significant after both Benjamini-Hochberg and Bonferroni correction. The strongest adjusted signal was COG_S (function unknown), highest in E (mean relative abundance 24.2%) and lowest in A (19.7%). COG_L was highest in A (6.70%) and lowest in D (5.76%). B2 showed maxima for COG_G and COG_Q, whereas A showed maxima for COG_C, COG_E, COG_F, COG_J, COG_L, and COG_P. These results support phylogroup-associated functional composition, but not a D-centered mobilome/defense enrichment model.",
    )

    replace_by_start(
        doc,
        "Figure 2. Revised phylogroup",
        "Figure 2. Phylogroup-stratified genome, GEM, secretion, membrane, and CRISPR spacer features. Boxplots show group distributions across genome size, GEM reactions/metabolites, SignalP-positive proteins per 1,000 CDS, TMHMM-positive proteins per 1,000 CDS (n = 999; one unresolved B1 genome treated as missing), and unique CRISPR spacer sequences. Supplementary Figure S1 summarizes phage- and plasmid-targeted CRISPR spacer matches.",
    )

    replace_by_start(
        doc,
        "The revised KEGG matrix",
        "The KEGG matrix contained 264 non-redundant ko pathways after removal of duplicated map columns. Hellinger-transformed pathway abundance PCA explained 19.1% and 9.0% of variance on PC1 and PC2, respectively, while binary pathway PCA explained 15.0% and 8.8%. Global PERMANOVA showed strong pathway structure by phylogroup: binary Jaccard PERMANOVA gave F = 288.60, R2 = 0.6997, p = 0.001, and Hellinger abundance PERMANOVA gave F = 142.58, R2 = 0.5351, p = 0.001. However, PERMDISP was also significant for both binary (F = 26.57, p = 0.001) and Hellinger analyses (F = 48.65, p = 0.001), indicating that the apparent phylogroup signal reflects both centroid separation and heterogeneous within-group dispersion. Therefore, the Hellinger R2 = 0.5351 should not be interpreted as a pure estimate of variance explained by discrete phylogroup identity.",
    )

    replace_exact(doc, "3.3.2 Corrected individual pathway-frequency differences", "3.3.2 Individual pathway-frequency differences")

    replace_contains(
        doc,
        "Colibactin and pks-island-related conclusions",
        "Colibactin and pks-island-related conclusions require direct gene/pathway annotation support beyond the high-value KEGG matrix; the global pathway analysis does not establish B2- or D-wide colibactin enrichment.",
    )
    replace_contains(
        doc,
        "Propionate, methylcitrate",
        "Propionate, methylcitrate, and mucus-sugar utilization pathways remain biologically plausible host-niche candidates, but they should be reported only when present among the significant KEGG pathways or verified by targeted annotation.",
    )
    replace_contains(
        doc,
        "Sugar phosphotransferase",
        "Sugar phosphotransferase and carbohydrate-utilization systems can be discussed as individual pathway-frequency differences when supported by adjusted p values, not as a broad phylogroup-wide pathway-space separation.",
    )
    replace_contains(
        doc,
        "putative central-metabolic retention",
        "Conversely, putative central-metabolic retention in A and B1 genomes should be interpreted as specific adjusted pathway-frequency contrasts, not as evidence that A/B1 and D/B2 occupy globally separated pathway landscapes.",
    )
    replace_contains(
        doc,
        "Together, the revised pathway-level",
        "Together, the pathway-level findings support a nuanced model: discrete phylogroups show selected pathway-frequency differences, but the high-value KEGG pathway composition as a whole is not strongly partitioned by phylogroup after accounting for dispersion and continuous genomic relatedness.",
    )

    replace_by_start(
        doc,
        "Figure 3. KEGG pathway landscape",
        "Figure 3. KEGG pathway landscape in the analyzed cohort. Hellinger-transformed pathway PCA shows strong global phylogroup-associated pathway structure (PERMANOVA R2 = 0.5351, p = 0.001), but PERMDISP is also significant (p = 0.001), indicating that both centroid separation and dispersion heterogeneity contribute. The heatmap shows the 15 pathways significant after BH and Bonferroni correction.",
    )
    replace_by_start(
        doc,
        "Figure 4. Corrected COG",
        "Figure 4. COG functional-category landscape. Heatmap shows mean relative abundance of the top BH-corrected COG categories by phylogroup; 22/26 COG categories remained significant after both BH and Bonferroni correction.",
    )

    replace_by_start(
        doc,
        "Boxplots stratified by phylogroup reveal",
        "Boxplots stratified by phylogroup reveal modest but structured differences in GEM topology. C had the highest mean reaction count (2,432.7), followed by B1 (2,424.3), A (2,407.4), D (2,406.7), F (2,402.8), G (2,382.8), B2 (2,380.2), U (2,361.5), and E (2,347.5). GEM metabolite counts were highest in C and B1 and lowest in E/U. These patterns indicate topology-level heterogeneity, but the magnitude is far smaller than the global pathway-space differences captured by the KEGG matrix.",
    )

    replace_by_start(
        doc,
        "The revised model-topology",
        "The model-topology results are best interpreted as compact, M9-constrained GEM variation superimposed on broader genome-content differences. They do not support a simple large-genome/pathogenic-lineage expansion model, because E has the largest mean genome size but the lowest mean GEM reaction and metabolite counts.",
    )

    replace_by_start(
        doc,
        "Our updated benchmark",
        "The computational benchmark uses a distribution-based assessment across 20 N50-stratified low-contiguity assemblies annotated with phylogroup and CDS metadata (A (n = 4), B1 (n = 7), B2 (n = 4), C (n = 1), E (n = 2), cryptic/U (n = 2); CDS_Count range 4,252-5,584). Optimized gapseq reduced median runtime from 19,821.5 s to 1,260.5 s, corresponding to a median 15.81-fold speedup and a narrow IQR of 15.40-16.35. Five reference model systems were organized for structure and phenotype benchmarking, with four paired runtime records showing a median 12.76-fold speedup and MG1655/iML1515 included in the reaction-set and SD01 FBA comparisons. This distributional benchmark is more useful for cohort-scale planning than a single point estimate.",
    )

    replace_by_start(
        doc,
        "The corrected analysis substantially changes",
        "The integrated analysis supports a conservative biological interpretation. KEGG pathway composition is strongly structured across the cohort, but FastANI analysis shows that most of this structure follows continuous core-genome relatedness rather than discrete phylogroup categories alone. COG composition, SignalP-positive protein burden, CRISPR spacer repertoire structure, CRISPR spacer target burden, and TMHMM-positive membrane-protein architecture retain measurable phylogroup-associated structure, with the strongest residual signal in the combined SignalP/TMHMM surface feature space and the strongest CRISPR target signal in plasmid-matched spacer burden.",
    )

    replace_by_start(
        doc,
        "The revised cohort also weakens",
        "The cohort does not support a D-centered virulence-enrichment narrative. D is not enriched for COG_L, COG_X, COG_V, secretome size, or membrane-protein burden in these outputs. Instead, B2 is notable for sparse but target-matched CRISPR spacer memory among spacer-bearing genomes and high membrane-protein counts, E for high genome size, COG_S abundance, and low spacer target recovery, and A/B1 for several COG categories and adjusted pathway-frequency maxima.",
    )

    replace_by_start(
        doc,
        "Several limitations remain.",
        "Several caveats should guide interpretation. First, pathway names are represented by KEGG identifiers only; pathway-name annotation would improve biological interpretability of the 15 adjusted signals. Second, phylogroup assignment was performed with EzClermont/ClermonTyping-compatible classification before cohort selection; a future MLST, Mash/ANI reference-placement, or core-gene tree concordance analysis would provide an additional independent phylogenomic check. Third, the residual phylogroup analysis uses ANI PCoA covariates and label-permutation sensitivity tests rather than a full core-gene tree with formal phylogenetic comparative modeling; it therefore distinguishes continuous genomic relatedness from categorical labels but does not replace tree-based evolutionary inference. Fourth, one B1 genome lacks detailed TMHMM output and was treated as missing (NA), not as a biological zero, in membrane-protein analyses. Fifth, CRISPR target matching depends on the completeness and taxonomic composition of the phage and plasmid databases used for spacer searches, and phage/plasmid labels can overlap when homologous mobile-element sequences are represented in both databases; therefore target absence should be interpreted as no detectable match in the queried databases rather than true absence of prior exposure. Sixth, the SD01 carbon-source benchmark uses the Monk et al. Table 5 in silico growth-screen matrix, not independent experimental Biolog data; 32/770 conditions remained unmapped automatically, iECs1327 required the iECSF_1327/SE15 proxy, and several gapseq reconstructions showed oxygen-condition-specific sensitivity losses. Reaction compaction should not be interpreted by itself as proof of functional equivalence; phenotype-level and substrate-level validation remain necessary in representative strains.",
    )

    replace_by_start(
        doc,
        "In conclusion, this revised analysis",
        "In conclusion, this study presents an optimized and systematically benchmarked E. coli GEM reconstruction framework applied to a 1,000-genome phylogroup-matched cohort. Stratified benchmarking across 20 assemblies showed a median 15.81-fold runtime acceleration, paired XML reaction matrices identified 111.0 optimized-only reactions per genome at the median, and five reference model systems supported reaction-set and phenotype-benchmark organization, with four paired runtime records showing a median 12.76-fold acceleration. Automated SD01 carbon-source mapping translated 738/770 aerobic/anaerobic environments into model-readable exchanges and produced reproducible FBA prediction matrices for curated BiGG and gapseq models, providing a transparent in silico consistency benchmark while preserving the need for independent experimental validation. FastANI residual analysis showed that KEGG pathway composition is strongly tied to continuous core-genome distance, whereas discrete phylogroup labels retain smaller residual effects for KEGG, COG, GEM topology, CRISPR spacer composition, CRISPR target burden, and SignalP/TMHMM surface architecture. The main validation priority is true multi-strain experimental Biolog or carbon-source data for functional accuracy benchmarking.",
    )

    try:
        doc.save(DOCX)
        saved = DOCX
    except PermissionError:
        doc.save(OUT)
        saved = OUT
    print(f"saved {saved}")
    print(f"backup {BACKUP}")


if __name__ == "__main__":
    main()
