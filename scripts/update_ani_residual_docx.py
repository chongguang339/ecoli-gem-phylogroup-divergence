from pathlib import Path
from shutil import copy2

from docx import Document
from docx.shared import Inches


ROOT = Path(r"D:\pipeline\analysis_v5")
DOCX = ROOT / "02.论文" / "Evolutionary Divergence of Genomic Defense.docx"
BACKUP = ROOT / "02.论文" / "Evolutionary Divergence of Genomic Defense.before_ani_residual_update.docx"
FIG3 = ROOT / "01.figures" / "Figure3_pathway_landscape.png"
FIG4 = ROOT / "01.figures" / "Figure4_cog_corrected_heatmap.png"
FIG5 = ROOT / "01.figures" / "Figure5_ani_residual_phylogroup.png"
INDEX = ROOT / "02.论文" / "docx_text_index_after_ani_residual_update.txt"


REPLACEMENTS = {
    17: "Results: We developed and deployed an optimized high-throughput annotation-and-reconstruction pipeline on 10,951 high-quality E. coli genomes passing CheckM2 quality thresholds, from which 1,000 representative genomes were selected for multi-omic functional profiling and GEM reconstruction. The corrected phylogroup design substantially increased B2 representation (A = 298, B2 = 244, B1 = 200, D = 103, E = 73, C = 36, F = 31, G = 11, U = 4). In the updated benchmark, an N50-stratified 20-genome stress-test subset showed a median 15.81-fold acceleration of optimized gapseq relative to the default workflow (IQR 15.40-16.35; range 14.78-17.77), reducing median wall-clock time from 19,821.5 s to 1,260.5 s. Four available curated-reference genomes showed a median 12.76-fold acceleration (range 11.30-15.52). Core-genome FastANI distances revealed that pathway composition is strongly structured by continuous genomic relatedness (KEGG Mantel Spearman r = 0.760, p = 0.001), whereas discrete phylogroup labels retained only a small KEGG residual effect after ANI correction (unique R2 = 0.017, p = 0.001). Larger residual phylogroup effects persisted for COG composition (unique R2 = 0.062) and SignalP/TMHMM surface architecture (unique R2 = 0.081).",
    19: "Conclusions: This work establishes a corrected, phylogroup-resolved metabolic and functional atlas of E. coli and shows that conclusions drawn from unbalanced phylogroup designs can be materially altered when B2 is adequately represented. The revised data support lineage-associated differences in COG profiles, secretion signals, CRISPR spacer repertoires, membrane-protein architecture, and KEGG pathway composition; however, FastANI residual analysis indicates that much of the pathway signal tracks continuous core-genome relatedness rather than discrete phylogroup categories alone. These results provide a conservative and reproducible framework for SCI-level interpretation of E. coli functional divergence.",
    84: "KEGG pathway matrices were analyzed using the revised high-value pathway table (`genome_kegg_pathway_matrix.csv`). The table contains paired `ko` and `map` columns that are exact duplicates; therefore, only the 264 non-redundant `ko` columns were retained. Because entries are count/abundance values rather than strict 0/1 calls, two complementary analyses were performed: binary presence/absence (value > 0) with Jaccard distances and abundance-based Hellinger transformation with Euclidean distances. PCA was used for ordination, PERMANOVA tested global phylogroup association, and PERMDISP assessed multivariate dispersion. Individual pathway frequency differences were tested by chi-square tests with Benjamini-Hochberg and Bonferroni correction across all 264 pathways. To determine whether pathway structure reflected discrete phylogroup labels or continuous genomic relatedness, FastANI-derived core-genome distances were analyzed in parallel with functional distance matrices.",
    92: "All statistical analyses were performed in Python using reproducible outputs from the corrected cohort. The revised cohort includes A (n = 298), B2 (n = 244), B1 (n = 200), D (n = 103), E (n = 73), C (n = 36), F (n = 31), G (n = 11), and U (n = 4). Genome-level continuous metrics were assessed by permutation ANOVA (999 permutations), and effect sizes were reported as eta-squared and omega-squared to distinguish statistical significance from biological magnitude. COG category comparisons were tested across 26 categories with both Benjamini-Hochberg and Bonferroni correction. KEGG pathway frequency tests were corrected across 264 non-redundant pathways. Balanced PERMANOVA sensitivity analyses were performed for pathway matrices to verify that global beta-diversity conclusions were not driven by sample-size imbalance. Core-genome ANI distance associations were tested by Mantel correlations (999 label permutations). Residual phylogroup effects were estimated by distance-based redundancy-style models using the first 20 ANI PCoA axes as continuous genomic covariates; unique phylogroup R2 was defined as R2(ANI PCs + phylogroup) - R2(ANI PCs).",
    141: "3.3 KEGG Pathway Composition Tracks Core-Genome ANI Structure with Limited Residual Phylogroup Effects",
    145: "The revised KEGG matrix contained 264 non-redundant ko pathways after removal of duplicated map columns. Hellinger-transformed pathway abundance PCA explained 19.1% and 9.0% of variance on PC1 and PC2, respectively, while binary pathway PCA explained 15.0% and 8.8%. Correctly aligned global PERMANOVA showed strong pathway structure by phylogroup: binary Jaccard PERMANOVA gave F = 288.60, R2 = 0.6997, p = 0.001, and Hellinger abundance PERMANOVA gave F = 142.58, R2 = 0.5351, p = 0.001. However, PERMDISP was also significant for both binary (F = 26.57, p = 0.001) and Hellinger analyses (F = 48.65, p = 0.001), indicating that the apparent phylogroup signal reflects both centroid separation and heterogeneous within-group dispersion. Therefore, the Hellinger R2 = 0.5351 should not be interpreted as a pure estimate of variance explained by discrete phylogroup identity.",
    147: "Balanced PERMANOVA sensitivity analyses confirmed that the global pathway signal was not an artifact of unequal group sizes. In repeated n = 31 subsampling across A, B2, B1, D, E, C, and F, mean Hellinger PERMANOVA R2 was 0.5784; in repeated n = 50 subsampling across A, B2, B1, D, and E, mean R2 was 0.5939. These balanced results support a genuine pathway-composition structure associated with lineage background, but PERMDISP and ANI analyses below show that this structure is better interpreted as continuous phylogenomic relatedness with smaller residual categorical phylogroup effects.",
    149: "3.3.2 Corrected individual pathway-frequency differences",
    170: "Figure 3. KEGG pathway landscape in the corrected cohort. Hellinger-transformed pathway PCA shows strong global phylogroup-associated pathway structure (PERMANOVA R2 = 0.5351, p = 0.001), but PERMDISP is also significant (p = 0.001), indicating that both centroid separation and dispersion heterogeneity contribute. The heatmap shows the 15 pathways significant after BH and Bonferroni correction.",
    185: "The corrected analysis substantially changes the biological interpretation of the study. KEGG pathway composition is strongly structured across the cohort, but FastANI analysis shows that most of this structure follows continuous core-genome relatedness rather than discrete phylogroup categories alone. COG composition, SignalP-positive protein burden, CRISPR spacer repertoire structure, and TMHMM-positive membrane-protein architecture retain measurable residual phylogroup effects after accounting for ANI PCs, with the strongest residual signal in the combined SignalP/TMHMM surface feature space.",
    188: "Several limitations remain. First, pathway names are represented by KEGG identifiers only; pathway-name annotation would improve biological interpretability of the 15 corrected signals. Second, the residual phylogroup analysis uses ANI PCoA covariates and label-permutation sensitivity tests rather than a full core-gene tree with formal phylogenetic comparative modeling; it therefore distinguishes continuous genomic relatedness from categorical labels but does not replace tree-based evolutionary inference. Third, although TMHMM coverage was improved from 846 to 999 genomes after reprocessing, one B1 genome still lacks detailed TMHMM output and was treated as missing in membrane-protein analyses. Fourth, CRISPR spacer sharing is sequence-level but not yet matched to phage/protospacer databases, so spacer targets and ecological exposures remain unresolved. Finally, all GEM predictions remain in silico and should be experimentally validated in representative strains.",
    194: "In conclusion, this revised analysis presents an optimized and more carefully validated E. coli GEM reconstruction framework applied to a corrected 1,000-genome phylogroup cohort. Stratified benchmarking across 20 assemblies showed a median 15.81-fold runtime acceleration, and four available curated-reference genomes showed a median 12.76-fold acceleration with mostly non-reference-supported reactions removed by optimization. The rebalanced cohort eliminated the earlier B2 sampling defect and substantially weakened the claim of D-specific virulence enrichment. FastANI residual analysis further showed that KEGG pathway composition is strongly tied to continuous core-genome distance, whereas discrete phylogroup labels retain smaller residual effects for KEGG, COG, GEM topology, CRISPR spacer composition, and SignalP/TMHMM surface architecture. Remaining publication-critical work is concentrated in adding true multi-strain experimental Biolog or carbon-source data for functional accuracy validation.",
}


INSERT_BEFORE_INDEX = 162


def main():
    copy2(DOCX, BACKUP)
    doc = Document(DOCX)
    for idx, text in REPLACEMENTS.items():
        doc.paragraphs[idx].text = text

    # Replace regenerated Figure 3 and Figure 4 images.
    if len(doc.inline_shapes) >= 4:
        for shape, fig, height in [
            (doc.inline_shapes[2], FIG3, 4.19),
            (doc.inline_shapes[3], FIG4, 3.63),
        ]:
            blip = shape._inline.graphic.graphicData.pic.blipFill.blip
            doc.part.related_parts[blip.embed]._blob = fig.read_bytes()
            shape.width = Inches(6.7)
            shape.height = Inches(height)

    target = doc.paragraphs[INSERT_BEFORE_INDEX]
    target.insert_paragraph_before(
        "3.3.3 Core-genome ANI distance and residual phylogroup effects"
    )
    target.insert_paragraph_before(
        "The FastANI distance matrix contained all 1,000 genomes, had no missing values, and was symmetrized by averaging reciprocal distances. Directional asymmetry was small (mean absolute asymmetry = 0.00047; maximum = 0.00317), with a median pairwise ANI distance of 0.02768 and a maximum of 0.04460. ANI PCoA showed strong clustering consistent with known E. coli phylogenomic structure; the first 20 positive ANI PCoA axes explained 76.0% of ANI distance variation and were used as continuous genomic covariates in residual models."
    )
    target.insert_paragraph_before(
        "Mantel analyses showed that KEGG pathway distance was strongly coupled to core-genome ANI distance (Spearman r = 0.760, p = 0.001), whereas COG composition (r = 0.118), GEM topology (r = 0.261), CRISPR spacer composition (r = 0.353), and combined SignalP/TMHMM surface features (r = 0.114) showed weaker but significant associations. After accounting for the 20 ANI PCoA axes, discrete phylogroup labels retained only a small residual effect for KEGG pathways (unique R2 = 0.017, p = 0.001), GEM topology (unique R2 = 0.020), and CRISPR spacer composition (unique R2 = 0.018), but larger residual effects for COG composition (unique R2 = 0.062) and SignalP/TMHMM surface architecture (unique R2 = 0.081). Thus, pathway composition primarily tracks continuous core-genome relatedness, while envelope/surface features and COG composition preserve more categorical phylogroup signal beyond ANI structure."
    )
    p_img = target.insert_paragraph_before()
    p_img.alignment = 1
    p_img.add_run().add_picture(str(FIG5), width=Inches(6.7))
    target.insert_paragraph_before(
        "Figure 5. Core-genome ANI structure and residual phylogroup effects. (A) ANI PCoA of the symmetrized FastANI distance matrix. (B) Mantel Spearman correlations between ANI distance and functional distance spaces. (C) Distance-based variance partitioning showing variance associated with ANI PCoA axes and the additional residual contribution of discrete phylogroup labels."
    )

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
