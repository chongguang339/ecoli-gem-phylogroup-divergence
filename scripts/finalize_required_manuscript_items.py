from __future__ import annotations

from pathlib import Path
import csv
import shutil
import zipfile

from docx import Document


ROOT = Path(r"D:\pipeline\analysis_v5")
DOCX = ROOT / "02.论文" / "Evolutionary Divergence of Genomic Defense.docx"
TEXT_INDEX = ROOT / "02.论文" / "docx_text_index_after_required_finalization.txt"
SUPP = ROOT / "05.supplementary_data"


def set_para_text(para, text: str) -> None:
    for run in para.runs:
        run.text = ""
    if para.runs:
        para.runs[0].text = text
    else:
        para.add_run(text)


def replace_starting(doc: Document, prefix: str, new_text: str) -> int:
    n = 0
    for para in doc.paragraphs:
        if para.text.strip().startswith(prefix):
            set_para_text(para, new_text)
            n += 1
    return n


def copy_if_exists(src: Path, dst_dir: Path, manifest_rows: list[dict], data_id: str, description: str) -> None:
    dst_dir.mkdir(parents=True, exist_ok=True)
    if src.exists():
        dst = dst_dir / src.name
        if not dst.exists() or dst.stat().st_size != src.stat().st_size:
            shutil.copy2(src, dst)
        status = "copied"
        size = src.stat().st_size
    else:
        dst = dst_dir / src.name
        status = "missing_needs_user_or_regeneration"
        size = ""
    manifest_rows.append(
        {
            "Supplementary_Data": data_id,
            "Description": description,
            "Source_Path": str(src),
            "Packaged_Path": str(dst.relative_to(SUPP)),
            "Status": status,
            "Bytes": size,
        }
    )


def write_readme(folder: Path, title: str, lines: list[str]) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "README.txt").write_text(title + "\n" + "=" * len(title) + "\n\n" + "\n".join(lines) + "\n", encoding="utf-8")


def package_supplementary_data() -> None:
    SUPP.mkdir(parents=True, exist_ok=True)
    manifest_rows: list[dict] = []

    s1 = SUPP / "Data_S1_cohort_and_phylogroup"
    for src, desc in [
        (SUPP / "Data_S1_cohort_and_phylogroup" / "selected_1000_genomes.tsv", "Selected 1,000-genome cohort table reconstructed from the final Genome_ID, phylogroup, and assembly metric matrix."),
        (ROOT / "00.data" / "genome_assembly_and_cog_metrics.csv", "Genome IDs, phylogroups, assembly metrics, and COG counts for the 1,000-genome cohort."),
        (ROOT / "00.data" / "recalculated_annotations.csv", "Per-genome GEM, CRISPR, SignalP, and TMHMM summary annotations."),
        (ROOT / "03.analysis" / "data_qc_summary.csv", "Quality-control summary for corrected analysis tables."),
    ]:
        copy_if_exists(src, s1, manifest_rows, "Data S1", desc)
    write_readme(s1, "Data S1: cohort and phylogroup tables", [
        "Contains the main 1,000-genome cohort matrix and per-genome annotation summaries used throughout the manuscript.",
        "Note: the manuscript still needs a final accession-level selected_1000_genomes.tsv if this file is maintained outside analysis_v5.",
    ])

    s2 = SUPP / "Data_S2_functional_annotation_matrices"
    for src, desc in [
        (ROOT / "00.data" / "genome_kegg_pathway_matrix.csv", "KEGG pathway abundance matrix."),
        (ROOT / "00.data" / "all_genomes_crispr_detailed.csv", "CRISPR array and subtype detail table."),
        (ROOT / "00.data" / "all_genomes_crispr_spacers_detailed.csv", "Spacer-level CRISPR detail table."),
        (ROOT / "00.data" / "all_genomes_signalp_positives_detailed.csv", "Protein-level SignalP-positive predictions."),
        (ROOT / "00.data" / "all_genomes_tmhmm_positives_detailed.csv", "Protein-level TMHMM-positive predictions after reprocessing missing genomes."),
    ]:
        copy_if_exists(src, s2, manifest_rows, "Data S2", desc)
    write_readme(s2, "Data S2: functional annotation matrices", [
        "Contains COG/KEGG, CRISPR, SignalP, and TMHMM inputs underlying the genome-level comparative analyses.",
        "Large protein-level SignalP/TMHMM tables are included as CSV files for reproducibility.",
    ])

    s3 = SUPP / "Data_S3_statistics_and_ani"
    for src, desc in [
        (ROOT / "03.analysis" / "continuous_phylogroup_effect_sizes_eta2.csv", "Effect sizes for continuous genome-level metrics."),
        (ROOT / "03.analysis" / "detail_metric_permutation_anova.csv", "Permutation ANOVA outputs for detailed metrics."),
        (ROOT / "03.analysis" / "cog_chisq_multipletest.csv", "COG multiple-testing corrected comparisons."),
        (ROOT / "03.analysis" / "pathway_chisq_multipletest.csv", "KEGG pathway multiple-testing corrected comparisons."),
        (ROOT / "03.analysis" / "pathway_hellinger_permanova_full.csv", "Full Hellinger PERMANOVA result."),
        (ROOT / "03.analysis" / "pathway_hellinger_permdisp_full.csv", "Full Hellinger PERMDISP result."),
        (ROOT / "03.analysis" / "pathway_hellinger_balanced_permanova_summary.csv", "Balanced PERMANOVA sensitivity summary."),
        (ROOT / "03.analysis" / "ani_mantel_distance_correlations.csv", "FastANI Mantel distance correlations."),
        (ROOT / "03.analysis" / "ani_residual_phylogroup_effects.csv", "Residual phylogroup effects after ANI PCoA covariates."),
        (ROOT / "00.data" / "core_genome_ani_distance_matrix.csv", "FastANI-derived core-genome distance matrix."),
    ]:
        copy_if_exists(src, s3, manifest_rows, "Data S3", desc)
    write_readme(s3, "Data S3: statistics and ANI residual analysis", [
        "Contains effect sizes, multiple-testing outputs, PERMANOVA/PERMDISP summaries, and ANI residual analyses.",
        "The very large raw FastANI text output is not duplicated here unless required by the target repository; it remains at 00.data/fastani_raw_output.txt.",
    ])

    s4 = SUPP / "Data_S4_runtime_and_reaction_benchmark"
    for src, desc in [
        (ROOT / "04.基准测试" / "modle_compare" / "benchmark_metrics.csv", "Curated-reference model runtime benchmark metrics."),
        (ROOT / "04.基准测试" / "stratified20_derived_metrics.csv", "N50-stratified 20-genome derived runtime benchmark metrics."),
        (ROOT / "00.data" / "all_20_genomes_metadata.csv", "Phylogroup and CDS metadata for the 20 benchmark genomes."),
        (ROOT / "03.analysis" / "benchmark_20_with_phylogroup_cds_reaction_delta.csv", "Benchmark table with phylogroup, CDS, and reaction delta annotations."),
        (ROOT / "04.基准测试" / "modle_compare" / "reaction_set_compare" / "reaction_set_compare_summary.csv", "Reaction-set comparison summary for default/optimized/reference models."),
        (ROOT / "04.基准测试" / "modle_compare" / "reaction_set_compare" / "reaction_annotation_coverage.csv", "Reaction annotation coverage for set comparisons."),
        (ROOT / "04.基准测试" / "modle_compare" / "reaction_set_compare" / "all_dropped_reactions_A_minus_B.csv", "All reactions present in default but absent from optimized reconstructions."),
    ]:
        copy_if_exists(src, s4, manifest_rows, "Data S4", desc)
    write_readme(s4, "Data S4: runtime and reaction-set benchmark", [
        "Contains the N50-stratified speed benchmark and reaction-set comparison outputs.",
        "Default/optimized model XML files remain in 04.基准测试/modle_compare and should be deposited in the external repository if file size permits.",
    ])

    s5 = SUPP / "Data_S5_sd01_exchange_mapping"
    for src, desc in [
        (ROOT / "03.analysis" / "monk2013_sd01_exchange_mapping" / "sd01_770_carbon_conditions_exchange_id_mapping.csv", "Automated mapping from SD01 carbon-source names to BiGG and gapseq exchange IDs."),
        (ROOT / "03.analysis" / "monk2013_sd01_exchange_mapping" / "sd01_770_carbon_conditions_unmapped.csv", "The 32 SD01 carbon-source conditions not mapped automatically."),
        (ROOT / "03.analysis" / "monk2013_sd01_exchange_mapping" / "sd01_exchange_mapping_summary.csv", "Mapping summary counts."),
        (ROOT / "03.analysis" / "monk2013_sd01_table5_truth" / "selected_models_table5_carbon_sources_long.csv", "Extracted SD01 Table 5 carbon-source in silico growth values for selected models."),
        (ROOT / "03.analysis" / "monk2013_sd01_table5_truth" / "requested_model_to_sd01_strain_mapping.csv", "Mapping between requested model labels and SD01 strain/model rows."),
    ]:
        copy_if_exists(src, s5, manifest_rows, "Data S5", desc)
    write_readme(s5, "Data S5: SD01 exchange mapping", [
        "Contains the automated 770-condition carbon-source mapping used for COBRApy simulations.",
        "SD01 Table 5 is treated as an in silico growth-screen matrix, not independent experimental Biolog truth.",
    ])

    s6 = SUPP / "Data_S6_sd01_fba_predictions"
    for src, desc in [
        (ROOT / "03.analysis" / "monk2013_sd01_cobra_fba" / "selected_models_770_carbon_fba_predictions_long.csv", "Long-format curated BiGG model FBA predictions for 770 SD01 carbon conditions."),
        (ROOT / "03.analysis" / "monk2013_sd01_cobra_fba" / "selected_models_770_carbon_fba_prediction_binary_wide.csv", "Wide binary prediction matrix for curated BiGG models."),
        (ROOT / "03.analysis" / "monk2013_sd01_cobra_fba" / "selected_models_770_carbon_fba_prediction_growth_rates_wide.csv", "Wide growth-rate prediction matrix for curated BiGG models."),
        (ROOT / "03.analysis" / "monk2013_sd01_cobra_fba" / "selected_models_770_carbon_fba_accuracy_metrics.csv", "Curated BiGG model confusion-matrix metrics."),
        (ROOT / "03.analysis" / "monk2013_sd01_cobra_fba" / "cobra_model_load_and_media_status.csv", "Curated model load and media coverage status."),
        (ROOT / "03.analysis" / "monk2013_sd01_available_benchmark_fba" / "available_benchmark_770_fba_predictions_long.csv", "Long-format available gapseq/reference benchmark FBA predictions."),
        (ROOT / "03.analysis" / "monk2013_sd01_available_benchmark_fba" / "available_benchmark_770_fba_accuracy_metrics.csv", "Available paired gapseq benchmark confusion-matrix metrics."),
        (ROOT / "03.analysis" / "monk2013_sd01_available_benchmark_fba" / "available_benchmark_model_load_status.csv", "Available benchmark model load status."),
    ]:
        copy_if_exists(src, s6, manifest_rows, "Data S6", desc)
    write_readme(s6, "Data S6: SD01 FBA predictions", [
        "Contains all 770-condition FBA prediction, growth-rate, load-status, and accuracy outputs.",
        "COBRApy simulations used GLPK, biomass > 1e-6 as the binary growth threshold, carbon uptake lower bound -20, and condition-specific oxygen bounds.",
    ])

    s7 = SUPP / "Data_S7_crispr_spacer_targets"
    target_dir = ROOT / "03.analysis" / "crispr_target_matching_analysis"
    for src, desc in [
        (target_dir / "crispr_target_match_summary.csv", "CRISPR spacer target-matching summary table."),
        (target_dir / "crispr_target_matching_summary.txt", "Plain-text CRISPR spacer target-matching analysis summary."),
        (target_dir / "crispr_target_matching_genome_summary.csv", "Genome-level CRISPR target burden table for spacer-bearing genomes."),
        (target_dir / "crispr_target_matching_all_genome_summary.csv", "Genome-level CRISPR target burden table with all genomes represented."),
        (target_dir / "crispr_target_matching_phylogroup_summary.csv", "Phylogroup-level CRISPR target burden summary."),
        (target_dir / "crispr_target_matching_phylogroup_summary_all_genomes.csv", "Phylogroup-level CRISPR target burden summary with all genomes represented."),
        (target_dir / "crispr_target_matching_phylogroup_effect_sizes.csv", "Effect sizes for CRISPR target burden comparisons."),
        (target_dir / "crispr_target_matching_phylogroup_effect_sizes_all_genomes.csv", "Effect sizes for all-genome CRISPR target burden comparisons."),
        (target_dir / "crispr_target_hits_best_mismatch_le2.csv", "High-confidence best CRISPR target hits using mismatch <= 2."),
        (target_dir / "crispr_target_hits_best_mismatch_le3.csv", "Sensitivity-threshold best CRISPR target hits using mismatch <= 3."),
        (target_dir / "crispr_target_hits_all_raw.csv", "Parsed raw CRISPR spacer target hits."),
        (target_dir / "crispr_spacer_occurrences_with_target_flags.csv", "Spacer occurrence table with phage/plasmid target flags."),
        (target_dir / "crispr_target_matching_top_targets.csv", "Top matched phage/plasmid target accessions or descriptions."),
        (ROOT / "03.analysis" / "crispr_target_matching_inputs" / "spacer_vs_phage.txt", "Raw spacer-vs-phage BLAST-like matches provided by user."),
        (ROOT / "03.analysis" / "crispr_target_matching_inputs" / "spacer_vs_plasmid.txt", "Raw spacer-vs-plasmid BLAST-like matches provided by user."),
    ]:
        copy_if_exists(src, s7, manifest_rows, "Data S7", desc)
    write_readme(s7, "Data S7: CRISPR spacer target matching", [
        "Contains high-confidence phage/plasmid spacer target-matching outputs and raw match tables where available.",
        "Target classes are database-match categories and should not be interpreted as mutually exclusive biological labels.",
    ])

    with (SUPP / "supplementary_data_manifest.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=["Supplementary_Data", "Description", "Source_Path", "Packaged_Path", "Status", "Bytes"])
        writer.writeheader()
        writer.writerows(manifest_rows)

    with zipfile.ZipFile(SUPP / "supplementary_data_package_index_only.zip", "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in SUPP.rglob("README.txt"):
            zf.write(path, path.relative_to(SUPP))
        zf.write(SUPP / "supplementary_data_manifest.csv", "supplementary_data_manifest.csv")


def update_docx() -> None:
    doc = Document(DOCX)
    replace_starting(
        doc,
        "A total of 10,951 complete Escherichia coli genome assemblies",
        "A total of 10,951 complete Escherichia coli genome assemblies were retrieved from the NCBI Assembly database after applying the download-time completeness filter. Genome quality was then assessed using CheckM2 v1.0.2 with the neural network-specific completeness model [12], and genomes passing strict quality thresholds (completeness >= 95%, contamination <= 5%, coding density >= 0.85) were retained for downstream sampling. This source set had a median genome size of 5.13 Mb, median GC content of 50.7%, and a median of 4,863 coding sequences (Table 1). Representative cohort and assembly metrics are provided in Supplementary Data S1.",
    )
    replace_starting(doc, "Table 2. Cohort composition", "Table 1. Cohort composition, benchmark datasets, and main output files.")
    replace_starting(doc, "A complete parameter specification is provided in Table 1.", "A complete parameter specification is provided in Table 2.")
    replace_starting(doc, "Table 1. Optimized gapseq", "Table 2. Optimized gapseq reconstruction parameters and rationale.")
    replace_starting(
        doc,
        "Structural reference-model benchmarking was performed for the four curated-reference genomes available",
        "Structural reference-model benchmarking was performed for the four curated-reference genomes available in the current `modle_compare` directory: iB21_1397, iEC1364_W, iECDH10B_1368, and iJN1463. For each genome, default and optimized gapseq reconstructions were compared with the corresponding BiGG reference model by reaction counts and by the fraction of reactions removed by optimization that were nevertheless present in the curated reference. We further implemented a reproducible SD01 carbon-source FBA benchmark. SD01 Table 5 condition names were parsed from 770 aerobic/anaerobic carbon-source environments and translated into exchange reactions using the SD01 Table 3 metabolite dictionary. The resulting map contained both BiGG-style exchange identifiers (e.g., EX_glc__D_e) for curated BiGG models and gapseq/ModelSEED-style exchange identifiers (e.g., EX_cpd00027_e0) for gapseq reconstructions. FBA was run with COBRApy v0.31.1 using the GLPK solver under an M9-like base medium derived from SD01 Table 6. Glucose was excluded from the base medium and one candidate carbon source was opened per condition with a lower bound of -20 mmol gDW-1 h-1. Oxygen uptake was set to -20 mmol gDW-1 h-1 for C_AERO conditions and 0 for C_ANAEROBIC conditions. Binary growth was defined as biomass flux > 1e-6, with a 1e-4 sensitivity threshold also exported. Solver timeouts were capped at 5 s per condition, and simulation/load status fields were retained in the output tables. Prediction matrices, growth-rate matrices, and confusion-matrix metrics were exported for all successfully loaded models. Because SD01 Table 5 is explicitly an in silico growth-screen matrix from Monk et al. rather than an independent experimental Biolog table, the analysis was interpreted as model-to-literature consistency benchmarking, not as experimental phenotype validation.",
    )
    replace_starting(
        doc,
        "Figure 1. Stratified runtime, reference-model, and carbon-source consistency benchmarking.",
        "Figure 1. Stratified runtime and reference-model benchmarking. (A) Speedup distribution across 20 assemblies spanning Actual_N50 = 22,535-675,909 bp and CDS_Count = 4,252-5,584; the subset includes A (n = 4), B1 (n = 7), B2 (n = 4), C (n = 1), E (n = 2), cryptic/U (n = 2). (B) Default versus optimized wall-clock runtime distributions. (C) Reaction-count compaction and reaction retention in the 20-genome stress-test subset; paired XML reaction matrices further identify optimized-only and default-only reactions for every genome. (D) Four available curated-reference genomes show that most reactions removed by the optimized workflow were not present in the matched BiGG reference model. SD01 Table 5 carbon-source mapping and 770-condition FBA outputs are provided in Supplementary Data S5-S6 as an in silico model-consistency benchmark, not an experimental Biolog validation dataset.",
    )
    replace_starting(
        doc,
        "Availability of data and materials:",
        "Availability of data and materials: All genome accession numbers, phylogroup assignments, annotation outputs, GEM topology metrics, benchmark metadata, reaction-set matrices, SD01 carbon-source exchange mappings, and FBA prediction/accuracy matrices are organized in `05.supplementary_data` and will be deposited with the manuscript as Supplementary Data Files. Genome sequences were obtained from NCBI Assembly (https://www.ncbi.nlm.nih.gov/assembly/). Modified gapseq scripts and pipeline configuration files are available at [repository URL to be supplied before submission].",
    )
    doc.save(DOCX)
    refreshed = Document(DOCX)
    with TEXT_INDEX.open("w", encoding="utf-8") as handle:
        for i, para in enumerate(refreshed.paragraphs):
            if para.text.strip():
                handle.write(f"[{i}] {para.text.strip()}\n")


def main() -> None:
    package_supplementary_data()
    update_docx()
    print(f"Updated DOCX: {DOCX}")
    print(f"Supplementary data package: {SUPP}")
    print(f"Text index: {TEXT_INDEX}")


if __name__ == "__main__":
    main()
