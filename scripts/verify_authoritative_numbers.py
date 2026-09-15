from __future__ import annotations

from pathlib import Path
import json
import os

import numpy as np
import pandas as pd


ROOT = Path(os.environ.get("ECOLI_GEM_ROOT", Path(__file__).resolve().parents[1]))
FINAL = ROOT / "18.mSystems_revision3_final_20260911"
OUT = FINAL / "07_internal_audit"


def close(actual, expected, tolerance=5e-5):
    return bool(abs(float(actual) - float(expected)) <= tolerance)


def record(rows, metric, actual, expected, source, tolerance=5e-5):
    if isinstance(expected, dict):
        passed = actual == expected
    else:
        passed = close(actual, expected, tolerance)
    rows.append({"metric": metric, "actual": actual, "expected_in_manuscript": expected, "pass": passed, "source": source})


def main() -> None:
    rows = []
    cohort_path = ROOT / "09.mSystems_revision_data_freeze" / "01_master_manifest" / "frozen_1000_master_manifest.csv"
    cohort = pd.read_csv(cohort_path)
    counts = cohort["Final_Phylogroup"].value_counts().to_dict()
    expected_counts = {"A": 298, "B2": 244, "B1": 200, "D": 103, "E": 73, "C": 36, "F": 31, "G": 11, "U": 4}
    record(rows, "Final cohort size", len(cohort), 1000, cohort_path)
    record(rows, "Final phylogroup counts", counts, expected_counts, cohort_path)
    record(rows, "Unique Genome_UID", cohort["Genome_UID"].nunique(), 1000, cohort_path)

    ablation_path = FINAL / "02_analysis_evidence" / "01_ablation" / "benchmark_20x6_authoritative_version_summary.csv"
    ablation = pd.read_csv(ablation_path).set_index("Version_ID")
    for version, runtime, reactions, gapfill in [
        ("V1", 19821.5, 3057.0, 46.5), ("V2", 1288.2642915, 2708.5, 95.0),
        ("V3", 1282.6394691, 2708.5, 95.0), ("V4", 1302.2974771, 2708.5, 95.0),
        ("V5", 1051.3987052, 2404.0, 66.5), ("V6", 1260.5, 2403.5, 67.0),
    ]:
        record(rows, f"{version} median runtime", ablation.loc[version, "Wallclock_Seconds_Median"], runtime, ablation_path)
        record(rows, f"{version} median final reactions", ablation.loc[version, "Final_Reactions_SBML_Median"], reactions, ablation_path)
        record(rows, f"{version} median gapfill", ablation.loc[version, "Gapfilled_Reactions_DraftFinal_Median"], gapfill, ablation_path)

    runs_path = FINAL / "02_analysis_evidence" / "01_ablation" / "benchmark_20x6_authoritative_run_table.csv"
    runs = pd.read_csv(runs_path)
    record(rows, "Ablation completed runs", len(runs), 120, runs_path)

    db_path = FINAL / "02_analysis_evidence" / "02_database_comparison" / "v4_v5_reaction_linked_database_summary.csv"
    db = pd.read_csv(db_path).set_index("Metric")["Value"]
    for metric, expected in {
        "V4_only_reaction_occurrences": 6012,
        "V4_only_unique_reactions": 584,
        "V4_only_reaction_occurrences_with_direct_protein": 5516,
        "V4_only_unique_direct_proteins": 6817,
        "V4_only_unique_EC_numbers": 294,
        "V4_only_unique_MetaCyc_pathways": 539,
        "reaction_linked_proteins_both_high_confidence_percent": 78.71497726,
        "reaction_linked_proteins_broad_only_high_confidence_percent": 19.46604078,
        "reaction_linked_proteins_neither_high_confidence_percent": 1.81898196,
    }.items():
        record(rows, metric, db.loc[metric], expected, db_path)

    qc_path = ROOT / "00.data" / "cobrapy_qc_remote_verified_models.csv"
    qc = pd.read_csv(qc_path)
    record(rows, "GEM QC model count", len(qc), 1000, qc_path)
    for column, metric, expected in [
        ("Total_Reactions", "Cohort median reactions", 2403),
        ("Total_Metabolites", "Cohort median metabolites", 1874),
        ("Total_Genes", "Cohort median genes", 1214),
        ("Biomass_Flux", "Cohort median biomass flux", 0.583),
    ]:
        record(rows, metric, qc[column].median(), expected, qc_path, tolerance=0.001)
    for column, metric, expected in [
        ("Mass_Balanced_Ratio(%)", "Mean mass-balanced ratio", 91.64),
        ("Blocked_Reactions_Ratio(%)", "Mean blocked-reaction ratio", 32.98),
        ("Deadend_Metabolites_Ratio(%)", "Mean dead-end ratio", 54.09),
    ]:
        record(rows, metric, qc[column].mean(), expected, qc_path, tolerance=0.01)
    record(rows, "Positive biomass models", int((qc["Biomass_Flux"] > 0).sum()), 1000, qc_path)

    prov_path = ROOT / "10.mSystems_revision2_model_evidence" / "02_model_evidence" / "reaction_evidence_counts_by_genome.csv"
    prov = pd.read_csv(prov_path)
    record(rows, "Reaction-provenance model count", len(prov), 1000, prov_path)
    record(rows, "Total gap-filled reaction occurrences", int(prov["Gap_Filled_Reactions"].sum()), 72296, prov_path)
    record(rows, "Median gap-filled reactions", prov["Gap_Filled_Reactions"].median(), 71, prov_path)
    record(rows, "Median gap-filled fraction percent", prov["Gap_Filled_Fraction"].median() * 100, 2.95, prov_path, tolerance=0.01)
    record(rows, "Median gene-evidenced fraction percent", prov["Gene_Evidenced_Fraction"].median() * 100, 81.83, prov_path, tolerance=0.01)

    filter_path = ROOT / "15.kegg_ontology_filtering_revision3" / "tables" / "kegg_ontology_filter_decision_table.csv"
    filt = pd.read_csv(filter_path)
    status_col = next(c for c in filt.columns if c.lower() in {"decision", "filter_decision", "status", "retain", "retain_for_bacterial_profile"})
    values = filt[status_col].astype(str).str.lower()
    retained = int(values.str.contains("retain|keep|true").sum())
    removed = len(filt) - retained
    record(rows, "KEGG original features", len(filt), 264, filter_path)
    record(rows, "KEGG retained features", retained, 115, filter_path)
    record(rows, "KEGG removed features", removed, 149, filter_path)

    kegg_path = ROOT / "15.kegg_ontology_filtering_revision3" / "tables" / "kegg_filtered_global_statistics.csv"
    kegg = pd.read_csv(kegg_path).set_index("test")
    record(rows, "Filtered KEGG PERMANOVA R2", kegg.loc["PERMANOVA_oneway_phylogroup", "R2"], 0.6558906102, kegg_path)
    record(rows, "Filtered KEGG PERMDISP F", kegg.loc["PERMDISP_phylogroup", "pseudo_F"], 1078.0950294, kegg_path)
    mantel_row = next(i for i in kegg.index if str(i).startswith("Mantel"))
    residual_row = next(i for i in kegg.index if str(i).startswith("ANI20"))
    record(rows, "Filtered KEGG Mantel r", kegg.loc[mantel_row, "mantel_r"], 0.6949735950, kegg_path)
    record(rows, "Filtered KEGG residual R2", kegg.loc[residual_row, "unique_R2"], 0.0041869015, kegg_path)

    factorial_path = ROOT / "16.factorial_medium_thread_benchmark" / "remote_downloaded" / "20260911_fact20_gapseq" / "extracted" / "tables" / "factorial_benchmark_enriched_condition_summary.csv"
    factorial = pd.read_csv(factorial_path)
    for medium, reactions, gapfill in [("default", 2404, 62), ("m9_glucose", 2403.5, 63), ("defined_complete", 2401, 62.5)]:
        subset = factorial[factorial["Medium_ID"] == medium]
        record(rows, f"{medium} median final reactions", subset["Final_Reactions_Median"].median(), reactions, factorial_path)
        record(rows, f"{medium} median gapfill", subset["Gapfill_Added_Sum_Median"].median(), gapfill, factorial_path)

    memote_path = ROOT / "14.mSystems_revision3_submission_candidate" / "02_analysis_evidence" / "memote_18_scores.csv"
    memote = pd.read_csv(memote_path)
    record(rows, "MEMOTE report count", len(memote), 18, memote_path)
    record(rows, "Median MEMOTE total score percent", memote["total_score"].median() * 100, 88.03, memote_path, tolerance=0.01)
    record(rows, "Mean MEMOTE consistency percent", memote["consistency"].mean() * 100, 95.90, memote_path, tolerance=0.01)
    record(rows, "MEMOTE gene annotation percent", memote["annotation_gene"].median() * 100, 33.33, memote_path, tolerance=0.01)

    result = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    result.to_csv(OUT / "authoritative_number_crosscheck.csv", index=False)
    summary = {
        "checks": len(result),
        "passed": int(result["pass"].sum()),
        "failed": int((~result["pass"]).sum()),
        "failed_metrics": result.loc[~result["pass"], "metric"].tolist(),
    }
    (OUT / "authoritative_number_crosscheck.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    if summary["failed"]:
        print(result.loc[~result["pass"]].to_string(index=False))


if __name__ == "__main__":
    main()
