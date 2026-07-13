from pathlib import Path

import pandas as pd


ROOT = Path(r"D:\pipeline\analysis_v5")
OUT = ROOT / "03.analysis"
SUPP = ROOT / "05.supplementary_data" / "Data_S1_cohort_and_phylogroup"


def main():
    sampling = pd.DataFrame(
        [
            ("A", 3228, 296, 29.48),
            ("B2", 2661, 243, 24.30),
            ("B1", 2183, 199, 19.93),
            ("D", 1123, 103, 10.25),
            ("E", 791, 72, 7.22),
            ("C", 396, 36, 3.62),
            ("F", 335, 31, 3.06),
            ("G", 125, 11, 1.14),
            ("U", 48, 4, 0.44),
            ("CRYPTIC", 29, 3, 0.26),
            ("U/CRYPTIC", 16, 1, 0.15),
            ("EC_CONTROL_FAIL", 16, 1, 0.15),
        ],
        columns=["Source_Phylogroup_Label", "Source_Count", "Planned_Sample_Count", "Source_Percent"],
    )
    sampling["Sampling_Fraction"] = sampling["Planned_Sample_Count"] / sampling["Source_Count"]
    sampling["Allocation_Note"] = "Proportional allocation from 10,951 EzClermont-typed genomes; non-standard labels harmonized for final analytical categories."
    sampling_path = OUT / "phylogroup_proportional_sampling_plan.csv"
    sampling.to_csv(sampling_path, index=False)
    (SUPP / sampling_path.name).write_text(sampling_path.read_text(encoding="utf-8"), encoding="utf-8")

    recalc = pd.read_csv(ROOT / "00.data" / "recalculated_annotations.csv")
    qc = recalc[["Genome_ID", "Phylogroup", "GEM_Genes", "GEM_Reactions", "GEM_Metabolites"]].copy()
    qc["Viability_Status"] = "Viable_in_cohort_topology_table"
    qc["Blocked_Reaction_Fraction"] = "not_computed_from_current_outputs"
    qc["Dead_End_Metabolites"] = "not_computed_from_current_outputs"
    qc["Mass_Balance_Status"] = "not_computed_from_current_outputs"
    qc["Exchange_Reaction_Count"] = "not_computed_from_current_outputs"
    qc["Biomass_Reaction_ID"] = "not_computed_from_current_outputs"
    qc["Gapfilled_Reaction_Count"] = "not_computed_from_current_outputs"
    qc_path = OUT / "gem_model_qc_available_metrics.csv"
    qc.to_csv(qc_path, index=False)
    (SUPP / qc_path.name).write_text(qc_path.read_text(encoding="utf-8"), encoding="utf-8")
    print(sampling_path)
    print(qc_path)


if __name__ == "__main__":
    main()
