from __future__ import annotations

from pathlib import Path
import re

import pandas as pd


ROOT = Path(r"D:\pipeline\analysis_v5")
XLSX = Path(r"D:\pipeline\data\sd01.xlsx")
OUT = ROOT / "03.analysis" / "monk2013_sd01_table5_truth"


REQUESTED = ["iJO1366", "iECs1327", "ic_1306", "iUTI89_1310", "iML1515"]


def normalize_model(s: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def condition_class(cond: str) -> str:
    parts = str(cond).split("_", 2)
    if len(parts) >= 2:
        return parts[0] + "_" + parts[1]
    return str(cond)


def clean_condition(cond: str) -> str:
    s = str(cond)
    for prefix in [
        "C_AERO_",
        "C_ANAEROBIC_",
        "N_AERO_",
        "N_ANAEROBIC_",
        "P_AERO_",
        "P_ANAEROBIC_",
        "S_AERO_",
        "S_ANAEROBIC_",
    ]:
        if s.startswith(prefix):
            return s[len(prefix) :]
    return s


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    table1 = pd.read_excel(XLSX, sheet_name="Table 1", header=2)
    mapping_table = table1[["Strain", "Model Name"]].dropna().copy()
    mapping_table["Model_Normalized"] = mapping_table["Model Name"].map(normalize_model)
    mapping_table.to_csv(OUT / "sd01_table1_model_to_strain_mapping.csv", index=False, encoding="utf-8-sig")

    model_by_norm = dict(zip(mapping_table["Model_Normalized"], mapping_table["Model Name"]))
    strain_by_norm = dict(zip(mapping_table["Model_Normalized"], mapping_table["Strain"]))

    selected = []
    candidates = []
    for requested in REQUESTED:
        norm = normalize_model(requested)
        if norm in strain_by_norm:
            selected.append(
                {
                    "Requested_Model": requested,
                    "SD01_Model_Name": model_by_norm[norm],
                    "SD01_Strain": strain_by_norm[norm],
                    "Match_Status": "exact_model_match",
                    "Truth_Source_Note": "Direct SD01 Table 1 model-to-strain mapping.",
                }
            )
            continue

        if requested == "iML1515":
            mg_norm = normalize_model("iJO1366")
            selected.append(
                {
                    "Requested_Model": requested,
                    "SD01_Model_Name": model_by_norm[mg_norm],
                    "SD01_Strain": strain_by_norm[mg_norm],
                    "Match_Status": "same_strain_mg1655_proxy",
                    "Truth_Source_Note": "iML1515 is an E. coli K-12 MG1655 reconstruction; SD01 provides MG1655 truth through iJO1366.",
                }
            )
            continue

        # Fuzzy rescue for the user-provided iECs1327 spelling.
        cand = mapping_table[
            mapping_table["Model_Normalized"].str.contains(norm, regex=False)
            | mapping_table["Model Name"].astype(str).str.contains(requested, case=False, regex=False)
            | mapping_table["Model Name"].astype(str).str.contains("1327", case=False, regex=False)
            | mapping_table["Strain"].astype(str).str.contains(requested, case=False, regex=False)
        ].copy()
        if not cand.empty:
            cand.insert(0, "Requested_Model", requested)
            candidates.append(cand)
        if requested == "iECs1327":
            # Closest SD01 model with 1327 suffix. This is not an exact spelling match.
            row = mapping_table[mapping_table["Model Name"].eq("iECSF_1327")].iloc[0]
            selected.append(
                {
                    "Requested_Model": requested,
                    "SD01_Model_Name": row["Model Name"],
                    "SD01_Strain": row["Strain"],
                    "Match_Status": "fuzzy_1327_match_verify",
                    "Truth_Source_Note": "No exact iECs1327 model exists in SD01 Table 1; extracted closest 1327-suffix model iECSF_1327 (E. coli SE15). Verify that this is the intended strain.",
                }
            )
        else:
            selected.append(
                {
                    "Requested_Model": requested,
                    "SD01_Model_Name": "",
                    "SD01_Strain": "",
                    "Match_Status": "not_found",
                    "Truth_Source_Note": "No exact or supported proxy mapping found.",
                }
            )

    selected_df = pd.DataFrame(selected)
    selected_df.to_csv(OUT / "requested_model_to_sd01_strain_mapping.csv", index=False, encoding="utf-8-sig")
    if candidates:
        pd.concat(candidates, ignore_index=True).to_csv(OUT / "requested_model_fuzzy_candidates.csv", index=False, encoding="utf-8-sig")

    raw = pd.read_excel(XLSX, sheet_name="Table 5", header=None)
    headers = raw.iloc[2].tolist()
    table = raw.iloc[3:].copy()
    table.columns = headers
    table = table.rename(columns={"Strain:": "Strain"})
    table = table[table["Strain"].notna()].copy()
    table = table[~table["Strain"].astype(str).str.startswith("Total Growing Strains")]

    condition_cols = [
        c for c in table.columns
        if isinstance(c, str)
        and (
            c.startswith("C_AERO_")
            or c.startswith("C_ANAEROBIC_")
            or c.startswith("N_AERO_")
            or c.startswith("N_ANAEROBIC_")
            or c.startswith("P_AERO_")
            or c.startswith("P_ANAEROBIC_")
            or c.startswith("S_AERO_")
            or c.startswith("S_ANAEROBIC_")
        )
    ]

    table5_csv = OUT / "sd01_table5_all_strains_all_conditions_growth_rates.csv"
    table[["Strain"] + condition_cols].to_csv(table5_csv, index=False, encoding="utf-8-sig")

    valid_selected = selected_df[selected_df["SD01_Strain"].astype(str).str.len() > 0].copy()
    selected_rows = []
    for row in valid_selected.itertuples(index=False):
        sdf = table[table["Strain"].eq(row.SD01_Strain)]
        if sdf.empty:
            selected_rows.append({**row._asdict(), "Extraction_Status": "strain_not_found_in_table5"})
            continue
        record = sdf.iloc[0]
        for cond in condition_cols:
            value = pd.to_numeric(record[cond], errors="coerce")
            selected_rows.append(
                {
                    "Requested_Model": row.Requested_Model,
                    "SD01_Model_Name": row.SD01_Model_Name,
                    "SD01_Strain": row.SD01_Strain,
                    "Match_Status": row.Match_Status,
                    "Truth_Source_Note": row.Truth_Source_Note,
                    "Condition": cond,
                    "Condition_Class": condition_class(cond),
                    "Substrate_or_Exchange": clean_condition(cond),
                    "Growth_Rate": value,
                    "Growth_Binary": int(value > 0) if pd.notna(value) else pd.NA,
                }
            )

    long = pd.DataFrame(selected_rows)
    all_long = long.copy()
    carbon_long = long[long["Condition_Class"].isin(["C_AERO", "C_ANAEROBIC"])].copy()

    all_long.to_csv(OUT / "selected_models_table5_all_conditions_long.csv", index=False, encoding="utf-8-sig")
    carbon_long.to_csv(OUT / "selected_models_table5_carbon_sources_long.csv", index=False, encoding="utf-8-sig")

    for stem, df in [
        ("selected_models_table5_all_conditions", all_long),
        ("selected_models_table5_carbon_sources", carbon_long),
    ]:
        rates = df.pivot_table(index=["Condition", "Condition_Class", "Substrate_or_Exchange"], columns="Requested_Model", values="Growth_Rate", aggfunc="first").reset_index()
        binary = df.pivot_table(index=["Condition", "Condition_Class", "Substrate_or_Exchange"], columns="Requested_Model", values="Growth_Binary", aggfunc="first").reset_index()
        rates.to_csv(OUT / f"{stem}_growth_rates_wide.csv", index=False, encoding="utf-8-sig")
        binary.to_csv(OUT / f"{stem}_binary_truth_wide.csv", index=False, encoding="utf-8-sig")

    summary = (
        carbon_long.groupby(["Requested_Model", "SD01_Model_Name", "SD01_Strain", "Match_Status", "Condition_Class"])
        .agg(
            n_conditions=("Condition", "nunique"),
            growth_positive=("Growth_Binary", "sum"),
            growth_negative=("Growth_Binary", lambda x: int((x == 0).sum())),
            median_growth_rate=("Growth_Rate", "median"),
            max_growth_rate=("Growth_Rate", "max"),
        )
        .reset_index()
    )
    summary["growth_positive_fraction"] = summary["growth_positive"] / summary["n_conditions"]
    summary.to_csv(OUT / "selected_models_table5_carbon_source_summary.csv", index=False, encoding="utf-8-sig")

    readme = OUT / "README_selected_truth_tables.txt"
    readme.write_text(
        "Monk et al. 2013 PNAS SD01 Table 5 extraction\n"
        "==============================================\n\n"
        "Source workbook: D:\\pipeline\\data\\sd01.xlsx\n"
        "Source sheet: Table 5, 'In silico growth screens'.\n\n"
        "Important interpretation note:\n"
        "Table 5 contains in silico growth screens from Monk et al. 2013, not experimental Biolog measurements.\n"
        "The extracted Growth_Binary column is defined as Growth_Rate > 0.\n\n"
        "Requested model mapping:\n"
        + selected_df.to_string(index=False)
        + "\n\n"
        "The requested model name iECs1327 has no exact SD01 Table 1 match. The extracted row uses iECSF_1327 / E. coli SE15 as the closest 1327-suffix match and must be verified.\n"
        "iML1515 is not in SD01 because it was published later; the extracted truth uses the same MG1655 strain row as iJO1366.\n\n"
        "Files:\n"
        "- selected_models_table5_carbon_sources_long.csv: C_AERO and C_ANAEROBIC long-format truth table.\n"
        "- selected_models_table5_carbon_sources_binary_truth_wide.csv: 770 carbon conditions x requested models, binary truth.\n"
        "- selected_models_table5_carbon_sources_growth_rates_wide.csv: 770 carbon conditions x requested models, growth rates.\n"
        "- selected_models_table5_all_conditions_long.csv: all C/N/P/S aerobic/anaerobic conditions.\n"
        "- requested_model_to_sd01_strain_mapping.csv: model-to-strain mapping and match status.\n",
        encoding="utf-8",
    )

    print(f"Wrote outputs to {OUT}")
    print("Selected mapping:")
    print(selected_df.to_string(index=False))
    print("\nCarbon source summary:")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
