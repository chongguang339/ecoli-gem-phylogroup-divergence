from __future__ import annotations

from pathlib import Path
import re
import sys
import time

sys.path.insert(0, r"D:\pipeline\analysis_v5\.python_packages")

import pandas as pd
from cobra.io import read_sbml_model


ROOT = Path(r"D:\pipeline\analysis_v5")
SD01 = Path(r"D:\pipeline\data\sd01.xlsx")
TRUTH_DIR = ROOT / "03.analysis" / "monk2013_sd01_table5_truth"
MAP_DIR = ROOT / "03.analysis" / "monk2013_sd01_exchange_mapping"
MODEL_DIR = ROOT / "03.analysis" / "bigg_downloaded_models"
OUT = ROOT / "03.analysis" / "monk2013_sd01_cobra_fba"
SOLVER_TIMEOUT_SECONDS = 5


MODEL_FILES = {
    "iJO1366": MODEL_DIR / "iJO1366.xml",
    "iECs1327": MODEL_DIR / "iECSF_1327.xml",
    "ic_1306": MODEL_DIR / "ic_1306.xml",
    "iUTI89_1310": MODEL_DIR / "iUTI89_1310.xml",
    "iML1515": MODEL_DIR / "iML1515.xml",
}


def bigg_from_met(met: str) -> str:
    met = str(met)
    met = re.sub(r"-([A-Z](?:$|_))", r"__\1", met)
    met = met.replace("-", "__")
    return met


def reaction_from_formula(formula: str) -> str:
    if not isinstance(formula, str):
        return ""
    m = re.search(r"([A-Za-z0-9_\-\+]+)\[([a-z])\]", formula)
    if not m:
        return ""
    met = bigg_from_met(m.group(1))
    comp = m.group(2)
    return f"EX_{met}_{comp}"


def parse_table6_media() -> tuple[list[dict], dict[str, list[dict]]]:
    t6 = pd.read_excel(SD01, sheet_name="Table 6", header=None)
    base = []
    for _, row in t6.iloc[4:29, :5].iterrows():
        rxn = reaction_from_formula(row[2])
        if not rxn:
            continue
        base.append(
            {
                "Reaction_ID": rxn,
                "Reaction_Name": row[1],
                "Lower_Bound": float(row[3]),
                "Upper_Bound": float(row[4]),
                "Original_Abbreviation": row[0],
            }
        )
    aux = {}
    current = None
    for _, row in t6.iloc[4:, 7:13].iterrows():
        strain = row.iloc[0]
        if isinstance(strain, str) and strain.strip():
            current = strain.strip()
            aux.setdefault(current, [])
            continue
        if current and isinstance(row.iloc[1], str) and row.iloc[1].startswith("EX_"):
            rxn = reaction_from_formula(row.iloc[3])
            if rxn:
                aux[current].append(
                    {
                        "Reaction_ID": rxn,
                        "Reaction_Name": row.iloc[2],
                        "Lower_Bound": float(row.iloc[4]),
                        "Upper_Bound": float(row.iloc[5]),
                        "Original_Abbreviation": row.iloc[1],
                    }
                )
    return base, aux


def set_lb_if_present(model, rxn_id: str, lb: float, ub: float | None = None) -> bool:
    if rxn_id not in model.reactions:
        return False
    rxn = model.reactions.get_by_id(rxn_id)
    rxn.lower_bound = min(rxn.lower_bound, lb)
    if ub is not None:
        rxn.upper_bound = max(rxn.upper_bound, ub)
    return True


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    truth = pd.read_csv(TRUTH_DIR / "selected_models_table5_carbon_sources_long.csv")
    truth = truth.rename(columns={"Match_Status": "Model_Match_Status"})
    mapping = pd.read_csv(MAP_DIR / "sd01_770_carbon_conditions_exchange_id_mapping.csv")
    mapping = mapping.rename(columns={"Match_Status": "Exchange_Mapping_Status"})
    truth_map = truth.merge(
        mapping[["Condition", "COBRA_BiGG_Reaction_ID", "Exchange_Mapping_Status"]],
        on="Condition",
        how="left",
    )
    model_map = pd.read_csv(TRUTH_DIR / "requested_model_to_sd01_strain_mapping.csv")
    base_media, aux_media = parse_table6_media()
    pd.DataFrame(base_media).to_csv(OUT / "table6_base_m9_media_bigg_ids.csv", index=False, encoding="utf-8-sig")
    aux_rows = []
    for strain, vals in aux_media.items():
        for v in vals:
            aux_rows.append({"Strain": strain, **v})
    pd.DataFrame(aux_rows).to_csv(OUT / "table6_auxotrophy_supplements_bigg_ids.csv", index=False, encoding="utf-8-sig")

    results = []
    model_status = []
    t0 = time.time()
    for req_model, path in MODEL_FILES.items():
        map_row = model_map[model_map["Requested_Model"].eq(req_model)]
        sd01_strain = map_row["SD01_Strain"].iloc[0] if not map_row.empty else ""
        if not path.exists():
            model_status.append({"Requested_Model": req_model, "Model_Path": str(path), "Status": "missing"})
            continue
        print(f"Loading {req_model}: {path}")
        model = read_sbml_model(str(path))
        model.solver = "glpk"
        try:
            model.solver.configuration.timeout = SOLVER_TIMEOUT_SECONDS
        except Exception:
            pass
        exchange_ids = {r.id for r in model.exchanges}
        base_present = [m["Reaction_ID"] for m in base_media if m["Reaction_ID"] in model.reactions]
        aux_present = [m["Reaction_ID"] for m in aux_media.get(sd01_strain, []) if m["Reaction_ID"] in model.reactions]
        model_status.append(
            {
                "Requested_Model": req_model,
                "Model_Path": str(path),
                "Status": "loaded",
                "COBRA_Model_ID": model.id,
                "n_reactions": len(model.reactions),
                "n_metabolites": len(model.metabolites),
                "n_genes": len(model.genes),
                "n_exchange_reactions": len(exchange_ids),
                "Base_Media_Reactions_Present": len(base_present),
                "Auxotrophy_Reactions_Present": len(aux_present),
                "SD01_Strain": sd01_strain,
            }
        )
        # Initialize a closed uptake environment once per model; only oxygen and carbon-source
        # bounds are toggled inside the condition loop.
        for ex in model.exchanges:
            ex.lower_bound = 0.0
            if ex.upper_bound < 1000:
                ex.upper_bound = 1000.0
        for m in base_media:
            rid = m["Reaction_ID"]
            if rid in {"EX_glc__D_e", "EX_o2_e"}:
                continue
            set_lb_if_present(model, rid, m["Lower_Bound"], m["Upper_Bound"])
        for m in aux_media.get(sd01_strain, []):
            set_lb_if_present(model, m["Reaction_ID"], m["Lower_Bound"], m["Upper_Bound"])

        sub = truth_map[truth_map["Requested_Model"].eq(req_model)].copy()
        condition_counter = 0
        for row in sub.itertuples(index=False):
            condition_counter += 1
            rxn_id = row.COBRA_BiGG_Reaction_ID if isinstance(row.COBRA_BiGG_Reaction_ID, str) else ""

            present_base = sum(1 for m in base_media if m["Reaction_ID"] in model.reactions and m["Reaction_ID"] != "EX_glc__D_e")
            present_aux = sum(1 for m in aux_media.get(sd01_strain, []) if m["Reaction_ID"] in model.reactions)

            # Toggle oxygen per condition class.
            if "EX_o2_e" in model.reactions:
                model.reactions.get_by_id("EX_o2_e").lower_bound = -20.0 if row.Condition_Class == "C_AERO" else 0.0

            previous_lb = None
            if rxn_id and rxn_id in model.reactions:
                carbon_present = True
                carbon_rxn = model.reactions.get_by_id(rxn_id)
                previous_lb = carbon_rxn.lower_bound
                carbon_rxn.lower_bound = -20.0
                try:
                    growth = float(model.slim_optimize(error_value=0.0))
                    solver_status = str(getattr(model.solver, "status", ""))
                    status = solver_status if solver_status else "optimized"
                except Exception as e:
                    growth = 0.0
                    status = f"optimization_error:{type(e).__name__}"
                carbon_rxn.lower_bound = previous_lb
            else:
                carbon_present = False
                growth = 0.0
                status = "missing_carbon_exchange"

            results.append(
                {
                    "Requested_Model": req_model,
                    "SD01_Strain": sd01_strain,
                    "Condition": row.Condition,
                    "Condition_Class": row.Condition_Class,
                    "Substrate_or_Exchange": row.Substrate_or_Exchange,
                    "Truth_Growth_Rate": row.Growth_Rate,
                    "Truth_Binary": row.Growth_Binary,
                    "Mapped_COBRA_Exchange": rxn_id,
                    "Mapping_Status": row.Exchange_Mapping_Status,
                    "Carbon_Exchange_Present": carbon_present,
                    "Base_Media_Reactions_Present": present_base,
                    "Auxotrophy_Reactions_Present": present_aux,
                    "Predicted_Growth_Rate": growth,
                    "Predicted_Binary_1e_6": int(growth > 1e-6),
                    "Predicted_Binary_1e_4": int(growth > 1e-4),
                    "Simulation_Status": status,
                    "Model_Match_Status": row.Model_Match_Status,
                    "Exchange_Mapping_Status": row.Exchange_Mapping_Status,
                }
            )
            if condition_counter % 100 == 0:
                print(f"  {req_model}: {condition_counter}/{len(sub)} conditions", flush=True)
        print(f"Finished {req_model}", flush=True)
        pd.DataFrame(results).to_csv(OUT / "selected_models_770_carbon_fba_predictions_long.partial.csv", index=False, encoding="utf-8-sig")

    res = pd.DataFrame(results)
    status_df = pd.DataFrame(model_status)
    res.to_csv(OUT / "selected_models_770_carbon_fba_predictions_long.csv", index=False, encoding="utf-8-sig")
    status_df.to_csv(OUT / "cobra_model_load_and_media_status.csv", index=False, encoding="utf-8-sig")

    if not res.empty:
        pred_wide = res.pivot_table(index=["Condition", "Condition_Class", "Substrate_or_Exchange"], columns="Requested_Model", values="Predicted_Binary_1e_6", aggfunc="first").reset_index()
        rate_wide = res.pivot_table(index=["Condition", "Condition_Class", "Substrate_or_Exchange"], columns="Requested_Model", values="Predicted_Growth_Rate", aggfunc="first").reset_index()
        pred_wide.to_csv(OUT / "selected_models_770_carbon_fba_prediction_binary_wide.csv", index=False, encoding="utf-8-sig")
        rate_wide.to_csv(OUT / "selected_models_770_carbon_fba_prediction_growth_rates_wide.csv", index=False, encoding="utf-8-sig")

        metrics = []
        for model_name, g in res.groupby("Requested_Model"):
            for cond_class, h in g.groupby("Condition_Class"):
                truth_bin = h["Truth_Binary"].astype(int)
                pred = h["Predicted_Binary_1e_6"].astype(int)
                tp = int(((truth_bin == 1) & (pred == 1)).sum())
                tn = int(((truth_bin == 0) & (pred == 0)).sum())
                fp = int(((truth_bin == 0) & (pred == 1)).sum())
                fn = int(((truth_bin == 1) & (pred == 0)).sum())
                n = len(h)
                acc = (tp + tn) / n if n else 0
                sens = tp / (tp + fn) if tp + fn else 0
                spec = tn / (tn + fp) if tn + fp else 0
                bal = (sens + spec) / 2
                denom = ((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)) ** 0.5
                mcc = ((tp * tn - fp * fn) / denom) if denom else 0
                metrics.append(
                    {
                        "Requested_Model": model_name,
                        "Condition_Class": cond_class,
                        "n_conditions": n,
                        "mapped_conditions": int(h["Mapped_COBRA_Exchange"].astype(str).str.len().gt(0).sum()),
                        "carbon_exchange_present": int(h["Carbon_Exchange_Present"].sum()),
                        "TP": tp,
                        "TN": tn,
                        "FP": fp,
                        "FN": fn,
                        "Accuracy": acc,
                        "Sensitivity": sens,
                        "Specificity": spec,
                        "Balanced_Accuracy": bal,
                        "MCC": mcc,
                    }
                )
        metrics_df = pd.DataFrame(metrics)
        metrics_df.to_csv(OUT / "selected_models_770_carbon_fba_accuracy_metrics.csv", index=False, encoding="utf-8-sig")
    elapsed = time.time() - t0
    print(f"Wrote outputs to {OUT}")
    print(status_df.to_string(index=False))
    if not res.empty:
        print(metrics_df.to_string(index=False))
    print(f"Elapsed seconds: {elapsed:.1f}")


if __name__ == "__main__":
    main()
