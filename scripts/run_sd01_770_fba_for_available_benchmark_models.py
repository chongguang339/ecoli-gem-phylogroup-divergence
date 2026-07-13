from __future__ import annotations

from pathlib import Path
import math
import re
import sys
import time

sys.path.insert(0, r"D:\pipeline\analysis_v5\.python_packages")

import pandas as pd
from cobra.io import read_sbml_model


ROOT = Path(r"D:\pipeline\analysis_v5")
V4_MODEL = Path(r"D:\pipeline\analysis_v4\data\model")
SD01 = Path(r"D:\pipeline\data\sd01.xlsx")
MAP_FILE = ROOT / "03.analysis" / "monk2013_sd01_exchange_mapping" / "sd01_770_carbon_conditions_exchange_id_mapping.csv"
OUT = ROOT / "03.analysis" / "monk2013_sd01_available_benchmark_fba"
MODLE_COMPARE = ROOT / "04.基准测试" / "modle_compare"

SOLVER_TIMEOUT_SECONDS = 5
GROWTH_THRESHOLD = 1e-6

BENCHMARK_MODELS = {
    "iB21_1397": {
        "sd01_strain": "Escherichia coli BL21(DE3) AM946981",
        "literature": MODLE_COMPARE / "00.data" / "iB21_1397" / "iB21_1397.xml",
        "gapseq_default": MODLE_COMPARE / "iB21_1397" / "results_default" / "GCA_000022665.2_ASM2266v1_genomic.xml",
        "gapseq_optimized": MODLE_COMPARE / "iB21_1397" / "results_optimized" / "GCA_000022665.2_ASM2266v1_genomic.xml",
    },
    "iEC1364_W": {
        "sd01_strain": "Escherichia coli str. K-12 substr. W3110",
        "literature": MODLE_COMPARE / "00.data" / "iEC1364_W" / "iEC1364_W.xml",
        "gapseq_default": MODLE_COMPARE / "iEC1364_W" / "results_default" / "GCA_000299415.1_ASM29941v1_genomic.xml",
        "gapseq_optimized": MODLE_COMPARE / "iEC1364_W" / "results_optimized" / "GCA_000299415.1_ASM29941v1_genomic.xml",
    },
    "iECDH10B_1368": {
        "sd01_strain": "Escherichia coli str. K-12 substr. DH10B",
        "literature": MODLE_COMPARE / "00.data" / "iECDH10B_1368" / "iECDH10B_1368.xml",
        "gapseq_default": MODLE_COMPARE / "iECDH10B_1368" / "results_default" / "GCF_000019425.1_ASM1942v1_genomic.xml",
        "gapseq_optimized": MODLE_COMPARE / "iECDH10B_1368" / "results_optimized" / "GCF_000019425.1_ASM1942v1_genomic.xml",
    },
    "iJN1463": {
        "sd01_strain": "Escherichia coli Nissle 1917",
        "literature": MODLE_COMPARE / "00.data" / "iJN1463" / "iJN1463.xml",
        "gapseq_default": MODLE_COMPARE / "iJN1463" / "results_default" / "cds_from_genomic.xml",
        "gapseq_optimized": MODLE_COMPARE / "iJN1463" / "results_optimized" / "cds_from_genomic.xml",
    },
    "iML1515": {
        "sd01_strain": "Escherichia coli str. K-12 substr. MG1655",
        "literature": V4_MODEL / "iML1515.xml",
        "gapseq_default": V4_MODEL / "gapseq默认参数模型" / "GCA_000005845.2_ASM584v2_genomic.xml",
        "gapseq_optimized": V4_MODEL / "diamond_enterobacterales_fastdb" / "GCF_000005845.2_ASM584v2_genomic.xml",
    },
}

BIGG_BASE = {
    "EX_ca2_e",
    "EX_cl_e",
    "EX_co2_e",
    "EX_cobalt2_e",
    "EX_cu2_e",
    "EX_fe2_e",
    "EX_fe3_e",
    "EX_h_e",
    "EX_h2o_e",
    "EX_k_e",
    "EX_mg2_e",
    "EX_mn2_e",
    "EX_mobd_e",
    "EX_na1_e",
    "EX_ni2_e",
    "EX_pi_e",
    "EX_sel_e",
    "EX_slnt_e",
    "EX_so4_e",
    "EX_tungs_e",
    "EX_zn2_e",
}

# ModelSEED exchange ids used by gapseq for a minimal M9-like background.
GAPSEQ_BASE = {
    "EX_cpd00001_e0",  # H2O
    "EX_cpd00009_e0",  # phosphate
    "EX_cpd00011_e0",  # CO2
    "EX_cpd00013_e0",  # ammonia
    "EX_cpd00030_e0",  # Mn2+
    "EX_cpd00034_e0",  # Zn2+
    "EX_cpd00048_e0",  # sulfate
    "EX_cpd00058_e0",  # Cu2+
    "EX_cpd00063_e0",  # Ca2+
    "EX_cpd00067_e0",  # H+
    "EX_cpd00099_e0",  # chloride
    "EX_cpd00149_e0",  # Co2+
    "EX_cpd00205_e0",  # K+
    "EX_cpd00244_e0",  # Ni2+
    "EX_cpd00254_e0",  # Mg2+
    "EX_cpd00971_e0",  # Na+
    "EX_cpd10515_e0",  # Fe2+
    "EX_cpd10516_e0",  # Fe3+
    "EX_cpd11574_e0",  # molybdate
}


def sanitize_id(text: str) -> str:
    clean = re.sub(r"\W+", "_", text).strip("_")
    if not clean or clean[0].isdigit():
        clean = "model_" + clean
    return clean


def read_model_with_repair(path: Path, model_id: str):
    try:
        return read_sbml_model(str(path)), "loaded"
    except Exception as first_error:
        text = path.read_text(encoding="utf-8", errors="ignore")
        if re.search(r"<model\s+(?![^>]*\bid=)", text):
            repaired = re.sub(r"<model\s+", f'<model id="{sanitize_id(model_id)}" ', text, count=1)
            repaired_path = OUT / "tmp_repaired_sbml" / f"{sanitize_id(model_id)}__{path.name}"
            repaired_path.parent.mkdir(parents=True, exist_ok=True)
            repaired_path.write_text(repaired, encoding="utf-8")
            return read_sbml_model(str(repaired_path)), f"loaded_after_adding_model_id:{type(first_error).__name__}"
        raise


def load_truth() -> pd.DataFrame:
    raw = pd.read_excel(SD01, sheet_name="Table 5", header=None)
    header = raw.iloc[2].tolist()
    data = raw.iloc[3:].copy()
    data.columns = header
    data = data.rename(columns={"Strain:": "SD01_Strain"})
    condition_cols = [c for c in data.columns if isinstance(c, str) and (c.startswith("C_AERO_") or c.startswith("C_ANAEROBIC_"))]
    rows = []
    for model_id, spec in BENCHMARK_MODELS.items():
        match = data[data["SD01_Strain"].eq(spec["sd01_strain"])]
        if match.empty:
            for cond in condition_cols:
                rows.append({"Benchmark_Model": model_id, "SD01_Strain": spec["sd01_strain"], "Condition": cond, "Truth_Growth_Rate": math.nan, "Truth_Binary": math.nan, "Truth_Status": "sd01_strain_not_found"})
            continue
        values = pd.to_numeric(match.iloc[0][condition_cols], errors="coerce")
        for cond, val in values.items():
            rows.append({"Benchmark_Model": model_id, "SD01_Strain": spec["sd01_strain"], "Condition": cond, "Truth_Growth_Rate": val, "Truth_Binary": int((0 if pd.isna(val) else val) > GROWTH_THRESHOLD), "Truth_Status": "extracted_from_sd01_table5_in_silico"})
    return pd.DataFrame(rows)


def configure_medium(model, variant: str) -> None:
    for rxn in model.exchanges:
        rxn.lower_bound = 0.0
        if rxn.upper_bound < 1000:
            rxn.upper_bound = 1000.0
    base = GAPSEQ_BASE if variant.startswith("gapseq") else BIGG_BASE
    for rid in base:
        if rid in model.reactions:
            model.reactions.get_by_id(rid).lower_bound = -1000.0


def set_oxygen(model, variant: str, aerobic: bool) -> None:
    rid = "EX_cpd00007_e0" if variant.startswith("gapseq") else "EX_o2_e"
    if rid in model.reactions:
        model.reactions.get_by_id(rid).lower_bound = -20.0 if aerobic else 0.0


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    truth = load_truth()
    mapping = pd.read_csv(MAP_FILE)
    mapping = mapping[["Condition", "Condition_Class", "Clean_Substrate", "COBRA_BiGG_Reaction_ID", "COBRA_Gapseq_Reaction_ID", "Match_Status"]]
    task = truth.merge(mapping, on="Condition", how="left")
    task.to_csv(OUT / "available_benchmark_sd01_770_truth_with_exchange_mapping.csv", index=False, encoding="utf-8-sig")

    rows = []
    status_rows = []
    t0 = time.time()
    for model_id, spec in BENCHMARK_MODELS.items():
        for variant in ["literature", "gapseq_default", "gapseq_optimized"]:
            path = spec[variant]
            if not path.exists():
                status_rows.append({"Benchmark_Model": model_id, "Variant": variant, "Model_Path": str(path), "Load_Status": "missing"})
                continue
            print(f"Loading {model_id} {variant}: {path}", flush=True)
            try:
                model, load_status = read_model_with_repair(path, f"{model_id}_{variant}")
                model.solver = "glpk"
                try:
                    model.solver.configuration.timeout = SOLVER_TIMEOUT_SECONDS
                except Exception:
                    pass
            except Exception as exc:
                status_rows.append({"Benchmark_Model": model_id, "Variant": variant, "Model_Path": str(path), "Load_Status": f"load_error:{type(exc).__name__}"})
                continue
            configure_medium(model, variant)
            status_rows.append({
                "Benchmark_Model": model_id,
                "Variant": variant,
                "Model_Path": str(path),
                "Load_Status": load_status,
                "COBRA_Model_ID": model.id,
                "n_reactions": len(model.reactions),
                "n_metabolites": len(model.metabolites),
                "n_genes": len(model.genes),
                "n_exchange_reactions": len(model.exchanges),
            })
            sub = task[task["Benchmark_Model"].eq(model_id)]
            rxn_col = "COBRA_Gapseq_Reaction_ID" if variant.startswith("gapseq") else "COBRA_BiGG_Reaction_ID"
            for idx, row in enumerate(sub.itertuples(index=False), start=1):
                rid = getattr(row, rxn_col)
                rid = rid if isinstance(rid, str) else ""
                set_oxygen(model, variant, row.Condition_Class == "C_AERO")
                old_lb = None
                if rid and rid in model.reactions:
                    carbon_present = True
                    rxn = model.reactions.get_by_id(rid)
                    old_lb = rxn.lower_bound
                    rxn.lower_bound = -20.0
                    try:
                        growth = float(model.slim_optimize(error_value=0.0))
                        sim_status = str(getattr(model.solver, "status", "")) or "optimized"
                    except Exception as exc:
                        growth = 0.0
                        sim_status = f"optimization_error:{type(exc).__name__}"
                    rxn.lower_bound = old_lb
                else:
                    carbon_present = False
                    growth = 0.0
                    sim_status = "missing_carbon_exchange"
                rows.append({
                    "Benchmark_Model": model_id,
                    "Variant": variant,
                    "Condition": row.Condition,
                    "Condition_Class": row.Condition_Class,
                    "Clean_Substrate": row.Clean_Substrate,
                    "Truth_Growth_Rate": row.Truth_Growth_Rate,
                    "Truth_Binary": row.Truth_Binary,
                    "Truth_Status": row.Truth_Status,
                    "Mapped_Exchange": rid,
                    "Carbon_Exchange_Present": carbon_present,
                    "Predicted_Growth_Rate": growth,
                    "Predicted_Binary_1e_6": int(growth > GROWTH_THRESHOLD),
                    "Simulation_Status": sim_status,
                })
                if idx % 100 == 0:
                    print(f"  {model_id} {variant}: {idx}/{len(sub)}", flush=True)
            pd.DataFrame(rows).to_csv(OUT / "available_benchmark_770_fba_predictions_long.partial.csv", index=False, encoding="utf-8-sig")

    res = pd.DataFrame(rows)
    status = pd.DataFrame(status_rows)
    res.to_csv(OUT / "available_benchmark_770_fba_predictions_long.csv", index=False, encoding="utf-8-sig")
    status.to_csv(OUT / "available_benchmark_model_load_status.csv", index=False, encoding="utf-8-sig")
    if not res.empty:
        wide = res.pivot_table(index=["Benchmark_Model", "Variant"], columns=["Condition"], values="Predicted_Binary_1e_6", aggfunc="first")
        wide.to_csv(OUT / "available_benchmark_770_fba_prediction_binary_matrix.csv", encoding="utf-8-sig")
        metrics = []
        for (model_id, variant, cond_class), g in res.dropna(subset=["Truth_Binary"]).groupby(["Benchmark_Model", "Variant", "Condition_Class"]):
            truth_bin = g["Truth_Binary"].astype(int)
            pred = g["Predicted_Binary_1e_6"].astype(int)
            tp = int(((truth_bin == 1) & (pred == 1)).sum())
            tn = int(((truth_bin == 0) & (pred == 0)).sum())
            fp = int(((truth_bin == 0) & (pred == 1)).sum())
            fn = int(((truth_bin == 1) & (pred == 0)).sum())
            n = len(g)
            sens = tp / (tp + fn) if tp + fn else 0
            spec = tn / (tn + fp) if tn + fp else 0
            denom = ((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)) ** 0.5
            metrics.append({
                "Benchmark_Model": model_id,
                "Variant": variant,
                "Condition_Class": cond_class,
                "n_conditions": n,
                "carbon_exchange_present": int(g["Carbon_Exchange_Present"].sum()),
                "TP": tp,
                "TN": tn,
                "FP": fp,
                "FN": fn,
                "Accuracy": (tp + tn) / n if n else 0,
                "Sensitivity": sens,
                "Specificity": spec,
                "Balanced_Accuracy": (sens + spec) / 2,
                "MCC": ((tp * tn - fp * fn) / denom) if denom else 0,
            })
        metrics_df = pd.DataFrame(metrics)
        metrics_df.to_csv(OUT / "available_benchmark_770_fba_accuracy_metrics.csv", index=False, encoding="utf-8-sig")
        print(metrics_df.to_string(index=False))
    print(status.to_string(index=False))
    print(f"Wrote outputs to {OUT}")
    print(f"Elapsed seconds: {time.time() - t0:.1f}")


if __name__ == "__main__":
    main()
