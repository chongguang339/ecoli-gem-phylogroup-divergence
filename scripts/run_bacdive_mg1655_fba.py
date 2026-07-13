from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
from cobra.io import read_sbml_model


ROOT = Path(__file__).resolve().parents[1]
TRUTH_XLSX = ROOT / "04.基准测试" / "K12 MG1655真实碳源数据.xlsx"
OUT = ROOT / "03.analysis" / "bacdive_mg1655_fba"
SUPP = ROOT / "05.supplementary_data" / "Data_S8_bacdive_mg1655_fba"

MODEL_SPECS = [
    {
        "model_id": "iML1515",
        "variant": "curated_reference",
        "path": Path(r"D:\pipeline\analysis_v4\data\model\iML1515.xml"),
    },
    {
        "model_id": "MG1655_gapseq_default",
        "variant": "gapseq_default",
        "path": ROOT / "mg1655_ascii_models" / "MG1655_gapseq_default.xml",
    },
    {
        "model_id": "MG1655_gapseq_optimized",
        "variant": "gapseq_optimized",
        "path": Path(r"D:\pipeline\analysis_v4\data\model\diamond_enterobacterales_fastdb\GCF_000005845.2_ASM584v2_genomic.xml"),
    },
]

BASE_METABOLITES = [
    "ca2_e",
    "cl_e",
    "co2_e",
    "cobalt2_e",
    "cu2_e",
    "fe2_e",
    "fe3_e",
    "h_e",
    "h2o_e",
    "k_e",
    "mg2_e",
    "mn2_e",
    "mobd_e",
    "na1_e",
    "nh4_e",
    "ni2_e",
    "o2_e",
    "pi_e",
    "so4_e",
    "tungs_e",
    "zn2_e",
]

BASE_MEDIUM_NAME_HINTS = {
    "ca2_e": ["ca2", "calcium"],
    "cl_e": ["chloride", "cl-"],
    "co2_e": ["co2", "carbon dioxide"],
    "cobalt2_e": ["cobalt", "co2+"],
    "cu2_e": ["cu2", "copper"],
    "fe2_e": ["fe2", "fe2+"],
    "fe3_e": ["fe3", "fe3+"],
    "h_e": ["h+", "proton"],
    "h2o_e": ["h2o", "water"],
    "k_e": ["k+", "potassium"],
    "mg2_e": ["mg", "magnesium"],
    "mn2_e": ["mn2", "manganese"],
    "mobd_e": ["molybdate", "mobd"],
    "na1_e": ["na+", "sodium"],
    "nh4_e": ["nh4", "ammonium", "nh3"],
    "ni2_e": ["ni2", "nickel"],
    "o2_e": ["o2", "oxygen"],
    "pi_e": ["phosphate"],
    "so4_e": ["sulfate"],
    "tungs_e": ["tungstate", "tungs"],
    "zn2_e": ["zn2", "zinc"],
}

NATIVE_CARBON_LIKE_EXCHANGES = {
    "EX_cpd00794_e0",  # trehalose-related
    "EX_cpd00027_e0",  # D-glucose
    "EX_cpd00051_e0",  # L-arginine
    "EX_cpd00154_e0",  # xylose
    "EX_cpd00122_e0",  # N-acetyl-D-glucosamine
    "EX_cpd00232_e0",  # Neu5Ac
    "EX_cpd00082_e0",  # D-fructose
    "EX_cpd00208_e0",  # lactate
    "EX_cpd00382_e0",  # melibiose/melitose
    "EX_cpd00396_e0",  # L-rhamnose
    "EX_cpd00751_e0",  # L-fucose
    "EX_cpd22614_e0",  # sulfoquinovose
}

EXCHANGE_ALIASES = {
    "rham__L_e": ["EX_rmn_e", "EX_cpd00396_e0"],
    "inositol_e": ["EX_inost_e", "EX_cpd00121_e0"],
    "sbt__D_e": ["EX_sbt__D_e", "EX_cpd00588_e0"],
}

CARBON_ASSIMILATION_SUBSET = {
    "D-Glucose",
    "L-Arabinose",
    "Melibiose",
    "L-Rhamnose",
    "Sorbitol (D-Glucitol)",
    "D-Mannitol",
    "Sucrose",
    "Myo-Inositol",
    "Citrate",
    "Amygdalin",
}

GROWTH_THRESHOLD = 1e-6
UPTAKE = 10.0


def normalize(text: object) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(text).lower())


def read_truth() -> pd.DataFrame:
    truth = pd.read_excel(TRUTH_XLSX)
    rename = {}
    for col in truth.columns:
        text = str(col)
        if "Carbon_Source" in text:
            rename[col] = "carbon_source"
        elif "BiGG_Metabolite" in text:
            rename[col] = "bigg_metabolite_id"
        elif "Exchange_Reaction" in text:
            rename[col] = "exchange_reaction"
        elif "Experimental_Growth" in text:
            rename[col] = "truth"
        elif "注" in text or "note" in text.lower():
            rename[col] = "note"
    truth = truth.rename(columns=rename)
    required = ["carbon_source", "bigg_metabolite_id", "exchange_reaction", "truth"]
    missing = [c for c in required if c not in truth.columns]
    if missing:
        raise ValueError(f"Missing required columns in BacDive workbook: {missing}")
    truth = truth[required + [c for c in ["note"] if c in truth.columns]].copy()
    truth["truth"] = pd.to_numeric(truth["truth"], errors="coerce").astype("Int64")
    truth["panel"] = np.where(
        truth["carbon_source"].isin(CARBON_ASSIMILATION_SUBSET),
        "carbon_assimilation_subset",
        "all_bacdive_panel_only",
    )
    truth["source"] = "BacDive 4902 / API 20E phenotype table; DOI 10.13145/bacdive4902.20260601.11"
    return truth


def build_exchange_lookup(model) -> dict[str, str]:
    lookup: dict[str, str] = {}
    for rxn in model.exchanges:
        candidates = {rxn.id, rxn.name}
        for met in rxn.metabolites:
            candidates.update(
                {
                    met.id,
                    met.name,
                    met.id.replace("_e", ""),
                    met.id.replace("[e]", ""),
                }
            )
        for candidate in candidates:
            key = normalize(candidate)
            if key:
                lookup.setdefault(key, rxn.id)
    return lookup


def resolve_exchange(row: pd.Series, lookup: dict[str, str], model) -> str | None:
    alias_key = str(row.get("bigg_metabolite_id", ""))
    for alias in EXCHANGE_ALIASES.get(alias_key, []):
        if alias in model.reactions:
            return alias

    candidates = [
        row.get("exchange_reaction", ""),
        row.get("bigg_metabolite_id", ""),
        str(row.get("bigg_metabolite_id", "")).replace("_e", ""),
        row.get("carbon_source", ""),
    ]
    for candidate in candidates:
        if isinstance(candidate, str) and candidate in model.reactions:
            return candidate
        key = normalize(candidate)
        if key in lookup:
            return lookup[key]
        if len(key) >= 4:
            for lookup_key, rxn_id in lookup.items():
                if len(lookup_key) >= 4 and (key in lookup_key or lookup_key in key):
                    return rxn_id
    return None


def resolve_by_hints(lookup: dict[str, str], hints: list[str]) -> str | None:
    for hint in hints:
        key = normalize(hint)
        if key in lookup:
            return lookup[key]
        if len(key) >= 4:
            for lookup_key, rxn_id in lookup.items():
                if len(lookup_key) >= 4 and (key in lookup_key or lookup_key in key):
                    return rxn_id
    return None


def set_medium(model, lookup: dict[str, str], carbon_exchange: str | None, oxygen: bool) -> tuple[int, int]:
    medium = {}
    resolved_base = 0
    for met in BASE_METABOLITES:
        if met == "o2_e" and not oxygen:
            continue
        hints = [met, f"EX_{met}", *BASE_MEDIUM_NAME_HINTS.get(met, [])]
        rxn_id = resolve_by_hints(lookup, hints)
        if rxn_id:
            medium[rxn_id] = UPTAKE if met != "o2_e" else 20.0
            resolved_base += 1
    if carbon_exchange and carbon_exchange in model.reactions:
        medium[carbon_exchange] = UPTAKE
    model.medium = medium
    return resolved_base, len(medium)


def is_native_noncarbon_exchange(reaction) -> bool:
    text = normalize(f"{reaction.id} {reaction.name} " + " ".join(m.name for m in reaction.metabolites))
    noncarbon_hints = [
        "phosphate",
        "sulfate",
        "h2o",
        "water",
        "oxygen",
        "o2",
        "h2s",
        "proton",
        "hplus",
        "ca2",
        "calcium",
        "cl",
        "chloride",
        "co2",
        "cobalt",
        "cu2",
        "copper",
        "fe2",
        "fe3",
        "mn2",
        "zinc",
        "zn2",
        "mg",
        "magnesium",
        "na",
        "sodium",
        "nh3",
        "nh4",
        "ammonium",
        "ni2",
        "nickel",
        "potassium",
    ]
    return any(normalize(hint) in text for hint in noncarbon_hints)


def set_native_carbon_stripped_medium(model, carbon_exchange: str | None, oxygen: bool) -> tuple[int, int]:
    medium = {}
    for reaction_id, uptake in model.medium.items():
        reaction = model.reactions.get_by_id(reaction_id)
        if not oxygen and ("oxygen" in reaction.name.lower() or reaction.id.lower().find("o2") >= 0):
            continue
        if reaction_id not in NATIVE_CARBON_LIKE_EXCHANGES:
            medium[reaction_id] = uptake
    if carbon_exchange and carbon_exchange in model.reactions:
        medium[carbon_exchange] = UPTAKE
    model.medium = medium
    return len(medium) - int(bool(carbon_exchange)), len(medium)


def fba_predict(
    model,
    lookup: dict[str, str],
    exchange_id: str | None,
    oxygen: bool,
    medium_strategy: str,
) -> tuple[float, int, str, int, int]:
    if not exchange_id:
        return 0.0, 0, "exchange_not_resolved", 0, 0
    try:
        with model:
            if medium_strategy == "native_carbon_stripped":
                resolved_base, medium_size = set_native_carbon_stripped_medium(model, exchange_id, oxygen)
            else:
                resolved_base, medium_size = set_medium(model, lookup, exchange_id, oxygen)
            solution = model.optimize()
            growth = 0.0 if solution.status != "optimal" else solution.objective_value
            if growth is None or math.isnan(growth) or growth < 0:
                growth = 0.0
            growth = float(growth)
            return growth, int(growth > GROWTH_THRESHOLD), "ok", resolved_base, medium_size
    except Exception as exc:
        return 0.0, 0, f"fba_error: {exc}", 0, 0


def safe_div(num: float, den: float) -> float:
    return float("nan") if den == 0 else num / den


def mcc(tp: int, tn: int, fp: int, fn: int) -> float:
    den = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float("nan") if den == 0 else ((tp * tn) - (fp * fn)) / den


def summarize_metrics(df: pd.DataFrame, panel: str) -> dict[str, object]:
    sub = df.copy() if panel == "all_bacdive_panel" else df[df["panel"] == panel].copy()
    sub = sub[sub["status"].eq("ok") & sub["truth"].notna()]
    truth = sub["truth"].astype(int)
    pred = sub["prediction"].astype(int)
    tp = int(((truth == 1) & (pred == 1)).sum())
    tn = int(((truth == 0) & (pred == 0)).sum())
    fp = int(((truth == 0) & (pred == 1)).sum())
    fn = int(((truth == 1) & (pred == 0)).sum())
    sens = safe_div(tp, tp + fn)
    spec = safe_div(tn, tn + fp)
    return {
        "panel": panel,
        "n_tested": int(len(sub)),
        "positives": int((truth == 1).sum()),
        "negatives": int((truth == 0).sum()),
        "TP": tp,
        "TN": tn,
        "FP": fp,
        "FN": fn,
        "Accuracy": safe_div(tp + tn, tp + tn + fp + fn),
        "Sensitivity": sens,
        "Specificity": spec,
        "Balanced_Accuracy": np.nanmean([sens, spec]),
        "MCC": mcc(tp, tn, fp, fn),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    SUPP.mkdir(parents=True, exist_ok=True)
    truth = read_truth()
    truth.to_csv(OUT / "bacdive_mg1655_truth_table.csv", index=False)

    predictions = []
    load_rows = []
    for spec in MODEL_SPECS:
        model_path = spec["path"]
        row = {
            "model_id": spec["model_id"],
            "variant": spec["variant"],
            "model_path": str(model_path),
            "loaded": False,
            "n_reactions": np.nan,
            "n_metabolites": np.nan,
            "n_genes": np.nan,
            "message": "",
        }
        if not model_path.exists():
            row["message"] = "model_file_not_found"
            load_rows.append(row)
            continue
        try:
            model = read_sbml_model(str(model_path))
            lookup = build_exchange_lookup(model)
            row.update(
                {
                    "loaded": True,
                    "n_reactions": len(model.reactions),
                    "n_metabolites": len(model.metabolites),
                    "n_genes": len(model.genes),
                }
            )
        except Exception as exc:
            row["message"] = f"load_error: {exc}"
            load_rows.append(row)
            continue
        load_rows.append(row)

        for oxygen_label, oxygen in [("aerobic", True), ("anaerobic_sensitivity", False)]:
            for _, truth_row in truth.iterrows():
                exchange_id = resolve_exchange(truth_row, lookup, model)
                medium_strategy = "native_carbon_stripped" if spec["variant"].startswith("gapseq") else "m9_like"
                growth, pred, status, resolved_base, medium_size = fba_predict(
                    model, lookup, exchange_id, oxygen, medium_strategy
                )
                predictions.append(
                    {
                        "model_id": spec["model_id"],
                        "variant": spec["variant"],
                        "condition": oxygen_label,
                        "carbon_source": truth_row["carbon_source"],
                        "bigg_metabolite_id": truth_row["bigg_metabolite_id"],
                        "expected_exchange_reaction": truth_row["exchange_reaction"],
                        "resolved_exchange_reaction": exchange_id or "",
                        "truth": int(truth_row["truth"]),
                        "prediction": pred,
                        "growth_rate": growth,
                        "panel": truth_row["panel"],
                        "status": status,
                        "resolved_base_medium_exchanges": resolved_base,
                        "medium_exchange_count": medium_size,
                        "medium_strategy": medium_strategy,
                        "model_path": str(model_path),
                    }
                )

    pred_df = pd.DataFrame(predictions)
    load_df = pd.DataFrame(load_rows)
    metrics = []
    for (model_id, variant, condition), sub in pred_df.groupby(["model_id", "variant", "condition"], dropna=False):
        for panel in ["carbon_assimilation_subset", "all_bacdive_panel"]:
            metrics.append(
                {
                    "model_id": model_id,
                    "variant": variant,
                    "condition": condition,
                    **summarize_metrics(sub, panel),
                }
            )
    metrics_df = pd.DataFrame(metrics)

    pred_df.to_csv(OUT / "bacdive_mg1655_fba_predictions_long.csv", index=False)
    metrics_df.to_csv(OUT / "bacdive_mg1655_fba_accuracy_metrics.csv", index=False)
    load_df.to_csv(OUT / "bacdive_mg1655_model_load_status.csv", index=False)

    for path in OUT.glob("*.csv"):
        target = SUPP / path.name
        target.write_bytes(path.read_bytes())

    print(metrics_df.to_string(index=False))
    print(f"Wrote outputs to {OUT}")


if __name__ == "__main__":
    main()
