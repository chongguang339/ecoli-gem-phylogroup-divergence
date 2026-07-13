from __future__ import annotations

from pathlib import Path
import re

import pandas as pd


ROOT = Path(r"D:\pipeline\analysis_v5")
SD01 = Path(r"D:\pipeline\data\sd01.xlsx")
TRUTH_DIR = ROOT / "03.analysis" / "monk2013_sd01_table5_truth"
OUT = ROOT / "03.analysis" / "monk2013_sd01_exchange_mapping"


def norm_text(s: object) -> str:
    s = str(s).strip().lower()
    s = s.replace("α", "alpha").replace("β", "beta")
    s = re.sub(r"^m_", "", s)
    s = s.replace(" exchange", "")
    s = re.sub(r"[^a-z0-9]+", "", s)
    return s


def clean_substrate(condition: str) -> str:
    s = str(condition)
    for prefix in ["C_AERO_", "C_ANAEROBIC_"]:
        if s.startswith(prefix):
            s = s[len(prefix):]
    s = re.sub(r"\s+[Ee]xchange$", "", s).strip()
    return s


def bigg_from_met(met: str) -> str:
    met = str(met)
    # BiGG encodes stereochemistry separators as double underscores.
    met = re.sub(r"-([A-Z](?:$|_))", r"__\1", met)
    met = met.replace("-", "__")
    return met


def seed_id(cross_ref: object) -> str:
    if pd.isna(cross_ref):
        return ""
    m = re.search(r"seed:(cpd\d+)", str(cross_ref))
    return m.group(1) if m else ""


def collect_reactions(xml_path: Path) -> set[str]:
    import xml.etree.ElementTree as ET

    ids = set()
    if not xml_path.exists():
        return ids
    for _, el in ET.iterparse(xml_path, events=("start",)):
        if el.tag.split("}", 1)[-1] == "reaction":
            rid = el.attrib.get("id")
            if rid:
                ids.add(rid)
    return ids


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    truth = pd.read_csv(TRUTH_DIR / "selected_models_table5_carbon_sources_long.csv")
    conds = (
        truth[["Condition", "Condition_Class", "Substrate_or_Exchange"]]
        .drop_duplicates()
        .sort_values(["Condition_Class", "Condition"])
        .reset_index(drop=True)
    )
    conds["Clean_Substrate"] = conds["Condition"].map(clean_substrate)
    conds["Substrate_Norm"] = conds["Clean_Substrate"].map(norm_text)

    t3 = pd.read_excel(SD01, sheet_name="Table 3", header=1)
    t3 = t3[["Abbreviation", "Metabolite", "Name", "Compartment", "Cross Reference"]].dropna(subset=["Name"]).copy()
    t3["Name_Norm"] = t3["Name"].map(norm_text)
    t3["Seed_ID"] = t3["Cross Reference"].map(seed_id)

    mappings = []
    for row in conds.itertuples(index=False):
        candidates = t3[t3["Name_Norm"].eq(row.Substrate_Norm)].copy()
        match_type = "exact_name"
        if candidates.empty:
            # Fallback: some Table 3 names include M_ prefixes or formula suffixes.
            candidates = t3[t3["Name_Norm"].str.contains(row.Substrate_Norm, regex=False, na=False)].copy()
            match_type = "contains_name" if not candidates.empty else "unmapped"
        if candidates.empty:
            mappings.append(
                {
                    "Condition": row.Condition,
                    "Condition_Class": row.Condition_Class,
                    "Clean_Substrate": row.Clean_Substrate,
                    "Match_Status": "unmapped",
                    "Match_Type": "",
                    "Table3_Abbreviation": "",
                    "Table3_Metabolite": "",
                    "Table3_Name": "",
                    "Table3_Compartment": "",
                    "Seed_ID": "",
                    "BiGG_Metabolite_ID": "",
                    "BiGG_Reaction_ID": "",
                    "COBRA_BiGG_Reaction_ID": "",
                    "SBML_Reaction_ID": "",
                    "Gapseq_Reaction_ID": "",
                    "COBRA_Gapseq_Reaction_ID": "",
                    "Candidate_Count": 0,
                }
            )
            continue
        candidates["comp_rank"] = candidates["Compartment"].map({"e": 0, "p": 1, "c": 2}).fillna(9)
        candidates = candidates.sort_values(["comp_rank", "Abbreviation"])
        best = candidates.iloc[0]
        abbr = str(best["Abbreviation"])
        met = str(best["Metabolite"])
        comp = str(best["Compartment"])
        if comp == "e" and abbr.endswith("_e"):
            bigg_met = bigg_from_met(abbr[:-2])
        else:
            bigg_met = bigg_from_met(met)
        bigg_rxn = f"EX_{bigg_met}_e"
        seed = best["Seed_ID"]
        mappings.append(
            {
                "Condition": row.Condition,
                "Condition_Class": row.Condition_Class,
                "Clean_Substrate": row.Clean_Substrate,
                "Match_Status": "mapped",
                "Match_Type": match_type,
                "Table3_Abbreviation": abbr,
                "Table3_Metabolite": met,
                "Table3_Name": best["Name"],
                "Table3_Compartment": comp,
                "Seed_ID": seed,
                "BiGG_Metabolite_ID": f"{bigg_met}_e",
                "BiGG_Reaction_ID": bigg_rxn,
                "COBRA_BiGG_Reaction_ID": bigg_rxn,
                "SBML_Reaction_ID": f"R_{bigg_rxn}",
                "Gapseq_Reaction_ID": f"R_EX_{seed}_e0" if seed else "",
                "COBRA_Gapseq_Reaction_ID": f"EX_{seed}_e0" if seed else "",
                "Candidate_Count": int(len(candidates)),
            }
        )

    mapping = pd.DataFrame(mappings)

    # Validate against currently available example models.
    xmls = {
        "iML1515_literature": Path(r"D:\pipeline\analysis_v4\data\model\iML1515.xml"),
        "MG1655_gapseq_optimized": Path(r"D:\pipeline\analysis_v4\data\model\diamond_enterobacterales_fastdb\GCF_000005845.2_ASM584v2_genomic.xml"),
        "MG1655_gapseq_default": Path(r"D:\pipeline\analysis_v4\data\model\gapseq默认参数模型\GCA_000005845.2_ASM584v2_genomic.xml"),
    }
    for label, path in xmls.items():
        rxns = collect_reactions(path)
        mapping[f"Present_In_{label}_BiGG_SBML"] = mapping["SBML_Reaction_ID"].isin(rxns)
        mapping[f"Present_In_{label}_Gapseq"] = mapping["Gapseq_Reaction_ID"].isin(rxns)

    mapping.to_csv(OUT / "sd01_770_carbon_conditions_exchange_id_mapping.csv", index=False, encoding="utf-8-sig")
    mapping[mapping["Match_Status"].eq("unmapped")].to_csv(OUT / "sd01_770_carbon_conditions_unmapped.csv", index=False, encoding="utf-8-sig")

    summary = {
        "n_conditions": len(mapping),
        "n_mapped": int(mapping["Match_Status"].eq("mapped").sum()),
        "n_unmapped": int(mapping["Match_Status"].eq("unmapped").sum()),
        "n_unique_bigg_exchange_ids": int(mapping.loc[mapping["Match_Status"].eq("mapped"), "BiGG_Reaction_ID"].nunique()),
        "n_unique_gapseq_exchange_ids": int(mapping.loc[mapping["Gapseq_Reaction_ID"].astype(str).str.len() > 0, "Gapseq_Reaction_ID"].nunique()),
    }
    pd.DataFrame([summary]).to_csv(OUT / "sd01_exchange_mapping_summary.csv", index=False, encoding="utf-8-sig")
    (OUT / "README_exchange_mapping.txt").write_text(
        "SD01 carbon-source exchange ID mapping\n"
        "=====================================\n\n"
        "Input: selected_models_table5_carbon_sources_long.csv and SD01 Table 3 metabolite dictionary.\n"
        "Output maps 770 C_AERO/C_ANAEROBIC conditions to candidate BiGG and gapseq exchange reaction IDs.\n\n"
        "BiGG_Reaction_ID is intended for COBRApy models after loading, e.g. EX_glc__D_e.\n"
        "SBML_Reaction_ID is the raw SBML XML reaction id, e.g. R_EX_glc__D_e.\n"
        "Gapseq_Reaction_ID is generated from Table 3 seed:cpd IDs when available, e.g. R_EX_cpd00027_e0.\n\n"
        "Mapping status should be checked before FBA. Unmapped or absent exchange reactions cannot be tested without manual curation.\n",
        encoding="utf-8",
    )
    print(summary)
    print(f"Wrote {OUT}")
    print("Top unmapped:")
    print(mapping[mapping["Match_Status"].eq("unmapped")].head(20).to_string(index=False))


if __name__ == "__main__":
    main()
