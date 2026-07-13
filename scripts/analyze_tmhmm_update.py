from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(r"D:\pipeline\analysis_v5")
DATA = ROOT / "00.data"
OUT = ROOT / "03.analysis"
OUT.mkdir(exist_ok=True)


def eta_omega(df, metric, group_col="Phylogroup"):
    sub = df[[group_col, metric]].dropna()
    groups = list(sub.groupby(group_col))
    n = len(sub)
    k = len(groups)
    grand = sub[metric].mean()
    ssb = sum(len(g) * (g[metric].mean() - grand) ** 2 for _, g in groups)
    ssw = sum(((g[metric] - g[metric].mean()) ** 2).sum() for _, g in groups)
    sst = ssb + ssw
    msw = ssw / (n - k) if n > k else np.nan
    eta = ssb / sst if sst else np.nan
    omega = (ssb - (k - 1) * msw) / (sst + msw) if sst and not np.isnan(msw) else np.nan
    return n, k, eta, omega


def main():
    ann = pd.read_csv(DATA / "recalculated_annotations.csv")
    asm = pd.read_csv(DATA / "genome_assembly_and_cog_metrics.csv")
    tm = pd.read_csv(DATA / "all_genomes_tmhmm_positives_detailed.csv")

    counts = (
        tm.groupby("Genome_ID")
        .agg(
            TMHMM_Proteins_Detail=("Protein_ID", "nunique"),
            TMHMM_Total_Helices=("TMH_Count", "sum"),
            TMHMM_Mean_Helices=("TMH_Count", "mean"),
            TMHMM_Max_Helices=("TMH_Count", "max"),
        )
        .reset_index()
    )

    updated = ann.merge(counts, on="Genome_ID", how="left")
    updated["TMHMM_Proteins_Detail"] = updated["TMHMM_Proteins_Detail"].fillna(0).astype(int)
    updated["TMHMM_Total_Helices"] = updated["TMHMM_Total_Helices"].fillna(0).astype(int)
    updated["TMHMM_Mean_Helices"] = updated["TMHMM_Mean_Helices"].fillna(0)
    updated["TMHMM_Max_Helices"] = updated["TMHMM_Max_Helices"].fillna(0).astype(int)
    updated["Transmembrane_Proteins_Previous"] = updated["Transmembrane_Proteins"]
    updated["Transmembrane_Proteins"] = updated["TMHMM_Proteins_Detail"]
    updated["TMHMM_Proteins_per_1000_CDS"] = np.nan

    if "Total_CDS" in asm.columns:
        updated = updated.merge(asm[["Genome_ID", "Total_CDS"]], on="Genome_ID", how="left")
        updated["TMHMM_Proteins_per_1000_CDS"] = (
            updated["Transmembrane_Proteins"] / updated["Total_CDS"] * 1000
        )

    detail_ids = set(tm["Genome_ID"])
    missing = updated.loc[~updated["Genome_ID"].isin(detail_ids)].copy()
    zero = updated.loc[updated["Transmembrane_Proteins"].eq(0)].copy()
    changed = updated.loc[
        updated["Transmembrane_Proteins"].ne(updated["Transmembrane_Proteins_Previous"])
    ].copy()

    updated.to_csv(OUT / "recalculated_annotations_tmhmm_updated.csv", index=False)
    missing.to_csv(OUT / "tmhmm_missing_after_update.csv", index=False)
    zero.to_csv(OUT / "tmhmm_zero_after_update.csv", index=False)
    changed.to_csv(OUT / "tmhmm_changed_by_update.csv", index=False)

    metrics = [
        "Genome_Size_Mb",
        "GC_Percent",
        "Total_CDS",
        "GEM_Genes",
        "GEM_Reactions",
        "GEM_Metabolites",
        "CRISPR_Arrays",
        "Cas_Operons",
        "Secreted_Proteins",
        "Transmembrane_Proteins",
        "TMHMM_Proteins_per_1000_CDS",
    ]
    eff_df = asm[["Genome_ID", "Phylogroup", "Genome_Size_Mb", "GC_Percent", "Total_CDS"]].merge(
        updated[
            [
                "Genome_ID",
                "GEM_Genes",
                "GEM_Reactions",
                "GEM_Metabolites",
                "CRISPR_Arrays",
                "Cas_Operons",
                "Secreted_Proteins",
                "Transmembrane_Proteins",
                "TMHMM_Proteins_per_1000_CDS",
            ]
        ],
        on="Genome_ID",
        how="inner",
    )
    rows = []
    for metric in metrics:
        n, k, eta, omega = eta_omega(eff_df, metric)
        rows.append(
            {
                "metric": metric,
                "n": n,
                "groups": k,
                "eta_squared": eta,
                "omega_squared": omega,
            }
        )
    pd.DataFrame(rows).to_csv(OUT / "continuous_phylogroup_effect_sizes_eta2_tmhmm_updated.csv", index=False)

    group_summary = (
        updated.groupby("Phylogroup")
        .agg(
            n=("Genome_ID", "count"),
            tmhmm_median=("Transmembrane_Proteins", "median"),
            tmhmm_mean=("Transmembrane_Proteins", "mean"),
            tmhmm_min=("Transmembrane_Proteins", "min"),
            tmhmm_max=("Transmembrane_Proteins", "max"),
            tmhmm_per_1000_cds_median=("TMHMM_Proteins_per_1000_CDS", "median"),
            total_helices_median=("TMHMM_Total_Helices", "median"),
        )
        .reset_index()
        .sort_values("Phylogroup")
    )
    group_summary.to_csv(OUT / "tmhmm_group_summary_after_update.csv", index=False)

    print(f"annotation genomes: {len(ann)}")
    print(f"TMHMM detail rows: {len(tm)}")
    print(f"TMHMM detail genomes: {len(detail_ids)}")
    print(f"missing after update: {len(missing)}")
    print(f"zero after update: {len(zero)}")
    print(f"changed counts: {len(changed)}")
    print("\nTMHMM group summary")
    print(group_summary.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    print("\nEffect sizes")
    print(pd.read_csv(OUT / "continuous_phylogroup_effect_sizes_eta2_tmhmm_updated.csv").to_string(index=False, float_format=lambda x: f"{x:.4f}"))


if __name__ == "__main__":
    main()
