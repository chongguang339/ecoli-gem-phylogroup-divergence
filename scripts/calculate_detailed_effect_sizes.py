from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(r"D:\pipeline\analysis_v5")
ANALYSIS = ROOT / "03.analysis"


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
    signalp = pd.read_csv(ANALYSIS / "signalp_genome_counts.csv")
    crispr = pd.read_csv(ANALYSIS / "crispr_genome_counts.csv")
    tmhmm = pd.read_csv(ANALYSIS / "tmhmm_genome_counts.csv")
    merged = (
        signalp.merge(
            crispr[
                [
                    "Genome_ID",
                    "CRISPR_Arrays_Detail",
                    "Estimated_Spacers",
                    "Spacer_Sequences",
                    "Unique_Spacer_Sequences",
                    "Shared_Spacer_Fraction",
                ]
            ],
            on="Genome_ID",
        )
        .merge(
            tmhmm[
                [
                    "Genome_ID",
                    "TMHMM_Proteins",
                    "TMHMM_Helices",
                    "TMHMM_Proteins_per_1000_CDS",
                ]
            ],
            on="Genome_ID",
        )
    )
    metrics = [
        "SignalP_Total",
        "SignalP_per_1000_CDS",
        "CRISPR_Arrays_Detail",
        "Estimated_Spacers",
        "Spacer_Sequences",
        "Unique_Spacer_Sequences",
        "Shared_Spacer_Fraction",
        "TMHMM_Proteins",
        "TMHMM_Helices",
        "TMHMM_Proteins_per_1000_CDS",
    ]
    rows = []
    for metric in metrics:
        n, k, eta, omega = eta_omega(merged, metric)
        rows.append(
            {
                "metric": metric,
                "n": n,
                "groups": k,
                "eta_squared": eta,
                "omega_squared": omega,
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(ANALYSIS / "detail_metric_effect_sizes_eta2.csv", index=False)
    print(out.to_string(index=False, float_format=lambda x: f"{x:.4f}"))


if __name__ == "__main__":
    main()
