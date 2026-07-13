from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(r"D:\pipeline\analysis_v5")
IN_DIR = ROOT / "03.analysis" / "crispr_target_matching_inputs"
PHAGE = IN_DIR / "spacer_vs_phage.txt"
PLASMID = IN_DIR / "spacer_vs_plasmid.txt"
UNIQUE = IN_DIR / "unique_crispr_spacers_for_target_matching.csv"
OCC = IN_DIR / "crispr_spacer_occurrences_for_target_matching.csv"
COHORT = ROOT / "00.data" / "recalculated_annotations.csv"
OUT = ROOT / "03.analysis" / "crispr_target_matching_analysis"


BLAST_COLS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "stitle",
]


def parse_uid(qseqid: str) -> str:
    return str(qseqid).split("|", 1)[0]


def read_hits(path: Path, target_type: str, unique: pd.DataFrame) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame(columns=BLAST_COLS + ["Target_Type"])
    hits = pd.read_csv(path, sep="\t", names=BLAST_COLS, header=None)
    hits["Target_Type"] = target_type
    hits["Unique_Spacer_ID"] = hits["qseqid"].map(parse_uid)
    hits = hits.merge(
        unique[["Unique_Spacer_ID", "Sequence", "n_occurrences", "n_genomes", "n_phylogroups", "phylogroups"]],
        on="Unique_Spacer_ID",
        how="left",
    )
    hits["query_length"] = hits["Sequence"].astype(str).str.len()
    hits["qcov"] = hits["length"] / hits["query_length"]
    return hits


def best_per_spacer_target(hits: pd.DataFrame) -> pd.DataFrame:
    if hits.empty:
        return hits
    sort_cols = ["mismatch", "qcov", "pident", "bitscore", "evalue"]
    ascending = [True, False, False, False, True]
    return (
        hits.sort_values(sort_cols, ascending=ascending)
        .groupby(["Unique_Spacer_ID", "Target_Type"], as_index=False)
        .head(1)
        .reset_index(drop=True)
    )


def permutation_anova(values: np.ndarray, groups: np.ndarray, n_perm: int = 999, seed: int = 7) -> tuple[float, float, float, float]:
    mask = np.isfinite(values) & pd.notna(groups)
    values = values[mask]
    groups = groups[mask]
    overall = values.mean()
    ss_total = float(((values - overall) ** 2).sum())
    levels = pd.unique(groups)
    ss_between = 0.0
    ss_within = 0.0
    for g in levels:
        v = values[groups == g]
        ss_between += len(v) * float((v.mean() - overall) ** 2)
        ss_within += float(((v - v.mean()) ** 2).sum())
    df_between = len(levels) - 1
    df_within = len(values) - len(levels)
    ms_between = ss_between / df_between if df_between else np.nan
    ms_within = ss_within / df_within if df_within else np.nan
    f_obs = ms_between / ms_within if ms_within and ms_within > 0 else np.nan
    eta2 = ss_between / ss_total if ss_total > 0 else np.nan
    omega2 = (ss_between - df_between * ms_within) / (ss_total + ms_within) if ss_total > 0 and np.isfinite(ms_within) else np.nan
    rng = np.random.default_rng(seed)
    ge = 0
    for _ in range(n_perm):
        perm = rng.permutation(groups)
        ss_b = 0.0
        for g in levels:
            v = values[perm == g]
            ss_b += len(v) * float((v.mean() - overall) ** 2)
        ms_b = ss_b / df_between
        f_perm = ms_b / ms_within if ms_within and ms_within > 0 else np.nan
        if np.isfinite(f_perm) and f_perm >= f_obs:
            ge += 1
    p = (ge + 1) / (n_perm + 1)
    return float(f_obs), float(p), float(eta2), float(omega2)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    unique = pd.read_csv(UNIQUE)
    occ = pd.read_csv(OCC)
    phage = read_hits(PHAGE, "phage", unique)
    plasmid = read_hits(PLASMID, "plasmid", unique)
    all_hits = pd.concat([phage, plasmid], ignore_index=True)
    all_hits.to_csv(OUT / "crispr_target_hits_all_raw.csv", index=False, encoding="utf-8-sig")

    summary_rows = []
    filtered_sets = {}
    for label, mismatch_max in [("mismatch_le2", 2), ("mismatch_le3", 3)]:
        filt = all_hits[
            (all_hits["mismatch"] <= mismatch_max)
            & (all_hits["qcov"] >= 0.80)
            & (all_hits["pident"] >= 90.0)
            & (all_hits["evalue"] <= 1e-3)
        ].copy()
        best = best_per_spacer_target(filt)
        filtered_sets[label] = best
        best.to_csv(OUT / f"crispr_target_hits_best_{label}.csv", index=False, encoding="utf-8-sig")
        for target in ["phage", "plasmid"]:
            sub = best[best["Target_Type"] == target]
            summary_rows.append(
                {
                    "Filter": label,
                    "Target_Type": target,
                    "Raw_hit_rows": int((all_hits["Target_Type"] == target).sum()),
                    "Filtered_hit_rows": int((filt["Target_Type"] == target).sum()),
                    "Unique_targeted_spacers": int(sub["Unique_Spacer_ID"].nunique()),
                    "Spacer_occurrences_represented": int(sub.drop_duplicates("Unique_Spacer_ID")["n_occurrences"].sum()),
                    "Genomes_represented_by_targeted_spacers_sum": int(sub.drop_duplicates("Unique_Spacer_ID")["n_genomes"].sum()),
                    "Median_best_pident": float(sub["pident"].median()) if len(sub) else np.nan,
                    "Median_best_qcov": float(sub["qcov"].median()) if len(sub) else np.nan,
                }
            )

    target_summary = pd.DataFrame(summary_rows)
    target_summary.to_csv(OUT / "crispr_target_match_summary.csv", index=False, encoding="utf-8-sig")

    primary = filtered_sets["mismatch_le2"]
    phage_ids = set(primary.loc[primary["Target_Type"] == "phage", "Unique_Spacer_ID"])
    plasmid_ids = set(primary.loc[primary["Target_Type"] == "plasmid", "Unique_Spacer_ID"])
    targeted_ids = phage_ids | plasmid_ids

    occ2 = occ.copy()
    occ2["Is_Phage_Targeted"] = occ2["Sequence"].map(dict(zip(unique["Sequence"], unique["Unique_Spacer_ID"]))).isin(phage_ids)
    occ2["Is_Plasmid_Targeted"] = occ2["Sequence"].map(dict(zip(unique["Sequence"], unique["Unique_Spacer_ID"]))).isin(plasmid_ids)
    occ2["Is_Targeted"] = occ2["Is_Phage_Targeted"] | occ2["Is_Plasmid_Targeted"]
    occ2.to_csv(OUT / "crispr_spacer_occurrences_with_target_flags.csv", index=False, encoding="utf-8-sig")

    genome_rows = []
    for genome, gdf in occ2.groupby("Genome_ID"):
        uniq_total = gdf["Sequence"].nunique()
        total = len(gdf)
        phage_count = int(gdf["Is_Phage_Targeted"].sum())
        plasmid_count = int(gdf["Is_Plasmid_Targeted"].sum())
        targeted_count = int(gdf["Is_Targeted"].sum())
        genome_rows.append(
            {
                "Genome_ID": genome,
                "Phylogroup": gdf["Phylogroup"].dropna().iloc[0] if gdf["Phylogroup"].notna().any() else np.nan,
                "Spacer_occurrences": total,
                "Unique_spacers": uniq_total,
                "Phage_targeted_occurrences": phage_count,
                "Plasmid_targeted_occurrences": plasmid_count,
                "Any_targeted_occurrences": targeted_count,
                "Phage_targeted_fraction": phage_count / total if total else np.nan,
                "Plasmid_targeted_fraction": plasmid_count / total if total else np.nan,
                "Any_targeted_fraction": targeted_count / total if total else np.nan,
                "Has_phage_targeted_spacer": phage_count > 0,
                "Has_plasmid_targeted_spacer": plasmid_count > 0,
                "Has_any_targeted_spacer": targeted_count > 0,
            }
        )
    genome_df = pd.DataFrame(genome_rows)
    genome_df.to_csv(OUT / "crispr_target_matching_genome_summary.csv", index=False, encoding="utf-8-sig")

    group_summary = (
        genome_df.groupby("Phylogroup")
        .agg(
            n_genomes_with_spacers=("Genome_ID", "nunique"),
            median_spacer_occurrences=("Spacer_occurrences", "median"),
            median_phage_targeted_fraction=("Phage_targeted_fraction", "median"),
            median_plasmid_targeted_fraction=("Plasmid_targeted_fraction", "median"),
            median_any_targeted_fraction=("Any_targeted_fraction", "median"),
            genomes_with_phage_targets=("Has_phage_targeted_spacer", "sum"),
            genomes_with_plasmid_targets=("Has_plasmid_targeted_spacer", "sum"),
            genomes_with_any_targets=("Has_any_targeted_spacer", "sum"),
        )
        .reset_index()
    )
    for col in ["genomes_with_phage_targets", "genomes_with_plasmid_targets", "genomes_with_any_targets"]:
        group_summary[col] = group_summary[col].astype(int)
    group_summary["fraction_genomes_with_phage_targets"] = group_summary["genomes_with_phage_targets"] / group_summary["n_genomes_with_spacers"]
    group_summary["fraction_genomes_with_plasmid_targets"] = group_summary["genomes_with_plasmid_targets"] / group_summary["n_genomes_with_spacers"]
    group_summary["fraction_genomes_with_any_targets"] = group_summary["genomes_with_any_targets"] / group_summary["n_genomes_with_spacers"]
    group_summary.to_csv(OUT / "crispr_target_matching_phylogroup_summary.csv", index=False, encoding="utf-8-sig")

    stat_rows = []
    for metric in ["Phage_targeted_fraction", "Plasmid_targeted_fraction", "Any_targeted_fraction"]:
        f, p, eta2, omega2 = permutation_anova(genome_df[metric].to_numpy(float), genome_df["Phylogroup"].to_numpy(object))
        stat_rows.append({"Metric": metric, "F": f, "p_perm": p, "eta2": eta2, "omega2": omega2, "n_genomes_with_spacers": len(genome_df)})
    stats = pd.DataFrame(stat_rows)
    stats.to_csv(OUT / "crispr_target_matching_phylogroup_effect_sizes.csv", index=False, encoding="utf-8-sig")

    cohort = pd.read_csv(COHORT, usecols=["Genome_ID", "Phylogroup", "CRISPR_Arrays", "Cas_Operons"])
    all_genomes = cohort.merge(genome_df.drop(columns=["Phylogroup"]), on="Genome_ID", how="left")
    count_cols = [
        "Spacer_occurrences",
        "Unique_spacers",
        "Phage_targeted_occurrences",
        "Plasmid_targeted_occurrences",
        "Any_targeted_occurrences",
    ]
    bool_cols = ["Has_phage_targeted_spacer", "Has_plasmid_targeted_spacer", "Has_any_targeted_spacer"]
    for col in count_cols:
        all_genomes[col] = all_genomes[col].fillna(0).astype(int)
    for col in bool_cols:
        all_genomes[col] = all_genomes[col].fillna(False).astype(bool)
    # Fractions remain conditional on spacer-bearing genomes; no-spacer genomes are set to zero for burden analyses.
    for col in ["Phage_targeted_fraction", "Plasmid_targeted_fraction", "Any_targeted_fraction"]:
        all_genomes[col] = all_genomes[col].fillna(0.0)
    all_genomes["Has_any_spacer"] = all_genomes["Spacer_occurrences"] > 0
    all_genomes.to_csv(OUT / "crispr_target_matching_all_genome_summary.csv", index=False, encoding="utf-8-sig")

    all_group_summary = (
        all_genomes.groupby("Phylogroup")
        .agg(
            n_genomes=("Genome_ID", "nunique"),
            genomes_with_spacers=("Has_any_spacer", "sum"),
            median_spacer_occurrences_all=("Spacer_occurrences", "median"),
            median_phage_targets_all=("Phage_targeted_occurrences", "median"),
            median_plasmid_targets_all=("Plasmid_targeted_occurrences", "median"),
            median_any_targets_all=("Any_targeted_occurrences", "median"),
            genomes_with_phage_targets=("Has_phage_targeted_spacer", "sum"),
            genomes_with_plasmid_targets=("Has_plasmid_targeted_spacer", "sum"),
            genomes_with_any_targets=("Has_any_targeted_spacer", "sum"),
        )
        .reset_index()
    )
    for col in [
        "genomes_with_spacers",
        "genomes_with_phage_targets",
        "genomes_with_plasmid_targets",
        "genomes_with_any_targets",
    ]:
        all_group_summary[col] = all_group_summary[col].astype(int)
    all_group_summary["fraction_genomes_with_spacers"] = all_group_summary["genomes_with_spacers"] / all_group_summary["n_genomes"]
    all_group_summary["fraction_genomes_with_phage_targets"] = all_group_summary["genomes_with_phage_targets"] / all_group_summary["n_genomes"]
    all_group_summary["fraction_genomes_with_plasmid_targets"] = all_group_summary["genomes_with_plasmid_targets"] / all_group_summary["n_genomes"]
    all_group_summary["fraction_genomes_with_any_targets"] = all_group_summary["genomes_with_any_targets"] / all_group_summary["n_genomes"]
    all_group_summary.to_csv(OUT / "crispr_target_matching_phylogroup_summary_all_genomes.csv", index=False, encoding="utf-8-sig")

    all_stat_rows = []
    for metric in ["Phage_targeted_occurrences", "Plasmid_targeted_occurrences", "Any_targeted_occurrences"]:
        f, p, eta2, omega2 = permutation_anova(all_genomes[metric].to_numpy(float), all_genomes["Phylogroup"].to_numpy(object))
        all_stat_rows.append({"Metric": metric, "F": f, "p_perm": p, "eta2": eta2, "omega2": omega2, "n_genomes": len(all_genomes)})
    all_stats = pd.DataFrame(all_stat_rows)
    all_stats.to_csv(OUT / "crispr_target_matching_phylogroup_effect_sizes_all_genomes.csv", index=False, encoding="utf-8-sig")

    top_targets = (
        primary.groupby(["Target_Type", "sseqid", "stitle"])
        .agg(n_unique_spacers=("Unique_Spacer_ID", "nunique"), median_pident=("pident", "median"), median_qcov=("qcov", "median"))
        .reset_index()
        .sort_values(["Target_Type", "n_unique_spacers", "median_pident"], ascending=[True, False, False])
    )
    top_targets.to_csv(OUT / "crispr_target_matching_top_targets.csv", index=False, encoding="utf-8-sig")

    lines = [
        "CRISPR spacer target matching analysis",
        "======================================",
        "",
        "Primary high-confidence filter: mismatch <= 2, query coverage >= 80%, pident >= 90%, evalue <= 1e-3.",
        f"Raw phage hit rows: {(all_hits['Target_Type'] == 'phage').sum():,}",
        f"Raw plasmid hit rows: {(all_hits['Target_Type'] == 'plasmid').sum():,}",
    ]
    for row in target_summary[target_summary["Filter"] == "mismatch_le2"].itertuples(index=False):
        lines.append(
            f"{row.Target_Type}: {row.Unique_targeted_spacers:,} unique targeted spacers; "
            f"{row.Spacer_occurrences_represented:,} spacer occurrences represented; "
            f"median best pident {row.Median_best_pident:.2f}%, median qcov {row.Median_best_qcov:.2f}."
        )
    lines.extend(
        [
            "",
            "Genome-level target fractions are computed only among genomes with spacers in the spacer table.",
            f"Genomes with spacers: {len(genome_df):,}",
            f"All-cohort genomes analyzed for target burden: {len(all_genomes):,}",
            "",
            "Phylogroup effect sizes:",
        ]
    )
    for row in stats.itertuples(index=False):
        lines.append(f"{row.Metric}: F={row.F:.4f}, p={row.p_perm:.3f}, eta2={row.eta2:.4f}, omega2={row.omega2:.4f}")
    lines.append("")
    lines.append("All-cohort target-burden effect sizes (no-spacer genomes counted as zero target burden):")
    for row in all_stats.itertuples(index=False):
        lines.append(f"{row.Metric}: F={row.F:.4f}, p={row.p_perm:.3f}, eta2={row.eta2:.4f}, omega2={row.omega2:.4f}")
    (OUT / "crispr_target_matching_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("\n".join(lines))
    print(f"\nWrote outputs to {OUT}")


if __name__ == "__main__":
    main()
