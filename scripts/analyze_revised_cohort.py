from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "00.data"
OUT = ROOT / "03.analysis"
OUT.mkdir(exist_ok=True)

PHYLO_ORDER = ["A", "B2", "B1", "D", "E", "C", "F", "G", "U"]
RNG = np.random.default_rng(20260622)


def gammln(xx: float) -> float:
    return math.lgamma(xx)


def gammq(a: float, x: float) -> float:
    """Regularized upper incomplete gamma Q(a, x), adapted from Numerical Recipes."""
    if x < 0 or a <= 0:
        return float("nan")
    if x == 0:
        return 1.0
    eps = 3.0e-12
    itmax = 1000
    fpmin = 1.0e-300
    gln = gammln(a)
    if x < a + 1.0:
        ap = a
        total = 1.0 / a
        delta = total
        for _ in range(itmax):
            ap += 1.0
            delta *= x / ap
            total += delta
            if abs(delta) < abs(total) * eps:
                p = total * math.exp(-x + a * math.log(x) - gln)
                return max(0.0, min(1.0, 1.0 - p))
        return float("nan")
    b = x + 1.0 - a
    c = 1.0 / fpmin
    d = 1.0 / b
    h = d
    for i in range(1, itmax + 1):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < fpmin:
            d = fpmin
        c = b + an / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            return max(0.0, min(1.0, math.exp(-x + a * math.log(x) - gln) * h))
    return float("nan")


def chisq_pvalue(stat: float, dof: int) -> float:
    return gammq(dof / 2.0, stat / 2.0)


def bh_adjust(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    n = len(p)
    adj = np.empty(n, dtype=float)
    running = 1.0
    for i in range(n - 1, -1, -1):
        running = min(running, ranked[i] * n / (i + 1))
        adj[order[i]] = running
    return np.clip(adj, 0.0, 1.0)


def chi_square_table(obs: np.ndarray) -> tuple[float, int, float]:
    obs = np.asarray(obs, dtype=float)
    row = obs.sum(axis=1, keepdims=True)
    col = obs.sum(axis=0, keepdims=True)
    total = obs.sum()
    if total <= 0:
        return float("nan"), 0, float("nan")
    expected = row @ col / total
    valid = expected > 0
    stat = float(((obs[valid] - expected[valid]) ** 2 / expected[valid]).sum())
    dof = (obs.shape[0] - 1) * (obs.shape[1] - 1)
    return stat, dof, chisq_pvalue(stat, dof)


def parse_signalp(path: Path) -> tuple[pd.DataFrame, int]:
    rows = []
    repaired = 0
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        if header[:5] != ["Genome_ID", "Protein_ID", "Prediction_Type", "SP_Probability", "CS_Position"]:
            raise ValueError(f"Unexpected SignalP header: {header}")
        for row in reader:
            if not row:
                continue
            if len(row) != 5:
                repaired += 1
            genome = row[0]
            pred = row[-3]
            try:
                sp_prob = float(row[-2])
                cs_pos = float(row[-1])
            except ValueError:
                continue
            protein = ",".join(row[1:-3])
            rows.append((genome, protein, pred, sp_prob, cs_pos))
    return pd.DataFrame(rows, columns=header), repaired


def pca_svd(x: np.ndarray, n_components: int = 2) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray(x, dtype=float)
    x = x - x.mean(axis=0)
    sd = x.std(axis=0, ddof=1)
    sd[sd == 0] = 1.0
    z = x / sd
    u, s, _ = np.linalg.svd(z, full_matrices=False)
    scores = u[:, :n_components] * s[:n_components]
    explained = (s**2) / np.sum(s**2)
    return scores, explained[:n_components]


def jaccard_distance(binary: np.ndarray) -> np.ndarray:
    x = binary.astype(bool)
    intersection = x.astype(np.uint16) @ x.astype(np.uint16).T
    row_sum = x.sum(axis=1).astype(np.uint16)
    union = row_sum[:, None] + row_sum[None, :] - intersection
    dist = np.ones_like(intersection, dtype=float)
    np.divide(union - intersection, union, out=dist, where=union > 0)
    np.fill_diagonal(dist, 0.0)
    return dist


def euclidean_distance(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    sq = np.sum(x * x, axis=1)
    d2 = sq[:, None] + sq[None, :] - 2.0 * (x @ x.T)
    np.maximum(d2, 0.0, out=d2)
    return np.sqrt(d2)


def permanova_from_dist(d: np.ndarray, labels: np.ndarray, permutations: int = 999) -> dict:
    labels = np.asarray(labels)
    groups = pd.unique(labels)
    n = len(labels)
    g = len(groups)
    d2 = d**2
    total_ss = float(d2.sum() / n)

    def stats_for(lab):
        within = 0.0
        for group in groups:
            idx = np.where(lab == group)[0]
            if len(idx) > 1:
                within += float(d2[np.ix_(idx, idx)].sum() / len(idx))
        between = total_ss - within
        f_stat = (between / (g - 1)) / (within / (n - g))
        return between / total_ss, f_stat

    r2, f_stat = stats_for(labels)
    if permutations <= 0:
        return {
            "n": n,
            "groups": g,
            "permutations": 0,
            "F": f_stat,
            "R2": r2,
            "p": np.nan,
        }
    ge = 0
    for _ in range(permutations):
        _, pf = stats_for(RNG.permutation(labels))
        if pf >= f_stat - 1e-12:
            ge += 1
    return {
        "n": n,
        "groups": g,
        "permutations": permutations,
        "F": f_stat,
        "R2": r2,
        "p": (ge + 1) / (permutations + 1),
    }


def pcoa_coordinates(d: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n = d.shape[0]
    d2 = d**2
    j = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * j @ d2 @ j
    eigvals, eigvecs = np.linalg.eigh(b)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    pos = eigvals > 1e-10
    eigvals = eigvals[pos]
    eigvecs = eigvecs[:, pos]
    coords = eigvecs * np.sqrt(eigvals)
    return coords, eigvals


def permdisp_from_coords(coords: np.ndarray, labels: np.ndarray, permutations: int = 999) -> dict:
    labels = np.asarray(labels)
    groups = pd.unique(labels)

    def distances_to_centroid(lab):
        out = np.zeros(len(lab), dtype=float)
        for group in groups:
            idx = np.where(lab == group)[0]
            centroid = coords[idx].mean(axis=0)
            out[idx] = np.linalg.norm(coords[idx] - centroid, axis=1)
        return out

    def anova_f(values, lab):
        grand = values.mean()
        ss_between = 0.0
        ss_within = 0.0
        for group in groups:
            v = values[lab == group]
            ss_between += len(v) * (v.mean() - grand) ** 2
            ss_within += ((v - v.mean()) ** 2).sum()
        df_between = len(groups) - 1
        df_within = len(values) - len(groups)
        return float((ss_between / df_between) / (ss_within / df_within))

    values = distances_to_centroid(labels)
    f_stat = anova_f(values, labels)
    ge = 0
    for _ in range(permutations):
        perm = RNG.permutation(labels)
        if anova_f(values, perm) >= f_stat - 1e-12:
            ge += 1
    return {
        "n": len(labels),
        "groups": len(groups),
        "permutations": permutations,
        "F": f_stat,
        "p": (ge + 1) / (permutations + 1),
    }, values


def summarize_metric(df: pd.DataFrame, metric: str) -> pd.DataFrame:
    return (
        df.groupby("Phylogroup")[metric]
        .agg(count="count", mean="mean", median="median", std="std", min="min", max="max")
        .reindex(PHYLO_ORDER)
    )


def permutation_anova(values: np.ndarray, labels: np.ndarray, permutations: int = 999) -> tuple[float, float]:
    labels = np.asarray(labels)
    groups = pd.unique(labels)

    def f_stat(v, lab):
        grand = v.mean()
        ssb = 0.0
        ssw = 0.0
        for g in groups:
            vals = v[lab == g]
            ssb += len(vals) * (vals.mean() - grand) ** 2
            ssw += ((vals - vals.mean()) ** 2).sum()
        return float((ssb / (len(groups) - 1)) / (ssw / (len(v) - len(groups))))

    observed = f_stat(values, labels)
    ge = 0
    for _ in range(permutations):
        if f_stat(values, RNG.permutation(labels)) >= observed - 1e-12:
            ge += 1
    return observed, (ge + 1) / (permutations + 1)


def main():
    assembly = pd.read_csv(DATA / "genome_assembly_and_cog_metrics.csv")
    annot = pd.read_csv(DATA / "recalculated_annotations.csv")
    base = assembly[["Genome_ID", "Phylogroup", "Genome_Size_Mb", "Total_CDS"]].merge(
        annot.drop(columns=["Phylogroup"]), on="Genome_ID", validate="one_to_one"
    )

    pathway = pd.read_csv(DATA / "genome_kegg_pathway_matrix.csv")
    crispr = pd.read_csv(DATA / "all_genomes_crispr_detailed.csv")
    spacers = pd.read_csv(DATA / "all_genomes_crispr_spacers_detailed.csv")
    tmhmm = pd.read_csv(DATA / "all_genomes_tmhmm_positives_detailed.csv")
    signalp, repaired_signalp_rows = parse_signalp(DATA / "all_genomes_signalp_positives_detailed.csv")

    qc = []
    for name, df in [
        ("genome_assembly_and_cog_metrics.csv", assembly),
        ("recalculated_annotations.csv", annot),
        ("genome_kegg_pathway_matrix.csv", pathway),
        ("all_genomes_signalp_positives_detailed.csv", signalp),
        ("all_genomes_tmhmm_positives_detailed.csv", tmhmm),
        ("all_genomes_crispr_detailed.csv", crispr),
        ("all_genomes_crispr_spacers_detailed.csv", spacers),
    ]:
        qc.append(
            {
                "file": name,
                "rows": len(df),
                "columns": len(df.columns),
                "unique_genomes": df["Genome_ID"].nunique() if "Genome_ID" in df.columns else np.nan,
                "overlap_with_assembly": len(set(assembly["Genome_ID"]) & set(df["Genome_ID"]))
                if "Genome_ID" in df.columns
                else np.nan,
            }
        )
    qc.append(
        {
            "file": "all_genomes_signalp_positives_detailed.csv repaired_rows",
            "rows": repaired_signalp_rows,
            "columns": np.nan,
            "unique_genomes": np.nan,
            "overlap_with_assembly": np.nan,
        }
    )
    pd.DataFrame(qc).to_csv(OUT / "data_qc_summary.csv", index=False)

    # SignalP aggregation.
    signalp_counts = (
        signalp.groupby(["Genome_ID", "Prediction_Type"]).size().unstack(fill_value=0).reset_index()
    )
    for col in ["SP", "LIPO", "TAT", "TATLIPO", "PILIN"]:
        if col not in signalp_counts:
            signalp_counts[col] = 0
    signalp_counts["SignalP_Total"] = signalp_counts.drop(columns=["Genome_ID"]).sum(axis=1)
    signalp_counts = base[["Genome_ID", "Phylogroup", "Total_CDS", "Genome_Size_Mb"]].merge(
        signalp_counts, on="Genome_ID", how="left"
    ).fillna(0)
    signalp_counts["SignalP_per_1000_CDS"] = signalp_counts["SignalP_Total"] / signalp_counts["Total_CDS"] * 1000
    signalp_counts.to_csv(OUT / "signalp_genome_counts.csv", index=False)
    sig_summary = []
    for metric in ["SignalP_Total", "SP", "LIPO", "TAT", "SignalP_per_1000_CDS"]:
        s = summarize_metric(signalp_counts, metric)
        s.insert(0, "metric", metric)
        sig_summary.append(s.reset_index())
    pd.concat(sig_summary).to_csv(OUT / "signalp_phylogroup_summary.csv", index=False)

    # CRISPR array and spacer aggregation.
    crispr["Spacer_Count_Est"] = (crispr["N_repeats"].astype(float) - 1).clip(lower=0)
    crispr_counts = crispr.groupby("Genome_ID").agg(
        CRISPR_Arrays_Detail=("CRISPR", "count"),
        Trusted_Arrays=("Trusted", "sum"),
        Total_Repeats=("N_repeats", "sum"),
        Estimated_Spacers=("Spacer_Count_Est", "sum"),
        Mean_Repeat_Identity=("Repeat_identity", "mean"),
        Mean_Spacer_Identity=("Spacer_identity", "mean"),
        Mean_Subtype_Probability=("Subtype_probability", "mean"),
    )
    subtype = crispr.groupby(["Genome_ID", "Subtype"]).size().unstack(fill_value=0)
    crispr_counts = crispr_counts.join(subtype, how="left").reset_index()
    spacer_counts = spacers.groupby("Genome_ID").agg(
        Spacer_Sequences=("Sequence", "count"),
        Unique_Spacer_Sequences=("Sequence", "nunique"),
        Mean_Spacer_Length=("Sequence", lambda s: s.astype(str).str.len().mean()),
    )
    shared_spacers = spacers.groupby("Sequence")["Genome_ID"].nunique()
    shared_sequence_set = set(shared_spacers[shared_spacers > 1].index)
    spacers["Shared_Spacer"] = spacers["Sequence"].isin(shared_sequence_set)
    shared_counts = spacers.groupby("Genome_ID").agg(
        Shared_Spacer_Sequences=("Shared_Spacer", "sum"),
    )
    spacer_counts = spacer_counts.join(shared_counts, how="left")
    spacer_counts["Shared_Spacer_Fraction"] = (
        spacer_counts["Shared_Spacer_Sequences"] / spacer_counts["Spacer_Sequences"]
    )
    crispr_counts = crispr_counts.merge(spacer_counts.reset_index(), on="Genome_ID", how="left")
    crispr_counts = base[["Genome_ID", "Phylogroup"]].merge(crispr_counts, on="Genome_ID", how="left").fillna(0)
    crispr_counts.to_csv(OUT / "crispr_genome_counts.csv", index=False)
    crispr_summary = []
    for metric in [
        "CRISPR_Arrays_Detail",
        "Trusted_Arrays",
        "Total_Repeats",
        "Estimated_Spacers",
        "Spacer_Sequences",
        "Unique_Spacer_Sequences",
        "Shared_Spacer_Fraction",
        "I-E",
        "I-F",
    ]:
        if metric in crispr_counts:
            s = summarize_metric(crispr_counts, metric)
            s.insert(0, "metric", metric)
            crispr_summary.append(s.reset_index())
    pd.concat(crispr_summary).to_csv(OUT / "crispr_phylogroup_summary.csv", index=False)

    spacer_phylo = spacers.merge(base[["Genome_ID", "Phylogroup"]], on="Genome_ID", how="left")
    spacer_group_presence = (
        spacer_phylo.dropna(subset=["Phylogroup"])
        .groupby(["Sequence", "Phylogroup"])["Genome_ID"]
        .nunique()
        .unstack(fill_value=0)
    )
    spacer_group_presence["n_phylogroups"] = (spacer_group_presence[PHYLO_ORDER] > 0).sum(axis=1)
    spacer_group_presence["total_genomes"] = spacer_group_presence[PHYLO_ORDER].sum(axis=1)
    spacer_group_presence.sort_values(["total_genomes", "n_phylogroups"], ascending=False).to_csv(
        OUT / "crispr_spacer_sharing_by_sequence.csv"
    )

    # TMHMM aggregation from detailed rows.
    tmhmm_counts = tmhmm.groupby("Genome_ID").agg(
        TMHMM_Proteins=("Protein_ID", "count"),
        TMHMM_Helices=("TMH_Count", "sum"),
        TMHMM_Mean_Helices_Per_Protein=("TMH_Count", "mean"),
        TMHMM_Max_Helices=("TMH_Count", "max"),
    )
    tmhmm_counts = base[["Genome_ID", "Phylogroup", "Total_CDS"]].merge(
        tmhmm_counts.reset_index(), on="Genome_ID", how="left"
    )
    tmhmm_counts["TMHMM_Proteins_per_1000_CDS"] = tmhmm_counts["TMHMM_Proteins"] / tmhmm_counts["Total_CDS"] * 1000
    tmhmm_counts.to_csv(OUT / "tmhmm_genome_counts.csv", index=False)
    tmhmm_summary = []
    for metric in [
        "TMHMM_Proteins",
        "TMHMM_Helices",
        "TMHMM_Mean_Helices_Per_Protein",
        "TMHMM_Max_Helices",
        "TMHMM_Proteins_per_1000_CDS",
    ]:
        s = summarize_metric(tmhmm_counts, metric)
        s.insert(0, "metric", metric)
        tmhmm_summary.append(s.reset_index())
    pd.concat(tmhmm_summary).to_csv(OUT / "tmhmm_phylogroup_summary.csv", index=False)

    # KEGG pathway: ko columns only because map columns duplicate the same identifiers.
    ko_cols = [c for c in pathway.columns if c.startswith("ko")]
    map_cols = [c for c in pathway.columns if c.startswith("map")]
    map_duplicate = True
    for ko in ko_cols:
        mp = "map" + ko[2:]
        if mp in pathway.columns and not np.array_equal(pathway[ko].to_numpy(), pathway[mp].to_numpy()):
            map_duplicate = False
            break
    path_base = base[["Genome_ID", "Phylogroup"]].merge(pathway[["Genome_ID"] + ko_cols], on="Genome_ID")
    x_counts = path_base[ko_cols].to_numpy(dtype=float)
    x_binary = (x_counts > 0).astype(np.uint8)
    labels = path_base["Phylogroup"].to_numpy()

    scores, explained = pca_svd(x_binary, 4)
    pca_df = path_base[["Genome_ID", "Phylogroup"]].copy()
    for i in range(scores.shape[1]):
        pca_df[f"PC{i+1}"] = scores[:, i]
    pca_df.to_csv(OUT / "pathway_binary_pca_scores.csv", index=False)
    pd.DataFrame(
        {"component": [f"PC{i+1}" for i in range(len(explained))], "variance_explained": explained}
    ).to_csv(OUT / "pathway_binary_pca_variance.csv", index=False)

    row_sums = x_counts.sum(axis=1, keepdims=True)
    rel_counts = np.divide(x_counts, row_sums, out=np.zeros_like(x_counts), where=row_sums > 0)
    x_hellinger = np.sqrt(rel_counts)
    hscores, hexplained = pca_svd(x_hellinger, 4)
    hpca_df = path_base[["Genome_ID", "Phylogroup"]].copy()
    for i in range(hscores.shape[1]):
        hpca_df[f"PC{i+1}"] = hscores[:, i]
    hpca_df.to_csv(OUT / "pathway_hellinger_pca_scores.csv", index=False)
    pd.DataFrame(
        {"component": [f"PC{i+1}" for i in range(len(hexplained))], "variance_explained": hexplained}
    ).to_csv(OUT / "pathway_hellinger_pca_variance.csv", index=False)

    d = jaccard_distance(x_binary)
    full_perm = permanova_from_dist(d, labels, permutations=999)
    coords, eigvals = pcoa_coordinates(d)
    permdisp, dist_to_centroid = permdisp_from_coords(coords, labels, permutations=999)
    centroid_df = path_base[["Genome_ID", "Phylogroup"]].copy()
    centroid_df["Jaccard_PCoA_Distance_To_Group_Centroid"] = dist_to_centroid
    centroid_df.to_csv(OUT / "pathway_permdisp_distances.csv", index=False)
    pd.DataFrame([full_perm]).to_csv(OUT / "pathway_permanova_full.csv", index=False)
    pd.DataFrame([permdisp]).to_csv(OUT / "pathway_permdisp_full.csv", index=False)

    hd = euclidean_distance(x_hellinger)
    h_full_perm = permanova_from_dist(hd, labels, permutations=999)
    h_permdisp, h_dist_to_centroid = permdisp_from_coords(x_hellinger, labels, permutations=999)
    h_centroid_df = path_base[["Genome_ID", "Phylogroup"]].copy()
    h_centroid_df["Hellinger_Distance_To_Group_Centroid"] = h_dist_to_centroid
    h_centroid_df.to_csv(OUT / "pathway_hellinger_permdisp_distances.csv", index=False)
    pd.DataFrame([h_full_perm]).to_csv(OUT / "pathway_hellinger_permanova_full.csv", index=False)
    pd.DataFrame([h_permdisp]).to_csv(OUT / "pathway_hellinger_permdisp_full.csv", index=False)

    # Balanced sensitivity designs.
    balanced_rows = []
    designs = [
        ("n31_A_B2_B1_D_E_C_F", ["A", "B2", "B1", "D", "E", "C", "F"], 31),
        ("n50_A_B2_B1_D_E", ["A", "B2", "B1", "D", "E"], 50),
    ]
    for design, groups, n_each in designs:
        for rep in range(100):
            chosen = []
            for group in groups:
                idx = np.where(labels == group)[0]
                chosen.extend(RNG.choice(idx, size=n_each, replace=False).tolist())
            chosen = np.array(chosen, dtype=int)
            sub_d = d[np.ix_(chosen, chosen)]
            sub_labels = labels[chosen]
            stats = permanova_from_dist(sub_d, sub_labels, permutations=0)
            balanced_rows.append(
                {
                    "design": design,
                    "replicate": rep + 1,
                    "n_each": n_each,
                    "groups": len(groups),
                    "n_total": len(chosen),
                    "F": stats["F"],
                    "R2": stats["R2"],
                }
            )
    balanced = pd.DataFrame(balanced_rows)
    balanced.to_csv(OUT / "pathway_balanced_permanova_replicates.csv", index=False)
    balanced.groupby("design")[["F", "R2"]].agg(["mean", "std", "median", "min", "max"]).to_csv(
        OUT / "pathway_balanced_permanova_summary.csv"
    )

    h_balanced_rows = []
    for design, groups, n_each in designs:
        for rep in range(100):
            chosen = []
            for group in groups:
                idx = np.where(labels == group)[0]
                chosen.extend(RNG.choice(idx, size=n_each, replace=False).tolist())
            chosen = np.array(chosen, dtype=int)
            sub_d = hd[np.ix_(chosen, chosen)]
            sub_labels = labels[chosen]
            stats = permanova_from_dist(sub_d, sub_labels, permutations=0)
            h_balanced_rows.append(
                {
                    "design": design,
                    "replicate": rep + 1,
                    "n_each": n_each,
                    "groups": len(groups),
                    "n_total": len(chosen),
                    "F": stats["F"],
                    "R2": stats["R2"],
                }
            )
    h_balanced = pd.DataFrame(h_balanced_rows)
    h_balanced.to_csv(OUT / "pathway_hellinger_balanced_permanova_replicates.csv", index=False)
    h_balanced.groupby("design")[["F", "R2"]].agg(["mean", "std", "median", "min", "max"]).to_csv(
        OUT / "pathway_hellinger_balanced_permanova_summary.csv"
    )

    # Pathway global chi-square tests on binary presence.
    path_tests = []
    for col in ko_cols:
        present = path_base[col].to_numpy() > 0
        obs = []
        for group in PHYLO_ORDER:
            m = path_base["Phylogroup"].to_numpy() == group
            obs.append([int((present & m).sum()), int((~present & m).sum())])
        stat, dof, p = chi_square_table(np.array(obs))
        freq = path_base.assign(present=present).groupby("Phylogroup")["present"].mean().reindex(PHYLO_ORDER)
        path_tests.append(
            {
                "Pathway": col,
                "chi2": stat,
                "df": dof,
                "p": p,
                "overall_frequency": float(present.mean()),
                "max_group": freq.idxmax(),
                "max_frequency": float(freq.max()),
                "min_group": freq.idxmin(),
                "min_frequency": float(freq.min()),
            }
        )
    path_tests = pd.DataFrame(path_tests)
    path_tests["p_bh"] = bh_adjust(path_tests["p"].to_numpy())
    path_tests["p_bonferroni"] = np.minimum(1.0, path_tests["p"] * len(path_tests))
    path_tests.sort_values(["p_bh", "p"]).to_csv(OUT / "pathway_chisq_multipletest.csv", index=False)

    # COG global chi-square tests: group counts for one COG vs all other COG counts.
    cog_cols = [c for c in assembly.columns if c.startswith("COG_")]
    cog_tests = []
    cog_total = assembly[cog_cols].sum(axis=1)
    for col in cog_cols:
        obs = []
        for group in PHYLO_ORDER:
            m = assembly["Phylogroup"] == group
            c = int(assembly.loc[m, col].sum())
            rest = int((cog_total[m] - assembly.loc[m, col]).sum())
            obs.append([c, rest])
        stat, dof, p = chi_square_table(np.array(obs))
        rel = (
            assembly.assign(rel=assembly[col] / cog_total)
            .groupby("Phylogroup")["rel"]
            .mean()
            .reindex(PHYLO_ORDER)
            * 100
        )
        cog_tests.append(
            {
                "COG": col,
                "chi2": stat,
                "df": dof,
                "p": p,
                "max_group": rel.idxmax(),
                "max_relative_percent": float(rel.max()),
                "min_group": rel.idxmin(),
                "min_relative_percent": float(rel.min()),
            }
        )
    cog_tests = pd.DataFrame(cog_tests)
    cog_tests["p_bh"] = bh_adjust(cog_tests["p"].to_numpy())
    cog_tests["p_bonferroni"] = np.minimum(1.0, cog_tests["p"] * len(cog_tests))
    cog_tests.sort_values(["p_bh", "p"]).to_csv(OUT / "cog_chisq_multipletest.csv", index=False)

    # Metric permutation ANOVA summaries for newly aggregated detail metrics.
    metric_rows = []
    merged_metrics = base[["Genome_ID", "Phylogroup", "Genome_Size_Mb", "Total_CDS"]].merge(
        signalp_counts[["Genome_ID", "SignalP_Total", "SignalP_per_1000_CDS"]], on="Genome_ID", how="left"
    ).merge(
        crispr_counts[
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
        how="left",
    ).merge(
        tmhmm_counts[["Genome_ID", "TMHMM_Proteins", "TMHMM_Helices", "TMHMM_Proteins_per_1000_CDS"]],
        on="Genome_ID",
        how="left",
    )
    for metric in [
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
    ]:
        sub_metric = merged_metrics[["Phylogroup", metric]].dropna()
        f, p = permutation_anova(
            sub_metric[metric].to_numpy(float),
            sub_metric["Phylogroup"].to_numpy(),
            permutations=999,
        )
        metric_rows.append({"metric": metric, "F_permutation_anova": f, "p": p, "permutations": 999})
    pd.DataFrame(metric_rows).to_csv(OUT / "detail_metric_permutation_anova.csv", index=False)

    summary = [
        f"KEGG pathway columns: {len(ko_cols)} ko columns and {len(map_cols)} map columns; map duplicates ko: {map_duplicate}.",
        f"Pathway matrix values are counts/abundances, not strict 0/1; binary presence was defined as value > 0 for Jaccard/PERMANOVA.",
        f"SignalP rows parsed: {len(signalp)}; repaired rows with unescaped commas: {repaired_signalp_rows}.",
        f"TMHMM detail rows: {len(tmhmm)} across {tmhmm['Genome_ID'].nunique()} genomes.",
        f"CRISPR spacer rows: {len(spacers)}; unique spacer sequences: {spacers['Sequence'].nunique()}.",
        f"Full pathway PERMANOVA: F={full_perm['F']:.4f}, R2={full_perm['R2']:.4f}, p={full_perm['p']:.4g}.",
        f"Full pathway PERMDISP: F={permdisp['F']:.4f}, p={permdisp['p']:.4g}.",
        f"Pathway PCA variance: PC1={explained[0]*100:.2f}%, PC2={explained[1]*100:.2f}%.",
        f"Hellinger pathway PERMANOVA: F={h_full_perm['F']:.4f}, R2={h_full_perm['R2']:.4f}, p={h_full_perm['p']:.4g}.",
        f"Hellinger pathway PERMDISP: F={h_permdisp['F']:.4f}, p={h_permdisp['p']:.4g}.",
        f"Hellinger pathway PCA variance: PC1={hexplained[0]*100:.2f}%, PC2={hexplained[1]*100:.2f}%.",
    ]
    (OUT / "analysis_summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print("\n".join(summary))
    print(f"Outputs written to {OUT}")


if __name__ == "__main__":
    main()
