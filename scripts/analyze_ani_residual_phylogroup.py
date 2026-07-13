from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(r"D:\pipeline\analysis_v5")
DATA = ROOT / "00.data"
ANALYSIS = ROOT / "03.analysis"
ANALYSIS.mkdir(exist_ok=True)

PHYLO_ORDER = ["A", "B2", "B1", "D", "E", "C", "F", "G", "U"]


def norm_id(x: str) -> str:
    x = str(x).strip()
    for suffix in [".fna", ".fa", ".fasta"]:
        if x.endswith(suffix):
            x = x[: -len(suffix)]
    return x


def align_matrix(df: pd.DataFrame, ids: list[str]) -> tuple[pd.DataFrame, list[str]]:
    tmp = df.copy()
    tmp["norm_id"] = tmp["Genome_ID"].map(norm_id)
    tmp = tmp.drop_duplicates("norm_id", keep="first").set_index("norm_id")
    keep = [i for i in ids if i in tmp.index]
    return tmp.loc[keep].copy(), keep


def read_ani_matrix():
    raw = pd.read_csv(DATA / "core_genome_ani_distance_matrix.csv")
    row_ids = raw.iloc[:, 0].map(norm_id).tolist()
    col_ids = [norm_id(c) for c in raw.columns[1:]]
    mat = raw.iloc[:, 1:].to_numpy(float)
    row_index = {g: i for i, g in enumerate(row_ids)}
    col_index = {g: i for i, g in enumerate(col_ids)}
    common = [g for g in row_ids if g in col_index]
    row_pos = [row_index[g] for g in common]
    col_pos = [col_index[g] for g in common]
    mat = mat[np.ix_(row_pos, col_pos)]
    asym = np.abs(mat - mat.T)
    sym = (mat + mat.T) / 2.0
    np.fill_diagonal(sym, 0.0)
    qc = {
        "raw_rows": len(row_ids),
        "raw_cols": len(col_ids),
        "common_square_ids": len(common),
        "missing_values": int(np.isnan(mat).sum()),
        "diag_abs_max_before_sym": float(np.nanmax(np.abs(np.diag(mat)))),
        "asymmetry_abs_mean": float(np.nanmean(asym)),
        "asymmetry_abs_median": float(np.nanmedian(asym)),
        "asymmetry_abs_max": float(np.nanmax(asym)),
        "distance_min_sym": float(np.nanmin(sym)),
        "distance_median_upper_sym": float(np.nanmedian(sym[np.triu_indices_from(sym, 1)])),
        "distance_max_sym": float(np.nanmax(sym)),
    }
    return common, sym, qc


def euclidean_distance(X: np.ndarray) -> np.ndarray:
    X = np.asarray(X, dtype=float)
    sq = np.sum(X * X, axis=1)
    d2 = sq[:, None] + sq[None, :] - 2 * (X @ X.T)
    d2[d2 < 0] = 0
    return np.sqrt(d2)


def hellinger_matrix(X: np.ndarray) -> np.ndarray:
    X = np.asarray(X, dtype=float)
    row_sum = X.sum(axis=1)
    row_sum[row_sum == 0] = 1.0
    return np.sqrt(X / row_sum[:, None])


def zscore(X: np.ndarray) -> np.ndarray:
    X = np.asarray(X, dtype=float)
    mu = np.nanmean(X, axis=0)
    sd = np.nanstd(X, axis=0)
    sd[sd == 0] = 1.0
    X = np.where(np.isnan(X), mu, X)
    return (X - mu) / sd


def rank_vector(x: np.ndarray) -> np.ndarray:
    return pd.Series(x).rank(method="average").to_numpy(float)


def corr(x: np.ndarray, y: np.ndarray, method: str) -> float:
    if method == "spearman":
        x, y = rank_vector(x), rank_vector(y)
    x = x - x.mean()
    y = y - y.mean()
    denom = np.sqrt(np.sum(x * x) * np.sum(y * y))
    return float(np.sum(x * y) / denom) if denom else np.nan


def mantel(feature_d: np.ndarray, ani_d: np.ndarray, permutations: int = 999, seed: int = 7):
    rng = np.random.default_rng(seed)
    n = feature_d.shape[0]
    iu = np.triu_indices(n, 1)
    x = ani_d[iu]
    y = feature_d[iu]
    x_rank = rank_vector(x)
    y_rank = rank_vector(y)
    rank_mat = np.zeros_like(feature_d, dtype=float)
    rank_mat[iu] = y_rank
    rank_mat[(iu[1], iu[0])] = y_rank
    observed_pearson = corr(x, y, "pearson")
    observed_spearman = corr(x_rank, y_rank, "pearson")
    pearson_ge = 0
    spearman_ge = 0
    for _ in range(permutations):
        perm = rng.permutation(n)
        yp = feature_d[np.ix_(perm, perm)][iu]
        yrp = rank_mat[np.ix_(perm, perm)][iu]
        rp = corr(x, yp, "pearson")
        rs = corr(x_rank, yrp, "pearson")
        pearson_ge += abs(rp) >= abs(observed_pearson)
        spearman_ge += abs(rs) >= abs(observed_spearman)
    return {
        "mantel_pearson_r": observed_pearson,
        "mantel_pearson_p": (pearson_ge + 1) / (permutations + 1),
        "mantel_spearman_r": observed_spearman,
        "mantel_spearman_p": (spearman_ge + 1) / (permutations + 1),
        "permutations": permutations,
    }


def gower_center(D: np.ndarray) -> np.ndarray:
    n = D.shape[0]
    J = np.eye(n) - np.ones((n, n)) / n
    return -0.5 * J @ (D * D) @ J


def pcoa(D: np.ndarray):
    B = gower_center(D)
    vals, vecs = np.linalg.eigh(B)
    order = np.argsort(vals)[::-1]
    vals = vals[order]
    vecs = vecs[:, order]
    positive = vals > 1e-12
    vals_pos = vals[positive]
    vecs_pos = vecs[:, positive]
    scores = vecs_pos * np.sqrt(vals_pos)
    return vals_pos, scores


def design_matrix_ani(scores: np.ndarray, vals: np.ndarray, target_var: float = 0.80, max_axes: int = 20):
    explained = vals / vals.sum()
    k = int(np.searchsorted(np.cumsum(explained), target_var) + 1)
    k = min(max(k, 5), max_axes, scores.shape[1])
    return scores[:, :k], k, float(explained[:k].sum())


def dummy_matrix(labels: np.ndarray) -> np.ndarray:
    cats = [g for g in PHYLO_ORDER if g in set(labels)]
    # Drop the first category because the intercept is included.
    return np.column_stack([(labels == g).astype(float) for g in cats[1:]])


def r2_from_B(B: np.ndarray, X: np.ndarray) -> float:
    intercept = np.ones((B.shape[0], 1))
    X = np.column_stack([intercept, X])
    denom = np.trace(B)
    if not denom:
        return np.nan
    Q, _ = np.linalg.qr(X)
    return float(np.sum(Q * (B @ Q)) / denom)


def r2_db(D: np.ndarray, X: np.ndarray) -> float:
    return r2_from_B(gower_center(D), X)


def residual_phylogroup(feature_d: np.ndarray, labels: np.ndarray, ani_scores: np.ndarray, permutations: int = 999, seed: int = 11):
    rng = np.random.default_rng(seed)
    B = gower_center(feature_d)
    G = dummy_matrix(labels)
    r2_group = r2_from_B(B, G)
    r2_ani = r2_from_B(B, ani_scores)
    r2_full = r2_from_B(B, np.column_stack([ani_scores, G]))
    unique_group = r2_full - r2_ani
    ge = 0
    for _ in range(permutations):
        perm_labels = labels[rng.permutation(len(labels))]
        Gp = dummy_matrix(perm_labels)
        stat = r2_from_B(B, np.column_stack([ani_scores, Gp])) - r2_ani
        ge += stat >= unique_group
    return {
        "phylogroup_only_R2": r2_group,
        "ANI_PC_R2": r2_ani,
        "full_ANI_PC_plus_phylogroup_R2": r2_full,
        "unique_phylogroup_R2_after_ANI": unique_group,
        "unique_phylogroup_p_label_permutation": (ge + 1) / (permutations + 1),
        "permutations": permutations,
    }


def build_feature_distances(ids: list[str]):
    assembly = pd.read_csv(DATA / "genome_assembly_and_cog_metrics.csv")
    annot = pd.read_csv(DATA / "recalculated_annotations.csv")
    pathway = pd.read_csv(DATA / "genome_kegg_pathway_matrix.csv")
    signalp = pd.read_csv(ANALYSIS / "signalp_genome_counts.csv")
    tmhmm = pd.read_csv(ANALYSIS / "tmhmm_genome_counts.csv")
    crispr_spacers = pd.read_csv(DATA / "all_genomes_crispr_spacers_detailed.csv")

    assembly_a, keep = align_matrix(assembly, ids)
    annot_a, keep2 = align_matrix(annot, keep)
    pathway_a, keep3 = align_matrix(pathway, keep2)
    signalp_a, keep4 = align_matrix(signalp, keep3)
    tmhmm_a, keep5 = align_matrix(tmhmm, keep4)
    # Keep all data types aligned to the most restrictive set.
    keep = keep5
    assembly_a, _ = align_matrix(assembly, keep)
    annot_a, _ = align_matrix(annot, keep)
    pathway_a, _ = align_matrix(pathway, keep)
    signalp_a, _ = align_matrix(signalp, keep)
    tmhmm_a, _ = align_matrix(tmhmm, keep)

    labels = assembly_a["Phylogroup"].to_numpy(str)

    ko_cols = [c for c in pathway_a.columns if c.startswith("ko")]
    kegg_x = hellinger_matrix(pathway_a[ko_cols].to_numpy(float))
    kegg_d = euclidean_distance(kegg_x)

    cog_cols = [c for c in assembly_a.columns if c.startswith("COG_")]
    cog_x = hellinger_matrix(assembly_a[cog_cols].to_numpy(float))
    cog_d = euclidean_distance(cog_x)

    gem_cols = ["GEM_Genes", "GEM_Reactions", "GEM_Metabolites"]
    gem_x = zscore(annot_a[gem_cols].to_numpy(float))
    gem_d = euclidean_distance(gem_x)

    surf = signalp_a[["SignalP_per_1000_CDS", "SP", "LIPO", "TAT"]].copy()
    for c in ["SP", "LIPO", "TAT"]:
        surf[c] = surf[c] / signalp_a["Total_CDS"] * 1000
    surf["TMHMM_Proteins_per_1000_CDS"] = tmhmm_a["TMHMM_Proteins_per_1000_CDS"]
    surf["TMHMM_Helices_per_1000_CDS"] = tmhmm_a["TMHMM_Helices"] / tmhmm_a["Total_CDS"] * 1000
    surface_d = euclidean_distance(zscore(surf.to_numpy(float)))

    spacer_norm = crispr_spacers[["Genome_ID", "Sequence"]].copy()
    spacer_norm["norm_id"] = spacer_norm["Genome_ID"].map(norm_id)
    spacer_norm = spacer_norm[spacer_norm["norm_id"].isin(keep)]
    genome_index = {g: i for i, g in enumerate(keep)}
    spacer_codes = {s: i for i, s in enumerate(pd.unique(spacer_norm["Sequence"]))}
    X = np.zeros((len(keep), len(spacer_codes)), dtype=np.uint8)
    for row in spacer_norm.itertuples(index=False):
        X[genome_index[row.norm_id], spacer_codes[row.Sequence]] = 1
    shared = X @ X.T
    counts = X.sum(axis=1).astype(float)
    union = counts[:, None] + counts[None, :] - shared
    crispr_d = np.ones_like(union, dtype=float)
    mask = union > 0
    crispr_d[mask] = 1.0 - shared[mask] / union[mask]
    crispr_d[~mask] = 0.0
    np.fill_diagonal(crispr_d, 0.0)

    features = {
        "KEGG_Hellinger": kegg_d,
        "COG_Hellinger": cog_d,
        "GEM_topology_zEuclidean": gem_d,
        "CRISPR_spacer_Jaccard": crispr_d,
        "SignalP_TMHMM_zEuclidean": surface_d,
    }
    metadata = pd.DataFrame({"norm_id": keep, "Genome_ID": assembly_a["Genome_ID"].to_numpy(), "Phylogroup": labels})
    return keep, labels, features, metadata


def main():
    ani_ids, ani_d, qc = read_ani_matrix()
    feature_ids, labels, features, metadata = build_feature_distances(ani_ids)
    pos = {g: i for i, g in enumerate(ani_ids)}
    idx = [pos[g] for g in feature_ids]
    ani_d = ani_d[np.ix_(idx, idx)]

    vals, ani_scores_all = pcoa(ani_d)
    ani_scores, n_axes, var_explained = design_matrix_ani(ani_scores_all, vals)

    qc.update(
        {
            "analysis_genomes_after_feature_intersection": len(feature_ids),
            "ANI_PCoA_positive_axes": int(len(vals)),
            "ANI_PC_axes_used_for_residual_models": int(n_axes),
            "ANI_PC_variance_explained_used": var_explained,
        }
    )
    pd.DataFrame([qc]).to_csv(ANALYSIS / "ani_distance_qc.csv", index=False)
    metadata.to_csv(ANALYSIS / "ani_residual_analysis_genomes.csv", index=False)
    pd.DataFrame(
        ani_scores[:, : min(20, ani_scores.shape[1])],
        index=metadata["Genome_ID"],
        columns=[f"ANI_PC{i+1}" for i in range(min(20, ani_scores.shape[1]))],
    ).to_csv(ANALYSIS / "ani_pcoa_scores.csv")

    mantel_rows = []
    residual_rows = []
    for i, (name, D) in enumerate(features.items()):
        m = mantel(D, ani_d, permutations=999, seed=100 + i)
        m.insert if False else None
        mantel_rows.append({"feature_space": name, "n": len(feature_ids), **m})
        r = residual_phylogroup(D, labels, ani_scores, permutations=999, seed=200 + i)
        residual_rows.append(
            {
                "feature_space": name,
                "n": len(feature_ids),
                "ANI_PC_axes": n_axes,
                "ANI_PC_variance_explained": var_explained,
                **r,
            }
        )

    mantel_df = pd.DataFrame(mantel_rows)
    residual_df = pd.DataFrame(residual_rows)
    mantel_df.to_csv(ANALYSIS / "ani_mantel_distance_correlations.csv", index=False)
    residual_df.to_csv(ANALYSIS / "ani_residual_phylogroup_effects.csv", index=False)

    with (ANALYSIS / "ani_residual_analysis_summary.txt").open("w", encoding="utf-8") as f:
        f.write("Core-genome ANI residual phylogroup analysis\n")
        f.write("================================================\n")
        for k, v in qc.items():
            f.write(f"{k}: {v}\n")
        f.write("\nMantel correlations against ANI distance\n")
        f.write(mantel_df.to_string(index=False, float_format=lambda x: f"{x:.5f}"))
        f.write("\n\nResidual phylogroup effect after ANI PCs\n")
        f.write(residual_df.to_string(index=False, float_format=lambda x: f"{x:.5f}"))
        f.write("\n")

    print("ANI QC")
    for k, v in qc.items():
        print(f"{k}: {v}")
    print("\nMantel")
    print(mantel_df.to_string(index=False, float_format=lambda x: f"{x:.5f}"))
    print("\nResidual")
    print(residual_df.to_string(index=False, float_format=lambda x: f"{x:.5f}"))


if __name__ == "__main__":
    main()
