from __future__ import annotations

import json
import argparse
from pathlib import Path

import numpy as np
import pandas as pd


AXIS_COUNTS = (5, 10, 20, 40)
PERMUTATIONS = 999
SEED = 20260912
GROUP_ORDER = ("A", "B2", "B1", "D", "E", "C", "F", "G", "U")


def hellinger(x: np.ndarray) -> np.ndarray:
    row_sums = x.sum(axis=1)
    row_sums[row_sums == 0] = 1.0
    return np.sqrt(x / row_sums[:, None])


def pcoa(distance: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n = distance.shape[0]
    j = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * j @ (distance ** 2) @ j
    vals, vecs = np.linalg.eigh(b)
    order = vals.argsort()[::-1]
    vals = vals[order]
    vecs = vecs[:, order]
    positive = vals > 1e-12
    return vecs[:, positive] * np.sqrt(vals[positive]), vals


def dummies(groups: pd.Series) -> np.ndarray:
    return np.column_stack([(groups.to_numpy(str) == group).astype(float)
                            for group in GROUP_ORDER[1:]])


def design(parts: list[np.ndarray], n: int) -> np.ndarray:
    return np.column_stack([np.ones(n), *parts])


def fit_sse(x: np.ndarray, y: np.ndarray) -> float:
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    residual = y - x @ beta
    return float(np.sum(residual ** 2))


def r2(sse: float, sst: float) -> float:
    return 1.0 - sse / sst


def adjusted_r2(value: float, n: int, predictors: int) -> float:
    return 1.0 - (1.0 - value) * (n - 1) / (n - predictors - 1)


def evaluate(y: np.ndarray, pcs: np.ndarray, groups: pd.Series, permutations: int, seed: int) -> dict:
    n = len(groups)
    g = dummies(groups)
    x0 = design([], n)
    xa = design([pcs], n)
    xg = design([g], n)
    xf = design([pcs, g], n)
    centered = y - y.mean(axis=0, keepdims=True)
    sst = float(np.sum(centered ** 2))
    sse0 = fit_sse(x0, y)
    ssea = fit_sse(xa, y)
    sseg = fit_sse(xg, y)
    ssef = fit_sse(xf, y)
    r2_ani = r2(ssea, sst)
    r2_group = r2(sseg, sst)
    r2_full = r2(ssef, sst)
    unique = (ssea - ssef) / sst
    df_num = g.shape[1]
    df_den = n - xf.shape[1]
    f_obs = ((ssea - ssef) / df_num) / (ssef / df_den)

    beta, *_ = np.linalg.lstsq(xa, y, rcond=None)
    fitted = xa @ beta
    residual = y - fitted
    rng = np.random.default_rng(seed)
    ge = 0
    for _ in range(permutations):
        yp = fitted + residual[rng.permutation(n), :]
        sr = fit_sse(xa, yp)
        sf = fit_sse(xf, yp)
        fp = ((sr - sf) / df_num) / (sf / df_den)
        ge += fp >= f_obs

    return {
        "n_genomes": n,
        "n_features": y.shape[1],
        "phylogroup_only_R2": r2_group,
        "phylogroup_only_adjusted_R2": adjusted_r2(r2_group, n, g.shape[1]),
        "ANI_only_R2": r2_ani,
        "ANI_only_adjusted_R2": adjusted_r2(r2_ani, n, pcs.shape[1]),
        "combined_R2": r2_full,
        "combined_adjusted_R2": adjusted_r2(r2_full, n, pcs.shape[1] + g.shape[1]),
        "unique_phylogroup_R2_after_ANI": unique,
        "partial_F": f_obs,
        "p_value_freedman_lane": (ge + 1) / (permutations + 1),
        "permutations": permutations,
        "permutation_seed": seed,
        "df_num": df_num,
        "df_den": df_den,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="ANI-PCoA dimensional sensitivity for ontology-filtered KEGG profiles")
    parser.add_argument("--ani-distance", required=True, type=Path, help="Square ANI distance CSV/CSV.GZ indexed by Genome_UID")
    parser.add_argument("--kegg-binary", required=True, type=Path, help="Binary KEGG feature matrix with Genome_UID and Phylogroup")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--permutations", type=int, default=999)
    parser.add_argument("--seed", type=int, default=20260912)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    kegg = pd.read_csv(args.kegg_binary)
    feature_cols = [c for c in kegg.columns if c.startswith("kegg_pathway_")]
    ids = kegg["Genome_UID"].astype(str).tolist()
    groups = kegg["Phylogroup"].astype(str)
    y = hellinger(kegg[feature_cols].to_numpy(float))

    ani = pd.read_csv(args.ani_distance).set_index("Genome_UID")
    ani = ani.loc[ids, ids].apply(pd.to_numeric, errors="coerce")
    distance = ani.to_numpy(float)
    if not np.isfinite(distance).all():
        raise ValueError("ANI distance matrix contains non-finite values after cohort alignment")
    asymmetry = float(np.max(np.abs(distance - distance.T)))
    distance = (distance + distance.T) / 2.0
    np.fill_diagonal(distance, 0.0)

    scores, eigenvalues = pcoa(distance)
    positive = eigenvalues[eigenvalues > 1e-12]
    negative = eigenvalues[eigenvalues < -1e-12]
    explained = positive / positive.sum()
    cumulative = np.cumsum(explained)

    eigen_df = pd.DataFrame({
        "axis": np.arange(1, len(positive) + 1),
        "eigenvalue": positive,
        "explained_variance_fraction_of_positive_eigenvalues": explained,
        "cumulative_explained_variance_fraction": cumulative,
    })
    eigen_df.to_csv(args.output_dir / "ANI_PCoA_positive_eigenvalue_spectrum.csv", index=False)

    rows = []
    for i, count in enumerate(AXIS_COUNTS):
        if count > scores.shape[1]:
            continue
        result = evaluate(y, scores[:, :count], groups, args.permutations, args.seed + i)
        result.update({
            "ANI_PCoA_axes": count,
            "ANI_PCoA_cumulative_variance_fraction": float(cumulative[count - 1]),
        })
        rows.append(result)
    results = pd.DataFrame(rows)
    first = ["ANI_PCoA_axes", "ANI_PCoA_cumulative_variance_fraction"]
    results = results[first + [c for c in results.columns if c not in first]]
    results.to_csv(args.output_dir / "KEGG_ANI_axis_sensitivity_5_10_20_40.csv", index=False)

    qc = {
        "analysis_date": "2026-09-12",
        "ANI_distance_source": str(args.ani_distance),
        "KEGG_feature_source": str(args.kegg_binary),
        "n_genomes": len(ids),
        "n_KEGG_features": len(feature_cols),
        "maximum_pre_symmetrization_asymmetry": asymmetry,
        "positive_eigenvalue_count": int(len(positive)),
        "negative_eigenvalue_count": int(len(negative)),
        "negative_eigenvalue_sum": float(negative.sum()) if len(negative) else 0.0,
        "positive_eigenvalue_sum": float(positive.sum()),
        "negative_to_positive_absolute_sum_ratio": float(abs(negative.sum()) / positive.sum()) if len(negative) else 0.0,
        "axis_counts_tested": list(AXIS_COUNTS),
        "permutations_per_model": args.permutations,
        "permutation_method": "Freedman-Lane permutation of reduced-model residual rows",
    }
    (args.output_dir / "ANI_PCoA_axis_sensitivity_QC.json").write_text(json.dumps(qc, indent=2), encoding="utf-8")

    lines = [
        "# ANI PCoA axis sensitivity analysis",
        "",
        f"The analysis included {len(ids):,} genomes and {len(feature_cols)} ontology-filtered KEGG pathway features.",
        f"The ANI PCoA contained {len(positive)} positive and {len(negative)} negative eigenvalues. The absolute negative-to-positive eigenvalue-sum ratio was {qc['negative_to_positive_absolute_sum_ratio']:.6f}.",
        "",
        "| ANI axes | Cumulative positive-eigenvalue variance | ANI-only R2 | Phylogroup-only R2 | Combined R2 | Unique phylogroup R2 | Permutation p |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in results.itertuples(index=False):
        lines.append(
            f"| {row.ANI_PCoA_axes} | {row.ANI_PCoA_cumulative_variance_fraction:.4f} | "
            f"{row.ANI_only_R2:.4f} | {row.phylogroup_only_R2:.4f} | {row.combined_R2:.4f} | "
            f"{row.unique_phylogroup_R2_after_ANI:.4f} | {row.p_value_freedman_lane:.3f} |"
        )
    (args.output_dir / "ANI_AXIS_SENSITIVITY_REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(results.to_string(index=False))
    print(json.dumps(qc, indent=2))


if __name__ == "__main__":
    main()
