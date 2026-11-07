import csv
import json
from pathlib import Path

import numpy as np


def squareform(v):
    n = int((1 + np.sqrt(1 + 8 * len(v))) / 2)
    a = np.zeros((n, n), dtype=float)
    i = np.triu_indices(n, 1)
    a[i] = v
    return a + a.T

ROOT = Path(r"D:\pipeline\analysis_v5\19.mSystems_revision4_20261028")
DATA = ROOT / "01_analysis" / "server_pair_summary_20261101" / "summary"
sample = list(csv.DictReader((DATA / "paired_reaction_counts_180.tsv").open(), delimiter="\t"))
labels = np.array([r["Phylogroup"] for r in sample])
z = np.load(DATA / "paired_network_distances_180.npz")
n = len(labels)
ani = squareform(z["ani"])
h = np.eye(n) - np.ones((n, n)) / n
eig, vec = np.linalg.eigh(-0.5 * h @ (ani**2) @ h)
order = np.argsort(eig)[::-1]
pcs = vec[:, order[:20]] * np.sqrt(eig[order[:20]])
pair = np.tril_indices(n, -1)
x = np.column_stack([np.ones(len(pair[0]))] + [np.abs(pcs[pair[0], k] - pcs[pair[1], k]) for k in range(20)])
residual = {}
for name in ("broad_all", "focused_all", "broad_gene_evidenced", "focused_gene_evidenced"):
    target = squareform(z[name])[pair]
    beta = np.linalg.lstsq(x, target, rcond=None)[0]
    a = np.zeros((n, n))
    a[pair] = target - x @ beta
    residual[name] = a + a.T

def effect(mat, ix):
    lo = np.tril_indices(len(ix), -1)
    keep = ix[lo[0]] != ix[lo[1]]
    a, b = ix[lo[0]][keep], ix[lo[1]][keep]
    y = mat[a, b]
    same = labels[a] == labels[b]
    mean = y.mean()
    ss = ((y - mean) ** 2).sum()
    return (same.sum() * (y[same].mean() - mean) ** 2 + (~same).sum() * (y[~same].mean() - mean) ** 2) / ss

rng = np.random.default_rng(20261101)
groups = [np.where(labels == g)[0] for g in sorted(set(labels))]
draws = {"all": [], "gene_evidenced": []}
for _ in range(1000):
    ix = np.concatenate([rng.choice(g, size=len(g), replace=True) for g in groups])
    for key in draws:
        draws[key].append(effect(residual["broad_" + key], ix) - effect(residual["focused_" + key], ix))
result = {}
for key, vals in draws.items():
    v = np.array(vals)
    result[key] = {
        "observed_difference": effect(residual["broad_" + key], np.arange(n)) - effect(residual["focused_" + key], np.arange(n)),
        "bootstrap_median": float(np.median(v)),
        "bootstrap_percentile_2_5": float(np.quantile(v, 0.025)),
        "bootstrap_percentile_97_5": float(np.quantile(v, 0.975)),
        "bootstrap_draws": len(v),
    }
out = DATA / "paired_lineage_effect_bootstrap.json"
out.write_text(json.dumps(result, indent=2))
print(json.dumps(result, indent=2))
