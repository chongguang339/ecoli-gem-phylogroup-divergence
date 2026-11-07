"""Recompute KEGG de-cloning effects after retaining map01053.

GEM rows are copied unchanged because the GEM inputs have not changed.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

ROOT = Path(r"D:\pipeline\analysis_v5")
sys.path.insert(0, str(ROOT / ".python_libs"))
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import spearmanr

OLD = ROOT / "19.mSystems_revision4_20261028"
NEW = OLD / "12_KEGG_RESPONSE_CORRECTION_20261106" / "01_analysis" / "tables"
spec = importlib.util.spec_from_file_location("declone_original", OLD / "scripts" / "analyze_declone_and_kegg.py")
mod = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(mod)

kegg = pd.read_csv(NEW / "kegg_pathway_genome_matrix_bacterial_filtered_binary.csv").set_index("Genome_UID")
features = [c for c in kegg if c.startswith("kegg_pathway_")]
assert len(features) == 116
ani = pd.read_csv(mod.ANI_FILE).set_index("Genome_UID").loc[kegg.index, kegg.index]
ani_array = ani.to_numpy(float)
ani_array = (ani_array + ani_array.T) / 2
np.fill_diagonal(ani_array, 0)
clones = pd.read_csv(mod.CLONE_FILE)
old = pd.read_csv(OLD / "01_analysis" / "decloned_core_effects.csv")
out = old.copy()
genome_to_uid = pd.Series(kegg.index, index=kegg.Genome_ID)

for threshold in (99.0, 99.5, 99.9):
    genomes = clones.loc[
        (clones.ani_threshold_percent == threshold)
        & (clones.is_representative.astype(str).str.lower() == "true"),
        "Genome_ID",
    ].tolist()
    uids = genome_to_uid.loc[genomes].tolist()
    positions = kegg.index.get_indexer(uids)
    labels = kegg.loc[uids, "Phylogroup"].to_numpy(str)
    y = kegg.loc[uids, features].to_numpy(float)
    z = np.sqrt(y / y.sum(axis=1, keepdims=True))
    d = squareform(pdist(z, metric="euclidean"))
    da = ani_array[np.ix_(positions, positions)]
    pairs = np.triu_indices(len(uids), 1)
    r = float(spearmanr(da[pairs], d[pairs]).statistic)
    pcs, eig = mod.pcoa(da)
    res = mod.residual_effect(z, pcs, labels, mod.SEED + int(threshold * 10) + 3)
    mask = (out.ANI_threshold == threshold) & out.profile.eq("KEGG")
    assert mask.sum() == 1
    out.loc[mask, "Mantel_Spearman_r"] = r
    out.loc[mask, "unique_R2"] = res["unique_R2"]
    out.loc[mask, "partial_F"] = res["partial_F"]
    out.loc[mask, "p_value_Freedman_Lane"] = res["p_value_Freedman_Lane"]
    out.loc[mask, "n_variable_features"] = int((y.var(axis=0) > 0).sum())
    print(threshold, r, res["unique_R2"], res["p_value_Freedman_Lane"], flush=True)

out.to_csv(NEW / "decloned_core_effects_corrected.csv", index=False)
