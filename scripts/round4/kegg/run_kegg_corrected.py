from __future__ import annotations

import json
import math
import os
import sys
from pathlib import Path

ROOT = Path(r"D:\pipeline\analysis_v5")
sys.path.insert(0, str(ROOT / ".python_libs"))
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "19.mSystems_revision4_20261028" / "12_KEGG_RESPONSE_CORRECTION_20261106" / ".matplotlib"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import chi2_contingency, spearmanr
from scipy.spatial.distance import pdist, squareform


OUT = ROOT / "19.mSystems_revision4_20261028" / "12_KEGG_RESPONSE_CORRECTION_20261106" / "01_analysis"
TABLES = OUT / "tables"
FIGURES = OUT / "figures"
REPORTS = OUT / "reports"
for directory in [TABLES, FIGURES, REPORTS]:
    directory.mkdir(parents=True, exist_ok=True)

KEGG_MATRIX = ROOT / "09.mSystems_revision_data_freeze" / "02_annotation_matrices" / "kegg_pathway_genome_matrix_binary.csv"
KEGG_DICTIONARY = ROOT / "09.mSystems_revision_data_freeze" / "02_annotation_matrices" / "kegg_pathway_feature_dictionary.csv"
GAPSEQ_KEGG_TBL = ROOT / "gapseq" / "gapseq" / "dat" / "kegg_pwy.tbl"
ANI_DISTANCE = ROOT / "09.mSystems_revision_data_freeze" / "01_master_manifest" / "harmonized_results" / "core_genome_ani_distance_matrix_by_genome_uid.csv.gz"
ANI_PCOA = ROOT / "09.mSystems_revision_data_freeze" / "01_master_manifest" / "harmonized_results" / "03.analysis__ani_pcoa_scores.csv"
DECLONE = ROOT / "13.revision3_biosample_clone_u_memote" / "02_clone_sensitivity" / "ani_declone_cluster_membership.csv"

GROUP_ORDER = ["A", "B2", "B1", "D", "E", "C", "F", "G", "U"]
GROUP_COLORS = {
    "A": "#4E79A7",
    "B2": "#E15759",
    "B1": "#59A14F",
    "D": "#F28E2B",
    "E": "#B07AA1",
    "C": "#76B7B2",
    "F": "#EDC948",
    "G": "#9C755F",
    "U": "#BAB0AC",
}
RNG_SEED = 20260909
PERMUTATIONS = 999

plt.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 7,
        "axes.titlesize": 8,
        "axes.labelsize": 7,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "legend.fontsize": 6,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)
sns.set_style("whitegrid", {"grid.color": "#E7EBEF", "grid.linewidth": 0.45})


def canonical_id(value: object) -> str:
    text = str(value).strip()
    if text.startswith("map") or text.startswith("ko"):
        text = text[3:] if text.startswith("map") else text[2:]
    if text.endswith(".0"):
        text = text[:-2]
    return text.zfill(5)


def read_gapseq_kegg_metadata() -> pd.DataFrame:
    meta = pd.read_csv(GAPSEQ_KEGG_TBL, sep="\t", dtype=str)
    meta["Canonical_Pathway_ID"] = meta["id"].map(canonical_id)
    hierarchy_parts = meta["hierarchy"].fillna("").str.split(";", expand=True)
    for idx, col in enumerate(["kegg_root", "level1", "level2", "level3"]):
        if idx < hierarchy_parts.shape[1]:
            meta[col] = hierarchy_parts[idx].fillna("")
        else:
            meta[col] = ""
    keep_cols = ["Canonical_Pathway_ID", "id", "name", "hierarchy", "level1", "level2", "level3", "taxrange", "reaNr", "ecNr"]
    return meta[keep_cols].drop_duplicates("Canonical_Pathway_ID")


def classify_feature(row: pd.Series) -> tuple[bool, str, str]:
    level1 = str(row.get("level1", ""))
    level2 = str(row.get("level2", ""))
    name = str(row.get("Description", "") or row.get("name", ""))
    hierarchy = str(row.get("hierarchy", ""))
    map_id = "map" + canonical_id(row.get("Canonical_Pathway_ID", ""))
    text = f"{name} {hierarchy}".lower()

    explicit_bacterial_allowlist = {
        "map01053",  # siderophore biosynthesis; includes enterobactin-relevant reactions
        "map00130",  # ubiquinone and other terpenoid-quinone biosynthesis
        "map00900",  # terpenoid backbone biosynthesis
        "map02010",  # ABC transporters
        "map02020",  # two-component system
        "map02024",  # quorum sensing
        "map02025",  # biofilm formation - Pseudomonas aeruginosa
        "map02026",  # biofilm formation - Escherichia coli
        "map02030",  # bacterial chemotaxis
        "map02040",  # flagellar assembly
        "map02060",  # phosphotransferase system
        "map03010",  # ribosome
        "map03018",  # RNA degradation
        "map03020",  # RNA polymerase
        "map03030",  # DNA replication
        "map03060",  # protein export
        "map03070",  # bacterial secretion system
        "map03410",  # base excision repair
        "map03420",  # nucleotide excision repair
        "map03430",  # mismatch repair
        "map03440",  # homologous recombination
        "map04122",  # sulfur relay system
        "map05111",  # biofilm formation - Vibrio cholerae
    }
    if map_id in explicit_bacterial_allowlist:
        return True, "retained_explicit_bacterial_allowlist", "Explicitly retained bacterial/prokaryotic KEGG feature."

    if level1 in {"Human Diseases", "Organismal Systems"}:
        return False, "excluded_host_or_organismal_system", f"KEGG category {level1} is outside bacterial pathway profiling."

    euk_or_host_terms = [
        "cancer",
        "carcinoma",
        "leukemia",
        "diabetes",
        "addiction",
        "alzheimer",
        "parkinson",
        "huntington",
        "cardiac",
        "myocardial",
        "adrenergic",
        "dopaminergic",
        "cholinergic",
        "serotonergic",
        "synaptic",
        "axon",
        "neuro",
        "immune",
        "endocrine",
        "hormone",
        "oocyte",
        "meiosis",
        "fanconi",
        "apoptosis",
        "platelet",
        "osteoclast",
        "insulin",
        "estrogen",
        "progesterone",
        "steroid hormone",
        "plant hormone",
        "brassinosteroid",
        "steroid biosynthesis",
        "primary bile acid",
        "secondary bile acid",
        "arachidonic acid",
        "linoleic acid",
        "alpha-linolenic acid",
        "sphingolipid",
        "glycosphingolipid",
        "glycosaminoglycan",
        "n-glycan",
        "retinol",
        "zeatin",
        "carotenoid biosynthesis",
        "flavonoid biosynthesis",
        "flavone",
        "anthocyanin",
        "stilbenoid",
        "isoquinoline alkaloid",
        "tropane",
        "glucosinolate",
        "phenylpropanoid biosynthesis",
        "betalain",
        "photosynthesis",
        "porphyrin and chlorophyll",
        "carbon fixation in photosynthetic organisms",
        "tetracycline biosynthesis",
        "biosynthesis of antibiotics",
        "nonribosomal peptide",
        "enediyne antibiotics",
        "vancomycin group antibiotics",
        "ansamycins",
        "mapk signaling",
        "nf-kappa b",
        "hif-1",
        "foxo",
        "phosphatidylinositol signaling",
        "phospholipase d signaling",
        "mtor",
        "pi3k-akt",
        "ampk signaling",
        "calcium signaling",
        "ribosome biogenesis in eukaryotes",
        "nucleocytoplasmic transport",
        "spliceosome",
        "endoplasmic reticulum",
    ]
    for term in euk_or_host_terms:
        if term in text:
            return False, "excluded_non_bacterial_or_host_specific_keyword", f"Feature name/hierarchy contains '{term}'."

    excluded_metabolic_subclasses = {
        "Biosynthesis of other secondary metabolites",
    }
    if level2 in excluded_metabolic_subclasses:
        return False, "excluded_secondary_metabolite_biosynthesis_without_ecoli_specificity", f"KEGG metabolic subcategory {level2} was not used for E. coli pathway profiling."

    if level2 == "Metabolism of terpenoids and polyketides":
        return False, "excluded_specialized_terpenoid_polyketide_pathway", "Specialized terpenoid/polyketide features were excluded unless explicitly retained as core quinone/isoprenoid metabolism."

    if level2 == "Signal transduction":
        return False, "excluded_eukaryote_or_host_signal_transduction", "Only explicitly bacterial signal-transduction features were retained."

    if level1 in {"Metabolism", "Genetic Information Processing", "Environmental Information Processing"}:
        return True, "retained_bacterial_relevant_core_category", f"Retained KEGG category: {level1}; subcategory: {level2}."

    if level1 == "Cellular Processes":
        bacterial_terms = ["prokaryotic", "bacterial", "biofilm", "flagellar", "quorum", "two-component", "secretion"]
        if any(term in text for term in bacterial_terms):
            return True, "retained_bacterial_cellular_process", "Cellular-process term has explicit bacterial/prokaryotic relevance."
        return False, "excluded_cellular_process_without_bacterial_specificity", "Cellular-process category lacks explicit bacterial/prokaryotic specificity."

    if not level1:
        return False, "excluded_missing_ontology_metadata", "No KEGG hierarchy metadata available in local gapseq pathway table."
    return False, "excluded_other_kegg_category", f"KEGG category {level1} not used for bacterial pathway profiling."


def build_filter_table() -> pd.DataFrame:
    dictionary = pd.read_csv(KEGG_DICTIONARY, dtype=str)
    dictionary["Canonical_Pathway_ID"] = dictionary["Canonical_Pathway_ID"].map(canonical_id)
    meta = read_gapseq_kegg_metadata()
    merged = dictionary.merge(meta, on="Canonical_Pathway_ID", how="left", suffixes=("", "_gapseq"))
    if "name" in merged.columns:
        merged["Description"] = merged["Description"].fillna("")
        merged.loc[merged["Description"].eq(""), "Description"] = merged.loc[merged["Description"].eq(""), "name"].fillna("")
    decisions = merged.apply(classify_feature, axis=1, result_type="expand")
    decisions.columns = ["Retain_For_Bacterial_Profile", "Ontology_Filter_Class", "Ontology_Filter_Reason"]
    merged = pd.concat([merged, decisions], axis=1)
    merged["KEGG_Map_ID"] = "map" + merged["Canonical_Pathway_ID"]
    columns = [
        "Matrix_Feature",
        "Canonical_Pathway_ID",
        "KEGG_Map_ID",
        "Description",
        "hierarchy",
        "level1",
        "level2",
        "level3",
        "Retain_For_Bacterial_Profile",
        "Ontology_Filter_Class",
        "Ontology_Filter_Reason",
    ]
    merged[columns].to_csv(TABLES / "kegg_ontology_filter_decision_table.csv", index=False, encoding="utf-8-sig")
    return merged


def hellinger_transform(x: np.ndarray) -> np.ndarray:
    x = x.astype(float)
    row_sums = x.sum(axis=1)
    row_sums[row_sums == 0] = 1.0
    return np.sqrt(x / row_sums[:, None])


def pairwise_euclidean(x: np.ndarray) -> np.ndarray:
    return squareform(pdist(x, metric="euclidean"))


def pcoa(distance: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n = distance.shape[0]
    d2 = distance ** 2
    j = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * j @ d2 @ j
    vals, vecs = np.linalg.eigh(b)
    order = vals.argsort()[::-1]
    vals = vals[order]
    vecs = vecs[:, order]
    positive = vals > 1e-12
    coords = vecs[:, positive] * np.sqrt(vals[positive])
    return coords, vals


def pca_scores(x: np.ndarray, n_components: int = 10) -> tuple[np.ndarray, np.ndarray]:
    centered = x - x.mean(axis=0, keepdims=True)
    u, s, vt = np.linalg.svd(centered, full_matrices=False)
    scores = u[:, :n_components] * s[:n_components]
    explained = (s ** 2) / np.sum(s ** 2)
    return scores, explained[:n_components]


def permanova_oneway(distance: np.ndarray, labels: np.ndarray, permutations: int = PERMUTATIONS, seed: int = RNG_SEED) -> dict:
    labels = np.asarray(labels)
    unique = pd.unique(labels)
    n = len(labels)
    k = len(unique)
    d2 = distance ** 2

    def ss_values(lbls: np.ndarray) -> tuple[float, float, float, float]:
        total = d2[np.triu_indices(n, 1)].sum() / n
        within = 0.0
        for group in pd.unique(lbls):
            idx = np.where(lbls == group)[0]
            if len(idx) > 1:
                sub = d2[np.ix_(idx, idx)]
                within += sub[np.triu_indices(len(idx), 1)].sum() / len(idx)
        between = total - within
        f_value = (between / (k - 1)) / (within / (n - k))
        r2 = between / total if total else math.nan
        return total, within, between, f_value

    total, within, between, f_obs = ss_values(labels)
    rng = np.random.default_rng(seed)
    perm_f = np.empty(permutations)
    for i in range(permutations):
        perm_f[i] = ss_values(rng.permutation(labels))[3]
    p_value = (np.sum(perm_f >= f_obs) + 1) / (permutations + 1)
    return {
        "test": "PERMANOVA_oneway_phylogroup",
        "n_genomes": n,
        "n_groups": k,
        "df_between": k - 1,
        "df_within": n - k,
        "SS_total": total,
        "SS_within": within,
        "SS_between": between,
        "pseudo_F": f_obs,
        "R2": between / total,
        "p_value": p_value,
        "permutations": permutations,
    }


def permdisp(distance: np.ndarray, labels: np.ndarray, ids: list[str], permutations: int = PERMUTATIONS, seed: int = RNG_SEED) -> tuple[pd.DataFrame, dict]:
    coords, _ = pcoa(distance)
    labels = np.asarray(labels)

    def distances_to_centroid(lbls: np.ndarray) -> np.ndarray:
        out = np.zeros(len(lbls))
        for group in pd.unique(lbls):
            idx = np.where(lbls == group)[0]
            centroid = coords[idx].mean(axis=0)
            out[idx] = np.linalg.norm(coords[idx] - centroid, axis=1)
        return out

    def f_stat(values: np.ndarray, lbls: np.ndarray) -> float:
        grand = values.mean()
        ss_between = sum(((values[lbls == group].mean() - grand) ** 2) * (lbls == group).sum() for group in pd.unique(lbls))
        ss_within = sum(((values[lbls == group] - values[lbls == group].mean()) ** 2).sum() for group in pd.unique(lbls))
        df_between = len(pd.unique(lbls)) - 1
        df_within = len(lbls) - len(pd.unique(lbls))
        return (ss_between / df_between) / (ss_within / df_within)

    obs_dist = distances_to_centroid(labels)
    f_obs = f_stat(obs_dist, labels)
    rng = np.random.default_rng(seed)
    perm_f = np.empty(permutations)
    for i in range(permutations):
        perm_labels = rng.permutation(labels)
        perm_f[i] = f_stat(distances_to_centroid(perm_labels), perm_labels)
    p_value = (np.sum(perm_f >= f_obs) + 1) / (permutations + 1)
    dist_df = pd.DataFrame({"Genome_UID": ids, "Phylogroup": labels, "distance_to_group_centroid": obs_dist})
    summary = {
        "test": "PERMDISP_phylogroup",
        "n_genomes": len(ids),
        "n_groups": len(pd.unique(labels)),
        "pseudo_F": f_obs,
        "p_value": p_value,
        "permutations": permutations,
    }
    return dist_df, summary


def mantel_spearman(distance_a: np.ndarray, distance_b: np.ndarray, permutations: int = PERMUTATIONS, seed: int = RNG_SEED) -> dict:
    iu = np.triu_indices(distance_a.shape[0], 1)
    a = distance_a[iu]
    b = distance_b[iu]
    obs = spearmanr(a, b).statistic
    rng = np.random.default_rng(seed)
    perm_stats = np.empty(permutations)
    for i in range(permutations):
        perm = rng.permutation(distance_b.shape[0])
        perm_b = distance_b[np.ix_(perm, perm)]
        perm_stats[i] = spearmanr(a, perm_b[iu]).statistic
    p_value = (np.sum(np.abs(perm_stats) >= abs(obs)) + 1) / (permutations + 1)
    return {
        "test": "Mantel_Spearman_KEGG_Hellinger_vs_FastANI_distance",
        "n_genomes": distance_a.shape[0],
        "mantel_r": obs,
        "p_value_two_sided": p_value,
        "permutations": permutations,
    }


def design_matrix(ani_pcoa: pd.DataFrame, groups: pd.Series, ids: list[str], include_group: bool) -> np.ndarray:
    pcoa_cols = [f"ANI_PC{i}" for i in range(1, 21)]
    parts = [np.ones((len(ids), 1)), ani_pcoa.set_index("Genome_UID").loc[ids, pcoa_cols].to_numpy(float)]
    if include_group:
        dummies = pd.get_dummies(groups.loc[ids].astype(str))
        for group in GROUP_ORDER[1:]:
            if group in dummies.columns:
                parts.append(dummies[[group]].to_numpy(float))
            else:
                parts.append(np.zeros((len(ids), 1)))
    return np.hstack(parts)


def residual_phylogroup_effect(y: np.ndarray, ani_pcoa: pd.DataFrame, groups: pd.Series, ids: list[str], permutations: int = PERMUTATIONS, seed: int = RNG_SEED) -> dict:
    x_reduced = design_matrix(ani_pcoa, groups, ids, include_group=False)
    x_full = design_matrix(ani_pcoa, groups, ids, include_group=True)

    def sse(x: np.ndarray, target: np.ndarray) -> float:
        beta, *_ = np.linalg.lstsq(x, target, rcond=None)
        resid = target - x @ beta
        return float(np.sum(resid ** 2))

    y_centered = y - y.mean(axis=0, keepdims=True)
    sst = float(np.sum(y_centered ** 2))
    sse_reduced = sse(x_reduced, y)
    sse_full = sse(x_full, y)
    df_num = x_full.shape[1] - x_reduced.shape[1]
    df_den = len(ids) - x_full.shape[1]
    f_obs = ((sse_reduced - sse_full) / df_num) / (sse_full / df_den)
    unique_r2 = (sse_reduced - sse_full) / sst

    beta_red, *_ = np.linalg.lstsq(x_reduced, y, rcond=None)
    fitted_red = x_reduced @ beta_red
    resid_red = y - fitted_red
    rng = np.random.default_rng(seed)
    perm_f = np.empty(permutations)
    for i in range(permutations):
        y_perm = fitted_red + resid_red[rng.permutation(len(ids)), :]
        sr = sse(x_reduced, y_perm)
        sf = sse(x_full, y_perm)
        perm_f[i] = ((sr - sf) / df_num) / (sf / df_den)
    p_value = (np.sum(perm_f >= f_obs) + 1) / (permutations + 1)
    return {
        "test": "ANI20_conditioned_unique_phylogroup_effect_multivariate_lm",
        "n_genomes": len(ids),
        "n_features": y.shape[1],
        "df_num": df_num,
        "df_den": df_den,
        "unique_R2": unique_r2,
        "partial_F": f_obs,
        "p_value": p_value,
        "permutations": permutations,
    }


def differential_features(matrix: pd.DataFrame, feature_cols: list[str], filter_table: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for feature in feature_cols:
        table = []
        for group in GROUP_ORDER:
            vals = matrix.loc[matrix["Phylogroup"] == group, feature].astype(int)
            table.append([int(vals.sum()), int((1 - vals).sum())])
        contingency = np.asarray(table)
        if contingency[:, 0].sum() == 0 or contingency[:, 1].sum() == 0:
            chi2, p_value = 0.0, 1.0
            test_status = "not_tested_invariant_feature"
        else:
            try:
                chi2, p_value, _, _ = chi2_contingency(contingency, correction=False)
                test_status = "tested_chi_square"
            except ValueError:
                chi2, p_value = np.nan, 1.0
                test_status = "not_tested_zero_expected_frequency"
        n = sum(sum(row) for row in table)
        r = len(table)
        c = 2
        cramer = math.sqrt(chi2 / (n * (min(r - 1, c - 1)))) if n and np.isfinite(chi2) else math.nan
        freqs = {group: table[i][0] / sum(table[i]) if sum(table[i]) else math.nan for i, group in enumerate(GROUP_ORDER)}
        best_group = max(freqs, key=freqs.get)
        worst_group = min(freqs, key=freqs.get)
        rows.append(
            {
                "Matrix_Feature": feature,
                "test_status": test_status,
                "chi_square": chi2,
                "p_value": p_value,
                "Cramers_V": cramer,
                "max_frequency_group": best_group,
                "max_frequency": freqs[best_group],
                "min_frequency_group": worst_group,
                "min_frequency": freqs[worst_group],
                "frequency_range": freqs[best_group] - freqs[worst_group],
                **{f"frequency_{group}": freqs[group] for group in GROUP_ORDER},
            }
        )
    out = pd.DataFrame(rows).sort_values("p_value").reset_index(drop=True)
    m = len(out)
    ranks = np.arange(1, m + 1)
    bh = out["p_value"].to_numpy(float) * m / ranks
    out["BH_q_value"] = np.minimum.accumulate(bh[::-1])[::-1].clip(max=1.0)
    out["Bonferroni_p"] = (out["p_value"] * m).clip(upper=1.0)
    annot = filter_table[["Matrix_Feature", "KEGG_Map_ID", "Description", "hierarchy", "Ontology_Filter_Class"]]
    return out.merge(annot, on="Matrix_Feature", how="left")


def balanced_subset(ids: list[str], groups: pd.Series, n_per_group: int = 4) -> list[str]:
    rng = np.random.default_rng(RNG_SEED)
    selected = []
    for group in GROUP_ORDER:
        group_ids = [x for x in ids if groups.loc[x] == group]
        if len(group_ids) < n_per_group:
            continue
        chosen = sorted(rng.choice(group_ids, size=n_per_group, replace=False).tolist())
        selected.extend(chosen)
    return selected


def plot_filtered_outputs(matrix: pd.DataFrame, feature_cols: list[str], scores: pd.DataFrame, explained: np.ndarray, permdisp_df: pd.DataFrame, diff: pd.DataFrame, permanova_summary: dict, mantel_summary: dict) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(8.2, 6.2))
    fig.subplots_adjust(left=0.09, right=0.98, bottom=0.10, top=0.93, hspace=0.46, wspace=0.58)

    ax = axes[0, 0]
    for group in GROUP_ORDER:
        sub = scores[scores["Phylogroup"] == group]
        ax.scatter(sub["PC1"], sub["PC2"], s=10, alpha=0.70, color=GROUP_COLORS[group], label=group, linewidths=0)
    ax.set_title("A  Filtered KEGG Hellinger PCA", loc="left", fontweight="bold")
    ax.set_xlabel(f"PC1 ({explained[0] * 100:.1f}%)")
    ax.set_ylabel(f"PC2 ({explained[1] * 100:.1f}%)")
    ax.legend(ncol=3, frameon=False, loc="best", handletextpad=0.2, columnspacing=0.5)
    sns.despine(ax=ax)

    ax = axes[0, 1]
    top = (
        diff[diff["test_status"].eq("tested_chi_square")]
        .sort_values(["BH_q_value", "Cramers_V"], ascending=[True, False])
        .head(8)
    )
    heat = []
    labels = []
    for _, row in top.iterrows():
        heat.append([row[f"frequency_{group}"] for group in GROUP_ORDER])
        desc = str(row["Description"])
        if len(desc) > 26:
            desc = desc[:25] + "..."
        labels.append(f"{row['KEGG_Map_ID']} {desc}")
    sns.heatmap(
        pd.DataFrame(heat, index=labels, columns=GROUP_ORDER),
        cmap=sns.color_palette("Blues", as_cmap=True),
        vmin=0,
        vmax=1,
        linewidths=0.25,
        linecolor="white",
        cbar_kws={"label": "Genome frequency", "shrink": 0.75},
        ax=ax,
    )
    ax.set_title("B  Top retained bacterial-relevant features", loc="left", fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="y", rotation=0, labelsize=5.5)
    ax.tick_params(axis="x", rotation=0)

    ax = axes[1, 0]
    sns.boxplot(data=permdisp_df, x="Phylogroup", y="distance_to_group_centroid", order=GROUP_ORDER, color="#8CB7CF", fliersize=1.4, linewidth=0.7, ax=ax)
    sns.stripplot(data=permdisp_df, x="Phylogroup", y="distance_to_group_centroid", order=GROUP_ORDER, color="#1F3A4D", size=1.2, jitter=0.18, alpha=0.45, ax=ax)
    ax.set_title("C  PERMDISP after ontology filtering", loc="left", fontweight="bold")
    ax.set_xlabel("Phylogroup")
    ax.set_ylabel("Distance to group centroid")
    sns.despine(ax=ax)

    ax = axes[1, 1]
    residual = pd.read_csv(TABLES / "kegg_filtered_global_statistics.csv")
    residual = residual.loc[residual["test"].eq("ANI20_conditioned_unique_phylogroup_effect_multivariate_lm"), "unique_R2"].iloc[0]
    bars = pd.DataFrame(
        [
            {"Metric": "PERMANOVA\nphylogroup R2", "Value": permanova_summary["R2"]},
            {"Metric": "Mantel\nANI r", "Value": mantel_summary["mantel_r"]},
            {"Metric": "ANI-corrected\nunique R2", "Value": residual},
        ]
    )
    ax.bar(bars["Metric"], bars["Value"], color=["#7AA6C2", "#76A889", "#D6A85B"], width=0.55)
    for i, value in enumerate(bars["Value"]):
        ax.text(i, value + 0.018, f"{value:.3f}", ha="center", va="bottom", fontsize=7)
    ax.set_ylim(0, max(0.85, bars["Value"].max() + 0.12))
    ax.set_ylabel("Statistic")
    ax.set_title("D  Global structure retained after filtering", loc="left", fontweight="bold")
    sns.despine(ax=ax)

    fig.savefig(FIGURES / "Figure3_KEGG_filtered_revision3.png", dpi=600, bbox_inches="tight", pad_inches=0.05)
    fig.savefig(FIGURES / "Figure3_KEGG_filtered_revision3.pdf", bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)


def main() -> None:
    filter_table = build_filter_table()
    matrix = pd.read_csv(KEGG_MATRIX)
    feature_cols_all = [c for c in matrix.columns if c.startswith("kegg_pathway_")]
    retained_features = filter_table.loc[filter_table["Retain_For_Bacterial_Profile"], "Matrix_Feature"].tolist()
    removed_features = [c for c in feature_cols_all if c not in set(retained_features)]
    filtered = matrix[["Genome_UID", "Genome_ID", "Phylogroup"] + retained_features].copy()
    filtered.to_csv(TABLES / "kegg_pathway_genome_matrix_bacterial_filtered_binary.csv", index=False, encoding="utf-8-sig")
    matrix[["Genome_UID", "Genome_ID", "Phylogroup"] + removed_features].to_csv(TABLES / "kegg_pathway_genome_matrix_removed_features_binary.csv", index=False, encoding="utf-8-sig")

    ids = filtered["Genome_UID"].astype(str).tolist()
    groups = filtered.set_index("Genome_UID")["Phylogroup"].astype(str)
    x = filtered[retained_features].to_numpy(float)
    xh = hellinger_transform(x)
    scores, explained = pca_scores(xh, n_components=10)
    score_df = pd.DataFrame(scores[:, :5], columns=[f"PC{i}" for i in range(1, 6)])
    score_df.insert(0, "Phylogroup", filtered["Phylogroup"].values)
    score_df.insert(0, "Genome_ID", filtered["Genome_ID"].values)
    score_df.insert(0, "Genome_UID", filtered["Genome_UID"].values)
    for i, ev in enumerate(explained[:5], start=1):
        score_df[f"PC{i}_explained_variance"] = ev
    score_df.to_csv(TABLES / "kegg_filtered_hellinger_pca_scores.csv", index=False, encoding="utf-8-sig")

    kegg_dist = pairwise_euclidean(xh)
    permanova_summary = permanova_oneway(kegg_dist, filtered["Phylogroup"].astype(str).to_numpy())
    permdisp_df, permdisp_summary = permdisp(kegg_dist, filtered["Phylogroup"].astype(str).to_numpy(), ids)
    permdisp_df.to_csv(TABLES / "kegg_filtered_permdisp_distances.csv", index=False, encoding="utf-8-sig")

    ani = pd.read_csv(ANI_DISTANCE).rename(columns={"Genome_UID": "Genome_UID"}).set_index("Genome_UID")
    ani = ani.loc[ids, ids].apply(pd.to_numeric, errors="coerce")
    ani_arr = ani.to_numpy(float)
    ani_arr = (ani_arr + ani_arr.T) / 2
    np.fill_diagonal(ani_arr, 0.0)
    mantel_summary = mantel_spearman(kegg_dist, ani_arr)

    ani_pcoa = pd.read_csv(ANI_PCOA)
    residual_summary = residual_phylogroup_effect(xh, ani_pcoa, groups, ids)

    diff = differential_features(filtered, retained_features, filter_table)
    diff.to_csv(TABLES / "kegg_filtered_phylogroup_differential_features.csv", index=False, encoding="utf-8-sig")

    rows = [permanova_summary, permdisp_summary, mantel_summary, residual_summary]
    pd.DataFrame(rows).to_csv(TABLES / "kegg_filtered_global_statistics.csv", index=False, encoding="utf-8-sig")

    balanced_ids = balanced_subset(ids, groups, n_per_group=4)
    if balanced_ids:
        idx = [ids.index(gid) for gid in balanced_ids]
        bal_dist = kegg_dist[np.ix_(idx, idx)]
        bal_groups = filtered.set_index("Genome_UID").loc[balanced_ids, "Phylogroup"].astype(str).to_numpy()
        bal_permanova = permanova_oneway(bal_dist, bal_groups, seed=RNG_SEED + 1)
        bal_permdisp_df, bal_permdisp = permdisp(bal_dist, bal_groups, balanced_ids, seed=RNG_SEED + 1)
        bal_permanova["scenario"] = "balanced_n4_per_group"
        bal_permdisp["scenario"] = "balanced_n4_per_group"
        pd.DataFrame([bal_permanova, bal_permdisp]).to_csv(TABLES / "kegg_filtered_balanced_n4_statistics.csv", index=False, encoding="utf-8-sig")
        bal_permdisp_df.to_csv(TABLES / "kegg_filtered_balanced_n4_permdisp_distances.csv", index=False, encoding="utf-8-sig")

    if DECLONE.exists():
        declone = pd.read_csv(DECLONE)
        sensitivity_rows = []
        distance_rows = []
        for threshold in [99.0, 99.5, 99.9]:
            sub = declone[(declone["ani_threshold_percent"] == threshold) & (declone["is_representative"] == True)]
            for include_u in [True, False]:
                sub_ids = [gid for gid in sub["Genome_ID"].astype(str) if gid in set(filtered["Genome_ID"])]
                genome_to_uid = filtered.set_index("Genome_ID")["Genome_UID"].to_dict()
                sub_uids = [genome_to_uid[gid] for gid in sub_ids]
                if not include_u:
                    sub_uids = [uid for uid in sub_uids if groups.loc[uid] != "U"]
                if len(sub_uids) < 10:
                    continue
                idx = [ids.index(uid) for uid in sub_uids]
                sub_dist = kegg_dist[np.ix_(idx, idx)]
                sub_groups = groups.loc[sub_uids].to_numpy()
                pm = permanova_oneway(sub_dist, sub_groups, seed=RNG_SEED + int(threshold * 10) + (0 if include_u else 1))
                pdist, pdsum = permdisp(sub_dist, sub_groups, sub_uids, seed=RNG_SEED + int(threshold * 10) + (10 if include_u else 11))
                pm.update({"analysis": "PERMANOVA", "declone_ani_threshold_percent": threshold, "include_U": include_u})
                pdsum.update({"analysis": "PERMDISP", "declone_ani_threshold_percent": threshold, "include_U": include_u})
                sensitivity_rows.extend([pm, pdsum])
                pdist["declone_ani_threshold_percent"] = threshold
                pdist["include_U"] = include_u
                distance_rows.append(pdist)
        pd.DataFrame(sensitivity_rows).to_csv(TABLES / "kegg_filtered_declone_and_U_sensitivity_statistics.csv", index=False, encoding="utf-8-sig")
        if distance_rows:
            pd.concat(distance_rows, ignore_index=True).to_csv(TABLES / "kegg_filtered_declone_and_U_permdisp_distances.csv", index=False, encoding="utf-8-sig")

    plot_filtered_outputs(filtered, retained_features, score_df, explained, permdisp_df, diff, permanova_summary, mantel_summary)

    class_counts = filter_table.groupby(["Retain_For_Bacterial_Profile", "Ontology_Filter_Class"]).size().reset_index(name="n_features")
    class_counts.to_csv(TABLES / "kegg_ontology_filter_summary_by_class.csv", index=False, encoding="utf-8-sig")
    summary = {
        "input_features": len(feature_cols_all),
        "retained_bacterial_relevant_features": len(retained_features),
        "removed_features": len(removed_features),
        "n_genomes": len(ids),
        "permutations": PERMUTATIONS,
        "pca_explained_pc1": float(explained[0]),
        "pca_explained_pc2": float(explained[1]),
        "permanova_R2": float(permanova_summary["R2"]),
        "permanova_p": float(permanova_summary["p_value"]),
        "permdisp_F": float(permdisp_summary["pseudo_F"]),
        "permdisp_p": float(permdisp_summary["p_value"]),
        "mantel_spearman_r": float(mantel_summary["mantel_r"]),
        "mantel_p": float(mantel_summary["p_value_two_sided"]),
        "ANI20_unique_phylogroup_R2": float(residual_summary["unique_R2"]),
        "ANI20_unique_phylogroup_p": float(residual_summary["p_value"]),
        "significant_features_BH_0_05": int((diff["BH_q_value"] < 0.05).sum()),
        "significant_features_bonferroni_0_05": int((diff["Bonferroni_p"] < 0.05).sum()),
    }
    (REPORTS / "kegg_filtered_analysis_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    readme = [
        "# KEGG Ontology Filtering and Reanalysis",
        "",
        "This folder contains the revision-specific KEGG ontology filtering and complete rerun of KEGG statistics.",
        "",
        f"- Input features: {summary['input_features']}",
        f"- Retained bacterial-relevant features: {summary['retained_bacterial_relevant_features']}",
        f"- Removed non-bacterial, eukaryote-specific, host/organismal-system, or unsupported features: {summary['removed_features']}",
        f"- PERMANOVA R2 after filtering: {summary['permanova_R2']:.4f}, p = {summary['permanova_p']:.3f}",
        f"- PERMDISP F after filtering: {summary['permdisp_F']:.4f}, p = {summary['permdisp_p']:.3f}",
        f"- FastANI Mantel Spearman r after filtering: {summary['mantel_spearman_r']:.4f}, p = {summary['mantel_p']:.3f}",
        f"- ANI20-conditioned unique phylogroup R2 after filtering: {summary['ANI20_unique_phylogroup_R2']:.4f}, p = {summary['ANI20_unique_phylogroup_p']:.3f}",
        "",
        "Filtering decisions are in `tables/kegg_ontology_filter_decision_table.csv`; this table is intended for supplementary inclusion.",
    ]
    (REPORTS / "README_KEGG_FILTERED_REANALYSIS.md").write_text("\n".join(readme), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
