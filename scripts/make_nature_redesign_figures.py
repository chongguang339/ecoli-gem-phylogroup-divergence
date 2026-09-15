from __future__ import annotations

import math
import os
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(os.environ.get("ECOLI_GEM_ROOT", Path(__file__).resolve().parents[1]))
HERE = Path(os.environ.get("ECOLI_FIGURE_WORKDIR", ROOT / "figure_workdir"))
sys.path.insert(0, str(HERE))

import fitz
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm
from matplotlib.ticker import FuncFormatter
from PIL import Image, ImageChops
from scipy.spatial.distance import pdist, squareform
from scipy.stats import spearmanr

from audit_panel_alignment import require_matplotlib_panel_alignment

PANELS = HERE / "single_panels"
ASSEMBLED = HERE / "assembled_submission"
QA = HERE / "qa"
ASSETS = HERE / "assets"
SOURCE_DATA = HERE / "source_data"
for p in (PANELS, ASSEMBLED, QA, ASSETS, SOURCE_DATA):
    p.mkdir(parents=True, exist_ok=True)

ORDER = ["A", "B2", "B1", "D", "E", "C", "F", "G", "U"]
COLORS = {"A":"#4477AA", "B2":"#EE6677", "B1":"#228833", "D":"#CCBB44",
          "E":"#66CCEE", "C":"#AA3377", "F":"#BBBBBB", "G":"#EE7733", "U":"#111111"}
NEUTRAL = "#6F7882"
BLUE = "#0072B2"
ORANGE = "#D55E00"
GREEN = "#009E73"

plt.rcParams.update({
    "font.family": "Arial", "font.size": 7, "axes.titlesize": 8,
    "axes.labelsize": 7.5, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "legend.fontsize": 7, "pdf.fonttype": 42, "ps.fonttype": 42,
    "svg.fonttype": "none", "axes.linewidth": 0.65,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.8, "ytick.major.size": 2.8,
    "axes.titleweight": "bold", "figure.facecolor": "white", "savefig.facecolor": "white",
})


def read(rel):
    return pd.read_csv(ROOT / rel)


RUNS = read("18.mSystems_revision3_final_20260911/02_analysis_evidence/01_ablation/benchmark_20x6_authoritative_run_table.csv")
ASSEMBLY_DF = read("00.data/genome_assembly_and_cog_metrics.csv")
QC_DF = read("00.data/cobrapy_qc_remote_verified_models.csv")
EVID = read("10.mSystems_revision2_model_evidence/02_model_evidence/reaction_evidence_counts_by_genome.csv")
KEGG = read("15.kegg_ontology_filtering_revision3/tables/kegg_pathway_genome_matrix_bacterial_filtered_binary.csv")
KEGG_PCA = read("15.kegg_ontology_filtering_revision3/tables/kegg_filtered_hellinger_pca_scores.csv")
KEGG_DIFF = read("15.kegg_ontology_filtering_revision3/tables/kegg_filtered_phylogroup_differential_features.csv")
KEGG_DISP = read("15.kegg_ontology_filtering_revision3/tables/kegg_filtered_permdisp_distances.csv")
KEGG_STAT = read("15.kegg_ontology_filtering_revision3/tables/kegg_filtered_global_statistics.csv")
GEM_STAT = read("10.mSystems_revision2_model_evidence/03_distance_stats/gem_distance_permanova_permdisp_ani_residual_summary.csv")
CRISPR = read("09.mSystems_revision_data_freeze/01_master_manifest/harmonized_results/03.analysis__crispr_genome_counts.csv")
CRISPR_TARGET = read("03.analysis/crispr_target_matching_analysis/crispr_target_matching_phylogroup_summary_all_genomes.csv")
SIGNALP = read("09.mSystems_revision_data_freeze/01_master_manifest/harmonized_results/03.analysis__signalp_genome_counts.csv")
TMHMM = read("09.mSystems_revision_data_freeze/01_master_manifest/harmonized_results/03.analysis__recalculated_annotations_tmhmm_updated.csv")
MONK = read("03.analysis/monk2013_sd01_available_benchmark_fba/available_benchmark_770_fba_accuracy_metrics.csv")
BAC_TRUTH = read("03.analysis/bacdive_mg1655_fba/bacdive_mg1655_truth_table.csv")
BAC_PRED = read("03.analysis/bacdive_mg1655_fba/bacdive_mg1655_fba_predictions_long.csv")
GAL_SENS = read("12.galardini_fba_benchmark/04_metrics/galardini_threshold_sensitivity_metrics.csv")
GAL_CORR = read("18.mSystems_revision3_final_20260911/12_nature_reviewer_presubmission/corrected_galardini/galardini_corrected_metrics_by_variant.csv")
MEM = read("14.mSystems_revision3_submission_candidate/02_analysis_evidence/memote_18_scores.csv")


def clean(ax, axis="y"):
    ax.spines[["top", "right"]].set_visible(False)
    if axis:
        ax.grid(axis=axis, color="#D9DEE3", lw=0.45, alpha=0.8)
        ax.set_axisbelow(True)


def letter(fig, value):
    fig.text(0.008, 0.988, value, va="top", ha="left", fontsize=9, fontweight="bold")


def save_panel(name, fig, panel_letter, align=False):
    letter(fig, panel_letter)
    fig.canvas.draw()
    if align and len(fig.axes) > 1:
        try:
            require_matplotlib_panel_alignment(
                fig, json_out=QA / f"{name}_alignment.json", strict=False,
                require_panel_labels=False,
            )
            (QA / f"{name}_alignment_error.txt").unlink(missing_ok=True)
        except Exception as exc:
            (QA / f"{name}_alignment_error.txt").write_text(str(exc), encoding="utf-8")
    fig.savefig(PANELS / f"{name}.pdf")
    fig.savefig(PANELS / f"{name}.svg")
    fig.savefig(PANELS / f"{name}.png", dpi=300)
    plt.close(fig)


def box_points(ax, df, col, ylabel, seed=20260911):
    vals = [pd.to_numeric(df.loc[df.Phylogroup.eq(g), col], errors="coerce").dropna().to_numpy() for g in ORDER]
    bp = ax.boxplot(vals, positions=np.arange(len(ORDER)), widths=0.58, patch_artist=True,
                    showfliers=False, medianprops={"color":"#1F252B", "lw":0.9},
                    whiskerprops={"color":"#59636D", "lw":0.6}, capprops={"color":"#59636D", "lw":0.6},
                    boxprops={"lw":0.6, "edgecolor":"#59636D"})
    for g, b in zip(ORDER, bp["boxes"]):
        b.set_facecolor(COLORS[g]); b.set_alpha(0.72)
    rng = np.random.default_rng(seed)
    for i, (g, a) in enumerate(zip(ORDER, vals)):
        jitter = rng.uniform(-0.22, 0.22, len(a))
        ax.scatter(i + jitter, a, s=3.5, color=COLORS[g], alpha=0.20, linewidth=0, rasterized=True)
    ax.set_xticks(range(len(ORDER)), ORDER)
    ax.set_xlabel("Phylogroup")
    ax.set_ylabel(ylabel)
    clean(ax)


def sbml_ids(path):
    out = set()
    for _, elem in ET.iterparse(path, events=("end",)):
        if elem.tag.rsplit("}", 1)[-1] == "reaction" and elem.get("id"):
            out.add(elem.get("id"))
        elem.clear()
    return out


def crop_workflow():
    src = ROOT / "绘图/v3.drawio.png"
    im = Image.open(src).convert("RGB")
    bg = Image.new("RGB", im.size, "white")
    diff = ImageChops.difference(im, bg).convert("L")
    bbox = diff.point(lambda x: 0 if x < 12 else 255).getbbox()
    if bbox:
        pad = 10
        bbox = (max(0, bbox[0]-pad), max(0, bbox[1]-pad), min(im.width, bbox[2]+pad), min(im.height, bbox[3]+pad))
        im = im.crop(bbox)
    dst = ASSETS / "Figure1A_workflow_cropped.png"
    im.save(dst, dpi=(600, 600))
    return dst


def panel_workflow():
    path = crop_workflow()
    im = Image.open(path)
    fig, ax = plt.subplots(figsize=(7.0, 2.55), constrained_layout=True)
    ax.imshow(im)
    ax.axis("off")
    save_panel("Figure1A_workflow", fig, "A")


def panel_1b():
    d = RUNS[RUNS.Version_ID.isin(["V1", "V6"])].pivot(index="Accession", columns="Version_ID", values="Wallclock_Seconds").dropna()
    speedup = (d.V1 / d.V6).sort_values()
    d.assign(Speedup_Fold=d.V1 / d.V6).reset_index().rename(
        columns={"V1":"Default_Wallclock_Seconds", "V6":"Optimized_Wallclock_Seconds"}
    ).to_csv(SOURCE_DATA / "Figure1B_runtime_and_speedup.csv", index=False)
    fig = plt.figure(figsize=(3.42, 2.45))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.45, 0.72], hspace=0.48,
                          left=.19, right=.98, bottom=.18, top=.93)
    ax = fig.add_subplot(gs[0])
    vals = [d.V1.to_numpy(), d.V6.to_numpy()]
    vp = ax.violinplot(vals, positions=[0, 1], widths=.72, showextrema=False)
    for body, color in zip(vp["bodies"], [NEUTRAL, BLUE]):
        body.set_facecolor(color); body.set_edgecolor("none"); body.set_alpha(.25)
    rng = np.random.default_rng(20260912)
    for x, (arr, color) in enumerate(zip(vals, [NEUTRAL, BLUE])):
        ax.scatter(x + rng.uniform(-.16, .16, len(arr)), arr, s=10, color=color,
                   alpha=.68, edgecolor="white", lw=.25, zorder=3)
        q1, med, q3 = np.percentile(arr, [25, 50, 75])
        ax.vlines(x, q1, q3, color="#20262C", lw=2.3, zorder=4)
        ax.scatter(x, med, s=22, color="white", edgecolor="#20262C", lw=.7, zorder=5)
    ax.set_yscale("log")
    ax.set_xticks([0, 1], ["Default", "Optimized"])
    ax.set_ylabel("Wall-clock time (s)")
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}" if v >= 1 else f"{v:g}"))
    clean(ax)
    ax2 = fig.add_subplot(gs[1])
    rank = np.arange(1, len(speedup) + 1)
    ax2.plot(rank, speedup.to_numpy(), color="#9BA3AA", lw=.65, zorder=1)
    ax2.scatter(rank, speedup.to_numpy(), s=13, color=ORANGE, edgecolor="white", lw=.3, zorder=2)
    ax2.axhline(speedup.median(), color="#20262C", lw=.8, ls="--")
    ax2.text(.03, .90, f"median {speedup.median():.2f}x", transform=ax2.transAxes,
             va="top", ha="left", fontsize=7,
             bbox={"facecolor":"white", "edgecolor":"none", "pad":1.2, "alpha":.9})
    ax2.set_xlim(.5, len(speedup) + .5)
    ax2.set_xticks([1, 5, 10, 15, 20])
    ax2.set_xlabel("Genome rank by speedup")
    ax2.set_ylabel("Speedup")
    clean(ax2)
    save_panel("Figure1B_paired_runtime", fig, "B")


def panel_1c():
    d = RUNS[RUNS.Version_ID.isin(["V1", "V6"])].pivot(index="Accession", columns="Version_ID", values="Wallclock_Seconds").dropna()
    n50 = read("04.基准测试/test_20/reaction_set_delta/reaction_set_delta_summary.csv")[["Accession", "Actual_N50"]].drop_duplicates()
    meta = read("00.data/all_20_genomes_metadata.csv")[["Accession", "Phylogroup"]]
    q = (d.V1 / d.V6).rename("Speedup").reset_index().merge(n50, on="Accession").merge(meta, on="Accession")
    fig, ax = plt.subplots(figsize=(3.42, 2.45), constrained_layout=True)
    for g in ORDER:
        z = q[q.Phylogroup.eq(g)]
        if len(z):
            ax.scatter(z.Actual_N50, z.Speedup, s=25, color=COLORS[g], edgecolor="white", lw=0.45, label=g)
    rho, p = spearmanr(q.Actual_N50, q.Speedup)
    coef = np.polyfit(q.Actual_N50, q.Speedup, 1)
    xx = np.linspace(q.Actual_N50.min(), q.Actual_N50.max(), 200)
    ax.plot(xx, np.polyval(coef, xx), color=ORANGE, lw=1.0)
    ax.set_title(f"Spearman r = {rho:.2f}; p = {p:.3g}", pad=3)
    ax.set_xlabel("Assembly N50 (bp)")
    ax.set_ylabel("Runtime speedup (fold)")
    clean(ax)
    ax.legend(ncol=5, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.19), handletextpad=0.15, columnspacing=0.45)
    save_panel("Figure1C_N50_speedup", fig, "C")


def authoritative_reactions():
    rows = []
    members = []
    for acc, sub in RUNS[RUNS.Version_ID.isin(["V1", "V6"])].groupby("Accession"):
        paths = {r.Version_ID: ROOT / Path(str(r.Final_XML_Relative_Path)) for _, r in sub.iterrows()}
        if set(paths) != {"V1", "V6"} or not all(p.exists() for p in paths.values()):
            raise FileNotFoundError(f"Missing authoritative V1/V6 pair: {acc}")
        a, b = sbml_ids(paths["V1"]), sbml_ids(paths["V6"])
        vals = {"Shared": a & b, "Default only": a - b, "Optimized only": b - a}
        rows.append({"Accession": acc, **{k: len(v) for k, v in vals.items()}})
        for kind, ids in vals.items():
            members.extend({"Accession":acc, "Reaction_set":kind, "Reaction_ID":rid} for rid in sorted(ids))
    pd.DataFrame(members).to_csv(HERE / "Figure1D_authoritative_SBML_reaction_membership.csv.gz", index=False, compression="gzip")
    out = pd.DataFrame(rows)
    out.to_csv(HERE / "Figure1D_authoritative_SBML_reaction_statistics.csv", index=False)
    return out


def panel_1d():
    rx = authoritative_reactions()
    labels = ["Shared", "Default only", "Optimized only"]
    colors = ["#56B4E9", "#999999", "#E69F00"]
    rx = rx.assign(default_fraction=rx["Default only"] / rx[labels].sum(axis=1)).sort_values("default_fraction", ascending=False)
    rx.to_csv(SOURCE_DATA / "Figure1D_authoritative_reaction_set_counts.csv", index=False)
    fig, ax = plt.subplots(figsize=(7.0, 2.25))
    fig.subplots_adjust(left=.09, right=.995, bottom=.24, top=.82)
    y = np.arange(len(rx))
    left = np.zeros(len(rx))
    for col, color in zip(labels, colors):
        vals = rx[col].to_numpy()
        med = np.median(vals)
        ax.barh(y, vals, left=left, height=.74, color=color, edgecolor="white", lw=.18,
                label=f"{col} (median {med:,.0f})")
        left += vals
    ax.invert_yaxis()
    ax.set_yticks([0, len(rx)-1], ["1", str(len(rx))])
    ax.set_ylabel("Genome rank")
    ax.set_xlabel("Reaction count in union of paired SBML models")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax.legend(frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(.5, 1.01),
              handlelength=1.3, columnspacing=1.1)
    clean(ax, axis="x")
    save_panel("Figure1D_authoritative_reaction_sets", fig, "D")


def triptych(name, letter_value, df, specs, size=(7.0, 1.85)):
    fig, axs = plt.subplots(1, 3, figsize=size)
    for ax, (col, lab) in zip(axs, specs):
        box_points(ax, df, col, lab)
    fig.subplots_adjust(left=.075, right=.992, bottom=.28, top=.90, wspace=.43)
    save_panel(name, fig, letter_value, align=True)


def violin_triptych(name, letter_value, df, specs, size=(7.0, 2.55)):
    """Dense raw-unit distributions with one horizontal violin per phylogroup."""
    id_cols = [c for c in ["Genome_ID", "Accession", "Phylogroup"] if c in df.columns]
    df[id_cols + [c for c, _ in specs]].to_csv(SOURCE_DATA / f"{name}.csv", index=False)
    fig, axs = plt.subplots(1, 3, figsize=size, sharey=True)
    y = np.arange(len(ORDER))
    for j, (ax, (col, lab)) in enumerate(zip(axs, specs)):
        arrays = [pd.to_numeric(df.loc[df.Phylogroup.eq(g), col], errors="coerce").dropna().to_numpy() for g in ORDER]
        valid = [(i, a) for i, a in enumerate(arrays) if len(a) >= 2 and np.nanstd(a) > 0]
        if valid:
            vp = ax.violinplot([a for _, a in valid], positions=[i for i, _ in valid], orientation="horizontal",
                               widths=.78, showextrema=False, points=80, bw_method=.35)
            for body, (i, _) in zip(vp["bodies"], valid):
                body.set_facecolor(COLORS[ORDER[i]]); body.set_edgecolor("none"); body.set_alpha(.42)
        for i, a in enumerate(arrays):
            if len(a) == 0:
                continue
            q1, med, q3 = np.percentile(a, [25, 50, 75])
            ax.hlines(i, q1, q3, color="#283038", lw=1.3, zorder=3)
            ax.scatter(med, i, s=14, color=COLORS[ORDER[i]], edgecolor="#20262C", lw=.45, zorder=4)
        ax.axvline(pd.to_numeric(df[col], errors="coerce").median(), color="#7C858D", lw=.65, ls="--", zorder=0)
        ax.set_xlabel(lab)
        ax.set_ylim(len(ORDER)-.55, -.55)
        ax.set_yticks(y)
        if j == 0:
            ax.set_yticklabels(ORDER)
        else:
            ax.tick_params(axis="y", labelleft=False)
        clean(ax, axis="x")
    fig.text(.012, .53, "Phylogroup", rotation=90, rotation_mode="anchor", va="center", fontsize=7.5)
    fig.subplots_adjust(left=.085, right=.975, bottom=.20, top=.93, wspace=.24)
    save_panel(name, fig, letter_value, align=True)


def panel_3a():
    # Exact overlapping coordinates are aggregated; colored diamonds show phylogroup centroids.
    grouped = KEGG_PCA.groupby(["PC1", "PC2"], as_index=False).size()
    fig, ax = plt.subplots(figsize=(3.42, 2.58))
    sizes = 5 + 9 * np.sqrt(grouped["size"])
    ax.scatter(grouped.PC1, grouped.PC2, s=sizes, color="#B8C0C7", alpha=0.55, linewidth=0, rasterized=True)
    handles = []
    for g in ORDER:
        q = KEGG_PCA[KEGG_PCA.Phylogroup.eq(g)]
        if len(q):
            x, y = q.PC1.mean(), q.PC2.mean()
            h = ax.scatter(x, y, marker="D", s=30, color=COLORS[g], edgecolor="white", lw=0.6, zorder=3, label=g)
            handles.append(h)
    ev1 = KEGG_PCA.PC1_explained_variance.iloc[0] * 100
    ev2 = KEGG_PCA.PC2_explained_variance.iloc[0] * 100
    ax.set_xlabel(f"PC1 ({ev1:.1f}%)")
    ax.set_ylabel(f"PC2 ({ev2:.1f}%)")
    clean(ax)
    ax.legend(handles=handles, ncol=5, frameon=False, loc="upper center", bbox_to_anchor=(.5, -.18),
              handletextpad=.15, columnspacing=.45)
    fig.subplots_adjust(left=.17, right=.97, bottom=.27, top=.88)
    save_panel("Figure3A_filtered_KEGG_PCA", fig, "A")


def panel_3b():
    q = KEGG_DIFF.query("test_status == 'tested_chi_square'").sort_values("Cramers_V", ascending=False).copy()
    short = {"map00624":"PAH degradation", "map00622":"Xylene degradation",
             "map00621":"Dioxin degradation", "map00365":"Furfural degradation"}
    labels = []
    for _, r in q.iterrows():
        sig = "q<0.001" if r.BH_q_value < 0.001 else "q=1.00 (ns)"
        labels.append(f"{r.KEGG_Map_ID}  {short.get(r.KEGG_Map_ID, r.Description)}  "
                      f"[{r.min_frequency_group} -> {r.max_frequency_group}]\n"
                      f"V={r.Cramers_V:.2f}; {sig}")
    q["label"] = labels
    fig, ax = plt.subplots(figsize=(3.42, 2.58), constrained_layout=True)
    y = np.arange(len(q))
    ax.hlines(y, q.min_frequency, q.max_frequency, color="#B6BEC6", lw=2)
    ax.scatter(q.min_frequency, y, s=25, color=NEUTRAL, label="Lowest group", zorder=3)
    ax.scatter(q.max_frequency, y, s=25, color=BLUE, label="Highest group", zorder=3)
    ax.set_yticks(y, q.label)
    ax.invert_yaxis()
    ax.set_xlim(-0.03, 1.04)
    ax.set_xlabel("Genome frequency")
    ax.tick_params(axis="y", length=0)
    clean(ax, axis="x")
    save_panel("Figure3B_KEGG_pathway_effects", fig, "B")


def panel_3d():
    per = KEGG_STAT.loc[KEGG_STAT.test.eq("PERMANOVA_oneway_phylogroup")].iloc[0]
    disp = KEGG_STAT.loc[KEGG_STAT.test.eq("PERMDISP_phylogroup")].iloc[0]
    man = KEGG_STAT.loc[KEGG_STAT.test.str.startswith("Mantel")].iloc[0]
    res = KEGG_STAT.loc[KEGG_STAT.test.str.startswith("ANI20")].iloc[0]
    rows = [("Centroid structure", "PERMANOVA R²", float(per.R2), float(per.p_value)),
            ("Dispersion heterogeneity", "PERMDISP F", float(disp.pseudo_F), float(disp.p_value)),
            ("Continuous relatedness", "Mantel r", float(man.mantel_r), float(man.p_value_two_sided)),
            ("Residual phylogroup", "ANI-adjusted R²", float(res.unique_R2), float(res.p_value))]
    fig, ax = plt.subplots(figsize=(3.42, 2.58), constrained_layout=True)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    xs = [0.02, 0.50, 0.78, 0.91]
    for x, h in zip(xs, ["Evidence", "Statistic", "Value", "P"]): ax.text(x, .92, h, fontweight="bold")
    for i, (meaning, stat, value, p) in enumerate(rows):
        yy = .73 - i*.19
        ax.text(xs[0], yy, meaning); ax.text(xs[1], yy, stat)
        ax.text(xs[2], yy, f"{value:.4f}" if value < 1 else f"{value:.1f}")
        ax.text(xs[3], yy, f"{p:.3f}")
        ax.plot([.02,.98], [yy-.075, yy-.075], color="#D9DEE3", lw=.6)
    save_panel("Figure3D_KEGG_effect_summary", fig, "D")


def pair_data():
    uids = KEGG.Genome_UID.tolist()
    ani = read("09.mSystems_revision_data_freeze/01_master_manifest/harmonized_results/core_genome_ani_distance_matrix_by_genome_uid.csv.gz").set_index("Genome_UID").loc[uids, uids]
    gem = read("10.mSystems_revision2_model_evidence/03_distance_stats/gem_gene_evidenced_jaccard_distance_by_genome_uid.csv.gz").set_index("Genome_UID").loc[uids, uids]
    feats = [c for c in KEGG if c.startswith("kegg_pathway_")]
    x = KEGG.set_index("Genome_UID").loc[uids, feats].astype(float).to_numpy()
    x = np.sqrt(x / np.maximum(x.sum(axis=1, keepdims=True), 1))
    i, j = np.triu_indices(len(uids), 1)
    return ani.to_numpy()[i, j], np.sqrt(((x[i]-x[j])**2).sum(axis=1)), gem.to_numpy()[i, j]


def hex_panel(x, y, ylabel, r, p, name, panel_letter):
    fig, ax = plt.subplots(figsize=(3.42, 2.62), constrained_layout=True)
    hb = ax.hexbin(x, y, gridsize=60, mincnt=1, norm=LogNorm(), cmap="viridis", linewidths=0, rasterized=True)
    ax.text(.04, .96, f"Mantel r = {r:.3f}\np = {p:.3f}", transform=ax.transAxes, va="top")
    ax.set_xlabel("FastANI distance")
    ax.set_ylabel(ylabel)
    clean(ax, axis=None)
    cb = fig.colorbar(hb, ax=ax, fraction=.045, pad=.025)
    cb.set_label("Pair count")
    ticks = [v for v in [1, 10, 100, 1000, 10000] if hb.norm.vmin <= v <= hb.norm.vmax]
    cb.set_ticks(ticks); cb.set_ticklabels([f"{v:,}" for v in ticks]); cb.ax.tick_params(labelsize=7)
    save_panel(name, fig, panel_letter)


def horizontal_effect(name, panel_letter, values, xlabel, colors, log=False, decimals=3):
    labels = ["KEGG", "Full GEM", "Gene-evidenced GEM", "Gap-filled reactions"]
    fig, ax = plt.subplots(figsize=(3.42, 2.62), constrained_layout=True)
    y = np.arange(4)
    left = min(values) / 1.4 if log else 0
    ax.hlines(y, left, values, color="#B6BEC6", lw=2)
    ax.scatter(values, y, s=34, color=colors, zorder=3)
    for yi, v in zip(y, values):
        if log:
            ax.annotate(f"{v:.{decimals}f}", (v, yi), xytext=(4, 7), textcoords="offset points", va="bottom")
        else:
            ax.text(max(values)*1.07, yi, f"{v:.{decimals}f}", va="center", ha="left")
    ax.set_yticks(y, labels); ax.invert_yaxis(); ax.set_xlabel(xlabel)
    if log:
        ax.set_xscale("log"); ax.set_xlim(7e-5, 8e-3)
        ax.set_xticks([1e-4, 1e-3, 1e-2]); ax.set_xticklabels(["0.0001", "0.001", "0.01"])
    else:
        ax.set_xlim(0, max(values)*1.20)
    clean(ax, axis="x")
    save_panel(name, fig, panel_letter)


def cog_panel():
    cogs = [f"COG_{c}" for c in "ABCDEFGHIJKLMNOPQRSTUVWXYZ" if f"COG_{c}" in ASSEMBLY_DF.columns and ASSEMBLY_DF[f"COG_{c}"].sum() > 0]
    cg = ASSEMBLY_DF.groupby("Phylogroup")[cogs].mean().reindex(ORDER).T
    z = cg.sub(cg.mean(axis=1), axis=0).div(cg.std(axis=1).replace(0, np.nan), axis=0).fillna(0)
    fig, ax = plt.subplots(figsize=(7.0, 3.15), constrained_layout=True)
    im = ax.imshow(z, aspect="auto", cmap="RdBu_r", vmin=-2.4, vmax=2.4)
    ax.set_xticks(range(9), ORDER); ax.set_yticks(range(len(cogs)), [c[-1] for c in cogs])
    ax.set_xlabel("Phylogroup"); ax.set_ylabel("COG functional class")
    cb = fig.colorbar(im, ax=ax, fraction=.022, pad=.018); cb.set_label("Row z score")
    save_panel("Figure5A_COG_profile", fig, "A")


def wilson(k, n, z=1.96):
    p = k/n
    den = 1 + z*z/n
    center = (p + z*z/(2*n))/den
    half = z*math.sqrt(p*(1-p)/n + z*z/(4*n*n))/den
    return center-half, center+half


def crispr_target_panel():
    q = CRISPR_TARGET.set_index("Phylogroup").reindex(ORDER)
    fig, ax = plt.subplots(figsize=(7.0, 2.25), constrained_layout=True)
    y = np.arange(len(ORDER)); off=.13
    for label, kcol, color, marker, delta in [
        ("Phage target", "genomes_with_phage_targets", BLUE, "o", -off),
        ("Plasmid target", "genomes_with_plasmid_targets", ORANGE, "s", off),
    ]:
        k = q[kcol].astype(int).to_numpy(); n = q.n_genomes.astype(int).to_numpy(); p = k/n
        ci = np.array([wilson(ki, ni) for ki, ni in zip(k, n)])
        ax.errorbar(p*100, y+delta, xerr=np.vstack([(p-ci[:,0])*100, (ci[:,1]-p)*100]),
                    fmt=marker, ms=4.5, capsize=2.2, lw=.8, color=color, label=label)
    ax.set_yticks(y, [f"{g}  (n={int(q.loc[g,'n_genomes'])})" for g in ORDER]); ax.invert_yaxis()
    ax.set_xlabel("Genomes with a target match (%, Wilson 95% CI)")
    ax.set_xlim(0, 104); clean(ax, axis="x")
    ax.legend(frameon=False, ncol=2, loc="lower center", bbox_to_anchor=(.5,1.01))
    save_panel("Figure5B_CRISPR_target_recovery", fig, "B")


def scatter_supp(name, letter_value, x, y, xlabel, ylabel):
    fig, ax = plt.subplots(figsize=(3.42, 2.65), constrained_layout=True)
    groups = QC_DF.Phylogroup if len(QC_DF) == len(x) else EVID.Phylogroup
    for g in ORDER:
        mask = groups.eq(g).to_numpy()
        ax.scatter(np.asarray(x)[mask], np.asarray(y)[mask], s=8, alpha=.48, color=COLORS[g], linewidth=0, rasterized=True)
    rho, p = spearmanr(x, y, nan_policy="omit")
    ax.set_title(f"Spearman r = {rho:.2f}; p = {p:.3g}", pad=3)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel); clean(ax)
    save_panel(name, fig, letter_value)


def validation_panels():
    mm = MONK[MONK.Variant.isin(["gapseq_default", "gapseq_optimized"])].groupby(["Benchmark_Model", "Variant"], as_index=False).Accuracy.mean()
    pv = mm.pivot(index="Benchmark_Model", columns="Variant", values="Accuracy").dropna()
    fig, ax = plt.subplots(figsize=(3.42, 2.65), constrained_layout=True)
    y=np.arange(len(pv)); ax.hlines(y,pv.gapseq_default,pv.gapseq_optimized,color="#AAB1B8",lw=1.5)
    ax.scatter(pv.gapseq_default,y,color=NEUTRAL,s=25,label="Default"); ax.scatter(pv.gapseq_optimized,y,color=BLUE,s=25,label="Optimized")
    ax.set_yticks(y,pv.index); ax.invert_yaxis(); ax.set_xlabel("Agreement with literature in silico matrix")
    ax.set_xlim(.5,.8); ax.legend(frameon=False,ncol=2,loc="lower center",bbox_to_anchor=(.5,1.01)); clean(ax, axis="x")
    save_panel("SupplementaryFigureS1A_Monk_consistency", fig, "A")

    bp = BAC_PRED[(BAC_PRED.condition.eq("aerobic")) &
                  (BAC_PRED.panel.eq("carbon_assimilation_subset")) &
                  (BAC_PRED.status.eq("ok"))].copy()
    # Only conditions with a resolved exchange reaction are evaluable. Amygdalin
    # is retained in the mapping audit as unresolved but must not enter the
    # agreement denominator or be treated as a predicted negative phenotype.
    resolved_sets = [set(z.carbon_source) for _, z in bp.groupby("model_id")]
    resolved_sources = set.intersection(*resolved_sets)
    sources = [s for s in BAC_TRUTH[BAC_TRUTH.panel.eq("carbon_assimilation_subset")].carbon_source.tolist()
               if s in resolved_sources]
    models = ["MG1655_gapseq_default", "MG1655_gapseq_optimized", "iML1515"]
    truth = bp.drop_duplicates("carbon_source").set_index("carbon_source").reindex(sources).truth.astype(float)
    pred = bp.pivot(index="model_id",columns="carbon_source",values="prediction").reindex(index=models,columns=sources).astype(float)
    mat = np.column_stack([truth.to_numpy(), pred.T.to_numpy()])
    accuracies = np.mean(mat[:, 1:] == mat[:, [0]], axis=0)
    xlabels = ["Observed", f"Default\n({accuracies[0]*100:.0f}% match)",
               f"Optimized\n({accuracies[1]*100:.0f}% match)", f"iML1515\n({accuracies[2]*100:.0f}% match)"]
    pd.DataFrame(mat, index=sources, columns=["Observed", "Default", "Optimized", "iML1515"]).rename_axis(
        "Carbon_Source"
    ).reset_index().to_csv(SOURCE_DATA / "SupplementaryFigureS1B_BacDive_predictions.csv", index=False)
    fig, ax = plt.subplots(figsize=(5.35, 2.85))
    fig.subplots_adjust(left=.25, right=.985, bottom=.21, top=.87)
    for i in range(len(sources)):
        if i % 2 == 0:
            ax.axhspan(i-.5, i+.5, color="#F4F6F7", zorder=0)
        for j in range(4):
            val = mat[i, j]
            mismatch = j > 0 and val != mat[i, 0]
            ax.scatter(j, i, s=54, facecolor=BLUE if val == 1 else "white",
                       edgecolor=ORANGE if mismatch else "#35404A",
                       linewidth=1.25 if mismatch else .65, zorder=3)
    ax.set_xlim(-.45, 3.45); ax.set_ylim(len(sources)-.5, -.5)
    ax.set_yticks(range(len(sources)), sources); ax.set_xticks(range(4), xlabels)
    ax.xaxis.tick_top(); ax.tick_params(axis="x", length=0, pad=3)
    ax.set_ylabel("Carbon source")
    for spine in ax.spines.values(): spine.set_visible(False)
    growth = ax.scatter([], [], s=42, facecolor=BLUE, edgecolor="#35404A", lw=.65, label="Growth")
    no_growth = ax.scatter([], [], s=42, facecolor="white", edgecolor="#35404A", lw=.65, label="No growth")
    mismatch = ax.scatter([], [], s=42, facecolor="white", edgecolor=ORANGE, lw=1.25, label="Prediction mismatch")
    ax.legend(handles=[growth, no_growth, mismatch], frameon=False, ncol=3, loc="upper center",
              bbox_to_anchor=(.5, -.11), handletextpad=.35, columnspacing=.9)
    save_panel("SupplementaryFigureS1B_BacDive_predictions", fig, "B")

    gm = GAL_CORR.copy()
    fig, ax = plt.subplots(figsize=(3.42, 2.65))
    fig.subplots_adjust(left=.31, right=.98, bottom=.22, top=.88)
    labels={"gapseq_default":"Default gapseq","gapseq_optimized":"Optimized gapseq","curated_reference":"Curated reference"}
    gm["display"] = gm.variant.map(labels)
    gm = gm.set_index("display").reindex(["Default gapseq", "Optimized gapseq", "Curated reference"]).reset_index()
    gm.to_csv(SOURCE_DATA / "SupplementaryFigureS1C_Galardini_validation.csv", index=False)
    outcomes = [("TN", "#56B4E9"), ("FP", "#E69F00"), ("FN", "#CC6677"), ("TP", "#009E73")]
    y = np.arange(len(gm))
    for i, row in gm.iterrows():
        for j, (col, color) in enumerate(outcomes):
            val = int(row[col])
            size = 22 + 8.5 * val
            ax.scatter(j, i, s=size, facecolor=color, edgecolor="none",
                       alpha=.92 if val else .18, zorder=2)
            ax.text(j, i, str(val), ha="center", va="center", fontsize=7,
                    color="white" if val >= 5 else "#20262C", fontweight="bold", zorder=3)
        ax.text(3.55, i, f"n={int(row['n_model_observations'])}", va="center", ha="left", fontsize=7)
    ax.set_xlim(-.5, 4.15); ax.set_ylim(len(gm)-.55, -.55)
    ax.set_xticks(range(4), [x for x, _ in outcomes]); ax.set_yticks(y, gm.display)
    ax.set_xlabel("Prediction outcome")
    ax.set_title("Experiment-matched outcomes", pad=3)
    for spine in ax.spines.values(): spine.set_visible(False)
    ax.tick_params(axis="both", length=0)
    save_panel("SupplementaryFigureS1C_Galardini_validation", fig, "C")


def memote_panel():
    q=MEM.copy(); q["Phylogroup"]=q.report_file.str.extract(r"ECGEM_\d+_([^_]+)_")[0]
    q["Selection"]=np.where(q.report_file.str.contains("highest_qc"),"Highest assembly QC","ANI medoid")
    q["order"]=q.Phylogroup.map({g:i for i,g in enumerate(ORDER)}); q=q.sort_values(["order","Selection"])
    fig, axs=plt.subplots(1,3,figsize=(7.0,3.25),sharey=True)
    y=np.arange(len(q)); markers={"ANI medoid":"o","Highest assembly QC":"s"}
    for ax,col,lab in zip(axs,["total_score","consistency","annotation_gene"],["Total score (%)","Consistency (%)","Gene annotation (%)"]):
        for pos,(_,r) in enumerate(q.iterrows()):
            ax.scatter(float(r[col])*100,pos,s=25,marker=markers[r.Selection],color=COLORS[r.Phylogroup],edgecolor="white",lw=.4)
        ax.set_xlabel(lab); ax.set_yticks(y); clean(ax,axis="x")
    axs[0].set_yticklabels([f"{r.Phylogroup}  {r.report_file.split('_')[1]}" for _,r in q.iterrows()]); axs[0].invert_yaxis()
    for sel,m in markers.items(): axs[2].scatter([],[],marker=m,color="#555",label=sel)
    handles, labels = axs[2].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, loc="lower center", bbox_to_anchor=(.62,.01), ncol=2)
    fig.subplots_adjust(left=.13,right=.992,bottom=.22,top=.94,wspace=.32)
    save_panel("SupplementaryFigureS2C_MEMOTE_18_scores",fig,"C",align=True)


def assemble(output_name, page_size, placements):
    # placements: (panel pdf, x, y, width, height), units are inches from top-left.
    doc=fitz.open(); page=doc.new_page(width=page_size[0]*72,height=page_size[1]*72)
    opened=[]
    for filename,x,y,w,h in placements:
        src=fitz.open(PANELS/filename); opened.append(src)
        rect=fitz.Rect(x*72,y*72,(x+w)*72,(y+h)*72)
        page.show_pdf_page(rect,src,0,keep_proportion=True)
    pdf=ASSEMBLED/f"{output_name}.pdf"; doc.save(pdf,garbage=4,deflate=True); doc.close()
    for src in opened: src.close()
    d=fitz.open(pdf); pix=d[0].get_pixmap(dpi=300,alpha=False); png=ASSEMBLED/f"{output_name}.png"; pix.save(png); d.close()
    im=Image.open(png).convert("RGB"); im.save(ASSEMBLED/f"{output_name}.tiff",dpi=(300,300),compression="tiff_lzw")


def main():
    panel_workflow(); panel_1b(); panel_1c(); panel_1d()
    violin_triptych("Figure2A_genome_features","A",ASSEMBLY_DF,[("Genome_Size_Mb","Genome size (Mb)"),("Total_CDS","CDS count"),("GC_Percent","GC content (%)")])
    triptych("Figure2B_GEM_topology","B",QC_DF,[("Total_Genes","Model genes"),("Total_Reactions","Reactions"),("Total_Metabolites","Metabolites")])
    violin_triptych("Figure2C_GEM_quality","C",QC_DF,[("Mass_Balanced_Ratio(%)","Mass balanced (%)"),("Blocked_Reactions_Ratio(%)","Blocked reactions (%)"),("Deadend_Metabolites_Ratio(%)","Dead-end metabolites (%)")])
    triptych("Figure2D_gapfill_gene_evidence","D",EVID,[("Gap_Filled_Reactions","Gap-filled reactions"),("Gap_Filled_Fraction","Gap-filled fraction"),("Gene_Evidenced_Fraction","Gene-evidenced fraction")])
    panel_3a(); panel_3b();
    fig,ax=plt.subplots(figsize=(3.42,2.58),constrained_layout=True); box_points(ax,KEGG_DISP,"distance_to_group_centroid","Distance to group centroid"); save_panel("Figure3C_KEGG_PERMDISP",fig,"C")
    panel_3d()
    ani,kegg_d,gem_d=pair_data(); man=KEGG_STAT.loc[KEGG_STAT.test.str.startswith("Mantel")].iloc[0]
    gene=GEM_STAT.loc[GEM_STAT.Network.eq("gem_gene_evidenced")].iloc[0]
    hex_panel(ani,kegg_d,"KEGG Hellinger distance",float(man.mantel_r),float(man.p_value_two_sided),"Figure4A_ANI_KEGG_distance","A")
    hex_panel(ani,gem_d,"Gene-evidenced GEM Jaccard distance",float(gene.ANI_Mantel_Spearman_r),float(gene.ANI_Mantel_p),"Figure4B_ANI_gene_evidenced_GEM_distance","B")
    mantel=[float(man.mantel_r)]+GEM_STAT.set_index("Network").loc[["gem_full_final","gem_gene_evidenced","gem_gapfilled"],"ANI_Mantel_Spearman_r"].astype(float).tolist()
    horizontal_effect("Figure4C_ANI_Mantel_summary","C",mantel,"Mantel Spearman r with FastANI distance",[GREEN]*4)
    res=KEGG_STAT.loc[KEGG_STAT.test.str.startswith("ANI20")].iloc[0]
    residual=[float(res.unique_R2)]+GEM_STAT.set_index("Network").loc[["gem_full_final","gem_gene_evidenced","gem_gapfilled"],"Phylogroup_Unique_Residual_R2_after_ANI20"].astype(float).tolist()
    horizontal_effect("Figure4D_ANI_adjusted_residuals","D",residual,"Unique phylogroup residual R² after ANI adjustment",["#56B4E9","#999999",BLUE,"#E69F00"],log=True,decimals=4)
    cog_panel(); crispr_target_panel()
    fig,axs=plt.subplots(1,2,figsize=(3.42,2.45),constrained_layout=True); box_points(axs[0],CRISPR,"Unique_Spacer_Sequences","Unique spacers"); box_points(axs[1],CRISPR,"Trusted_Arrays","Trusted arrays"); save_panel("Figure5C_CRISPR_burden",fig,"C",align=True)
    fig,axs=plt.subplots(1,2,figsize=(3.42,2.45),constrained_layout=True); box_points(axs[0],SIGNALP,"SignalP_per_1000_CDS","SignalP per 1,000 CDS"); box_points(axs[1],TMHMM,"TMHMM_Proteins_per_1000_CDS","TMHMM per 1,000 CDS"); save_panel("Figure5D_surface_features",fig,"D",align=True)
    validation_panels()
    scatter_supp("SupplementaryFigureS2A_QC_relationship","A",QC_DF["Deadend_Metabolites_Ratio(%)"],QC_DF["Blocked_Reactions_Ratio(%)"],"Dead-end metabolites (%)","Blocked reactions (%)")
    scatter_supp("SupplementaryFigureS2B_evidence_relationship","B",EVID.Gap_Filled_Fraction,EVID.Gene_Evidenced_Fraction,"Gap-filled fraction","Gene-evidenced fraction")
    memote_panel()

    assemble("Figure_1",(7,8.20),[("Figure1A_workflow.pdf",0,0,7,2.55),("Figure1B_paired_runtime.pdf",0,2.65,3.42,2.45),("Figure1C_N50_speedup.pdf",3.58,2.65,3.42,2.45),("Figure1D_authoritative_reaction_sets.pdf",0,5.36,7,2.25)])
    assemble("Figure_2",(7,9.15),[("Figure2A_genome_features.pdf",0,0,7,2.55),("Figure2B_GEM_topology.pdf",0,2.68,7,1.85),("Figure2C_GEM_quality.pdf",0,4.66,7,2.55),("Figure2D_gapfill_gene_evidence.pdf",0,7.34,7,1.85)])
    assemble("Figure_3",(7,5.32),[("Figure3A_filtered_KEGG_PCA.pdf",0,0,3.42,2.58),("Figure3B_KEGG_pathway_effects.pdf",3.58,0,3.42,2.58),("Figure3C_KEGG_PERMDISP.pdf",0,2.74,3.42,2.58),("Figure3D_KEGG_effect_summary.pdf",3.58,2.74,3.42,2.58)])
    assemble("Figure_4",(7,5.40),[("Figure4A_ANI_KEGG_distance.pdf",0,0,3.42,2.62),("Figure4B_ANI_gene_evidenced_GEM_distance.pdf",3.58,0,3.42,2.62),("Figure4C_ANI_Mantel_summary.pdf",0,2.78,3.42,2.62),("Figure4D_ANI_adjusted_residuals.pdf",3.58,2.78,3.42,2.62)])
    assemble("Figure_5",(7,8.02),[("Figure5A_COG_profile.pdf",0,0,7,3.15),("Figure5B_CRISPR_target_recovery.pdf",0,3.30,7,2.25),("Figure5C_CRISPR_burden.pdf",0,5.57,3.42,2.45),("Figure5D_surface_features.pdf",3.58,5.57,3.42,2.45)])
    assemble("Supplementary_Figure_S1",(7,5.66),[("SupplementaryFigureS1A_Monk_consistency.pdf",0,0,3.42,2.65),("SupplementaryFigureS1C_Galardini_validation.pdf",3.58,0,3.42,2.65),("SupplementaryFigureS1B_BacDive_predictions.pdf",.825,2.81,5.35,2.85)])
    assemble("Supplementary_Figure_S2",(7,6.06),[("SupplementaryFigureS2A_QC_relationship.pdf",0,0,3.42,2.65),("SupplementaryFigureS2B_evidence_relationship.pdf",3.58,0,3.42,2.65),("SupplementaryFigureS2C_MEMOTE_18_scores.pdf",0,2.81,7,3.25)])

    files=[]
    for p in sorted(PANELS.glob("*.pdf")):
        files.append({"panel":p.stem,"pdf":p.name,"svg":p.with_suffix('.svg').name,"png":p.with_suffix('.png').name,"font":"Arial","minimum_font_pt":7})
    pd.DataFrame(files).to_csv(HERE/"FIGURE_FILE_MANIFEST.csv",index=False)
    print(f"Generated {len(files)} panels and 7 assembled figures in {HERE}")


if __name__ == "__main__":
    main()
