from __future__ import annotations

import html
import math
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "00.data"
OUT_DIR = ROOT / "01.figures"

ASSEMBLY_CSV = DATA_DIR / "genome_assembly_and_cog_metrics.csv"
ANNOT_CSV = DATA_DIR / "recalculated_annotations.csv"

PHYLO_ORDER = ["A", "B1", "B2", "C", "D", "E", "F", "G", "U"]
COG_COLS = [f"COG_{chr(ord('A') + i)}" for i in range(26)]

COLORS = {
    "A": "#4C78A8",
    "B1": "#F58518",
    "B2": "#54A24B",
    "C": "#E45756",
    "D": "#72B7B2",
    "E": "#B279A2",
    "F": "#FF9DA6",
    "G": "#9D755D",
    "U": "#BAB0AC",
}


def esc(value: object) -> str:
    return html.escape(str(value), quote=True)


def svg_text(x, y, text, size=12, fill="#27313f", anchor="middle", weight="400", rotate=None):
    transform = f' transform="rotate({rotate} {x:.2f} {y:.2f})"' if rotate else ""
    return (
        f'<text x="{x:.2f}" y="{y:.2f}" font-family="Arial, Helvetica, sans-serif" '
        f'font-size="{size}" fill="{fill}" text-anchor="{anchor}" font-weight="{weight}"{transform}>'
        f"{esc(text)}</text>"
    )


def line(x1, y1, x2, y2, stroke="#5b6472", width=1):
    return f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" stroke="{stroke}" stroke-width="{width}"/>'


def rect(x, y, w, h, fill, stroke="none", width=1, opacity=1.0, rx=0):
    return (
        f'<rect x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" rx="{rx}" '
        f'fill="{fill}" stroke="{stroke}" stroke-width="{width}" opacity="{opacity}"/>'
    )


def circle(x, y, r, fill, stroke="none", width=1, opacity=1.0):
    return (
        f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{r:.2f}" fill="{fill}" '
        f'stroke="{stroke}" stroke-width="{width}" opacity="{opacity}"/>'
    )


def nice_ticks(vmin: float, vmax: float, n: int = 5):
    if math.isclose(vmin, vmax):
        pad = abs(vmin) * 0.1 + 1
        vmin -= pad
        vmax += pad
    span = vmax - vmin
    raw = span / max(n - 1, 1)
    mag = 10 ** math.floor(math.log10(raw))
    norm = raw / mag
    if norm <= 1:
        step = mag
    elif norm <= 2:
        step = 2 * mag
    elif norm <= 5:
        step = 5 * mag
    else:
        step = 10 * mag
    start = math.floor(vmin / step) * step
    stop = math.ceil(vmax / step) * step
    ticks = []
    value = start
    while value <= stop + step * 0.5:
        ticks.append(value)
        value += step
    return ticks


def fmt_tick(value: float) -> str:
    if abs(value) >= 100:
        return f"{value:.0f}"
    if abs(value) >= 10:
        return f"{value:.1f}".rstrip("0").rstrip(".")
    return f"{value:.2f}".rstrip("0").rstrip(".")


def box_stats(values):
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    q1, med, q3 = np.percentile(arr, [25, 50, 75])
    iqr = q3 - q1
    low_fence = q1 - 1.5 * iqr
    high_fence = q3 + 1.5 * iqr
    inliers = arr[(arr >= low_fence) & (arr <= high_fence)]
    low = float(np.min(inliers)) if len(inliers) else float(np.min(arr))
    high = float(np.max(inliers)) if len(inliers) else float(np.max(arr))
    outliers = arr[(arr < low) | (arr > high)]
    return float(q1), float(med), float(q3), low, high, outliers


def draw_boxplot_panel(parts, df, metric, title, ylabel, x, y, w, h):
    plot_x, plot_y = x + 62, y + 38
    plot_w, plot_h = w - 82, h - 86
    values = df[metric].astype(float).to_numpy()
    pad = (np.nanmax(values) - np.nanmin(values)) * 0.08 + 1e-9
    ymin, ymax = float(np.nanmin(values) - pad), float(np.nanmax(values) + pad)
    ticks = nice_ticks(ymin, ymax, 5)
    ymin, ymax = min(ticks), max(ticks)

    def yy(v):
        return plot_y + plot_h - (v - ymin) / (ymax - ymin) * plot_h

    parts.append(svg_text(x + w / 2, y + 20, title, 15, weight="700"))
    parts.append(line(plot_x, plot_y, plot_x, plot_y + plot_h, "#3b4350", 1.2))
    parts.append(line(plot_x, plot_y + plot_h, plot_x + plot_w, plot_y + plot_h, "#3b4350", 1.2))
    for t in ticks:
        ty = yy(t)
        parts.append(line(plot_x - 5, ty, plot_x, ty, "#3b4350", 1))
        parts.append(line(plot_x, ty, plot_x + plot_w, ty, "#e3e7ed", 0.8))
        parts.append(svg_text(plot_x - 9, ty + 4, fmt_tick(t), 10, "#4b5563", "end"))
    parts.append(svg_text(x + 16, plot_y + plot_h / 2, ylabel, 11, "#374151", "middle", rotate=-90))

    step = plot_w / len(PHYLO_ORDER)
    rng = np.random.default_rng(7)
    for i, group in enumerate(PHYLO_ORDER):
        sub = df.loc[df["Phylogroup"] == group, metric].dropna().astype(float)
        cx = plot_x + step * (i + 0.5)
        bw = min(34, step * 0.52)
        if sub.empty:
            continue
        q1, med, q3, low, high, outliers = box_stats(sub)
        color = COLORS[group]
        parts.append(line(cx, yy(low), cx, yy(q1), color, 1.2))
        parts.append(line(cx, yy(q3), cx, yy(high), color, 1.2))
        parts.append(line(cx - bw * 0.35, yy(low), cx + bw * 0.35, yy(low), color, 1.2))
        parts.append(line(cx - bw * 0.35, yy(high), cx + bw * 0.35, yy(high), color, 1.2))
        parts.append(rect(cx - bw / 2, yy(q3), bw, max(1.2, yy(q1) - yy(q3)), color, "#303846", 1, 0.72, 2))
        parts.append(line(cx - bw / 2, yy(med), cx + bw / 2, yy(med), "#172033", 1.5))

        sample = sub.to_numpy()
        if len(sample) > 80:
            sample = rng.choice(sample, size=80, replace=False)
        for v in sample:
            jitter = rng.uniform(-bw * 0.42, bw * 0.42)
            parts.append(circle(cx + jitter, yy(float(v)), 1.45, color, opacity=0.20))

        for v in outliers[:25]:
            jitter = rng.uniform(-bw * 0.25, bw * 0.25)
            parts.append(circle(cx + jitter, yy(float(v)), 1.8, "#1f2937", opacity=0.40))
        parts.append(svg_text(cx, plot_y + plot_h + 18, group, 11, COLORS[group]))
    parts.append(svg_text(plot_x + plot_w / 2, y + h - 10, "Phylogroup", 11, "#374151"))


def write_svg(path: Path, width: int, height: int, body_parts):
    content = "\n".join(body_parts)
    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">\n'
        f'<rect width="{width}" height="{height}" fill="#ffffff"/>\n{content}\n</svg>\n'
    )
    path.write_text(svg, encoding="utf-8")


def make_fig1_boxplots(merged):
    parts = [
        svg_text(600, 28, "Genome Features and GEM Topology Across Phylogroups", 21, weight="700"),
        svg_text(600, 50, "Boxplots show medians, IQR, 1.5xIQR whiskers, sampled observations and outliers", 12, "#5b6472"),
    ]
    panels = [
        ("Genome_Size_Mb", "Genome size", "Mb"),
        ("Total_CDS", "Total CDS", "Count"),
        ("GEM_Reactions", "GEM reactions", "Count"),
        ("GEM_Metabolites", "GEM metabolites", "Count"),
    ]
    for idx, (metric, title, ylabel) in enumerate(panels):
        px = 38 + (idx % 2) * 580
        py = 76 + (idx // 2) * 352
        draw_boxplot_panel(parts, merged, metric, title, ylabel, px, py, 535, 318)
    write_svg(OUT_DIR / "fig1_genome_gem_boxplots.svg", 1200, 800, parts)


def pca_svd(matrix):
    x = np.asarray(matrix, dtype=float)
    x = x - x.mean(axis=0)
    sd = x.std(axis=0, ddof=1)
    sd[sd == 0] = 1
    z = x / sd
    u, s, vt = np.linalg.svd(z, full_matrices=False)
    scores = u[:, :2] * s[:2]
    explained = (s**2) / np.sum(s**2)
    return scores, explained[:2]


def make_fig2_pca(assembly):
    cog = assembly[COG_COLS].astype(float)
    rel = cog.div(cog.sum(axis=1), axis=0).fillna(0.0)
    scores, explained = pca_svd(rel.to_numpy())
    xvals, yvals = scores[:, 0], scores[:, 1]
    width, height = 980, 720
    margin = dict(left=86, right=190, top=82, bottom=76)
    plot_x, plot_y = margin["left"], margin["top"]
    plot_w = width - margin["left"] - margin["right"]
    plot_h = height - margin["top"] - margin["bottom"]
    xpad = (xvals.max() - xvals.min()) * 0.08
    ypad = (yvals.max() - yvals.min()) * 0.08
    xmin, xmax = xvals.min() - xpad, xvals.max() + xpad
    ymin, ymax = yvals.min() - ypad, yvals.max() + ypad

    def sx(v):
        return plot_x + (v - xmin) / (xmax - xmin) * plot_w

    def sy(v):
        return plot_y + plot_h - (v - ymin) / (ymax - ymin) * plot_h

    parts = [
        svg_text(width / 2, 32, "PCA of Relative COG Functional Profiles", 22, weight="700"),
        svg_text(width / 2, 54, "Each point is one genome; COG counts were normalized within genome and z-scored by COG class", 12, "#5b6472"),
        line(plot_x, plot_y, plot_x, plot_y + plot_h, "#3b4350", 1.2),
        line(plot_x, plot_y + plot_h, plot_x + plot_w, plot_y + plot_h, "#3b4350", 1.2),
    ]
    for t in nice_ticks(float(xmin), float(xmax), 6):
        tx = sx(t)
        parts.append(line(tx, plot_y, tx, plot_y + plot_h, "#e5e7eb", 0.8))
        parts.append(svg_text(tx, plot_y + plot_h + 19, fmt_tick(t), 10, "#4b5563"))
    for t in nice_ticks(float(ymin), float(ymax), 6):
        ty = sy(t)
        parts.append(line(plot_x - 5, ty, plot_x, ty, "#3b4350", 1))
        parts.append(line(plot_x, ty, plot_x + plot_w, ty, "#e5e7eb", 0.8))
        parts.append(svg_text(plot_x - 10, ty + 4, fmt_tick(t), 10, "#4b5563", "end"))
    parts.append(svg_text(plot_x + plot_w / 2, height - 28, f"PC1 ({explained[0] * 100:.1f}% variance)", 13, "#374151"))
    parts.append(svg_text(28, plot_y + plot_h / 2, f"PC2 ({explained[1] * 100:.1f}% variance)", 13, "#374151", rotate=-90))

    for group in PHYLO_ORDER:
        idx = assembly["Phylogroup"].to_numpy() == group
        for px, py in zip(xvals[idx], yvals[idx]):
            parts.append(circle(sx(px), sy(py), 3.2, COLORS[group], "#ffffff", 0.35, 0.72))
    lx, ly = plot_x + plot_w + 34, plot_y + 28
    parts.append(svg_text(lx, ly - 15, "Phylogroup", 13, "#27313f", "start", "700"))
    for i, group in enumerate(PHYLO_ORDER):
        yy = ly + i * 26
        parts.append(circle(lx + 8, yy - 4, 6, COLORS[group]))
        parts.append(svg_text(lx + 23, yy, f"{group} (n={int((assembly['Phylogroup'] == group).sum())})", 12, "#374151", "start"))
    write_svg(OUT_DIR / "fig2_cog_pca.svg", width, height, parts)
    return explained


def pairwise_euclidean(x):
    n = x.shape[0]
    dist = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            d = float(np.linalg.norm(x[i] - x[j]))
            dist[i, j] = dist[j, i] = d
    return dist


def average_linkage_order(x):
    n = x.shape[0]
    if n <= 2:
        return list(range(n))
    base_dist = pairwise_euclidean(x)
    clusters = [[i] for i in range(n)]
    while len(clusters) > 1:
        best = None
        best_d = float("inf")
        for i in range(len(clusters)):
            for j in range(i + 1, len(clusters)):
                vals = [base_dist[a, b] for a in clusters[i] for b in clusters[j]]
                d = float(np.mean(vals))
                if d < best_d:
                    best_d = d
                    best = (i, j)
        i, j = best
        clusters[i] = clusters[i] + clusters[j]
        del clusters[j]
    return clusters[0]


def heat_color(v, vmin, vmax):
    if vmax <= vmin:
        t = 0.5
    else:
        t = (v - vmin) / (vmax - vmin)
    t = max(0, min(1, t))
    stops = [
        (0.00, (49, 54, 149)),
        (0.25, (69, 117, 180)),
        (0.50, (255, 255, 191)),
        (0.75, (253, 174, 97)),
        (1.00, (165, 0, 38)),
    ]
    for (a, ca), (b, cb) in zip(stops, stops[1:]):
        if a <= t <= b:
            f = (t - a) / (b - a)
            rgb = tuple(round(ca[k] + f * (cb[k] - ca[k])) for k in range(3))
            return f"#{rgb[0]:02x}{rgb[1]:02x}{rgb[2]:02x}"
    return "#a50026"


def make_fig2_heatmap(assembly):
    rel = assembly[COG_COLS].astype(float)
    rel = rel.div(rel.sum(axis=1), axis=0).fillna(0.0) * 100
    means = rel.join(assembly["Phylogroup"]).groupby("Phylogroup")[COG_COLS].mean().loc[PHYLO_ORDER]
    row_order = [means.index[i] for i in average_linkage_order(means.to_numpy())]
    col_order = [means.columns[i] for i in average_linkage_order(means.to_numpy().T)]
    data = means.loc[row_order, col_order]

    width, height = 1180, 570
    left, top = 115, 96
    cell_w, cell_h = 34, 34
    parts = [
        svg_text(width / 2, 32, "Clustered Heatmap of Mean Relative COG Abundance by Phylogroup", 21, weight="700"),
        svg_text(width / 2, 54, "Values are group means after within-genome COG normalization (% of assigned COG genes)", 12, "#5b6472"),
    ]
    vmin, vmax = float(np.percentile(data.to_numpy(), 2)), float(np.percentile(data.to_numpy(), 98))
    for r, group in enumerate(row_order):
        parts.append(svg_text(left - 18, top + r * cell_h + cell_h * 0.65, group, 12, COLORS[group], "end", "700"))
        for c, cog in enumerate(col_order):
            val = float(data.loc[group, cog])
            parts.append(rect(left + c * cell_w, top + r * cell_h, cell_w, cell_h, heat_color(val, vmin, vmax), "#ffffff", 0.75))
    for c, cog in enumerate(col_order):
        parts.append(svg_text(left + c * cell_w + cell_w / 2, top - 9, cog.replace("COG_", ""), 11, "#374151", rotate=-45))
    parts.append(svg_text(left + len(col_order) * cell_w / 2, top + len(row_order) * cell_h + 44, "COG category", 12, "#374151"))
    parts.append(svg_text(32, top + len(row_order) * cell_h / 2, "Phylogroup", 12, "#374151", rotate=-90))

    lx = left + len(col_order) * cell_w + 58
    ly = top + 6
    grad_h = 220
    for i in range(100):
        val = vmin + (vmax - vmin) * (1 - i / 99)
        parts.append(rect(lx, ly + i * grad_h / 100, 24, grad_h / 100 + 0.5, heat_color(val, vmin, vmax)))
    parts.append(rect(lx, ly, 24, grad_h, "none", "#3b4350", 0.8))
    parts.append(svg_text(lx + 36, ly + 4, f"{vmax:.2f}%", 10, "#4b5563", "start"))
    parts.append(svg_text(lx + 36, ly + grad_h, f"{vmin:.2f}%", 10, "#4b5563", "start"))
    parts.append(svg_text(lx - 6, ly + grad_h + 30, "Relative abundance", 11, "#374151", "start"))
    write_svg(OUT_DIR / "fig2_cog_group_heatmap.svg", width, height, parts)
    return data


def make_fig3_defense_surface(annot):
    parts = [
        svg_text(600, 28, "CRISPR/Cas Defense and Surface-Associated Proteins Across Phylogroups", 20, weight="700"),
        svg_text(600, 50, "Distributions compare phage-defense markers with secretion and membrane-associated interaction capacity", 12, "#5b6472"),
    ]
    panels = [
        ("CRISPR_Arrays", "CRISPR arrays", "Count"),
        ("Cas_Operons", "Cas operons", "Count"),
        ("Secreted_Proteins", "Secreted proteins", "Count"),
        ("Transmembrane_Proteins", "Transmembrane proteins", "Count"),
    ]
    for idx, (metric, title, ylabel) in enumerate(panels):
        px = 38 + (idx % 2) * 580
        py = 76 + (idx // 2) * 352
        draw_boxplot_panel(parts, annot, metric, title, ylabel, px, py, 535, 318)
    write_svg(OUT_DIR / "fig3_defense_surface_boxplots.svg", 1200, 800, parts)


def write_summary(assembly, annot, explained, heat_data):
    merged = assembly.merge(annot, on=["Genome_ID", "Phylogroup"], validate="one_to_one")
    summary_cols = [
        "Genome_Size_Mb",
        "Total_CDS",
        "GEM_Reactions",
        "GEM_Metabolites",
        "CRISPR_Arrays",
        "Cas_Operons",
        "Secreted_Proteins",
        "Transmembrane_Proteins",
    ]
    group_summary = merged.groupby("Phylogroup")[summary_cols].agg(["count", "mean", "median", "std"]).loc[PHYLO_ORDER]
    group_summary.to_csv(OUT_DIR / "phylogroup_metric_summary.csv")
    heat_data.to_csv(OUT_DIR / "cog_relative_abundance_group_means_clustered.csv")
    (OUT_DIR / "pca_variance_summary.txt").write_text(
        f"PC1 variance: {explained[0] * 100:.4f}%\nPC2 variance: {explained[1] * 100:.4f}%\n",
        encoding="utf-8",
    )


def main():
    OUT_DIR.mkdir(exist_ok=True)
    assembly = pd.read_csv(ASSEMBLY_CSV)
    annot = pd.read_csv(ANNOT_CSV)
    merged = assembly.merge(annot, on=["Genome_ID", "Phylogroup"], validate="one_to_one")

    make_fig1_boxplots(merged)
    explained = make_fig2_pca(assembly)
    heat_data = make_fig2_heatmap(assembly)
    make_fig3_defense_surface(annot)
    write_summary(assembly, annot, explained, heat_data)

    print("Generated:")
    for path in sorted(OUT_DIR.iterdir()):
        print(f"  {path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
