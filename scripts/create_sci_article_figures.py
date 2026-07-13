from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "00.data"
ANALYSIS = ROOT / "03.analysis"
BENCH = ROOT / "04.基准测试"
SD01 = ROOT / "05.supplementary_data" / "Data_S6_sd01_fba_predictions"
OUT = ROOT / "01.figures" / "sci_article_set"
OUT.mkdir(parents=True, exist_ok=True)

PHYLO = ["A", "B2", "B1", "D", "E", "C", "F", "G", "U"]
COLORS = {
    "A": "#3B6EA8",
    "B2": "#2F9E44",
    "B1": "#E67700",
    "D": "#0CA678",
    "E": "#9C36B5",
    "C": "#C92A2A",
    "F": "#D6336C",
    "G": "#7950F2",
    "U": "#868E96",
}


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    candidates = [
        "C:/Windows/Fonts/arialbd.ttf" if bold else "C:/Windows/Fonts/arial.ttf",
        "C:/Windows/Fonts/calibrib.ttf" if bold else "C:/Windows/Fonts/calibri.ttf",
        "C:/Windows/Fonts/msyhbd.ttc" if bold else "C:/Windows/Fonts/msyh.ttc",
    ]
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, size)
        except OSError:
            continue
    return ImageFont.load_default()


F10 = font(16)
F11 = font(18)
F12 = font(20)
F13 = font(22)
F14 = font(24)
F16 = font(28)
F18 = font(32)
F20 = font(36)
F24 = font(44)
F12B = font(20, True)
F14B = font(24, True)
F16B = font(28, True)
F18B = font(32, True)
F20B = font(36, True)
F24B = font(44, True)


def rgb(hex_value: str) -> tuple[int, int, int]:
    value = hex_value.lstrip("#")
    return tuple(int(value[i : i + 2], 16) for i in (0, 2, 4))


def draw_text(draw: ImageDraw.ImageDraw, xy, text, fnt=F12, fill=(33, 37, 41), anchor=None):
    draw.text(xy, str(text), font=fnt, fill=fill, anchor=anchor)


def save(img: Image.Image, name: str) -> None:
    png = OUT / f"{name}.png"
    img.save(png, dpi=(300, 300))
    img.convert("RGB").save(png.with_suffix(".pdf"), "PDF", resolution=300.0)


def fmt(value: float) -> str:
    if not np.isfinite(value):
        return "NA"
    if abs(value) >= 100:
        return f"{value:.0f}"
    if abs(value) >= 10:
        return f"{value:.1f}".rstrip("0").rstrip(".")
    return f"{value:.2f}".rstrip("0").rstrip(".")


def compact_fmt(value: float) -> str:
    if not np.isfinite(value):
        return "NA"
    sign = "-" if value < 0 else ""
    value = abs(value)
    if value >= 1_000_000:
        return f"{sign}{value / 1_000_000:.1f}M".replace(".0M", "M")
    if value >= 1_000:
        return f"{sign}{value / 1_000:.0f}K"
    return f"{sign}{fmt(value)}"


def short_p(p: float) -> str:
    if p < 0.001:
        return "<0.001"
    return f"{p:.3f}"


def nice_ticks(vmin: float, vmax: float, n: int = 5) -> list[float]:
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
        return [0, 1]
    span = vmax - vmin
    raw = span / max(1, n - 1)
    mag = 10 ** math.floor(math.log10(raw))
    norm = raw / mag
    step = mag if norm <= 1 else 2 * mag if norm <= 2 else 5 * mag if norm <= 5 else 10 * mag
    start = math.floor(vmin / step) * step
    stop = math.ceil(vmax / step) * step
    ticks = []
    current = start
    while current <= stop + step * 0.1:
        ticks.append(current)
        current += step
    return ticks


def box_stats(values: pd.Series | np.ndarray):
    x = pd.to_numeric(pd.Series(values), errors="coerce").dropna().to_numpy(float)
    if len(x) == 0:
        return None
    q1, med, q3 = np.percentile(x, [25, 50, 75])
    iqr = q3 - q1
    low = x[x >= q1 - 1.5 * iqr].min()
    high = x[x <= q3 + 1.5 * iqr].max()
    return q1, med, q3, low, high


def draw_panel_label(draw, x: int, y: int, label: str) -> None:
    draw_text(draw, (x, y), label, F20B, fill=(0, 0, 0), anchor="lt")


def draw_boxplot(
    draw: ImageDraw.ImageDraw,
    df: pd.DataFrame,
    metric: str,
    title: str,
    ylabel: str,
    box: tuple[int, int, int, int],
    groups: list[str] = PHYLO,
) -> None:
    x0, y0, x1, y1 = box
    draw_text(draw, ((x0 + x1) / 2, y0), title, F16B, anchor="ma")
    px0, py0, px1, py1 = x0 + 100, y0 + 58, x1 - 32, y1 - 76
    draw_text(draw, (px0, py0 - 18), ylabel, F10, fill=(73, 80, 87), anchor="ls")
    values = pd.to_numeric(df[metric], errors="coerce").dropna()
    ticks = nice_ticks(values.min(), values.max(), 5)
    ymin, ymax = min(ticks), max(ticks)

    def sx(i: int) -> float:
        return px0 + (i + 0.5) * (px1 - px0) / len(groups)

    def sy(v: float) -> float:
        return py1 - (v - ymin) / (ymax - ymin) * (py1 - py0)

    draw.line((px0, py0, px0, py1, px1, py1), fill=(73, 80, 87), width=2)
    for tick in ticks:
        yy = sy(tick)
        draw.line((px0, yy, px1, yy), fill=(222, 226, 230), width=1)
        draw_text(draw, (px0 - 10, yy), fmt(tick), F10, fill=(73, 80, 87), anchor="rm")
    rng = np.random.default_rng(42)
    for i, group in enumerate(groups):
        sub = pd.to_numeric(df.loc[df["Phylogroup"] == group, metric], errors="coerce").dropna()
        stats = box_stats(sub)
        if stats is None:
            continue
        q1, med, q3, low, high = stats
        cx = sx(i)
        bw = min(62, (px1 - px0) / len(groups) * 0.58)
        color = rgb(COLORS.get(group, "#868E96"))
        draw.line((cx, sy(low), cx, sy(q1)), fill=color, width=4)
        draw.line((cx, sy(q3), cx, sy(high)), fill=color, width=4)
        draw.line((cx - bw * 0.34, sy(low), cx + bw * 0.34, sy(low)), fill=color, width=4)
        draw.line((cx - bw * 0.34, sy(high), cx + bw * 0.34, sy(high)), fill=color, width=4)
        draw.rectangle((cx - bw / 2, sy(q3), cx + bw / 2, sy(q1)), fill=color + (210,), outline=(33, 37, 41), width=2)
        draw.line((cx - bw / 2, sy(med), cx + bw / 2, sy(med)), fill=(0, 0, 0), width=4)
        sample = sub.to_numpy()
        if len(sample) > 45:
            sample = rng.choice(sample, 45, replace=False)
        for val in sample:
            jitter = rng.uniform(-bw * 0.38, bw * 0.38)
            draw.ellipse((cx + jitter - 2.5, sy(val) - 2.5, cx + jitter + 2.5, sy(val) + 2.5), fill=color + (85,))
        draw_text(draw, (cx, py1 + 28), group, F11, fill=color, anchor="ma")


def draw_bar_chart(
    draw: ImageDraw.ImageDraw,
    labels: list[str],
    values: list[float],
    title: str,
    ylabel: str,
    box: tuple[int, int, int, int],
    colors: list[str] | None = None,
    value_suffix: str = "",
    y_max: float | None = None,
) -> None:
    x0, y0, x1, y1 = box
    draw_text(draw, ((x0 + x1) / 2, y0), title, F16B, anchor="ma")
    px0, py0, px1, py1 = x0 + 96, y0 + 58, x1 - 36, y1 - 88
    draw_text(draw, (px0, py0 - 18), ylabel, F10, fill=(73, 80, 87), anchor="ls")
    vmax = y_max if y_max is not None else max(values) * 1.16 if values else 1
    ticks = nice_ticks(0, vmax, 5)
    ymax = y_max if y_max is not None else max(ticks)

    def sy(v: float) -> float:
        return py1 - v / ymax * (py1 - py0)

    draw.line((px0, py0, px0, py1, px1, py1), fill=(73, 80, 87), width=2)
    for tick in ticks:
        yy = sy(tick)
        draw.line((px0, yy, px1, yy), fill=(222, 226, 230), width=1)
        draw_text(draw, (px0 - 8, yy), fmt(tick), F10, fill=(73, 80, 87), anchor="rm")
    step = (px1 - px0) / len(labels)
    bw = step * 0.62
    for i, (label, value) in enumerate(zip(labels, values)):
        cx = px0 + (i + 0.5) * step
        color = rgb(colors[i] if colors else "#3B6EA8")
        draw.rectangle((cx - bw / 2, sy(value), cx + bw / 2, py1), fill=color, outline=(33, 37, 41), width=1)
        draw_text(draw, (cx, sy(value) - 10), f"{fmt(value)}{value_suffix}", F10, fill=(33, 37, 41), anchor="ms")
        draw_text(draw, (cx, py1 + 30), label, F10, fill=(73, 80, 87), anchor="ma")


def heat_color(t: float) -> tuple[int, int, int]:
    t = max(0.0, min(1.0, float(t)))
    stops = [(49, 54, 149), (69, 117, 180), (224, 243, 248), (254, 224, 144), (215, 48, 39)]
    pos = t * (len(stops) - 1)
    idx = min(len(stops) - 2, int(pos))
    frac = pos - idx
    return tuple(round(stops[idx][k] + frac * (stops[idx + 1][k] - stops[idx][k])) for k in range(3))


def draw_heatmap(
    draw: ImageDraw.ImageDraw,
    matrix: pd.DataFrame,
    title: str,
    box: tuple[int, int, int, int],
    value_fmt: str = "{:.2f}",
    label_font=F10,
) -> None:
    x0, y0, x1, y1 = box
    draw_text(draw, ((x0 + x1) / 2, y0), title, F16B, anchor="ma")
    row_labels = list(matrix.index)
    col_labels = list(matrix.columns)
    left = x0 + 150
    top = y0 + 68
    right = x1 - 36
    bottom = y1 - 42
    cw = (right - left) / len(col_labels)
    ch = (bottom - top) / len(row_labels)
    values = matrix.to_numpy(float)
    vmin, vmax = np.nanpercentile(values, [2, 98])
    if vmax <= vmin:
        vmax = vmin + 1

    for r, row in enumerate(row_labels):
        yy = top + r * ch
        draw_text(draw, (left - 12, yy + ch / 2), row, label_font, fill=(52, 58, 64), anchor="rm")
        for c, col in enumerate(col_labels):
            val = float(matrix.loc[row, col])
            color = heat_color((val - vmin) / (vmax - vmin))
            xx = left + c * cw
            draw.rectangle((xx, yy, xx + cw, yy + ch), fill=color, outline=(255, 255, 255), width=1)
            if cw >= 54 and ch >= 25:
                draw_text(draw, (xx + cw / 2, yy + ch / 2), value_fmt.format(val), F10, fill=(20, 20, 20), anchor="mm")
    for c, col in enumerate(col_labels):
        draw_text(draw, (left + c * cw + cw / 2, top - 10), col, F10, fill=(52, 58, 64), anchor="ms")


def draw_scatter(
    draw: ImageDraw.ImageDraw,
    df: pd.DataFrame,
    xcol: str,
    ycol: str,
    title: str,
    xlabel: str,
    ylabel: str,
    box: tuple[int, int, int, int],
) -> None:
    x0, y0, x1, y1 = box
    draw_text(draw, ((x0 + x1) / 2, y0), title, F16B, anchor="ma")
    px0, py0, px1, py1 = x0 + 96, y0 + 58, x1 - 34, y1 - 82
    draw_text(draw, (px0, py0 - 18), ylabel, F10, fill=(73, 80, 87), anchor="ls")
    xs = pd.to_numeric(df[xcol], errors="coerce")
    ys = pd.to_numeric(df[ycol], errors="coerce")
    xmin, xmax = xs.quantile(0.01), xs.quantile(0.99)
    ymin, ymax = ys.quantile(0.01), ys.quantile(0.99)
    xpad = (xmax - xmin) * 0.08
    ypad = (ymax - ymin) * 0.08
    xmin, xmax = xmin - xpad, xmax + xpad
    ymin, ymax = ymin - ypad, ymax + ypad

    def sx(v: float) -> float:
        return px0 + (v - xmin) / (xmax - xmin) * (px1 - px0)

    def sy(v: float) -> float:
        return py1 - (v - ymin) / (ymax - ymin) * (py1 - py0)

    draw.rectangle((px0, py0, px1, py1), outline=(73, 80, 87), width=2)
    for tick in nice_ticks(float(xmin), float(xmax), 5):
        xx = sx(tick)
        if px0 <= xx <= px1:
            draw.line((xx, py0, xx, py1), fill=(233, 236, 239), width=1)
            draw.line((xx, py1, xx, py1 + 8), fill=(73, 80, 87), width=2)
            draw_text(draw, (xx, py1 + 24), compact_fmt(tick), F10, fill=(73, 80, 87), anchor="ma")
    for tick in nice_ticks(float(ymin), float(ymax), 5):
        yy = sy(tick)
        if py0 <= yy <= py1:
            draw.line((px0, yy, px1, yy), fill=(233, 236, 239), width=1)
            draw.line((px0 - 8, yy, px0, yy), fill=(73, 80, 87), width=2)
            draw_text(draw, (px0 - 12, yy), compact_fmt(tick), F10, fill=(73, 80, 87), anchor="rm")
    if "Phylogroup" in df.columns:
        draw_groups = [(group, df[df["Phylogroup"] == group], rgb(COLORS[group])) for group in PHYLO]
    else:
        draw_groups = [("", df, rgb("#3B6EA8"))]
    for _, sub, color in draw_groups:
        for _, row in sub.iterrows():
            x, y = float(row[xcol]), float(row[ycol])
            if not np.isfinite(x) or not np.isfinite(y):
                continue
            xx, yy = sx(x), sy(y)
            if px0 <= xx <= px1 and py0 <= yy <= py1:
                draw.ellipse((xx - 4, yy - 4, xx + 4, yy + 4), fill=color + (110,), outline=color + (180,))
    draw_text(draw, ((px0 + px1) / 2, y1 - 28), xlabel, F11, fill=(73, 80, 87), anchor="ma")


def draw_phylo_legend(draw, x: int, y: int) -> None:
    for i, group in enumerate(PHYLO):
        xx = x + i * 82
        color = rgb(COLORS[group])
        draw.ellipse((xx, y, xx + 18, y + 18), fill=color)
        draw_text(draw, (xx + 26, y + 9), group, F10, fill=(52, 58, 64), anchor="lm")


def figure1() -> None:
    sampling = pd.read_csv(ANALYSIS / "phylogroup_proportional_sampling_plan.csv")
    sampling = sampling[sampling["Source_Phylogroup_Label"].isin(PHYLO)].set_index("Source_Phylogroup_Label").reindex(PHYLO)
    benchmark = pd.read_csv(BENCH / "test_20" / "benchmark_metrics.csv")
    benchmark["Speedup"] = benchmark["Default_Time_Sec"] / benchmark["Optimized_Time_Sec"]
    benchmark["Reaction_Retention"] = benchmark["Optimized_Reactions"] / benchmark["Default_Reactions"] * 100
    ref = pd.read_csv(BENCH / "modle_compare" / "reaction_set_compare" / "reaction_set_compare_summary.csv")
    ref["Reaction_Retention"] = ref["B_optimized_reactions"] / ref["A_default_reactions"] * 100

    img = Image.new("RGB", (2600, 1700), "white")
    draw = ImageDraw.Draw(img, "RGBA")
    draw_text(draw, (1300, 34), "SCI Figure 1. Cohort design and workflow performance", F24B, anchor="ma")

    draw_panel_label(draw, 58, 92, "A")
    draw_bar_chart(
        draw,
        PHYLO,
        sampling["Planned_Sample_Count"].astype(float).tolist(),
        "Proportional phylogroup-stratified cohort",
        "Sampled genomes",
        (80, 105, 1220, 760),
        [COLORS[g] for g in PHYLO],
    )

    draw_panel_label(draw, 1310, 92, "B")
    draw_scatter(
        draw,
        benchmark.rename(columns={"Actual_N50": "N50"}),
        "N50",
        "Speedup",
        "Assembly quality versus acceleration",
        "Contig N50",
        "Speedup",
        (1330, 105, 2500, 760),
    )
    median_speed = benchmark["Speedup"].median()
    draw_text(draw, (2420, 180), f"median {median_speed:.1f}x", F18B, fill=(194, 65, 12), anchor="ra")

    draw_panel_label(draw, 58, 850, "C")
    runtime_df = pd.DataFrame(
        {
            "Variant": ["Default"] * len(benchmark) + ["Optimized"] * len(benchmark),
            "Runtime_h": list(benchmark["Default_Time_Sec"] / 3600) + list(benchmark["Optimized_Time_Sec"] / 3600),
        }
    )
    draw_bar_chart(
        draw,
        ["Default", "Optimized"],
        runtime_df.groupby("Variant")["Runtime_h"].median().reindex(["Default", "Optimized"]).tolist(),
        "Median runtime across 20 stress-test genomes",
        "Hours per genome",
        (80, 865, 1220, 1530),
        ["#868E96", "#2F9E44"],
    )
    draw_text(draw, (145, 1588), "Bars show medians; full distribution is reported in the benchmark table.", F11, fill=(73, 80, 87))

    draw_panel_label(draw, 1310, 850, "D")
    draw_bar_chart(
        draw,
        ["20 genomes", "Reference models"],
        [benchmark["Reaction_Retention"].median(), ref["Reaction_Retention"].median()],
        "Reaction-set retention after optimization",
        "Retained reactions (%)",
        (1330, 865, 2500, 1530),
        ["#3B6EA8", "#7950F2"],
        "%",
    )
    save(img, "SCI_Figure1_cohort_and_performance")


def figure2() -> None:
    assembly = pd.read_csv(DATA / "genome_assembly_and_cog_metrics.csv")
    annot = pd.read_csv(DATA / "recalculated_annotations.csv")
    df = assembly.merge(annot.drop(columns=["Phylogroup"]), on="Genome_ID")
    img = Image.new("RGB", (2600, 1700), "white")
    draw = ImageDraw.Draw(img, "RGBA")
    draw_text(draw, (1300, 34), "SCI Figure 2. Genome architecture and GEM topology across phylogroups", F24B, anchor="ma")
    panels = [
        ("Genome_Size_Mb", "Genome size", "Mb"),
        ("Total_CDS", "Protein-coding capacity", "CDS"),
        ("GEM_Reactions", "GEM reaction space", "Reactions"),
        ("GEM_Metabolites", "GEM metabolite space", "Metabolites"),
    ]
    boxes = [(80, 120, 1260, 760), (1340, 120, 2520, 760), (80, 890, 1260, 1530), (1340, 890, 2520, 1530)]
    for idx, ((metric, title, ylabel), box) in enumerate(zip(panels, boxes)):
        draw_panel_label(draw, box[0] - 22, box[1] - 28, chr(65 + idx))
        draw_boxplot(draw, df, metric, title, ylabel, box)
    draw_phylo_legend(draw, 880, 1608)
    save(img, "SCI_Figure2_genome_and_gem_topology")


def figure3() -> None:
    pca = pd.read_csv(ANALYSIS / "pathway_hellinger_pca_scores.csv")
    var = pd.read_csv(ANALYSIS / "pathway_hellinger_pca_variance.csv").set_index("component")["variance_explained"]
    perma = pd.read_csv(ANALYSIS / "pathway_hellinger_permanova_full.csv").iloc[0]
    disp_test = pd.read_csv(ANALYSIS / "pathway_hellinger_permdisp_full.csv").iloc[0]
    disp = pd.read_csv(ANALYSIS / "pathway_hellinger_permdisp_distances.csv")
    pathway_tests = pd.read_csv(ANALYSIS / "pathway_chisq_multipletest.csv").sort_values("p_bh").head(15)
    pathway = pd.read_csv(DATA / "genome_kegg_pathway_matrix.csv")
    phy = pd.read_csv(DATA / "genome_assembly_and_cog_metrics.csv")[["Genome_ID", "Phylogroup"]]
    heat = phy.merge(pathway[["Genome_ID"] + pathway_tests["Pathway"].tolist()], on="Genome_ID")
    heat = heat.groupby("Phylogroup")[pathway_tests["Pathway"].tolist()].apply(lambda x: (x > 0).mean()).reindex(PHYLO).T
    residual = pd.read_csv(ANALYSIS / "ani_residual_phylogroup_effects.csv")
    residual = residual.set_index("feature_space")

    img = Image.new("RGB", (2800, 1900), "white")
    draw = ImageDraw.Draw(img, "RGBA")
    draw_text(draw, (1400, 34), "SCI Figure 3. KEGG pathway structure with dispersion-aware interpretation", F24B, anchor="ma")

    draw_panel_label(draw, 58, 100, "A")
    draw_scatter(
        draw,
        pca,
        "PC1",
        "PC2",
        "Hellinger PCA of KEGG pathway profiles",
        f"PC1 ({var['PC1'] * 100:.1f}%)",
        f"PC2 ({var['PC2'] * 100:.1f}%)",
        (80, 120, 1320, 940),
    )
    draw_text(
        draw,
        (1120, 180),
        f"PERMANOVA R2={float(perma['R2']):.3f}, p={short_p(float(perma['p']))}",
        F14B,
        fill=(33, 37, 41),
        anchor="ra",
    )
    draw_text(draw, (1120, 220), f"PERMDISP F={float(disp_test['F']):.1f}, p={short_p(float(disp_test['p']))}", F14B, fill=(194, 65, 12), anchor="ra")

    draw_panel_label(draw, 1448, 100, "B")
    draw_boxplot(
        draw,
        disp.rename(columns={"Hellinger_Distance_To_Group_Centroid": "Distance"}),
        "Distance",
        "Within-group dispersion",
        "Distance to centroid",
        (1470, 120, 2700, 940),
    )

    draw_panel_label(draw, 58, 1040, "C")
    draw_heatmap(draw, heat, "Top BH-corrected KEGG pathway frequencies", (80, 1060, 1700, 1810), "{:.2f}", F10)

    draw_panel_label(draw, 1798, 1040, "D")
    r2_labels = ["KEGG", "COG", "GEM", "CRISPR", "Signal/TM"]
    r2_index = ["KEGG_Hellinger", "COG_Hellinger", "GEM_topology_zEuclidean", "CRISPR_spacer_Jaccard", "SignalP_TMHMM_zEuclidean"]
    values = [float(residual.loc[i, "unique_phylogroup_R2_after_ANI"]) * 100 for i in r2_index if i in residual.index]
    draw_bar_chart(
        draw,
        r2_labels[: len(values)],
        values,
        "Residual phylogroup effect after ANI",
        "Unique R2 (%)",
        (1820, 1060, 2700, 1810),
        ["#3B6EA8", "#E67700", "#2F9E44", "#C92A2A", "#7950F2"],
        "%",
    )
    save(img, "SCI_Figure3_kegg_ani_dispersion")


def figure4() -> None:
    assembly = pd.read_csv(DATA / "genome_assembly_and_cog_metrics.csv")
    cog_cols = [c for c in assembly.columns if c.startswith("COG_")]
    cog_tests = pd.read_csv(ANALYSIS / "cog_chisq_multipletest.csv").sort_values("p_bh").head(16)
    cog_rel = assembly[["Phylogroup"] + cog_cols].copy()
    totals = cog_rel[cog_cols].sum(axis=1).replace(0, np.nan)
    cog_rel[cog_cols] = cog_rel[cog_cols].div(totals, axis=0) * 100
    cog_heat = cog_rel.groupby("Phylogroup")[cog_tests["COG"].tolist()].mean().reindex(PHYLO).T
    crispr = pd.read_csv(ANALYSIS / "crispr_genome_counts.csv")
    signalp = pd.read_csv(ANALYSIS / "signalp_genome_counts.csv")
    tmhmm = pd.read_csv(ANALYSIS / "tmhmm_genome_counts.csv")
    aux = signalp[["Genome_ID", "Phylogroup", "SignalP_per_1000_CDS"]].merge(
        tmhmm[["Genome_ID", "TMHMM_Proteins_per_1000_CDS"]], on="Genome_ID", how="left"
    )

    img = Image.new("RGB", (2800, 1850), "white")
    draw = ImageDraw.Draw(img, "RGBA")
    draw_text(draw, (1400, 34), "SCI Figure 4. COG categories and auxiliary host-interaction layers", F24B, anchor="ma")

    draw_panel_label(draw, 58, 100, "A")
    draw_heatmap(draw, cog_heat, "BH-corrected COG category relative abundance", (80, 120, 1740, 1730), "{:.1f}", F10)

    draw_panel_label(draw, 1818, 100, "B")
    draw_boxplot(draw, crispr, "Unique_Spacer_Sequences", "CRISPR spacer repertoire", "Unique spacers", (1840, 120, 2720, 640))

    draw_panel_label(draw, 1818, 720, "C")
    draw_boxplot(draw, aux, "SignalP_per_1000_CDS", "Secreted-protein burden", "per 1000 CDS", (1840, 740, 2720, 1250))

    draw_panel_label(draw, 1818, 1320, "D")
    draw_boxplot(draw, aux, "TMHMM_Proteins_per_1000_CDS", "Transmembrane-protein burden", "per 1000 CDS", (1840, 1340, 2720, 1730))
    save(img, "SCI_Figure4_cog_crispr_secretome")


def figure5() -> None:
    selected = pd.read_csv(SD01 / "selected_models_770_carbon_fba_accuracy_metrics.csv")
    available = pd.read_csv(SD01 / "available_benchmark_770_fba_accuracy_metrics.csv")
    selected["Model"] = selected["Requested_Model"] + " " + selected["Condition_Class"].str.replace("C_", "", regex=False)
    selected_labels = selected["Requested_Model"].unique().tolist()
    selected_values = selected.groupby("Requested_Model")["Balanced_Accuracy"].mean().reindex(selected_labels) * 100
    available["Variant_Label"] = available["Variant"].str.replace("gapseq_", "", regex=False).str.title()
    avail_summary = available.groupby("Variant_Label")[["Balanced_Accuracy", "MCC"]].mean().reindex(["Default", "Optimized"])

    img = Image.new("RGB", (2400, 1500), "white")
    draw = ImageDraw.Draw(img, "RGBA")
    draw_text(draw, (1200, 34), "SCI Figure 5. SD01 model-consistency benchmarking", F24B, anchor="ma")

    draw_panel_label(draw, 58, 100, "A")
    draw_bar_chart(
        draw,
        selected_labels,
        selected_values.tolist(),
        "Curated reference models against SD01 matrix",
        "Balanced accuracy (%)",
        (80, 120, 1160, 1320),
        ["#3B6EA8", "#2F9E44", "#E67700", "#9C36B5", "#C92A2A"][: len(selected_labels)],
        "%",
        y_max=100,
    )

    draw_panel_label(draw, 1268, 100, "B")
    draw_bar_chart(
        draw,
        avail_summary.index.tolist(),
        (avail_summary["Balanced_Accuracy"] * 100).tolist(),
        "gapseq model consistency",
        "Balanced accuracy (%)",
        (1290, 120, 2320, 720),
        ["#868E96", "#2F9E44"],
        "%",
    )

    draw_panel_label(draw, 1268, 820, "C")
    draw_bar_chart(
        draw,
        avail_summary.index.tolist(),
        avail_summary["MCC"].tolist(),
        "Matthews correlation coefficient",
        "MCC",
        (1290, 840, 2320, 1320),
        ["#868E96", "#2F9E44"],
    )
    draw_text(
        draw,
        (1350, 1390),
        "Benchmarking is interpreted as model-consistency against the Monk et al. SD01\ncomputational growth matrix, not as independent Biolog validation.",
        F13,
        fill=(73, 80, 87),
    )
    save(img, "SCI_Figure5_sd01_model_consistency")


def main() -> None:
    figure1()
    figure2()
    figure3()
    figure4()
    figure5()
    print(f"Saved SCI figure set to {OUT}")


if __name__ == "__main__":
    main()
