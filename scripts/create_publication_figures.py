from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "00.data"
ANALYSIS = ROOT / "03.analysis"
FIG = ROOT / "01.figures"
FIG.mkdir(exist_ok=True)

PHYLO = ["A", "B2", "B1", "D", "E", "C", "F", "G", "U"]
COLORS = {
    "A": "#4C78A8",
    "B2": "#54A24B",
    "B1": "#F58518",
    "D": "#72B7B2",
    "E": "#B279A2",
    "C": "#E45756",
    "F": "#FF9DA6",
    "G": "#9D755D",
    "U": "#BAB0AC",
}


def font(size: int, bold: bool = False):
    paths = [
        "C:/Windows/Fonts/arialbd.ttf" if bold else "C:/Windows/Fonts/arial.ttf",
        "C:/Windows/Fonts/calibrib.ttf" if bold else "C:/Windows/Fonts/calibri.ttf",
    ]
    for p in paths:
        try:
            return ImageFont.truetype(p, size)
        except Exception:
            pass
    return ImageFont.load_default()


F10, F12, F14, F16, F20, F24 = [font(s) for s in [12, 14, 16, 19, 23, 28]]
F14B, F16B, F20B, F24B = [font(s, True) for s in [16, 19, 23, 28]]


def save_png_pdf(img: Image.Image, png_path: Path) -> None:
    img.save(png_path, dpi=(300, 300))
    img.convert("RGB").save(png_path.with_suffix(".pdf"), "PDF", resolution=300.0)


def hex_color(h: str):
    h = h.lstrip("#")
    return tuple(int(h[i : i + 2], 16) for i in (0, 2, 4))


def draw_text(draw, xy, text, fill=(30, 41, 59), fnt=F12, anchor=None):
    draw.text(xy, str(text), fill=fill, font=fnt, anchor=anchor)


def box_stats(vals):
    x = np.asarray(vals, dtype=float)
    x = x[np.isfinite(x)]
    q1, med, q3 = np.percentile(x, [25, 50, 75])
    iqr = q3 - q1
    low = x[x >= q1 - 1.5 * iqr].min()
    high = x[x <= q3 + 1.5 * iqr].max()
    return q1, med, q3, low, high


def nice_ticks(vmin, vmax, n=5):
    if vmax <= vmin:
        vmax = vmin + 1
    span = vmax - vmin
    raw = span / (n - 1)
    mag = 10 ** math.floor(math.log10(raw))
    norm = raw / mag
    step = mag if norm <= 1 else 2 * mag if norm <= 2 else 5 * mag if norm <= 5 else 10 * mag
    start = math.floor(vmin / step) * step
    stop = math.ceil(vmax / step) * step
    ticks = []
    t = start
    while t <= stop + step * 0.1:
        ticks.append(t)
        t += step
    return ticks


def fmt(x):
    if abs(x) >= 100:
        return f"{x:.0f}"
    if abs(x) >= 10:
        return f"{x:.1f}".rstrip("0").rstrip(".")
    return f"{x:.2f}".rstrip("0").rstrip(".")


def draw_boxplot(draw, df, metric, title, ylabel, box):
    x0, y0, x1, y1 = box
    draw_text(draw, ((x0 + x1) / 2, y0), title, fnt=F16B, anchor="ma")
    px0, py0, px1, py1 = x0 + 78, y0 + 42, x1 - 16, y1 - 54
    vals = df[metric].astype(float)
    pad = (vals.max() - vals.min()) * 0.07 + 1e-9
    ticks = nice_ticks(vals.min() - pad, vals.max() + pad)
    ymin, ymax = min(ticks), max(ticks)

    def sx(i):
        return px0 + (i + 0.5) * (px1 - px0) / len(PHYLO)

    def sy(v):
        return py1 - (v - ymin) / (ymax - ymin) * (py1 - py0)

    draw.line((px0, py0, px0, py1, px1, py1), fill=(55, 65, 81), width=2)
    for t in ticks:
        yy = sy(t)
        draw.line((px0, yy, px1, yy), fill=(229, 231, 235), width=1)
        draw_text(draw, (px0 - 8, yy), fmt(t), fill=(75, 85, 99), fnt=F10, anchor="rm")
    draw_text(draw, (x0 + 18, (py0 + py1) / 2), ylabel, fill=(55, 65, 81), fnt=F12, anchor="mm")

    rng = np.random.default_rng(5)
    for i, g in enumerate(PHYLO):
        sub = df.loc[df["Phylogroup"] == g, metric].dropna().astype(float).to_numpy()
        if len(sub) == 0:
            continue
        q1, med, q3, low, high = box_stats(sub)
        cx = sx(i)
        bw = min(48, (px1 - px0) / len(PHYLO) * 0.55)
        color = hex_color(COLORS[g])
        draw.line((cx, sy(low), cx, sy(q1)), fill=color, width=3)
        draw.line((cx, sy(q3), cx, sy(high)), fill=color, width=3)
        draw.line((cx - bw * 0.35, sy(low), cx + bw * 0.35, sy(low)), fill=color, width=3)
        draw.line((cx - bw * 0.35, sy(high), cx + bw * 0.35, sy(high)), fill=color, width=3)
        draw.rectangle((cx - bw / 2, sy(q3), cx + bw / 2, sy(q1)), outline=(31, 41, 55), fill=color, width=2)
        draw.line((cx - bw / 2, sy(med), cx + bw / 2, sy(med)), fill=(17, 24, 39), width=3)
        sample = sub if len(sub) <= 50 else rng.choice(sub, 50, replace=False)
        for v in sample:
            jitter = rng.uniform(-bw * 0.38, bw * 0.38)
            draw.ellipse((cx + jitter - 2, sy(v) - 2, cx + jitter + 2, sy(v) + 2), fill=color + (90,))
        draw_text(draw, (cx, py1 + 21), g, fill=color, fnt=F12, anchor="ma")


def heat_color(t):
    t = max(0.0, min(1.0, float(t)))
    stops = [(49, 54, 149), (69, 117, 180), (255, 255, 191), (253, 174, 97), (165, 0, 38)]
    pos = t * (len(stops) - 1)
    i = min(len(stops) - 2, int(pos))
    f = pos - i
    return tuple(round(stops[i][k] + f * (stops[i + 1][k] - stops[i][k])) for k in range(3))


def draw_heatmap(draw, data, row_labels, col_labels, box, title, value_fmt="{:.2f}"):
    x0, y0, x1, y1 = box
    draw_text(draw, ((x0 + x1) / 2, y0), title, fnt=F16B, anchor="ma")
    left, top = x0 + 110, y0 + 58
    w = (x1 - left - 20) / len(col_labels)
    h = (y1 - top - 40) / len(row_labels)
    arr = data.to_numpy(float)
    vmin, vmax = np.nanpercentile(arr, [3, 97])
    for r, row in enumerate(row_labels):
        draw_text(draw, (left - 12, top + r * h + h / 2), row, fnt=F11 if False else F12, anchor="rm")
        for c, col in enumerate(col_labels):
            val = data.loc[row, col]
            t = 0.5 if vmax <= vmin else (val - vmin) / (vmax - vmin)
            xx, yy = left + c * w, top + r * h
            draw.rectangle((xx, yy, xx + w, yy + h), fill=heat_color(t), outline=(255, 255, 255), width=1)
            if w > 44 and h > 24:
                draw_text(draw, (xx + w / 2, yy + h / 2), value_fmt.format(val), fill=(20, 24, 33), fnt=F10, anchor="mm")
    for c, col in enumerate(col_labels):
        draw_text(draw, (left + c * w + w / 2, top - 8), col, fnt=F10, anchor="ms")


def figure2():
    assembly = pd.read_csv(DATA / "genome_assembly_and_cog_metrics.csv")
    annot = pd.read_csv(DATA / "recalculated_annotations.csv")
    base = assembly.merge(annot.drop(columns=["Phylogroup"]), on="Genome_ID")
    sig = pd.read_csv(ANALYSIS / "signalp_genome_counts.csv")
    tm = pd.read_csv(ANALYSIS / "tmhmm_genome_counts.csv")
    cr = pd.read_csv(ANALYSIS / "crispr_genome_counts.csv")
    df = (
        base.merge(sig[["Genome_ID", "SignalP_per_1000_CDS"]], on="Genome_ID")
        .merge(tm[["Genome_ID", "TMHMM_Proteins_per_1000_CDS"]], on="Genome_ID")
        .merge(cr[["Genome_ID", "Unique_Spacer_Sequences"]], on="Genome_ID")
    )
    img = Image.new("RGB", (2400, 1600), "white")
    draw = ImageDraw.Draw(img, "RGBA")
    draw_text(draw, (1200, 34), "Figure 2. Revised phylogroup-stratified genome and interaction features", fnt=F24B, anchor="ma")
    metrics = [
        ("Genome_Size_Mb", "Genome size", "Mb"),
        ("GEM_Reactions", "GEM reactions", "count"),
        ("SignalP_per_1000_CDS", "SignalP positives", "per 1000 CDS"),
        ("TMHMM_Proteins_per_1000_CDS", "TMHMM proteins", "per 1000 CDS"),
        ("Unique_Spacer_Sequences", "Unique CRISPR spacers", "count"),
        ("GEM_Metabolites", "GEM metabolites", "count"),
    ]
    for i, (m, t, y) in enumerate(metrics):
        row, col = divmod(i, 3)
        draw_boxplot(draw, df, m, t, y, (40 + col * 780, 90 + row * 720, 760 + col * 780, 720 + row * 720))
    save_png_pdf(img, FIG / "Figure2_revised_functional_landscape.png")


def figure3():
    pca = pd.read_csv(ANALYSIS / "pathway_hellinger_pca_scores.csv")
    var = pd.read_csv(ANALYSIS / "pathway_hellinger_pca_variance.csv")
    perm = pd.read_csv(ANALYSIS / "pathway_hellinger_permanova_full.csv").iloc[0]
    disp = pd.read_csv(ANALYSIS / "pathway_hellinger_permdisp_full.csv").iloc[0]
    tests = pd.read_csv(ANALYSIS / "pathway_chisq_multipletest.csv").sort_values("p_bh").head(15)
    pathway = pd.read_csv(DATA / "genome_kegg_pathway_matrix.csv")
    phy = pd.read_csv(DATA / "genome_assembly_and_cog_metrics.csv")[["Genome_ID", "Phylogroup"]]
    path_base = phy.merge(pathway[["Genome_ID"] + tests["Pathway"].tolist()], on="Genome_ID")
    freq = path_base.groupby("Phylogroup")[tests["Pathway"].tolist()].apply(lambda x: (x > 0).mean()).reindex(PHYLO)
    img = Image.new("RGB", (2400, 1500), "white")
    draw = ImageDraw.Draw(img, "RGBA")
    draw_text(draw, (1200, 34), "Figure 3. KEGG pathway composition shows weak global phylogroup separation", fnt=F24B, anchor="ma")
    # PCA scatter
    x0, y0, x1, y1 = 70, 110, 1120, 1360
    draw_text(
        draw,
        ((x0 + x1) / 2, y0 - 24),
        f"Hellinger PCA: PERMANOVA R2={perm['R2']:.3f}, p={perm['p']:.3f}; PERMDISP p={disp['p']:.3f}",
        fnt=F16B,
        anchor="ma",
    )
    xs, ys = pca["PC1"].to_numpy(), pca["PC2"].to_numpy()
    xmin, xmax = xs.min(), xs.max()
    ymin, ymax = ys.min(), ys.max()
    padx, pady = (xmax - xmin) * 0.08, (ymax - ymin) * 0.08
    xmin, xmax, ymin, ymax = xmin - padx, xmax + padx, ymin - pady, ymax + pady
    def sx(v): return x0 + (v - xmin) / (xmax - xmin) * (x1 - x0)
    def sy(v): return y1 - (v - ymin) / (ymax - ymin) * (y1 - y0)
    draw.rectangle((x0, y0, x1, y1), outline=(55, 65, 81), width=2)
    for g in PHYLO:
        sub = pca[pca.Phylogroup == g]
        color = hex_color(COLORS[g])
        for _, r in sub.iterrows():
            cx, cy = sx(r.PC1), sy(r.PC2)
            draw.ellipse((cx - 4, cy - 4, cx + 4, cy + 4), fill=color + (150,), outline=(255, 255, 255, 140))
    draw_text(draw, ((x0+x1)/2, y1+38), f"PC1 ({var.loc[0,'variance_explained']*100:.1f}%)", fnt=F14, anchor="ma")
    draw_text(draw, (x0-35, (y0+y1)/2), f"PC2 ({var.loc[1,'variance_explained']*100:.1f}%)", fnt=F14, anchor="mm")
    lx, ly = x1 - 125, y0 + 20
    for i, g in enumerate(PHYLO):
        color = hex_color(COLORS[g])
        draw.rectangle((lx, ly + i * 28, lx + 18, ly + i * 28 + 18), fill=color)
        draw_text(draw, (lx + 26, ly + i * 28 + 10), g, fnt=F12, anchor="lm")
    draw_heatmap(draw, freq, PHYLO, tests["Pathway"].tolist(), (1180, 110, 2360, 1360), "Top corrected pathway-frequency differences", "{:.2f}")
    save_png_pdf(img, FIG / "Figure3_pathway_landscape.png")


def figure4():
    assembly = pd.read_csv(DATA / "genome_assembly_and_cog_metrics.csv")
    cog_tests = pd.read_csv(ANALYSIS / "cog_chisq_multipletest.csv").sort_values("p_bh").head(18)
    cogs = cog_tests["COG"].tolist()
    totals = assembly[[c for c in assembly.columns if c.startswith("COG_")]].sum(axis=1)
    rel = assembly[["Genome_ID", "Phylogroup"]].copy()
    for c in cogs:
        rel[c] = assembly[c] / totals * 100
    means = rel.groupby("Phylogroup")[cogs].mean().reindex(PHYLO)
    means.columns = [c.replace("COG_", "") for c in cogs]
    img = Image.new("RGB", (2400, 1300), "white")
    draw = ImageDraw.Draw(img, "RGBA")
    draw_text(draw, (1200, 34), "Figure 4. COG categories retain corrected phylogroup-associated shifts", fnt=F24B, anchor="ma")
    draw_heatmap(draw, means, PHYLO, list(means.columns), (80, 110, 2320, 1120), "Mean relative abundance of top BH-corrected COG categories (%)", "{:.1f}")
    y = 1180
    note = "22/26 COG categories significant after BH and Bonferroni correction; strongest maxima: COG_S in E, COG_L/C/J/F/P/E in A, COG_G/Q in B2."
    draw_text(draw, (1200, y), note, fnt=F16, anchor="ma")
    save_png_pdf(img, FIG / "Figure4_cog_corrected_heatmap.png")


def main():
    figure2()
    figure3()
    figure4()
    print("Generated publication figures in", FIG)


if __name__ == "__main__":
    main()
