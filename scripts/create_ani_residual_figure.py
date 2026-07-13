from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont


ROOT = Path(r"D:\pipeline\analysis_v5")
ANALYSIS = ROOT / "03.analysis"
FIG = ROOT / "01.figures"
FIG.mkdir(exist_ok=True)
OUT = FIG / "Figure5_ani_residual_phylogroup.png"

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


def font(size, bold=False):
    for p in [
        "C:/Windows/Fonts/arialbd.ttf" if bold else "C:/Windows/Fonts/arial.ttf",
        "C:/Windows/Fonts/calibrib.ttf" if bold else "C:/Windows/Fonts/calibri.ttf",
    ]:
        try:
            return ImageFont.truetype(p, size)
        except Exception:
            pass
    return ImageFont.load_default()


F12, F14, F16, F20, F24 = [font(s) for s in [14, 16, 19, 23, 28]]
F12B, F16B, F20B, F24B = [font(s, True) for s in [14, 19, 23, 28]]


def save_png_pdf(img, png_path):
    img.save(png_path, dpi=(300, 300))
    img.convert("RGB").save(png_path.with_suffix(".pdf"), "PDF", resolution=300.0)


def hex_color(h):
    h = h.lstrip("#")
    return tuple(int(h[i : i + 2], 16) for i in (0, 2, 4))


def text(draw, xy, s, f=F14, fill=(31, 41, 55), anchor=None):
    draw.text(xy, str(s), font=f, fill=fill, anchor=anchor)


def main():
    meta = pd.read_csv(ANALYSIS / "ani_residual_analysis_genomes.csv")
    scores = pd.read_csv(ANALYSIS / "ani_pcoa_scores.csv")
    mantel = pd.read_csv(ANALYSIS / "ani_mantel_distance_correlations.csv")
    residual = pd.read_csv(ANALYSIS / "ani_residual_phylogroup_effects.csv")
    qc = pd.read_csv(ANALYSIS / "ani_distance_qc.csv").iloc[0]
    scores = scores.rename(columns={scores.columns[0]: "Genome_ID"}).merge(meta[["Genome_ID", "Phylogroup"]], on="Genome_ID")

    img = Image.new("RGB", (2400, 1500), "white")
    draw = ImageDraw.Draw(img, "RGBA")
    text(draw, (1200, 34), "Figure 5. Core-genome ANI structure and residual phylogroup effects", F24B, anchor="ma")
    text(
        draw,
        (1200, 68),
        f"ANI distance matrix: n=1000, median distance={qc['distance_median_upper_sym']:.4f}; residual models used 20 ANI PCoA axes ({qc['ANI_PC_variance_explained_used']*100:.0f}% variance).",
        F14,
        fill=(75, 85, 99),
        anchor="ma",
    )

    # Panel A: ANI PCoA.
    x0, y0, x1, y1 = 90, 150, 930, 1320
    text(draw, (x0, y0 - 45), "A  Core-genome ANI PCoA", F20B)
    xs, ys = scores["ANI_PC1"].to_numpy(), scores["ANI_PC2"].to_numpy()
    xmin, xmax = xs.min(), xs.max()
    ymin, ymax = ys.min(), ys.max()
    padx, pady = (xmax - xmin) * 0.08, (ymax - ymin) * 0.08
    xmin, xmax, ymin, ymax = xmin - padx, xmax + padx, ymin - pady, ymax + pady
    sx = lambda v: x0 + (v - xmin) / (xmax - xmin) * (x1 - x0)
    sy = lambda v: y1 - (v - ymin) / (ymax - ymin) * (y1 - y0)
    draw.rectangle((x0, y0, x1, y1), outline=(55, 65, 81), width=2)
    for g in PHYLO:
        sub = scores[scores["Phylogroup"] == g]
        c = hex_color(COLORS[g])
        for r in sub.itertuples(index=False):
            cx, cy = sx(getattr(r, "ANI_PC1")), sy(getattr(r, "ANI_PC2"))
            draw.ellipse((cx - 4, cy - 4, cx + 4, cy + 4), fill=c + (145,), outline=(255, 255, 255, 110))
    text(draw, ((x0 + x1) / 2, y1 + 38), "ANI PC1", F14, anchor="ma")
    lx, ly = x1 - 95, y0 + 20
    for i, g in enumerate(PHYLO):
        c = hex_color(COLORS[g])
        draw.rectangle((lx, ly + i * 28, lx + 18, ly + i * 28 + 18), fill=c)
        text(draw, (lx + 26, ly + i * 28 + 10), g, F12, anchor="lm")

    # Panel B: Mantel correlations.
    x0, y0, x1, y1 = 1050, 150, 2300, 670
    text(draw, (x0, y0 - 45), "B  Functional distance correlation with ANI distance", F20B)
    draw.line((x0, y1, x1, y1), fill=(55, 65, 81), width=2)
    draw.line((x0, y0, x0, y1), fill=(55, 65, 81), width=2)
    labels = ["KEGG", "COG", "GEM", "CRISPR", "SignalP/TMHMM"]
    vals = mantel["mantel_spearman_r"].to_numpy()
    maxv = 0.8
    for t in np.linspace(0, maxv, 5):
        yy = y1 - t / maxv * (y1 - y0)
        draw.line((x0, yy, x1, yy), fill=(229, 231, 235), width=1)
        text(draw, (x0 - 8, yy), f"{t:.1f}", F12, fill=(75, 85, 99), anchor="rm")
    bw = (x1 - x0) / len(vals) * 0.55
    for i, (lab, val) in enumerate(zip(labels, vals)):
        cx = x0 + (i + 0.5) * (x1 - x0) / len(vals)
        yy = y1 - val / maxv * (y1 - y0)
        draw.rectangle((cx - bw / 2, yy, cx + bw / 2, y1), fill=(44, 123, 182), outline=(31, 41, 55))
        text(draw, (cx, yy - 10), f"{val:.2f}", F12B, anchor="mb")
        text(draw, (cx, y1 + 18), lab, F12, anchor="ma")

    # Panel C: variance partition.
    x0, y0, x1, y1 = 1050, 850, 2300, 1320
    text(draw, (x0, y0 - 45), "C  Residual phylogroup effect after ANI PCs", F20B)
    draw.line((x0, y1, x1, y1), fill=(55, 65, 81), width=2)
    draw.line((x0, y0, x0, y1), fill=(55, 65, 81), width=2)
    ani = residual["ANI_PC_R2"].to_numpy()
    uniq = residual["unique_phylogroup_R2_after_ANI"].to_numpy()
    maxv = 0.75
    for t in np.linspace(0, maxv, 6):
        yy = y1 - t / maxv * (y1 - y0)
        draw.line((x0, yy, x1, yy), fill=(229, 231, 235), width=1)
        text(draw, (x0 - 8, yy), f"{t:.2f}", F12, fill=(75, 85, 99), anchor="rm")
    group_w = (x1 - x0) / len(labels)
    bw = group_w * 0.25
    for i, lab in enumerate(labels):
        cx = x0 + (i + 0.5) * group_w
        y_ani = y1 - ani[i] / maxv * (y1 - y0)
        y_uni = y1 - uniq[i] / maxv * (y1 - y0)
        draw.rectangle((cx - bw - 5, y_ani, cx - 5, y1), fill=(117, 107, 177), outline=(31, 41, 55))
        draw.rectangle((cx + 5, y_uni, cx + bw + 5, y1), fill=(230, 126, 34), outline=(31, 41, 55))
        text(draw, (cx - bw / 2 - 5, y_ani - 8), f"{ani[i]:.2f}", F12, anchor="mb")
        text(draw, (cx + bw / 2 + 5, y_uni - 8), f"{uniq[i]:.2f}", F12, anchor="mb")
        text(draw, (cx, y1 + 18), lab, F12, anchor="ma")
    draw.rectangle((x0 + 830, y0 + 18, x0 + 850, y0 + 38), fill=(117, 107, 177))
    text(draw, (x0 + 860, y0 + 28), "ANI PCs", F12, anchor="lm")
    draw.rectangle((x0 + 960, y0 + 18, x0 + 980, y0 + 38), fill=(230, 126, 34))
    text(draw, (x0 + 990, y0 + 28), "Phylogroup residual", F12, anchor="lm")

    save_png_pdf(img, OUT)
    print(OUT)


if __name__ == "__main__":
    main()
