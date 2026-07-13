from pathlib import Path
import math

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont


ROOT = Path(r"D:\pipeline\analysis_v5")
BENCH = next(p.parent.parent for p in ROOT.rglob("benchmark_metrics.csv") if p.parent.name == "modle_compare")
FIG_DIR = ROOT / "01.figures"
OUT = FIG_DIR / "Figure1_benchmark_distribution.png"


def font(size, bold=False):
    candidates = [
        r"C:\Windows\Fonts\arialbd.ttf" if bold else r"C:\Windows\Fonts\arial.ttf",
        r"C:\Windows\Fonts\calibrib.ttf" if bold else r"C:\Windows\Fonts\calibri.ttf",
    ]
    for c in candidates:
        if Path(c).exists():
            return ImageFont.truetype(c, size)
    return ImageFont.load_default()


FONT = font(30)
FONT_SMALL = font(24)
FONT_TINY = font(21)
FONT_BOLD = font(35, True)
FONT_PANEL = font(28, True)


def save_png_pdf(img, png_path):
    img.save(png_path, dpi=(300, 300))
    img.convert("RGB").save(png_path.with_suffix(".pdf"), "PDF", resolution=300.0)


def percentile_text(values):
    arr = np.asarray(values, dtype=float)
    q1, med, q3 = np.percentile(arr, [25, 50, 75])
    return q1, med, q3


def scale(value, vmin, vmax, lo, hi, log=False):
    if log:
        value = math.log10(value)
        vmin = math.log10(vmin)
        vmax = math.log10(vmax)
    if vmax == vmin:
        return (lo + hi) / 2
    return lo + (value - vmin) * (hi - lo) / (vmax - vmin)


def axes(draw, box, xlabel, ylabel, xticks=None, yticks=None, xlog=False, yfmt="{:.0f}"):
    x0, y0, x1, y1 = box
    draw.line([(x0, y1), (x1, y1)], fill=(40, 45, 52), width=2)
    draw.line([(x0, y0), (x0, y1)], fill=(40, 45, 52), width=2)
    if xticks:
        for v, lab in xticks:
            x = scale(v, xticks[0][0], xticks[-1][0], x0, x1, log=xlog)
            draw.line([(x, y1), (x, y1 + 7)], fill=(40, 45, 52), width=2)
            draw.text((x, y1 + 10), lab, anchor="mt", font=FONT_TINY, fill=(55, 60, 68))
    if yticks:
        ymin, ymax = yticks[0], yticks[-1]
        for v in yticks:
            y = scale(v, ymin, ymax, y1, y0)
            draw.line([(x0 - 7, y), (x0, y)], fill=(40, 45, 52), width=2)
            draw.text((x0 - 12, y), yfmt.format(v), anchor="rm", font=FONT_TINY, fill=(55, 60, 68))
            draw.line([(x0, y), (x1, y)], fill=(231, 234, 239), width=1)
    draw.text(((x0 + x1) / 2, y1 + 48), xlabel, anchor="mm", font=FONT_SMALL, fill=(35, 40, 48))


def draw_boxplot(draw, values, x, y0, y1, vmin, vmax, color, label=None):
    arr = np.asarray(values, dtype=float)
    q1, med, q3 = np.percentile(arr, [25, 50, 75])
    mn, mx = arr.min(), arr.max()
    def sy(v):
        return scale(v, vmin, vmax, y1, y0)
    draw.line([(x, sy(mn)), (x, sy(mx))], fill=color, width=3)
    draw.line([(x - 18, sy(mn)), (x + 18, sy(mn))], fill=color, width=3)
    draw.line([(x - 18, sy(mx)), (x + 18, sy(mx))], fill=color, width=3)
    draw.rectangle([x - 42, sy(q3), x + 42, sy(q1)], fill=tuple(min(255, c + 45) for c in color), outline=color, width=3)
    draw.line([(x - 42, sy(med)), (x + 42, sy(med))], fill=(20, 25, 32), width=3)
    for i, v in enumerate(arr):
        jitter = ((i % 7) - 3) * 5
        draw.ellipse([x + jitter - 4, sy(v) - 4, x + jitter + 4, sy(v) + 4], fill=color)
    if label:
        draw.text((x, y1 + 12), label, anchor="mt", font=FONT_TINY, fill=(55, 60, 68))


def draw():
    FIG_DIR.mkdir(exist_ok=True)
    df = pd.read_csv(BENCH / "test_20" / "benchmark_metrics.csv")
    models = pd.read_csv(BENCH / "modle_compare" / "benchmark_metrics.csv")
    comp = pd.read_csv(BENCH / "modle_compare" / "reaction_set_compare" / "reaction_set_compare_summary.csv")

    for frame in (df, models):
        frame["speedup"] = frame["Default_Time_Sec"] / frame["Optimized_Time_Sec"]
        frame["reaction_retention"] = frame["Optimized_Reactions"] / frame["Default_Reactions"]

    df.to_csv(BENCH / "stratified20_derived_metrics.csv", index=False)
    models.to_csv(BENCH / "model4_derived_metrics.csv", index=False)

    img = Image.new("RGB", (2400, 1750), (255, 255, 255))
    draw = ImageDraw.Draw(img)
    draw.text((90, 50), "Figure 1. Stratified gapseq runtime and reference-model benchmarking", font=FONT_BOLD, fill=(25, 30, 38))
    draw.text((90, 92), "Benchmarking reports distributions across low-quality assemblies and available curated-reference genomes.", font=FONT_SMALL, fill=(90, 96, 106))

    accent = (32, 120, 142)
    orange = (220, 118, 45)
    green = (64, 142, 82)
    purple = (117, 95, 170)
    gray = (85, 92, 105)

    # Panel A
    box_a = (180, 210, 1040, 690)
    draw.text((90, 165), "A", font=FONT_PANEL, fill=(20, 25, 32))
    draw.text((130, 165), "N50-stratified speedup distribution (n=20)", font=FONT_PANEL, fill=(20, 25, 32))
    axes(draw, box_a, "Assembly N50 (bp, log scale)", "Speedup (default / optimized)",
         xticks=[(22535, "22k"), (100000, "100k"), (675909, "676k")],
         yticks=[14, 15, 16, 17, 18], xlog=True, yfmt="{:.0f}x")
    xmin, xmax = df["Actual_N50"].min(), df["Actual_N50"].max()
    ymin, ymax = 14, 18
    for _, r in df.iterrows():
        x = scale(r["Actual_N50"], xmin, xmax, box_a[0], box_a[2], log=True)
        y = scale(r["speedup"], ymin, ymax, box_a[3], box_a[1])
        draw.ellipse([x - 8, y - 8, x + 8, y + 8], fill=accent, outline=(20, 70, 88), width=1)
    _, med, _ = percentile_text(df["speedup"])
    ymed = scale(med, ymin, ymax, box_a[3], box_a[1])
    draw.line([(box_a[0], ymed), (box_a[2], ymed)], fill=orange, width=3)
    draw.text((box_a[2] - 6, ymed - 8), f"median {med:.2f}x", anchor="rs", font=FONT_TINY, fill=orange)

    # Panel B
    box_b = (1320, 210, 2220, 690)
    draw.text((1230, 165), "B", font=FONT_PANEL, fill=(20, 25, 32))
    draw.text((1270, 165), "Wall-clock runtime", font=FONT_PANEL, fill=(20, 25, 32))
    axes(draw, box_b, "Pipeline", "Runtime (s)",
         yticks=[0, 5000, 10000, 15000, 20000, 25000], yfmt="{:.0f}")
    draw_boxplot(draw, df["Default_Time_Sec"], 1580, box_b[1], box_b[3], 0, 25000, gray, "Default")
    draw_boxplot(draw, df["Optimized_Time_Sec"], 1960, box_b[1], box_b[3], 0, 25000, green, "Optimized")
    draw.text((box_b[0] + 20, box_b[1] + 20), "median: 19,821.5 s -> 1,260.5 s", font=FONT_TINY, fill=(70, 76, 86))

    # Panel C
    box_c = (180, 930, 1040, 1410)
    draw.text((90, 885), "C", font=FONT_PANEL, fill=(20, 25, 32))
    draw.text((130, 885), "Model compactness across stratified subset", font=FONT_PANEL, fill=(20, 25, 32))
    axes(draw, box_c, "Genome rank by optimized reaction count", "Reaction count",
         yticks=[2200, 2400, 2600, 2800, 3000, 3200], yfmt="{:.0f}")
    ordered = df.sort_values("Optimized_Reactions").reset_index(drop=True)
    for i, r in ordered.iterrows():
        x = scale(i, 0, len(ordered) - 1, box_c[0] + 25, box_c[2] - 25)
        yd = scale(r["Default_Reactions"], 2200, 3200, box_c[3], box_c[1])
        yo = scale(r["Optimized_Reactions"], 2200, 3200, box_c[3], box_c[1])
        draw.line([(x, yd), (x, yo)], fill=(205, 211, 219), width=2)
        draw.ellipse([x - 5, yd - 5, x + 5, yd + 5], fill=gray)
        draw.ellipse([x - 6, yo - 6, x + 6, yo + 6], fill=green)
    rr = df["reaction_retention"].median() * 100
    draw.text((box_c[0] + 20, box_c[1] + 20), f"median reaction retention: {rr:.1f}%", font=FONT_TINY, fill=(70, 76, 86))
    draw.rectangle([box_c[2] - 190, box_c[1] + 18, box_c[2] - 172, box_c[1] + 36], fill=gray)
    draw.text((box_c[2] - 160, box_c[1] + 16), "Default", font=FONT_TINY, fill=(70, 76, 86))
    draw.rectangle([box_c[2] - 190, box_c[1] + 48, box_c[2] - 172, box_c[1] + 66], fill=green)
    draw.text((box_c[2] - 160, box_c[1] + 46), "Optimized", font=FONT_TINY, fill=(70, 76, 86))

    # Panel D
    box_d = (1320, 930, 2220, 1410)
    draw.text((1230, 885), "D", font=FONT_PANEL, fill=(20, 25, 32))
    draw.text((1270, 885), "Reference-model benchmark (available n=4)", font=FONT_PANEL, fill=(20, 25, 32))
    axes(draw, box_d, "", "Fraction of dropped reactions in BiGG",
         yticks=[0, 5, 10, 15, 20], yfmt="{:.0f}%")
    labels = list(comp["model_id"])
    vals = comp["Dropped_present_in_C_fraction"].to_numpy() * 100
    for i, (lab, val) in enumerate(zip(labels, vals)):
        x = scale(i, 0, len(labels) - 1, box_d[0] + 70, box_d[2] - 70)
        y = scale(val, 0, 20, box_d[3], box_d[1])
        draw.rectangle([x - 36, y, x + 36, box_d[3]], fill=purple)
        label_lines = [lab]
        for j, line in enumerate(label_lines):
            draw.text((x, box_d[3] + 13 + j * 19), line, anchor="mt", font=FONT_TINY, fill=(55, 60, 68))
        draw.text((x, y - 10), f"{val:.1f}%", anchor="mb", font=FONT_TINY, fill=(55, 60, 68))
    draw.text((box_d[0] + 20, box_d[3] + 58), "Most removed reactions were not in the matched BiGG reference.", font=FONT_TINY, fill=(70, 76, 86))

    # Footer
    s1 = df["speedup"].to_numpy()
    s2 = models["speedup"].to_numpy()
    footer = (
        f"Stratified subset: speedup median {np.median(s1):.2f}x "
        f"(IQR {np.percentile(s1,25):.2f}-{np.percentile(s1,75):.2f}; range {s1.min():.2f}-{s1.max():.2f}). "
        f"Reference-model genomes: median {np.median(s2):.2f}x (range {s2.min():.2f}-{s2.max():.2f})."
    )
    draw.text((90, 1630), footer, font=FONT_SMALL, fill=(65, 72, 82))
    draw.text((90, 1665), "Note: SD01 carbon-source simulations are model-consistency benchmarks, not independent experimental Biolog validation.", font=FONT_TINY, fill=(112, 118, 128))

    save_png_pdf(img, OUT)
    print(OUT)


if __name__ == "__main__":
    draw()
