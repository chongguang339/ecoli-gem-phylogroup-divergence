from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
MEMOTE = ROOT / "memote_reports"
DATA = ROOT / "00.data"
ANALYSIS = ROOT / "03.analysis" / "gem_quality_control"
SUPP = ROOT / "05.supplementary_data" / "Data_S9_gem_quality_control"
FIG = ROOT / "01.figures"

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


def font(size: int, bold: bool = False):
    for path in [
        "C:/Windows/Fonts/arialbd.ttf" if bold else "C:/Windows/Fonts/arial.ttf",
        "C:/Windows/Fonts/calibrib.ttf" if bold else "C:/Windows/Fonts/calibri.ttf",
    ]:
        try:
            return ImageFont.truetype(path, size)
        except OSError:
            pass
    return ImageFont.load_default()


F10 = font(16)
F11 = font(18)
F12 = font(20)
F14 = font(24)
F16 = font(28)
F20 = font(36)
F24 = font(44)
F12B = font(20, True)
F14B = font(24, True)
F16B = font(28, True)
F20B = font(36, True)
F24B = font(44, True)


def rgb(hex_value: str) -> tuple[int, int, int]:
    h = hex_value.lstrip("#")
    return tuple(int(h[i : i + 2], 16) for i in (0, 2, 4))


def txt(draw, xy, value, fnt=F12, fill=(33, 37, 41), anchor=None):
    draw.text(xy, str(value), font=fnt, fill=fill, anchor=anchor)


def save(img: Image.Image, name: str) -> None:
    png = FIG / f"{name}.png"
    img.save(png, dpi=(300, 300))
    img.convert("RGB").save(png.with_suffix(".pdf"), "PDF", resolution=300)


def box_stats(values):
    x = pd.to_numeric(pd.Series(values), errors="coerce").dropna().to_numpy(float)
    q1, med, q3 = np.percentile(x, [25, 50, 75])
    iqr = q3 - q1
    lo = x[x >= q1 - 1.5 * iqr].min()
    hi = x[x <= q3 + 1.5 * iqr].max()
    return q1, med, q3, lo, hi


def draw_boxplot(draw, df, metric, title, ylabel, box):
    x0, y0, x1, y1 = box
    txt(draw, ((x0 + x1) / 2, y0), title, F16B, anchor="ma")
    px0, py0, px1, py1 = x0 + 105, y0 + 62, x1 - 35, y1 - 84
    values = pd.to_numeric(df[metric], errors="coerce")
    ymin = np.floor(values.min() * 10) / 10
    ymax = np.ceil(values.max() * 10) / 10
    if metric.startswith("Total"):
        ymin = np.floor(values.min() / 100) * 100
        ymax = np.ceil(values.max() / 100) * 100
    ticks = np.linspace(ymin, ymax, 5)

    def sx(i):
        return px0 + (i + 0.5) * (px1 - px0) / len(PHYLO)

    def sy(v):
        return py1 - (v - ymin) / (ymax - ymin) * (py1 - py0)

    draw.line((px0, py0, px0, py1, px1, py1), fill=(73, 80, 87), width=2)
    txt(draw, (px0, py0 - 18), ylabel, F10, fill=(73, 80, 87), anchor="ls")
    for t in ticks:
        yy = sy(t)
        draw.line((px0, yy, px1, yy), fill=(222, 226, 230), width=1)
        label = f"{t:.0f}" if metric.startswith("Total") else f"{t:.1f}"
        txt(draw, (px0 - 10, yy), label, F10, fill=(73, 80, 87), anchor="rm")
    rng = np.random.default_rng(7)
    for i, g in enumerate(PHYLO):
        sub = pd.to_numeric(df.loc[df["Phylogroup"] == g, metric], errors="coerce").dropna()
        if sub.empty:
            continue
        q1, med, q3, lo, hi = box_stats(sub)
        cx = sx(i)
        bw = min(56, (px1 - px0) / len(PHYLO) * 0.55)
        color = rgb(COLORS[g])
        draw.line((cx, sy(lo), cx, sy(q1)), fill=color, width=3)
        draw.line((cx, sy(q3), cx, sy(hi)), fill=color, width=3)
        draw.line((cx - bw * 0.35, sy(lo), cx + bw * 0.35, sy(lo)), fill=color, width=3)
        draw.line((cx - bw * 0.35, sy(hi), cx + bw * 0.35, sy(hi)), fill=color, width=3)
        draw.rectangle((cx - bw / 2, sy(q3), cx + bw / 2, sy(q1)), fill=color, outline=(33, 37, 41), width=2)
        draw.line((cx - bw / 2, sy(med), cx + bw / 2, sy(med)), fill=(0, 0, 0), width=3)
        sample = sub.to_numpy()
        if len(sample) > 35:
            sample = rng.choice(sample, 35, replace=False)
        for v in sample:
            jitter = rng.uniform(-bw * 0.35, bw * 0.35)
            draw.ellipse((cx + jitter - 2, sy(v) - 2, cx + jitter + 2, sy(v) + 2), fill=color + (75,))
        txt(draw, (cx, py1 + 28), g, F10, fill=color, anchor="ma")


def draw_summary_text(draw, summary):
    x0, y0 = 105, 1390
    lines = [
        f"n = {int(summary.loc['count', 'Total_Reactions'])} GEMs with COBRApy topology QC",
        f"Mass-balanced reactions: {summary.loc['mean', 'Mass_Balanced_Ratio(%)']:.2f}% +/- {summary.loc['std', 'Mass_Balanced_Ratio(%)']:.2f}%",
        f"Blocked reactions: {summary.loc['mean', 'Blocked_Reactions_Ratio(%)']:.2f}% +/- {summary.loc['std', 'Blocked_Reactions_Ratio(%)']:.2f}%",
        f"Dead-end metabolites: {summary.loc['mean', 'Deadend_Metabolites_Ratio(%)']:.2f}% +/- {summary.loc['std', 'Deadend_Metabolites_Ratio(%)']:.2f}%",
        "Representative MEMOTE HTML reports are provided as interactive supplementary files.",
    ]
    for i, line in enumerate(lines):
        txt(draw, (x0, y0 + i * 34), line, F14 if i == 0 else F12, fill=(52, 58, 64))


def main() -> None:
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    SUPP.mkdir(parents=True, exist_ok=True)
    qc = pd.read_csv(MEMOTE / "cobrapy_qc_results.csv")
    phy = pd.read_csv(DATA / "recalculated_annotations.csv")[["Genome_ID", "Phylogroup"]]
    merged = qc.merge(phy, on="Genome_ID", how="left")
    if merged["Phylogroup"].isna().any():
        raise ValueError("Some QC rows did not map to phylogroup labels.")

    metrics = [
        "Total_Reactions",
        "Total_Metabolites",
        "Mass_Balanced_Ratio(%)",
        "Blocked_Reactions_Ratio(%)",
        "Deadend_Metabolites_Ratio(%)",
    ]
    summary = merged[metrics].describe().loc[["count", "mean", "std", "min", "25%", "50%", "75%", "max"]]
    by_phylo = merged.groupby("Phylogroup")[metrics].agg(["count", "mean", "std", "median", "min", "max"]).reindex(PHYLO)
    html_inventory = pd.DataFrame(
        [
            {"report_file": p.name, "size_bytes": p.stat().st_size, "supplementary_role": "Representative MEMOTE HTML report"}
            for p in sorted(MEMOTE.glob("*.html"))
        ]
    )

    merged.to_csv(ANALYSIS / "cobrapy_qc_results_with_phylogroup.csv", index=False)
    summary.to_csv(ANALYSIS / "cobrapy_qc_summary_statistics.csv")
    by_phylo.to_csv(ANALYSIS / "cobrapy_qc_by_phylogroup.csv")
    html_inventory.to_csv(ANALYSIS / "memote_html_report_inventory.csv", index=False)

    for p in [ANALYSIS / "cobrapy_qc_results_with_phylogroup.csv", ANALYSIS / "cobrapy_qc_summary_statistics.csv", ANALYSIS / "cobrapy_qc_by_phylogroup.csv", ANALYSIS / "memote_html_report_inventory.csv"]:
        (SUPP / p.name).write_bytes(p.read_bytes())
    for p in sorted(MEMOTE.glob("*.html")):
        (SUPP / p.name).write_bytes(p.read_bytes())

    img = Image.new("RGB", (2600, 1600), "white")
    draw = ImageDraw.Draw(img, "RGBA")
    txt(draw, (1300, 34), "Supplementary Figure. Cohort-level GEM quality control", F24B, anchor="ma")
    txt(draw, (1300, 86), "COBRApy topology and stoichiometric checks across the 1,000 reconstructed GEMs", F14, fill=(73, 80, 87), anchor="ma")
    panels = [
        ("Total_Reactions", "Reaction count", "reactions"),
        ("Total_Metabolites", "Metabolite count", "metabolites"),
        ("Mass_Balanced_Ratio(%)", "Mass-balanced reactions", "%"),
        ("Blocked_Reactions_Ratio(%)", "Blocked reactions", "%"),
        ("Deadend_Metabolites_Ratio(%)", "Dead-end metabolites", "%"),
    ]
    boxes = [
        (70, 140, 860, 720),
        (910, 140, 1700, 720),
        (1750, 140, 2540, 720),
        (480, 800, 1270, 1360),
        (1320, 800, 2110, 1360),
    ]
    for (metric, title, ylabel), box in zip(panels, boxes):
        draw_boxplot(draw, merged, metric, title, ylabel, box)
    draw_summary_text(draw, summary)
    save(img, "FigureS_gem_quality_control")

    print(summary.round(3).to_string())
    print(f"Wrote {ANALYSIS}")


if __name__ == "__main__":
    main()
