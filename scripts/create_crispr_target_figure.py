from __future__ import annotations

from pathlib import Path

import pandas as pd
from PIL import Image, ImageDraw, ImageFont


ROOT = Path(r"D:\pipeline\analysis_v5")
IN = ROOT / "03.analysis" / "crispr_target_matching_analysis"
OUT = ROOT / "01.figures" / "FigureS1_crispr_target_matching.png"


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    candidates = [
        r"C:\Windows\Fonts\arialbd.ttf" if bold else r"C:\Windows\Fonts\arial.ttf",
        r"C:\Windows\Fonts\calibrib.ttf" if bold else r"C:\Windows\Fonts\calibri.ttf",
    ]
    for c in candidates:
        p = Path(c)
        if p.exists():
            return ImageFont.truetype(str(p), size)
    return ImageFont.load_default()


def draw_axes(draw: ImageDraw.ImageDraw, x0: int, y0: int, w: int, h: int, y_label: str, title: str) -> None:
    black = (30, 30, 30)
    grey = (210, 215, 220)
    draw.line((x0, y0 + h, x0 + w, y0 + h), fill=black, width=2)
    draw.line((x0, y0, x0, y0 + h), fill=black, width=2)
    for i in range(1, 5):
        y = y0 + h - int(h * i / 5)
        draw.line((x0, y, x0 + w, y), fill=grey, width=1)
        draw.text((x0 - 42, y - 9), f"{i/5:.1f}", fill=black, font=font(21))
    draw.text((x0, y0 - 46), title, fill=black, font=font(28, True))
    draw.text((x0 - 78, y0 - 28), y_label, fill=black, font=font(21))


def bar_group(draw: ImageDraw.ImageDraw, x0: int, y0: int, w: int, h: int, labels: list[str], series: list[tuple[str, list[float], tuple[int, int, int]]], title: str, y_label: str, ymax: float = 1.0) -> None:
    draw_axes(draw, x0, y0, w, h, y_label, title)
    n = len(labels)
    group_w = w / n
    bar_w = min(26, group_w / (len(series) + 1.3))
    for j, (_, vals, color) in enumerate(series):
        offset = (j - (len(series) - 1) / 2) * bar_w * 1.15
        for i, v in enumerate(vals):
            cx = x0 + group_w * (i + 0.5) + offset
            bh = int(h * min(v, ymax) / ymax)
            draw.rectangle((cx - bar_w / 2, y0 + h - bh, cx + bar_w / 2, y0 + h), fill=color)
    for i, lab in enumerate(labels):
        tx = x0 + group_w * (i + 0.5) - 8
        draw.text((tx, y0 + h + 12), lab, fill=(30, 30, 30), font=font(21))
    lx = x0 + 8
    ly = y0 + 10
    for name, _, color in series:
        draw.rectangle((lx, ly + 5, lx + 18, ly + 23), fill=color)
        draw.text((lx + 25, ly), name, fill=(30, 30, 30), font=font(20))
        ly += 26


def horizontal_bars(draw: ImageDraw.ImageDraw, x0: int, y0: int, w: int, h: int, data: pd.DataFrame, title: str) -> None:
    black = (30, 30, 30)
    draw.text((x0, y0 - 46), title, fill=black, font=font(28, True))
    maxv = max(1, int(data["n_unique_spacers"].max()))
    row_h = h / len(data)
    for idx, row in enumerate(data.itertuples(index=False)):
        y = y0 + idx * row_h
        color = (37, 99, 235) if row.Target_Type == "phage" else (220, 38, 38)
        bw = int(w * row.n_unique_spacers / maxv)
        draw.rectangle((x0 + 170, y + 4, x0 + 170 + bw, y + row_h - 4), fill=color)
        label = f"{str(row.Target_Type).title()} {row.sseqid}"
        draw.text((x0, y + 4), label[:22], fill=black, font=font(19))
        draw.text((x0 + 178 + bw, y + 4), str(int(row.n_unique_spacers)), fill=black, font=font(19))
    draw.line((x0 + 170, y0 + h + 4, x0 + 170 + w, y0 + h + 4), fill=black, width=2)
    draw.text((x0 + 170, y0 + h + 12), "Unique spacers matched", fill=black, font=font(19))


def main() -> None:
    all_group = pd.read_csv(IN / "crispr_target_matching_phylogroup_summary_all_genomes.csv")
    spacer_group = pd.read_csv(IN / "crispr_target_matching_phylogroup_summary.csv")
    top = pd.read_csv(IN / "crispr_target_matching_top_targets.csv")
    order = ["A", "B1", "B2", "C", "D", "E", "F", "G", "U"]
    all_group = all_group.set_index("Phylogroup").loc[order].reset_index()
    spacer_group = spacer_group.set_index("Phylogroup").loc[order].reset_index()

    img = Image.new("RGB", (2100, 720), "white")
    draw = ImageDraw.Draw(img)
    labels = all_group["Phylogroup"].tolist()
    bar_group(
        draw,
        95,
        110,
        520,
        390,
        labels,
        [
            ("Any spacer", all_group["fraction_genomes_with_spacers"].tolist(), (107, 114, 128)),
            ("Phage target", all_group["fraction_genomes_with_phage_targets"].tolist(), (37, 99, 235)),
            ("Plasmid target", all_group["fraction_genomes_with_plasmid_targets"].tolist(), (220, 38, 38)),
        ],
        "A. Target-bearing genomes",
        "Fraction",
    )
    bar_group(
        draw,
        720,
        110,
        420,
        390,
        labels,
        [
            ("Phage", spacer_group["median_phage_targeted_fraction"].tolist(), (37, 99, 235)),
            ("Plasmid", spacer_group["median_plasmid_targeted_fraction"].tolist(), (220, 38, 38)),
        ],
        "B. Spacer-bearing genomes",
        "Median fraction",
        ymax=0.9,
    )
    top2 = pd.concat([top[top["Target_Type"] == "phage"].head(6), top[top["Target_Type"] == "plasmid"].head(6)])
    top2 = top2.sort_values(["Target_Type", "n_unique_spacers"], ascending=[False, True])
    horizontal_bars(draw, 1240, 110, 620, 390, top2, "C. Top matched targets")
    draw.text(
        (95, 620),
        "High-confidence hits: mismatch <= 2, query coverage >= 80%, pident >= 90%, evalue <= 1e-3.",
        fill=(45, 45, 45),
        font=font(23),
    )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    img.save(OUT, dpi=(300, 300))
    img.convert("RGB").save(OUT.with_suffix(".pdf"), "PDF", resolution=300.0)
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
