from __future__ import annotations

from pathlib import Path

import pandas as pd
from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "03.analysis" / "bacdive_mg1655_fba"
FIG = ROOT / "01.figures"
FIG.mkdir(exist_ok=True)


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


F12, F14, F16, F20, F24 = [font(s) for s in [18, 22, 26, 32, 40]]
F12B, F14B, F16B, F20B, F24B = [font(s, True) for s in [18, 22, 26, 32, 40]]


def txt(draw, xy, value, fnt=F12, fill=(33, 37, 41), anchor=None):
    draw.text(xy, str(value), font=fnt, fill=fill, anchor=anchor)


def save(img: Image.Image, name: str):
    png = FIG / f"{name}.png"
    img.save(png, dpi=(300, 300))
    img.convert("RGB").save(png.with_suffix(".pdf"), "PDF", resolution=300)


def draw_bar_panel(draw, metrics: pd.DataFrame):
    x0, y0, x1, y1 = 90, 160, 1180, 900
    txt(draw, (x0, y0 - 58), "A", F24B)
    txt(draw, ((x0 + x1) / 2, y0 - 52), "BacDive MG1655 phenotype agreement", F20B, anchor="ma")
    plot = metrics[(metrics["condition"] == "aerobic") & metrics["panel"].isin(["carbon_assimilation_subset", "all_bacdive_panel"])].copy()
    order = {
        "MG1655_gapseq_default": 0,
        "MG1655_gapseq_optimized": 1,
        "iML1515": 2,
    }
    panel_order = {"carbon_assimilation_subset": 0, "all_bacdive_panel": 1}
    plot["sort_model"] = plot["model_id"].map(order)
    plot["sort_panel"] = plot["panel"].map(panel_order)
    plot = plot.sort_values(["sort_model", "sort_panel"])
    plot["label"] = plot["model_id"].replace({"MG1655_gapseq_default": "default", "MG1655_gapseq_optimized": "optimized"}) + "\n" + plot["panel"].replace(
        {"carbon_assimilation_subset": "carbon subset", "all_bacdive_panel": "all panel"}
    )
    labels = plot["label"].tolist()
    values = (plot["Balanced_Accuracy"] * 100).tolist()
    px0, py0, px1, py1 = x0 + 110, y0 + 60, x1 - 30, y1 - 110
    draw.line((px0, py0, px0, py1, px1, py1), fill=(73, 80, 87), width=2)
    for tick in [0, 25, 50, 75, 100]:
        yy = py1 - tick / 100 * (py1 - py0)
        draw.line((px0, yy, px1, yy), fill=(222, 226, 230), width=1)
        txt(draw, (px0 - 12, yy), f"{tick}", F12, fill=(73, 80, 87), anchor="rm")
    txt(draw, (px0, py0 - 20), "Balanced accuracy (%)", F12, fill=(73, 80, 87), anchor="ls")
    colors = [(134, 142, 150), (173, 181, 189), (47, 158, 68), (91, 183, 109), (59, 110, 168), (93, 135, 188)]
    step = (px1 - px0) / len(labels)
    for i, (label, val) in enumerate(zip(labels, values)):
        cx = px0 + (i + 0.5) * step
        bw = step * 0.55
        yy = py1 - val / 100 * (py1 - py0)
        draw.rectangle((cx - bw / 2, yy, cx + bw / 2, py1), fill=colors[i], outline=(33, 37, 41))
        txt(draw, (cx, yy - 10), f"{val:.1f}", F12B, anchor="ms")
        txt(draw, (cx, py1 + 28), label, F12, fill=(73, 80, 87), anchor="ma")


def draw_matrix_panel(draw, pred: pd.DataFrame):
    x0, y0, x1, y1 = 1280, 160, 2350, 900
    txt(draw, (x0, y0 - 58), "B", F24B)
    txt(draw, ((x0 + x1) / 2, y0 - 52), "Carbon-assimilation subset predictions", F20B, anchor="ma")
    sub = pred[(pred["condition"] == "aerobic") & (pred["panel"] == "carbon_assimilation_subset")].copy()
    sub = sub[sub["status"] == "ok"]
    sources = list(dict.fromkeys(sub["carbon_source"]))
    models = ["MG1655_gapseq_default", "MG1655_gapseq_optimized", "iML1515"]
    model_labels = ["Default", "Optimized", "iML1515"]
    left, top = x0 + 270, y0 + 80
    cw, ch = 190, 52
    for c, label in enumerate(model_labels):
        txt(draw, (left + c * cw + cw / 2, top - 18), label, F14B, anchor="ms")
    for r, source in enumerate(sources):
        yy = top + r * ch
        txt(draw, (left - 14, yy + ch / 2), source, F12, fill=(73, 80, 87), anchor="rm")
        for c, model in enumerate(models):
            row = sub[(sub["model_id"] == model) & (sub["carbon_source"] == source)]
            if row.empty:
                color, mark = (220, 220, 220), "NA"
            else:
                row = row.iloc[0]
                truth, pred_val = int(row["truth"]), int(row["prediction"])
                if truth == pred_val:
                    color, mark = (47, 158, 68), "OK"
                elif truth == 0 and pred_val == 1:
                    color, mark = (245, 159, 0), "FP"
                else:
                    color, mark = (201, 42, 42), "FN"
            xx = left + c * cw
            draw.rectangle((xx, yy, xx + cw - 8, yy + ch - 8), fill=color, outline=(255, 255, 255))
            txt(draw, (xx + cw / 2 - 4, yy + ch / 2 - 4), mark, F12B, fill=(255, 255, 255), anchor="mm")
    txt(draw, (left, y1 - 22), "OK, correct; FP, false positive growth prediction; FN, false negative.", F12, fill=(73, 80, 87))


def main():
    metrics = pd.read_csv(ANALYSIS / "bacdive_mg1655_fba_accuracy_metrics.csv")
    pred = pd.read_csv(ANALYSIS / "bacdive_mg1655_fba_predictions_long.csv")
    img = Image.new("RGB", (2400, 1050), "white")
    draw = ImageDraw.Draw(img, "RGBA")
    txt(draw, (1200, 42), "Supplementary Figure. Independent BacDive MG1655 phenotype check", F24B, anchor="ma")
    txt(
        draw,
        (1200, 92),
        "BacDive/API 20E entries are used as a small experimental phenotype panel; SD01 remains a model-consistency benchmark.",
        F14,
        fill=(73, 80, 87),
        anchor="ma",
    )
    draw_bar_panel(draw, metrics)
    draw_matrix_panel(draw, pred)
    save(img, "FigureS_bacdive_mg1655_fba")


if __name__ == "__main__":
    main()
