from pathlib import Path


SCRIPT = Path(r"D:\pipeline\analysis_v5\scripts\create_benchmark_figure.py")


def main():
    text = SCRIPT.read_text(encoding="utf-8")
    lines = []
    for line in text.splitlines():
        if line.startswith("BENCH = "):
            lines.append('BENCH = next(p.parent.parent for p in ROOT.rglob("benchmark_metrics.csv") if p.parent.name == "modle_compare")')
        else:
            lines.append(line)
    text = "\n".join(lines) + "\n"
    text = text.replace("FONT = font(26)", "FONT = font(30)")
    text = text.replace("FONT_SMALL = font(21)", "FONT_SMALL = font(24)")
    text = text.replace("FONT_TINY = font(18)", "FONT_TINY = font(21)")
    text = text.replace("FONT_BOLD = font(30, True)", "FONT_BOLD = font(35, True)")
    text = text.replace("FONT_PANEL = font(24, True)", "FONT_PANEL = font(28, True)")
    text = text.replace('pd.read_csv(BENCH / "benchmark_metrics.csv")', 'pd.read_csv(BENCH / "test_20" / "benchmark_metrics.csv")', 1)
    if "def save_png_pdf(img, png_path):" not in text:
        text = text.replace(
            "FONT_PANEL = font(28, True)\n\n\n",
            "FONT_PANEL = font(28, True)\n\n\n"
            "def save_png_pdf(img, png_path):\n"
            "    img.save(png_path, dpi=(300, 300))\n"
            "    img.convert(\"RGB\").save(png_path.with_suffix(\".pdf\"), \"PDF\", resolution=300.0)\n\n\n",
        )
    text = text.replace("img.save(OUT, dpi=(300, 300))", "save_png_pdf(img, OUT)")
    SCRIPT.write_text(text, encoding="utf-8")
    print(SCRIPT)


if __name__ == "__main__":
    main()
