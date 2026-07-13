from pathlib import Path

from docx import Document


DOCX = Path(r"D:\pipeline\analysis_v5\02.论文\Evolutionary Divergence of Genomic Defense.docx")

STARTS = (
    "Background:",
    "Results:",
    "Conclusions:",
    "Here we",
    "From this candidate",
    "Phylogroup assignment",
    "CRISPR-Cas immune",
    "Computational benchmarking",
    "Structural reference",
    "For the SD01",
    "KEGG pathway matrices",
    "All statistical",
    "The KEGG matrix",
    "Single-pathway",
    "The computational benchmark",
    "Several caveats",
    "In conclusion",
    "2.",
    "3.",
    "4.",
    "5.",
)


def main():
    doc = Document(DOCX)
    for i, para in enumerate(doc.paragraphs):
        text = para.text.strip()
        if i < 210 and text.startswith(STARTS):
            print(f"P{i}: {text}")


if __name__ == "__main__":
    main()
