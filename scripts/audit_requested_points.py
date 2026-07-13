from pathlib import Path

from docx import Document


DOCX = Path(r"D:\pipeline\analysis_v5\02.论文\Evolutionary Divergence of Genomic Defense.first_version_clean.docx")

TERMS = [
    "90",
    "ninety",
    "minute",
    "min",
    "representative",
    "quantile",
    "genome size",
    "GC",
    "N50",
    "iML1515",
    "four paired",
    "four curated",
    "median 12.76",
    "collection",
    "selected",
]


def main():
    doc = Document(DOCX)
    for i, para in enumerate(doc.paragraphs):
        text = para.text
        low = text.lower()
        if any(term.lower() in low for term in TERMS):
            print(f"P{i}: {text}")


if __name__ == "__main__":
    main()
