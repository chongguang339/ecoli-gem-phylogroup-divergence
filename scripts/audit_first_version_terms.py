import sys
from docx import Document
from pathlib import Path


DOCX = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(r"D:\pipeline\analysis_v5\02.论文\Evolutionary Divergence of Genomic Defense.docx")

TERMS = [
    "corrected",
    "revised",
    "updated",
    "previous",
    "earlier",
    "original",
    "initial",
    "former",
    "missing",
    "DeepTMHMM",
    "COBRApy-readable",
    "load-status",
    "single MG1655",
    "D, F, and G",
    "current limits",
    "B2 is adequately",
    "expanded from",
    "adequately represented",
]


def main():
    doc = Document(DOCX)
    for i, para in enumerate(doc.paragraphs):
        text = para.text
        low = text.lower()
        if any(term.lower() in low for term in TERMS):
            print(f"P{i}: {text}")
    for ti, table in enumerate(doc.tables):
        for ri, row in enumerate(table.rows):
            for ci, cell in enumerate(row.cells):
                text = cell.text
                low = text.lower()
                if any(term.lower() in low for term in TERMS):
                    print(f"T{ti}R{ri}C{ci}: {text}")


if __name__ == "__main__":
    main()
