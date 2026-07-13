from pathlib import Path

from docx import Document


DOCX = Path(r"D:\pipeline\analysis_v5\02.论文\Evolutionary Divergence of Genomic Defense.first_version_clean.docx")

PROHIBITED = [
    "DeepTMHMM",
    "COBRApy-readable",
    "load-status",
    "current limits",
    "corrected phylogroup",
    "revised cohort",
    "updated benchmark",
    "original single",
    "single MG1655",
    "D, F, and G",
    "previously missing",
    "adequately represented",
    "when B2",
    "after reprocessing",
    "former",
]


def main():
    doc = Document(DOCX)
    print(f"docx={DOCX}")
    print(f"paragraphs={len(doc.paragraphs)} tables={len(doc.tables)} inline_shapes={len(doc.inline_shapes)}")
    for term in PROHIBITED:
        hits = []
        for i, para in enumerate(doc.paragraphs):
            if term.lower() in para.text.lower():
                hits.append((f"P{i}", para.text))
        for ti, table in enumerate(doc.tables):
            for ri, row in enumerate(table.rows):
                for ci, cell in enumerate(row.cells):
                    if term.lower() in cell.text.lower():
                        hits.append((f"T{ti}R{ri}C{ci}", cell.text))
        print(f"{term}: {len(hits)}")
        for loc, text in hits[:3]:
            print(f"  {loc}: {text[:220]}")
    print("captions:")
    for i, para in enumerate(doc.paragraphs):
        text = para.text
        if text.startswith("Figure ") or text.startswith("Table "):
            print(f"  P{i}: {text[:220]}")


if __name__ == "__main__":
    main()
