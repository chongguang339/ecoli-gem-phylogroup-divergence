from pathlib import Path

from docx import Document


DOCX = Path(r"D:\pipeline\analysis_v5\02.论文\Evolutionary Divergence of Genomic Defense.docx")
TERMS = [
    "functional accuracy",
    "phenotype accuracy",
    "experimental phenotype",
    "four available",
    "4 available",
    "representative cohort",
    "quantile",
    "reference-strain",
    "GC-content",
    "model-load",
    "current re-run",
    "corrected",
    "revised",
    "updated benchmark",
    "single MG1655",
    "[Author",
    "[Funding",
    "[repository URL",
]


def main():
    doc = Document(DOCX)
    for term in TERMS:
        hits = []
        for i, para in enumerate(doc.paragraphs):
            if term.lower() in para.text.lower():
                hits.append((f"P{i}", para.text[:300]))
        for ti, table in enumerate(doc.tables):
            for ri, row in enumerate(table.rows):
                for ci, cell in enumerate(row.cells):
                    if term.lower() in cell.text.lower():
                        hits.append((f"T{ti}R{ri}C{ci}", cell.text[:300]))
        print(f"{term}: {len(hits)}")
        for loc, text in hits[:6]:
            print(f"  {loc}: {text}")


if __name__ == "__main__":
    main()
