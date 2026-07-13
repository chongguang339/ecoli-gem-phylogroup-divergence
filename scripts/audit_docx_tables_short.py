from pathlib import Path

from docx import Document


DOCX = Path(r"D:\pipeline\analysis_v5\02.论文\Evolutionary Divergence of Genomic Defense.docx")


def main():
    doc = Document(DOCX)
    print(f"paragraphs={len(doc.paragraphs)} tables={len(doc.tables)} figures={len(doc.inline_shapes)}")
    for ti, table in enumerate(doc.tables):
        print(f"TABLE {ti}")
        for ri, row in enumerate(table.rows[:8]):
            print(f"{ri}: " + " | ".join(cell.text.replace("\n", " ")[:180] for cell in row.cells))


if __name__ == "__main__":
    main()
