from pathlib import Path

from docx import Document


DOCX = Path(r"D:\pipeline\analysis_v5\02.论文\Evolutionary Divergence of Genomic Defense.docx")


def main():
    doc = Document(DOCX)
    print(f"docx={DOCX}")
    print(f"paragraphs={len(doc.paragraphs)} tables={len(doc.tables)} figures={len(doc.inline_shapes)}")
    print(f"title={doc.paragraphs[0].text}")
    for i, para in enumerate(doc.paragraphs):
        text = para.text.strip()
        if text.startswith("Figure ") or text.startswith("Table ") or text.startswith("3.1 "):
            print(f"P{i}: {text[:260]}")


if __name__ == "__main__":
    main()
