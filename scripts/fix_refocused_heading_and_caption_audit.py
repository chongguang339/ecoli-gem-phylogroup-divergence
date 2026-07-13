from pathlib import Path

from docx import Document


DOCX = Path(r"D:\pipeline\analysis_v5\02.论文\Evolutionary Divergence of Genomic Defense.docx")


def main():
    doc = Document(DOCX)
    for para in doc.paragraphs:
        if para.text.strip() == "3.1 GEM Workflow Benchmarking Shows":
            para.clear()
            para.add_run("3.1 GEM Workflow Benchmarking Shows 11-18-Fold Runtime Acceleration with Compact, Reference-Aware Models")
    doc.save(DOCX)
    doc = Document(DOCX)
    for i, para in enumerate(doc.paragraphs):
        text = para.text.strip()
        if text.startswith("Figure") or text.startswith("3.1 "):
            print(f"P{i}: {text}")


if __name__ == "__main__":
    main()
