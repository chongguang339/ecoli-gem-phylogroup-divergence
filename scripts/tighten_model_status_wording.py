from pathlib import Path

from docx import Document


DOCX = Path(r"D:\pipeline\analysis_v5\02.论文\Evolutionary Divergence of Genomic Defense.first_version_clean.docx")


def main():
    doc = Document(DOCX)
    for para in doc.paragraphs:
        if "model-load status" in para.text:
            text = para.text.replace("model-load status", "model inclusion and simulation status")
            para.clear()
            para.add_run(text)
    doc.save(DOCX)
    print(f"saved {DOCX}")


if __name__ == "__main__":
    main()
