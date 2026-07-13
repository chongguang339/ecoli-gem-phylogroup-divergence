from pathlib import Path

from docx import Document


DOCX = Path(r"D:\pipeline\analysis_v5\02.论文\Evolutionary Divergence of Genomic Defense.docx")


def main() -> None:
    doc = Document(DOCX)
    print("paragraphs", len(doc.paragraphs), "tables", len(doc.tables), "sections", len(doc.sections))
    print("--- headings/fig/table/ref ---")
    for i, para in enumerate(doc.paragraphs):
        text = para.text.strip()
        if not text:
            continue
        if (
            text[:2].isdigit()
            or text.startswith(
                (
                    "Abstract",
                    "Keywords",
                    "Figure",
                    "Table",
                    "Supplementary",
                    "References",
                    "Availability",
                    "Competing",
                    "Funding",
                    "Authors",
                )
            )
        ):
            print(f"{i}: {text[:240]}")

    print("--- table summaries ---")
    for ti, table in enumerate(doc.tables):
        print("table", ti + 1, "rows", len(table.rows), "cols", len(table.columns))
        for row in table.rows[:3]:
            print(" | ".join(cell.text.replace("\n", " / ")[:90] for cell in row.cells))
        print("...")

    print("--- refs ---")
    refs = []
    in_refs = False
    for para in doc.paragraphs:
        text = para.text.strip()
        if text == "References":
            in_refs = True
            continue
        if in_refs and text:
            refs.append(text)
    print("refs", len(refs))
    for ref in refs[:6] + refs[-6:]:
        print(ref[:360])

    print("--- risk term hits ---")
    full_text = "\n".join(para.text for para in doc.paragraphs)
    terms = [
        "updated",
        "re-run",
        "rerun",
        "previous",
        "current limitation",
        "DeepTMHMM",
        "functional accuracy",
        "gold standard",
        "excellent",
        "can explain",
        "56.1",
        "965",
        "B2群得到充分代表",
        "更新",
        "此前",
        "重跑",
        "the authors",
        "our previous",
        "unclear",
        "may be",
    ]
    for term in terms:
        print(term, full_text.lower().find(term.lower()))


if __name__ == "__main__":
    main()
