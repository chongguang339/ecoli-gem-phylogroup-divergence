from __future__ import annotations

from pathlib import Path
import re
import zipfile

from docx import Document


ROOT = Path(r"D:\pipeline\analysis_v5")
DOCX = ROOT / "02.论文" / "Evolutionary Divergence of Genomic Defense.docx"
TEXT_INDEX = ROOT / "02.论文" / "docx_text_index_after_submission_precheck.txt"
AUDIT = ROOT / "03.analysis" / "submission_precheck_audit.md"


def set_para_text(para, text: str) -> None:
    for run in para.runs:
        run.text = ""
    if para.runs:
        para.runs[0].text = text
    else:
        para.add_run(text)


def replace_contains(doc: Document, needle: str, repl: str) -> int:
    n = 0
    for para in doc.paragraphs:
        if needle in para.text:
            set_para_text(para, para.text.replace(needle, repl))
            n += 1
    return n


def replace_start(doc: Document, prefix: str, text: str) -> int:
    n = 0
    for para in doc.paragraphs:
        if para.text.strip().startswith(prefix):
            set_para_text(para, text)
            n += 1
    return n


def find_para(doc: Document, prefix: str):
    for para in doc.paragraphs:
        if para.text.strip().startswith(prefix):
            return para
    raise RuntimeError(f"Anchor not found: {prefix}")


def move_figure_block_before(doc: Document, caption_prefix: str, target_prefix: str) -> None:
    caption = find_para(doc, caption_prefix)
    prev = caption._p.getprevious()
    block = []
    if prev is not None:
        block.append(prev)
    block.append(caption._p)
    for el in block:
        parent = el.getparent()
        if parent is not None:
            parent.remove(el)
    target = find_para(doc, target_prefix)._p
    for el in reversed(block):
        target.addprevious(el)


def cite_and_reorder_refs(doc: Document) -> None:
    replace_contains(
        doc,
        "Because SD01 Table 5 is explicitly an in silico growth-screen matrix from Monk et al. rather than an independent experimental Biolog table",
        "Because SD01 Table 5 is explicitly an in silico growth-screen matrix from Monk et al. rather than an independent experimental Biolog table [19]",
    )
    replace_contains(doc, "FBA was run with COBRApy v0.31.1 using the GLPK solver", "FBA was run with COBRApy v0.31.1 [20] using the GLPK solver")
    replace_contains(
        doc,
        "Public BiGG SBML files were used for these curated models.",
        "Public BiGG SBML files were used for these curated models [21, 22].",
    )
    replace_contains(
        doc,
        "Second, the new SD01 770-condition benchmark tests whether curated BiGG models and paired gapseq reconstructions reproduce the Monk et al. in silico growth-screen matrix after automated exchange-ID mapping.",
        "Second, the new SD01 770-condition benchmark tests whether curated BiGG models and paired gapseq reconstructions reproduce the Monk et al. in silico growth-screen matrix after automated exchange-ID mapping [23].",
    )
    replace_contains(
        doc,
        "extension of this framework to other Enterobacteriaceae using the same Enterobacterales fast-pass database;",
        "extension of this framework to other Enterobacteriaceae and microbial ecosystem modeling using the same Enterobacterales fast-pass database [24];",
    )
    replace_contains(
        doc,
        "application of the metabolic distance landscape for antimicrobial susceptibility prediction, given that metabolic network topology is a known correlate of antibiotic sensitivity;",
        "application of the metabolic distance landscape for antimicrobial susceptibility prediction and microbial pathogenesis studies, given the known links between bacterial metabolism and pathogenicity [25];",
    )

    refs = {
        21: "[21] Norsigian CJ, Pusarla N, McConn JL, et al. BiGG Models 2020: multi-strain genome-scale models and expansion across the phylogenetic tree. Nucleic Acids Res. 2020;48(D1):D402-D406.",
        22: "[22] Orth JD, Conrad TM, Na J, Lerman JA, Nam H, Feist AM, Palsson BO. A comprehensive genome-scale reconstruction of Escherichia coli metabolism-2011. Mol Syst Biol. 2011;7:535.",
        23: "[23] Monk JM, Charusanti P, Aziz RK, Lerman JA, Premyodhin N, Orth JD, Feist AM, Palsson BO. Genome-scale metabolic reconstructions of multiple Escherichia coli strains highlight strain-specific adaptations to nutritional environments. Proc Natl Acad Sci USA. 2013;110(50):20338-20343.",
        24: "[24] Garza DR, Dutilh BE. From cultured to uncultured genome sequences: metagenomics and modeling microbial ecosystems. Cell Mol Life Sci. 2015;72(22):4287-4308.",
        25: "[25] Rohmer L, Hocquet D, Miller SI. Are pathogenic bacteria just looking for food? Metabolism and microbial pathogenesis. Trends Microbiol. 2011;19(7):341-348.",
    }
    for num, text in refs.items():
        replace_start(doc, f"[{num}]", text)


def audit_doc(doc: Document) -> str:
    body = "\n".join(p.text for p in doc.paragraphs[:216])
    refs = "\n".join(p.text for p in doc.paragraphs[216:])
    cited: list[int] = []
    for m in re.finditer(r"\[([0-9,\-\s]+)\]", body):
        for part in re.split(r",\s*", m.group(1)):
            if "-" in part:
                a, b = map(int, part.split("-"))
                cited.extend(range(a, b + 1))
            elif part.strip():
                cited.append(int(part))
    refnums = [int(m.group(1)) for m in re.finditer(r"^\[(\d+)\]", refs, re.M)]
    fig_caps = [(i, p.text.strip()) for i, p in enumerate(doc.paragraphs) if p.text.strip().startswith("Figure ")]
    table_caps = [(i, p.text.strip()) for i, p in enumerate(doc.paragraphs) if p.text.strip().startswith("Table ")]
    placeholders = [p.text.strip() for p in doc.paragraphs if "[" in p.text and ("Author" in p.text or "Funding" in p.text or "repository URL" in p.text)]
    return "\n".join(
        [
            "# Submission precheck audit",
            "",
            f"- DOCX zip test: {zipfile.ZipFile(DOCX).testzip()}",
            f"- Paragraphs: {len(doc.paragraphs)}",
            f"- Tables: {len(doc.tables)}",
            f"- Inline images: {len(doc.inline_shapes)}",
            f"- Figure captions in order: {fig_caps}",
            f"- Table captions in order: {table_caps}",
            f"- Cited references: {sorted(set(cited))}",
            f"- Reference list numbers: {refnums}",
            f"- Unused references: {sorted(set(refnums) - set(cited))}",
            f"- Missing references: {sorted(set(cited) - set(refnums))}",
            f"- Remaining required placeholders: {placeholders}",
        ]
    )


def main() -> None:
    doc = Document(DOCX)

    replace_start(doc, "Keywords:Escherichia", "Keywords: Escherichia coli, genome-scale metabolic model, phylogroup, gapseq, metabolic landscape, CRISPR, secretome, high-throughput genomics")
    replace_start(
        doc,
        "Background:Escherichia",
        "Background: Escherichia coli is characterized by an exceptionally open pan-genome, spanning clinical pathogens, commensal residents, and environmental isolates across eight recognized phylogroups. While individual genome-scale metabolic models (GEMs) have illuminated strain-specific physiology, few large-cohort studies have jointly resolved how genomic defense systems, secretome architecture, and metabolic network topology co-diverge across phylogroups at scale. A principal barrier has been computational: constructing even a single GEM with gap-filling can require >90 minutes per genome, rendering thousand-genome analyses practically infeasible with unoptimized pipelines.",
    )
    replace_contains(doc, "Guiyang，China", "Guiyang, China")

    # Place figure blocks in numerical order: Figure 3 and Figure 4 should precede Figure 5.
    move_figure_block_before(doc, "Figure 3.", "3.3.3 Core-genome ANI distance")
    move_figure_block_before(doc, "Figure 4.", "3.3.3 Core-genome ANI distance")

    cite_and_reorder_refs(doc)

    doc.save(DOCX)
    refreshed = Document(DOCX)
    with TEXT_INDEX.open("w", encoding="utf-8") as handle:
        for i, para in enumerate(refreshed.paragraphs):
            if para.text.strip():
                handle.write(f"[{i}] {para.text.strip()}\n")
    AUDIT.write_text(audit_doc(refreshed), encoding="utf-8")
    print(f"Updated {DOCX}")
    print(f"Audit: {AUDIT}")
    print(f"Text index: {TEXT_INDEX}")


if __name__ == "__main__":
    main()
