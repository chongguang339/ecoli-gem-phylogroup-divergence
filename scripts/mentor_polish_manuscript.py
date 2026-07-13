from __future__ import annotations

from pathlib import Path
import re
import zipfile

from docx import Document


ROOT = Path(r"D:\pipeline\analysis_v5")
DOCX = ROOT / "02.论文" / "Evolutionary Divergence of Genomic Defense.docx"
TEXT_INDEX = ROOT / "02.论文" / "docx_text_index_after_mentor_polish.txt"
AUDIT = ROOT / "03.analysis" / "mentor_polish_audit.md"


REPLACEMENTS = {
    "Evolutionary Divergence of Genomic Defense, Secretome, and Metabolic Networks in Escherichia coli: A High-Throughput GEM Cohort Analysis":
        "Phylogroup-Resolved Divergence of Defense, Secretome, and Metabolic Network Features in Escherichia coli: A High-Throughput GEM Cohort Analysis",
    "These results provide a conservative and reproducible framework for SCI-level interpretation of E. coli functional divergence.":
        "These results provide a conservative and reproducible framework for interpreting E. coli functional divergence in large genome cohorts.",
    "Beyond metabolic modeling, a complete picture of E. coli ecological diversification requires integrating additional functional dimensions.":
        "Beyond metabolic modeling, a broader view of E. coli ecological diversification requires integrating additional functional dimensions.",
    "CRISPR-Cas adaptive immune systems act as genomic shields that modulate horizontal gene transfer":
        "CRISPR-Cas adaptive immune systems modulate horizontal gene transfer",
    "Secreted proteins—identified via signal peptide prediction—represent the frontline molecular toolkit for host interaction, nutrient acquisition, and immune evasion;":
        "Secreted proteins—identified via signal peptide prediction—contribute to host interaction, nutrient acquisition, and immune evasion;",
    "No prior large-scale study has jointly quantified these four molecular dimensions—metabolic network topology, CRISPR immune architecture, secretome composition, and transmembrane repertoire—across a phylogroup-resolved thousand-genome E. coli cohort, leaving the co-evolutionary dependencies among these systems uncharacterized.":
        "Few large-scale studies have jointly quantified these four molecular dimensions—metabolic network topology, CRISPR immune architecture, secretome composition, and transmembrane repertoire—across a phylogroup-resolved thousand-genome E. coli cohort.",
    "3.4 Large-Scale Optimized GEMs Reveal Pan-Genome-Driven Metabolic Plasticity with Phylogroup-Characteristic Topological Fingerprints":
        "3.4 Cohort-Level GEM Topology Shows Compact but Phylogroup-Associated Network Variation",
    "At a mechanistic level, these data support a modular rather than monolithic model of E. coli diversification.":
        "At an interpretive level, these data support a modular rather than monolithic model of E. coli diversification.",
    "In conclusion, this revised analysis presents an optimized and more carefully validated E. coli GEM reconstruction framework applied to a corrected 1,000-genome phylogroup cohort.":
        "In conclusion, this revised analysis presents an optimized and systematically benchmarked E. coli GEM reconstruction framework applied to a corrected 1,000-genome phylogroup cohort.",
}


def set_para_text(para, text: str) -> None:
    for run in para.runs:
        run.text = ""
    if para.runs:
        para.runs[0].text = text
    else:
        para.add_run(text)


def main() -> None:
    doc = Document(DOCX)
    edits = []
    for para in doc.paragraphs:
        text = para.text
        new_text = text
        for old, new in REPLACEMENTS.items():
            if old in new_text:
                new_text = new_text.replace(old, new)
        if new_text != text:
            set_para_text(para, new_text)
            edits.append({"old": text, "new": new_text})

    doc.save(DOCX)
    refreshed = Document(DOCX)
    with TEXT_INDEX.open("w", encoding="utf-8") as handle:
        for i, para in enumerate(refreshed.paragraphs):
            if para.text.strip():
                handle.write(f"[{i}] {para.text.strip()}\n")

    refs_start = next(i for i, p in enumerate(refreshed.paragraphs) if p.text.strip() == "References")
    body = "\n".join(p.text for p in refreshed.paragraphs[:refs_start])
    refs = "\n".join(p.text for p in refreshed.paragraphs[refs_start:])
    cited = []
    for m in re.finditer(r"\[([0-9,\-\s]+)\]", body):
        for part in re.split(r",\s*", m.group(1)):
            if "-" in part:
                a, b = map(int, part.split("-"))
                cited.extend(range(a, b + 1))
            elif part.strip():
                cited.append(int(part))
    refnums = [int(m.group(1)) for m in re.finditer(r"^\[(\d+)\]", refs, re.M)]
    placeholders = [p.text.strip() for p in refreshed.paragraphs if "[Author" in p.text or "[Funding" in p.text or "[repository URL" in p.text]
    audit = "\n".join(
        [
            "# Mentor polish audit",
            "",
            f"- DOCX zip test: {zipfile.ZipFile(DOCX).testzip()}",
            f"- Paragraphs: {len(refreshed.paragraphs)}",
            f"- Tables: {len(refreshed.tables)}",
            f"- Inline images: {len(refreshed.inline_shapes)}",
            f"- Edits applied: {len(edits)}",
            f"- References: {len(refnums)} total; missing={sorted(set(cited)-set(refnums))}; unused={sorted(set(refnums)-set(cited))}",
            f"- Remaining required placeholders: {placeholders}",
        ]
    )
    AUDIT.write_text(audit, encoding="utf-8")
    print(audit)


if __name__ == "__main__":
    main()
