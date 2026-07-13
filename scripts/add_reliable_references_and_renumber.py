from __future__ import annotations

from pathlib import Path
import re

from docx import Document


ROOT = Path(r"D:\pipeline\analysis_v5")
DOCX = ROOT / "02.论文" / "Evolutionary Divergence of Genomic Defense.docx"
TEXT_INDEX = ROOT / "02.论文" / "docx_text_index_after_reference_expansion.txt"
AUDIT = ROOT / "03.analysis" / "reference_expansion_audit.md"


ADDED_REFS = {
    26: "Tenaillon O, Skurnik D, Picard B, Denamur E. The population genetics of commensal Escherichia coli. Nat Rev Microbiol. 2010;8(3):207-217.",
    27: "Lukjancenko O, Wassenaar TM, Ussery DW. Comparison of 61 sequenced Escherichia coli genomes. Microb Ecol. 2010;60(4):708-720.",
    28: "Page AJ, Cummins CA, Hunt M, et al. Roary: rapid large-scale prokaryote pan genome analysis. Bioinformatics. 2015;31(22):3691-3693.",
    29: "Jain C, Rodriguez-R LM, Phillippy AM, Konstantinidis KT, Aluru S. High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries. Nat Commun. 2018;9:5114.",
    30: "Ondov BD, Treangen TJ, Melsted P, et al. Mash: fast genome and metagenome distance estimation using MinHash. Genome Biol. 2016;17:132.",
    31: "Kanehisa M, Goto S. KEGG: Kyoto Encyclopedia of Genes and Genomes. Nucleic Acids Res. 2000;28(1):27-30.",
    32: "Anderson MJ. A new method for non-parametric multivariate analysis of variance. Austral Ecology. 2001;26(1):32-46.",
    33: "McArdle BH, Anderson MJ. Fitting multivariate models to community data: a comment on distance-based redundancy analysis. Ecology. 2001;82(1):290-297.",
    34: "Anderson MJ. Distance-based tests for homogeneity of multivariate dispersions. Biometrics. 2006;62(1):245-253.",
    35: "Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc Series B Stat Methodol. 1995;57(1):289-300.",
    36: "Orth JD, Thiele I, Palsson BO. What is flux balance analysis? Nat Biotechnol. 2010;28(3):245-248.",
    37: "Bordbar A, Monk JM, King ZA, Palsson BO. Constraint-based models predict metabolic and associated cellular functions. Nat Rev Genet. 2014;15(2):107-120.",
    38: "Barrangou R, Fremaux C, Deveau H, et al. CRISPR provides acquired resistance against viruses in prokaryotes. Science. 2007;315(5819):1709-1712.",
    39: "Makarova KS, Wolf YI, Alkhnbashi OS, et al. An updated evolutionary classification of CRISPR-Cas systems. Nat Rev Microbiol. 2015;13(11):722-736.",
    40: "Almagro Armenteros JJ, Tsirigos KD, Sonderby CK, et al. SignalP 5.0 improves signal peptide predictions using deep neural networks. Nat Biotechnol. 2019;37(4):420-423.",
}


def set_para_text(para, text: str) -> None:
    for run in para.runs:
        run.text = ""
    if para.runs:
        para.runs[0].text = text
    else:
        para.add_run(text)


def replace_once(doc: Document, old: str, new: str) -> int:
    n = 0
    for para in doc.paragraphs:
        if old in para.text:
            set_para_text(para, para.text.replace(old, new, 1))
            n += 1
    return n


def add_citation_insertions(doc: Document) -> None:
    replacements = [
        (
            "with virulent extraintestinal pathogenic E. coli (ExPEC) predominantly clustering in phylogroups B2 and D, commensal and environmental strains in A and B1, and Shiga toxin-producing strains frequently in E [5, 6].",
            "with virulent extraintestinal pathogenic E. coli (ExPEC) predominantly clustering in phylogroups B2 and D, commensal and environmental strains in A and B1, and Shiga toxin-producing strains frequently in E [5, 6, 26, 27].",
        ),
        (
            "one of the most open pan-genomes among bacteria—with a core genome of ~2,000 genes and a dispensable genome exceeding 18,000 accessory loci [3, 4].",
            "one of the most open pan-genomes among bacteria—with a core genome of ~2,000 genes and a dispensable genome exceeding 18,000 accessory loci [3, 4, 28].",
        ),
        (
            "enabling in silico growth simulation, essential gene identification, and pathway flux analysis [7, 8].",
            "enabling in silico growth simulation, essential gene identification, and pathway flux analysis [7, 8, 36, 37].",
        ),
        (
            "CRISPR-Cas adaptive immune systems act as genomic shields that modulate horizontal gene transfer, shaping which metabolic accessory elements a lineage can stably acquire and retain;",
            "CRISPR-Cas adaptive immune systems act as genomic shields that modulate horizontal gene transfer, shaping which metabolic accessory elements a lineage can stably acquire and retain [38, 39];",
        ),
        (
            "SignalP 6.0 [16] was used with the gram-negative organism model to identify signal peptides.",
            "SignalP 6.0 [16], building on deep-learning signal-peptide prediction frameworks [40], was used with the gram-negative organism model to identify signal peptides.",
        ),
        (
            "KEGG pathway matrices were analyzed using the revised high-value pathway table (`genome_kegg_pathway_matrix.csv`).",
            "KEGG pathway matrices were analyzed using the revised high-value pathway table (`genome_kegg_pathway_matrix.csv`) and the KEGG pathway framework [31].",
        ),
        (
            "Balanced PERMANOVA sensitivity analyses were performed for pathway matrices to verify that global beta-diversity conclusions were not driven by sample-size imbalance.",
            "Balanced PERMANOVA sensitivity analyses were performed for pathway matrices to verify that global beta-diversity conclusions were not driven by sample-size imbalance [32-34].",
        ),
        (
            "COG category comparisons were tested across 26 categories with both Benjamini-Hochberg and Bonferroni correction.",
            "COG category comparisons were tested across 26 categories with both Benjamini-Hochberg false-discovery-rate control [35] and Bonferroni correction.",
        ),
        (
            "Core-genome ANI distance associations were tested by Mantel correlations (999 label permutations).",
            "Core-genome ANI distance associations were tested by Mantel correlations (999 label permutations), using FastANI-derived genomic distances and Mash/ANI concepts for scalable genome comparison [29, 30].",
        ),
        (
            "PCA was used for ordination, PERMANOVA tested global phylogroup association, and PERMDISP assessed multivariate dispersion.",
            "PCA was used for ordination, PERMANOVA tested global phylogroup association [32, 33], and PERMDISP assessed multivariate dispersion [34].",
        ),
        (
            "FBA was run with COBRApy v0.31.1 [20] using the GLPK solver",
            "FBA was run with COBRApy v0.31.1 [20] using standard flux-balance-analysis logic [36] and the GLPK solver",
        ),
        (
            "The FastANI distance matrix contained all 1,000 genomes",
            "The FastANI distance matrix [29] contained all 1,000 genomes",
        ),
    ]
    for old, new in replacements:
        replace_once(doc, old, new)


def collect_references(doc: Document) -> tuple[int, dict[int, str]]:
    refs_start = None
    refs: dict[int, str] = {}
    for i, para in enumerate(doc.paragraphs):
        if para.text.strip() == "References":
            refs_start = i
            continue
        m = re.match(r"^\[(\d+)\]\s*(.+)", para.text.strip())
        if refs_start is not None and m:
            refs[int(m.group(1))] = m.group(2)
    if refs_start is None:
        raise RuntimeError("References heading not found")
    refs.update(ADDED_REFS)
    return refs_start, refs


def citation_values(text: str) -> list[int]:
    vals: list[int] = []
    for part in re.split(r",\s*", text):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a, b = map(int, part.split("-"))
            vals.extend(range(a, b + 1))
        else:
            vals.append(int(part))
    return vals


def format_citation(vals: list[int], mapping: dict[int, int]) -> str:
    return "[" + ", ".join(str(mapping[v]) for v in vals if v in mapping) + "]"


def renumber_by_first_appearance(doc: Document, refs_start: int, refs: dict[int, str]) -> None:
    body_paras = doc.paragraphs[:refs_start]
    first_order: list[int] = []
    for para in body_paras:
        for m in re.finditer(r"\[([0-9,\-\s]+)\]", para.text):
            for val in citation_values(m.group(1)):
                if val not in first_order:
                    first_order.append(val)
    mapping = {old: new for new, old in enumerate(first_order, start=1)}

    for para in body_paras:
        text = para.text
        if not text:
            continue

        def repl(match):
            return format_citation(citation_values(match.group(1)), mapping)

        new_text = re.sub(r"\[([0-9,\-\s]+)\]", repl, text)
        if new_text != text:
            set_para_text(para, new_text)

    # Remove old reference paragraphs and rebuild a clean ordered list after the heading.
    ref_heading = doc.paragraphs[refs_start]
    for para in list(doc.paragraphs[refs_start + 1 :]):
        if re.match(r"^\[\d+\]", para.text.strip()):
            para._element.getparent().remove(para._element)
    cursor = ref_heading._p
    for old_num in first_order:
        if old_num not in refs:
            continue
        new_num = mapping[old_num]
        p = doc.add_paragraph(f"[{new_num}] {refs[old_num]}")
        p._element.getparent().remove(p._element)
        cursor.addnext(p._element)
        cursor = p._element


def audit(doc: Document, refs_start: int) -> str:
    body = "\n".join(p.text for p in doc.paragraphs[:refs_start])
    refs_text = "\n".join(p.text for p in doc.paragraphs[refs_start:])
    cited: list[int] = []
    for m in re.finditer(r"\[([0-9,\-\s]+)\]", body):
        cited.extend(citation_values(m.group(1)))
    refnums = [int(m.group(1)) for m in re.finditer(r"^\[(\d+)\]", refs_text, re.M)]
    return "\n".join(
        [
            "# Reference expansion audit",
            "",
            f"- Cited unique references: {sorted(set(cited))}",
            f"- Reference list numbers: {refnums}",
            f"- Missing references: {sorted(set(cited) - set(refnums))}",
            f"- Unused references: {sorted(set(refnums) - set(cited))}",
            f"- Total references: {len(refnums)}",
        ]
    )


def main() -> None:
    doc = Document(DOCX)
    add_citation_insertions(doc)
    refs_start, refs = collect_references(doc)
    renumber_by_first_appearance(doc, refs_start, refs)
    doc.save(DOCX)
    refreshed = Document(DOCX)
    refs_start, _ = collect_references(refreshed)
    with TEXT_INDEX.open("w", encoding="utf-8") as handle:
        for i, para in enumerate(refreshed.paragraphs):
            if para.text.strip():
                handle.write(f"[{i}] {para.text.strip()}\n")
    AUDIT.write_text(audit(refreshed, refs_start), encoding="utf-8")
    print(f"Updated references in {DOCX}")
    print(f"Audit: {AUDIT}")


if __name__ == "__main__":
    main()
