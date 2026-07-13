from pathlib import Path
from shutil import copy2

from docx import Document
from docx.shared import Inches


ROOT = Path(r"D:\pipeline\analysis_v5")
DOCX = ROOT / "02.论文" / "Evolutionary Divergence of Genomic Defense.docx"
BACKUP = ROOT / "02.论文" / "Evolutionary Divergence of Genomic Defense.before_tmhmm_update.docx"
FIG2 = ROOT / "01.figures" / "Figure2_revised_functional_landscape.png"
INDEX = ROOT / "02.论文" / "docx_text_index_after_tmhmm_update.txt"


REPLACEMENTS = {
    56: "Transmembrane topology prediction. TMHMM v2.0 [17] was applied to all predicted protein sequences to enumerate proteins containing at least one transmembrane helix. After updating the 154 previously missing genomes, the revised TMHMM-positive protein table contains 1,125,213 protein-level records across 999 of 1,000 genomes. Per-genome analyses used TMHMM-positive protein count, total predicted helices, mean helices per membrane protein, maximum helices, and TMHMM proteins per 1,000 CDS. The single remaining genome without TMHMM detail output was retained as missing (NA), not as a biological zero.",
    92: "All statistical analyses were performed in Python using reproducible outputs from the corrected cohort. The revised cohort includes A (n = 298), B2 (n = 244), B1 (n = 200), D (n = 103), E (n = 73), C (n = 36), F (n = 31), G (n = 11), and U (n = 4). Genome-level continuous metrics were assessed by permutation ANOVA (999 permutations), and effect sizes were reported as eta-squared and omega-squared to distinguish statistical significance from biological magnitude. COG category comparisons were tested across 26 categories with both Benjamini-Hochberg and Bonferroni correction. KEGG pathway frequency tests were corrected across 264 non-redundant pathways. Balanced PERMANOVA sensitivity analyses were performed for pathway matrices to verify that global beta-diversity conclusions were not driven by sample-size imbalance.",
    133: "TMHMM-positive proteins were identified in 999 of 1,000 genomes after reprocessing the 154 previously missing genomes, yielding 1,125,213 protein-level records. Only one B1 genome (GCF_045037575.1_ASM4503757v1_genomic.fna) remained without TMHMM detail output and was treated as missing rather than as a true zero. TMHMM burden differed significantly across phylogroups (permutation ANOVA p = 0.001), with substantial effect sizes for absolute TMHMM-positive protein count (eta-squared = 0.415; omega-squared = 0.410) and total predicted helices (eta-squared = 0.525; omega-squared = 0.521), and a moderate effect after normalization by CDS (TMHMM proteins per 1,000 CDS: eta-squared = 0.289; omega-squared = 0.283). Median TMHMM-positive protein counts were highest in G (1,171), E (1,167), F (1,162), C (1,161), and B2 (1,153), whereas A and U were lower. After CDS normalization, the strongest medians were observed in C (238.7), A (237.6), B1 (236.4), B2 (236.3), and G (235.9) TMHMM-positive proteins per 1,000 CDS, indicating phylogroup-associated membrane-protein architecture without the earlier data-quality artifact caused by missing TMHMM files.",
    139: "Figure 2. Revised phylogroup-stratified genome, GEM, secretion, membrane, and CRISPR spacer features. Boxplots show group distributions across genome size, GEM reactions/metabolites, SignalP-positive proteins per 1,000 CDS, TMHMM-positive proteins per 1,000 CDS (n = 999; one unresolved B1 genome treated as missing), and unique CRISPR spacer sequences.",
    188: "Several limitations remain. First, pathway names are represented by KEGG identifiers only; pathway-name annotation would improve biological interpretability of the 15 corrected signals. Second, the current analysis does not include a core-genome ANI or phylogenetic distance matrix, so residual phylogroup effects after continuous phylogenetic relatedness cannot yet be quantified. Third, although TMHMM coverage was improved from 846 to 999 genomes after reprocessing, one B1 genome still lacks detailed TMHMM output and was treated as missing in membrane-protein analyses. Fourth, CRISPR spacer sharing is sequence-level but not yet matched to phage/protospacer databases, so spacer targets and ecological exposures remain unresolved. Finally, all GEM predictions remain in silico and should be experimentally validated in representative strains.",
}


def main():
    copy2(DOCX, BACKUP)
    doc = Document(DOCX)
    for idx, text in REPLACEMENTS.items():
        doc.paragraphs[idx].text = text

    if len(doc.inline_shapes) < 2:
        raise RuntimeError("Expected at least two inline images to replace Figure 2")
    shape = doc.inline_shapes[1]
    blip = shape._inline.graphic.graphicData.pic.blipFill.blip
    doc.part.related_parts[blip.embed]._blob = FIG2.read_bytes()
    shape.width = Inches(6.7)
    shape.height = Inches(4.47)

    doc.save(DOCX)
    reloaded = Document(DOCX)
    with INDEX.open("w", encoding="utf-8") as f:
        for i, paragraph in enumerate(reloaded.paragraphs):
            f.write(f"{i:04d}: {paragraph.text}\n")
    print(f"saved={DOCX}")
    print(f"backup={BACKUP}")
    print(f"index={INDEX}")
    print(f"inline_shapes={len(reloaded.inline_shapes)}")


if __name__ == "__main__":
    main()
