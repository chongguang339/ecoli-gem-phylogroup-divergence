from pathlib import Path

from docx import Document


DOCX = Path(r"D:\pipeline\analysis_v5\02.论文\Evolutionary Divergence of Genomic Defense.first_version_clean.docx")

REPLACEMENTS = {
    "Figure 4. COG functional-category landscape. Heatmap shows mean relative abundance of the top BH-corrected COG categories by phylogroup; 22/26 COG categories remained significant after both BH and Bonferroni correction.": "Figure 4. COG functional-category landscape. Heatmap shows mean relative abundance of the top BH-significant COG categories by phylogroup; 22/26 COG categories remained significant after both BH and Bonferroni correction.",
    "Future directions include: (i) extension of this framework to other Enterobacteriaceae and microbial ecosystem modeling using the same Enterobacterales fast-pass database [39]; (ii) integration of transcriptome-based flux context using resource allocation constraints; (iii) application of the metabolic distance landscape for antimicrobial susceptibility prediction and microbial pathogenesis studies, given the known links between bacterial metabolism and pathogenicity [40]; and (iv) development of a continuously updated web resource for browsing the 1,000-GEM metabolic atlas.": "Future directions include: (i) extension of this framework to other Enterobacteriaceae and microbial ecosystem modeling using the same Enterobacterales fast-pass database [39]; (ii) integration of transcriptome-based flux context using resource allocation constraints; (iii) application of the metabolic distance landscape for antimicrobial susceptibility prediction and microbial pathogenesis studies, given the known links between bacterial metabolism and pathogenicity [40]; and (iv) development of a continuously maintained web resource for browsing the 1,000-GEM metabolic atlas.",
}


def main():
    doc = Document(DOCX)
    count = 0
    for para in doc.paragraphs:
        old = para.text
        if old in REPLACEMENTS:
            para.clear()
            para.add_run(REPLACEMENTS[old])
            count += 1
    doc.save(DOCX)
    print(f"updated {count} paragraphs")


if __name__ == "__main__":
    main()
