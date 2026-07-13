from __future__ import annotations

import shutil
from pathlib import Path

from docx import Document


ROOT = Path(__file__).resolve().parents[1]
DOCX = ROOT / "02.论文" / "Evolutionary Divergence of Genomic Defense.docx"
BACKUP = ROOT / "02.论文" / "Evolutionary Divergence of Genomic Defense.before_final_tightening.docx"


def replace_para(paragraphs, prefix: str, text: str) -> bool:
    for para in paragraphs:
        if para.text.strip().startswith(prefix):
            para.text = text
            return True
    return False


def add_table1_bacdive_row(doc: Document) -> None:
    table = doc.tables[0]
    first_col = [row.cells[0].text.strip() for row in table.rows]
    if "BacDive MG1655 phenotype check" in first_col:
        return
    row = table.add_row()
    row.cells[0].text = "BacDive MG1655 phenotype check"
    row.cells[1].text = "15 API 20E/BacDive entries; 10-entry carbon-assimilation subset"
    row.cells[2].text = (
        "Carbon subset: default gapseq BA = 0.500; optimized gapseq BA = 0.667; "
        "iML1515 BA = 0.667"
    )
    row.cells[3].text = "Data S8; Figure S2"


def main() -> None:
    if not BACKUP.exists():
        shutil.copy2(DOCX, BACKUP)

    doc = Document(DOCX)
    ps = doc.paragraphs

    add_table1_bacdive_row(doc)

    replace_para(
        ps,
        "Structural reference-model benchmarking was organized across five model systems:",
        "Structural reference-model benchmarking was organized across five model systems: iB21_1397, iEC1364_W, iECDH10B_1368, iJN1463, and MG1655/iML1515. Curated SBML/XML files were parsed and compared with paired default and optimized gapseq reconstructions where corresponding local reconstructions were available. MG1655/iML1515 was retained as the fifth reference model system; for runtime summaries it used the benchmark-average runtime basis rather than an independently logged default-runtime measurement, and measured 20-genome runtime distributions are reported separately. The structural analysis focused on default-only, optimized-only, shared-reaction, and curated-reference-overlap metrics to distinguish reaction compaction from loss of reference-supported structure. SD01 Table 5 was analyzed as a model-consistency benchmark, and BacDive MG1655/API 20E data were analyzed as a separate small experimental phenotype panel.",
    )

    replace_para(
        ps,
        "For the independent BacDive MG1655 check,",
        "For the independent BacDive MG1655 check, iML1515 and paired MG1655 gapseq reconstructions were simulated under aerobic single-substrate FBA conditions after exchange-ID harmonization. iML1515 used a BiGG-compatible M9-like medium, whereas gapseq models used their native media after removal of default carbon-like exchanges before adding each BacDive substrate. Among resolved exchanges in the carbon-assimilation subset, the default reconstruction correctly classified 6/9 substrates (accuracy = 66.7%, sensitivity = 100.0%, specificity = 0.0%, balanced accuracy = 0.500), whereas both iML1515 and the optimized reconstruction correctly classified 7/9 substrates (accuracy = 77.8%, sensitivity = 100.0%, specificity = 33.3%, balanced accuracy = 0.667, MCC = 0.500). Across the full resolved BacDive panel, balanced accuracy was 0.500 for the default reconstruction (n = 13), 0.638 for the optimized reconstruction (n = 13), and 0.604 for iML1515 (n = 14). These results show that optimization reduced false-positive substrate predictions in this small MG1655 panel; the panel size and mixed API 20E assay types limit inference beyond these tested substrates.",
    )

    replace_para(
        ps,
        "For the available paired gapseq reconstructions, the SD01 consistency benchmark",
        "For available paired gapseq reconstructions, the SD01 consistency benchmark showed a modest average shift after optimization. Mean balanced accuracy increased from 0.562 for default reconstructions to 0.592 for optimized reconstructions, and mean MCC increased from 0.15 to 0.25. SD01 evaluates consistency with a published computational matrix, whereas the BacDive MG1655 panel provides a small independent experimental utilization check (Supplementary Figure S2 and Supplementary Data S8).",
    )

    replace_para(
        ps,
        "The computational benchmark anchors the manuscript's main contribution:",
        "The computational benchmark anchors the manuscript's main contribution: optimized gapseq performance was evaluated as a runtime distribution across 20 N50-stratified assemblies rather than as a single example run. The optimized workflow reduced reconstruction time and produced compact reaction sets, while reference-model comparisons and carbon-source benchmarks separated computational scalability from substrate-level phenotype agreement.",
    )

    replace_para(
        ps,
        "At an interpretive level, these data support",
        "Together, these data support a modular model of E. coli functional diversification. Immune memory, secretion signals, membrane topology, COG composition, and pathway frequency each vary by phylogroup, but they do not reduce to a single host-pathogen versus commensal axis. The analysis separates global multivariate structure from corrected feature-level contrasts and reports effect-size context alongside p values.",
    )

    replace_para(
        ps,
        "Several caveats should guide interpretation.",
        "The main interpretation is bounded by three constraints. First, KEGG identifiers define pathway-level signals and require pathway-specific annotation before mechanistic claims. Second, PERMANOVA results are interpreted together with PERMDISP and ANI residual models, so PERMANOVA R2 is not treated as pure categorical phylogroup-explained variance. Third, the BacDive MG1655 panel is a small API 20E-derived phenotype check and does not replace larger Biolog or phenotype-microarray validation.",
    )

    replace_para(
        ps,
        "Future directions include:",
        "Future work will require larger independent phenotype panels, model-quality audits with blocked-reaction and mass-balance metrics, and extension of the workflow to additional Enterobacteriaceae cohorts.",
    )

    replace_para(
        ps,
        "In conclusion, this study presents",
        "In conclusion, this study presents an optimized E. coli GEM reconstruction framework applied to a 1,000-genome phylogroup-stratified cohort. The benchmark demonstrates distribution-level runtime acceleration and compact reaction sets; SD01 and BacDive analyses separate model-consistency testing from a small independent phenotype check. FastANI residual analysis shows that KEGG pathway composition is strongly tied to continuous core-genome distance, while discrete phylogroup labels retain smaller residual effects. Auxiliary COG, CRISPR, SignalP, and TMHMM analyses provide context for lineage-associated functional variation without replacing substrate-level validation.",
    )

    doc.save(DOCX)
    print(f"Updated {DOCX}")
    print(f"Backup {BACKUP}")


if __name__ == "__main__":
    main()
