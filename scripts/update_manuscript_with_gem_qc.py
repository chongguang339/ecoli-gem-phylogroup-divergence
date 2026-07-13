from __future__ import annotations

import shutil
from pathlib import Path

from docx import Document


ROOT = Path(__file__).resolve().parents[1]
DOCX = ROOT / "02.论文" / "Evolutionary Divergence of Genomic Defense.docx"
BACKUP = ROOT / "02.论文" / "Evolutionary Divergence of Genomic Defense.before_gem_qc_update.docx"


def replace_para(paragraphs, prefix: str, text: str) -> bool:
    for para in paragraphs:
        if para.text.strip().startswith(prefix):
            para.text = text
            return True
    return False


def add_after(paragraphs, prefix: str, text: str) -> bool:
    for i, para in enumerate(paragraphs):
        if para.text.strip().startswith(prefix):
            new_para = para.insert_paragraph_before(text)
            # Move new paragraph after matched paragraph by swapping text with next paragraph pattern.
            old = para.text
            para.text = text
            new_para.text = old
            return True
    return False


def add_table1_row(doc: Document) -> None:
    table = doc.tables[0]
    first_col = [row.cells[0].text.strip() for row in table.rows]
    if "COBRApy/MEMOTE GEM quality control" in first_col:
        return
    row = table.add_row()
    row.cells[0].text = "COBRApy/MEMOTE GEM quality control"
    row.cells[1].text = "1,000 reconstructed GEMs; 2 representative MEMOTE HTML reports"
    row.cells[2].text = (
        "Mass-balanced reactions = 92.04% +/- 0.60%; blocked reactions = 41.39% +/- 0.79%; "
        "dead-end metabolites = 35.07% +/- 2.85%"
    )
    row.cells[3].text = "Data S9; Figure S3; MEMOTE HTML reports"


def main() -> None:
    if not BACKUP.exists():
        shutil.copy2(DOCX, BACKUP)

    doc = Document(DOCX)
    ps = doc.paragraphs
    add_table1_row(doc)

    replace_para(
        ps,
        "For each of the 1,000 phylogroup-stratified genomes, both the optimized",
        "For each of the 1,000 phylogroup-stratified genomes, optimized M9-constrained models were assessed for topology, solvability, and COBRApy-derived quality metrics. The primary topology fields were number of genes, reactions, and metabolites; biomass flux under rich medium was used as a solvability check. In addition, a cohort-level QC scan quantified total reactions, total metabolites, mass-balanced reaction ratio, blocked-reaction ratio, and dead-end metabolite ratio for all reconstructed GEMs. Two representative models were further inspected with MEMOTE v0.13.0 and the resulting interactive HTML reports were retained as supplementary files. MEMOTE annotation failures were interpreted as database-namespace compliance findings, particularly for gapseq/ModelSEED identifiers, rather than as direct evidence of metabolic inconsistency.",
    )

    add_after(
        ps,
        "Of the 1,000 attempted GEM reconstructions",
        "COBRApy quality control across all 1,000 reconstructed GEMs showed consistent cohort-level stoichiometric and topological properties (Supplementary Figure S3 and Supplementary Data S9). The models contained 2,285-3,100 reactions (median = 2,417) and 1,777-2,476 metabolites (median = 1,888). The mean mass-balanced reaction ratio was 92.04% +/- 0.60%, indicating high stoichiometric consistency across the automated reconstructions. Mean blocked-reaction ratio was 41.39% +/- 0.79%, and mean dead-end metabolite ratio was 35.07% +/- 2.85%, values consistent with automatically generated draft-scale reconstructions that have not undergone manual MEMOTE-guided curation.",
    )

    replace_para(
        ps,
        "Across the 1,000 viable models in the cohort-level topology tables",
        "Across the 1,000 viable models in the cohort-level topology tables, topology statistics were: median GEM_Genes = 8,007, median GEM_Reactions = 2,403, and median GEM_Metabolites = 12,354. These fields derive from the gapseq output summary table and are therefore interpreted as cohort-level reconstruction descriptors rather than one-to-one equivalents of manually curated iML1515 identifier counts. The independent COBRApy QC scan provides the SBML-parsed reaction/metabolite counts and mass-balance, blocked-reaction, and dead-end-metabolite ratios reported above.",
    )

    replace_para(
        ps,
        "Biomass flux values among viable models range",
        "Biomass flux values among viable models range from 0.526 to 0.584 mmol gDW-1 h-1 under standardized rich medium conditions, with median flux = 0.5826 mmol gDW-1 h-1. This confirms that the optimized reconstructions produced solvable biomass equations under the standardized rich-medium test condition. Solvability and mass-balance QC support model usability, while substrate-specific behavior is evaluated separately through SD01 model-consistency and BacDive MG1655 phenotype checks.",
    )

    replace_para(
        ps,
        "The main interpretation is bounded by three constraints.",
        "The main interpretation is bounded by four constraints. First, KEGG identifiers define pathway-level signals and require pathway-specific annotation before mechanistic claims. Second, PERMANOVA results are interpreted together with PERMDISP and ANI residual models, so PERMANOVA R2 is not treated as pure categorical phylogroup-explained variance. Third, the BacDive MG1655 panel is a small API 20E-derived phenotype check and does not replace larger Biolog or phenotype-microarray validation. Fourth, COBRApy/MEMOTE QC supports stoichiometric and topological consistency but does not replace manual model curation; blocked reactions and dead-end metabolites remain expected features of automated large-cohort reconstructions.",
    )

    replace_para(
        ps,
        "Future work will require larger independent phenotype panels",
        "Future work will require larger independent phenotype panels, targeted manual curation of representative GEMs with high blocked-reaction or dead-end-metabolite burden, and extension of the workflow to additional Enterobacteriaceae cohorts.",
    )

    replace_para(
        ps,
        "Availability of data and materials:",
        "Availability of data and materials: All genome accession numbers, phylogroup assignments, annotation outputs, GEM topology metrics, COBRApy/MEMOTE QC outputs, benchmark metadata, reaction-set matrices, SD01 carbon-source exchange mappings, SD01 FBA prediction/model-consistency matrices, and BacDive MG1655 phenotype-check FBA outputs are organized in the project supplementary data directory. Repository URL: pending insertion before submission.",
    )

    replace_para(
        ps,
        "Supplementary Data files include:",
        "Supplementary Data files include: Data S1, selected 1,000-genome accession and phylogroup table; Data S2, genome assembly, COG, KEGG, GEM topology, CRISPR, SignalP, and TMHMM matrices; Data S3, statistical-test outputs with effect sizes and multiple-testing correction; Data S4, benchmark runtime and reaction-set comparison matrices; Data S5, SD01 carbon-source exchange-ID mapping tables; Data S6, SD01 770-condition FBA prediction matrices, growth-rate matrices, model inclusion and simulation status tables; Data S7, CRISPR spacer target-matching outputs; Data S8, BacDive MG1655 phenotype truth table, model-load status, FBA predictions, and accuracy metrics; Data S9, COBRApy cohort-level GEM QC tables and representative MEMOTE HTML reports.",
    )

    doc.save(DOCX)
    print(f"Updated {DOCX}")
    print(f"Backup {BACKUP}")


if __name__ == "__main__":
    main()
