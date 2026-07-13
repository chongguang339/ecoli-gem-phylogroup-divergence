from __future__ import annotations

import shutil
from pathlib import Path

from docx import Document


ROOT = Path(__file__).resolve().parents[1]
DOCX = ROOT / "02.论文" / "Evolutionary Divergence of Genomic Defense.docx"
BACKUP = ROOT / "02.论文" / "Evolutionary Divergence of Genomic Defense.before_bacdive_mg1655_update.docx"


def replace_para(paragraphs, prefix: str, text: str) -> bool:
    for para in paragraphs:
        if para.text.strip().startswith(prefix):
            para.text = text
            return True
    return False


def append_after(paragraphs, prefix: str, text: str) -> bool:
    for para in paragraphs:
        if para.text.strip().startswith(prefix):
            para.insert_paragraph_before("")
            new_para = para.insert_paragraph_before(text)
            # Move insertion visually before matched para is awkward in python-docx; this is acceptable
            # for local manuscript insertion because the following replacement block controls final order.
            return True
    return False


def main() -> None:
    if not BACKUP.exists():
        shutil.copy2(DOCX, BACKUP)

    doc = Document(DOCX)
    ps = doc.paragraphs

    replace_para(
        ps,
        "For the SD01 770-condition BiGG-model benchmark",
        "For the SD01 770-condition BiGG-model benchmark, five curated models were evaluated: iJO1366, iML1515, ic_1306, iUTI89_1310, and the 1327-gene proxy iECSF_1327 for the requested iECs1327 comparison. Public BiGG SBML files were used for these curated models [29, 30]. Paired gapseq reconstructions were evaluated where local SBML files could be parsed by COBRApy. This SD01 analysis was treated as model-consistency benchmarking against the Monk et al. computational growth matrix, not as independent experimental phenotype validation.",
    )

    replace_para(
        ps,
        "The SD01 carbon-source benchmark separates three evidence levels.",
        "Carbon-source and phenotype evaluation was separated into two evidence classes. The SD01 770-condition benchmark tested whether curated BiGG models and paired gapseq reconstructions reproduced the Monk et al. computational growth matrix and was therefore interpreted as model-consistency benchmarking. In addition, an independent MG1655 phenotype check was constructed from the BacDive 4902/API 20E utilization table (DOI: 10.13145/bacdive4902.20260601.11), using the locally curated workbook K12_MG1655_BacDive_carbon_source_data.xlsx. Because the BacDive/API 20E panel combines fermentation, assimilation, degradation, and hydrolysis assays, the primary FBA summary used a stricter 10-entry carbon-assimilation subset, while the full 15-entry BacDive panel was retained in Supplementary Data S8.",
    )
    replace_para(
        ps,
        "Carbon-source and phenotype evaluation was separated into two evidence classes.",
        "Carbon-source and phenotype evaluation was separated into two evidence classes. The SD01 770-condition benchmark tested whether curated BiGG models and paired gapseq reconstructions reproduced the Monk et al. computational growth matrix and was therefore interpreted as model-consistency benchmarking. In addition, an independent MG1655 phenotype check was constructed from the BacDive 4902/API 20E utilization table (DOI: 10.13145/bacdive4902.20260601.11), using the locally curated workbook K12_MG1655_BacDive_carbon_source_data.xlsx. Because the BacDive/API 20E panel combines fermentation, assimilation, degradation, and hydrolysis assays, the primary FBA summary used a stricter 10-entry carbon-assimilation subset, while the full 15-entry BacDive panel was retained in Supplementary Data S8.",
    )

    replace_para(
        ps,
        "3.1.3 SD01 carbon-source consistency benchmarking",
        "3.1.3 SD01 model-consistency and BacDive MG1655 phenotype benchmarking",
    )

    replace_para(
        ps,
        "For MG1655, the M9-constrained optimized model achieved accuracy",
        "For the independent BacDive MG1655 check, iML1515 and the optimized MG1655 gapseq reconstruction were simulated under aerobic single-substrate FBA conditions after exchange-ID harmonization. iML1515 used a BiGG-compatible M9-like medium, whereas the optimized gapseq model used its native medium after removal of default carbon-like exchanges before adding each BacDive substrate. The default MG1655 gapseq SBML file did not pass COBRApy validation and was therefore reported in the model-load table rather than forced into the phenotype benchmark. Among resolved exchanges in the carbon-assimilation subset, both iML1515 and the optimized reconstruction correctly classified 7/9 substrates (accuracy = 77.8%). iML1515 achieved sensitivity = 100.0%, specificity = 33.3%, balanced accuracy = 0.667, and MCC = 0.500; the optimized reconstruction achieved the same balanced accuracy and MCC, with false-positive growth primarily involving negative BacDive substrates such as citrate and/or sucrose depending on model namespace. Across the full resolved BacDive panel, balanced accuracy was 0.604 for iML1515 (n = 14) and 0.638 for the optimized reconstruction (n = 13). These results provide a small independent phenotype sanity check, but their limited sample size and mixed API 20E assay types preclude strong claims of broad phenotype accuracy.",
    )

    replace_para(
        ps,
        "For the available paired gapseq reconstructions, the SD01 consistency benchmark",
        "For the available paired gapseq reconstructions, the SD01 consistency benchmark showed modest improvement on average but not uniform superiority of the optimized workflow. Mean balanced accuracy increased from 0.562 for default gapseq reconstructions to 0.592 for optimized reconstructions, and mean MCC increased from 0.15 to 0.25. These SD01 results should be read together with the BacDive MG1655 phenotype check (Supplementary Figure S2 and Supplementary Data S8): SD01 evaluates consistency with a published computational matrix, whereas BacDive provides a small independent experimental utilization panel.",
    )

    replace_para(
        ps,
        "Figure 1. Stratified runtime and reference-model and SD01 benchmarking.",
        "Figure 1. Stratified runtime and reference-model benchmarking. (A) Speedup distribution across 20 assemblies spanning Actual_N50 = 22,535-675,909 bp and CDS_Count = 4,252-5,584; the subset includes A (n = 4), B1 (n = 7), B2 (n = 4), C (n = 1), E (n = 2), cryptic/U (n = 2). (B) Default versus optimized wall-clock runtime distributions. (C) Reaction-count compaction and reaction retention in the 20-genome stress-test subset; paired XML reaction matrices further identify optimized-only and default-only reaction differences. (D) Reference-model comparison showing the fraction of removed reactions that overlapped matched BiGG/reference models.",
    )

    replace_para(
        ps,
        "Of the 1,000 attempted GEM reconstructions",
        "Of the 1,000 attempted GEM reconstructions, all 1,000 models achieved viable growth (positive biomass flux under rich medium conditions) following M9-constrained gap-filling. This outcome is an expected consequence of gap-filling by design-which adds reactions until biomass production is feasible-and is therefore not interpreted as an independent quality metric; the relevant measures of model fidelity are the SD01 model-consistency benchmark and the BacDive MG1655 phenotype check described in Section 3.1.3.",
    )

    replace_para(
        ps,
        "Several caveats should guide interpretation.",
        "Several caveats should guide interpretation. First, pathway names are represented primarily by KEGG identifiers; pathway-name annotation, representative KO/gene members, Cramer's V, and effect-size tables are needed for detailed biological interpretation of the 15 adjusted signals. Second, phylogroup-associated PERMANOVA results are interpreted together with significant PERMDISP and ANI residual models, not as pure categorical phylogroup-explained variance. Third, the BacDive MG1655 phenotype check is intentionally small and combines several API 20E assay types; it supports a limited independent sanity check but does not replace larger experimental Biolog or phenotype-microarray validation.",
    )

    replace_para(
        ps,
        "Availability of data and materials:",
        "Availability of data and materials: All genome accession numbers, phylogroup assignments, annotation outputs, GEM topology metrics, available GEM QC fields, benchmark metadata, reaction-set matrices, SD01 carbon-source exchange mappings, SD01 FBA prediction/model-consistency matrices, and BacDive MG1655 phenotype-check FBA outputs are organized in the project supplementary data directory. The final repository URL should be supplied at submission.",
    )

    replace_para(
        ps,
        "Supplementary Data files include:",
        "Supplementary Data files include: Data S1, selected 1,000-genome accession and phylogroup table; Data S2, genome assembly, COG, KEGG, GEM topology, CRISPR, SignalP, and TMHMM matrices; Data S3, statistical-test outputs with effect sizes and multiple-testing correction; Data S4, benchmark runtime and reaction-set comparison matrices; Data S5, SD01 carbon-source exchange-ID mapping tables; Data S6, SD01 770-condition FBA prediction matrices, growth-rate matrices, model inclusion and simulation status tables; Data S7, CRISPR spacer target-matching outputs; Data S8, BacDive MG1655 phenotype truth table, model-load status, FBA predictions, and accuracy metrics.",
    )

    doc.save(DOCX)
    print(f"Updated {DOCX}")
    print(f"Backup {BACKUP}")


if __name__ == "__main__":
    main()
