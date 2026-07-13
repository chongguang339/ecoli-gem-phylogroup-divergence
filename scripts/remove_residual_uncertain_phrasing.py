from __future__ import annotations

import shutil
from pathlib import Path

from docx import Document


ROOT = Path(__file__).resolve().parents[1]
DOCX = ROOT / "02.论文" / "Evolutionary Divergence of Genomic Defense.docx"
BACKUP = ROOT / "02.论文" / "Evolutionary Divergence of Genomic Defense.before_uncertain_phrase_cleanup.docx"


def replace_para(paragraphs, prefix: str, text: str) -> None:
    for para in paragraphs:
        if para.text.strip().startswith(prefix):
            para.text = text
            return


def main() -> None:
    if not BACKUP.exists():
        shutil.copy2(DOCX, BACKUP)
    doc = Document(DOCX)
    ps = doc.paragraphs

    replace_para(
        ps,
        "The reference-model structural benchmark comprised five model systems",
        "The reference-model structural benchmark comprised five model systems (Figure 1D). Four systems (iB21_1397, iEC1364_W, iECDH10B_1368, and iJN1463) had independently logged paired default/optimized runtime records and showed substantial acceleration, with a median speedup of 12.76-fold and a range of 11.30-15.52-fold. MG1655/iML1515 was retained as the fifth structural reference system and used the benchmark-average runtime basis for five-model summaries. Across the four independently timed reference systems, optimized reaction retention ranged from 77.0% to 87.5%. When reactions present in the default reconstruction but absent after optimization were compared against matched BiGG reference models, only 10.55-17.53% of removed reactions were present in the curated references. Thus, most removed reactions lacked matched reference-model support, while the reference-overlapping subset is reported explicitly in the reaction matrices.",
    )

    replace_para(
        ps,
        "Propionate, methylcitrate, and mucus-sugar utilization pathways remain",
        "Propionate, methylcitrate, and mucus-sugar utilization pathways are interpreted only as pathway-specific candidates when they appear among the corrected KEGG signals or targeted annotations.",
    )

    replace_para(
        ps,
        "Conversely, putative central-metabolic retention",
        "Conversely, putative central-metabolic retention in A and B1 genomes is interpreted as specific adjusted pathway-frequency contrasts, not as evidence that A/B1 and D/B2 occupy globally separated pathway landscapes.",
    )

    replace_para(
        ps,
        "Availability of data and materials:",
        "Availability of data and materials: All genome accession numbers, phylogroup assignments, annotation outputs, GEM topology metrics, available GEM QC fields, benchmark metadata, reaction-set matrices, SD01 carbon-source exchange mappings, SD01 FBA prediction/model-consistency matrices, and BacDive MG1655 phenotype-check FBA outputs are organized in the project supplementary data directory. Repository URL: pending insertion before submission.",
    )

    doc.save(DOCX)
    print(f"Updated {DOCX}")
    print(f"Backup {BACKUP}")


if __name__ == "__main__":
    main()
