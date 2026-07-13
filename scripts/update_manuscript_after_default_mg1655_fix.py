from __future__ import annotations

import shutil
from pathlib import Path

from docx import Document


ROOT = Path(__file__).resolve().parents[1]
DOCX = ROOT / "02.论文" / "Evolutionary Divergence of Genomic Defense.docx"
BACKUP = ROOT / "02.论文" / "Evolutionary Divergence of Genomic Defense.before_default_mg1655_path_fix.docx"


def replace_para(paragraphs, prefix: str, text: str) -> bool:
    for para in paragraphs:
        if para.text.strip().startswith(prefix):
            para.text = text
            return True
    return False


def main() -> None:
    if not BACKUP.exists():
        shutil.copy2(DOCX, BACKUP)
    doc = Document(DOCX)
    ps = doc.paragraphs

    replace_para(
        ps,
        "For the independent BacDive MG1655 check",
        "For the independent BacDive MG1655 check, iML1515 and paired MG1655 gapseq reconstructions were simulated under aerobic single-substrate FBA conditions after exchange-ID harmonization. iML1515 used a BiGG-compatible M9-like medium, whereas gapseq models used their native media after removal of default carbon-like exchanges before adding each BacDive substrate. The default MG1655 gapseq SBML initially failed COBRApy parsing because libSBML could not read the non-ASCII source path; copying the same SBML file to an ASCII-only project path resolved the issue, and the model was included in the benchmark. Among resolved exchanges in the carbon-assimilation subset, the default reconstruction correctly classified 6/9 substrates (accuracy = 66.7%, sensitivity = 100.0%, specificity = 0.0%, balanced accuracy = 0.500), whereas both iML1515 and the optimized reconstruction correctly classified 7/9 substrates (accuracy = 77.8%, sensitivity = 100.0%, specificity = 33.3%, balanced accuracy = 0.667, MCC = 0.500). Across the full resolved BacDive panel, balanced accuracy was 0.500 for the default reconstruction (n = 13), 0.638 for the optimized reconstruction (n = 13), and 0.604 for iML1515 (n = 14). These results provide a small independent phenotype sanity check and suggest that optimization reduced some false-positive substrate predictions, but their limited sample size and mixed API 20E assay types preclude strong claims of broad phenotype accuracy.",
    )

    replace_para(
        ps,
        "For the available paired gapseq reconstructions, the SD01 consistency benchmark",
        "For the available paired gapseq reconstructions, the SD01 consistency benchmark showed modest improvement on average but not uniform superiority of the optimized workflow. Mean balanced accuracy increased from 0.562 for default gapseq reconstructions to 0.592 for optimized reconstructions, and mean MCC increased from 0.15 to 0.25. These SD01 results should be read together with the BacDive MG1655 phenotype check (Supplementary Figure S2 and Supplementary Data S8): SD01 evaluates consistency with a published computational matrix, whereas BacDive provides a small independent experimental utilization panel in which the optimized MG1655 reconstruction improved specificity relative to the default gapseq reconstruction.",
    )

    doc.save(DOCX)
    print(f"Updated {DOCX}")
    print(f"Backup {BACKUP}")


if __name__ == "__main__":
    main()
