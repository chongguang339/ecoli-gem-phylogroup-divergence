from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd
from docx import Document


ROOT = Path(r"D:\pipeline\analysis_v5")
DOCX = ROOT / "02.论文" / "Evolutionary Divergence of Genomic Defense.docx"
BACKUP = ROOT / "02.论文" / "Evolutionary Divergence of Genomic Defense.before_ezclermont_crispr_update.docx"
CHECKLIST = ROOT / "03.analysis" / "reviewer_issue_resolution_checklist.csv"
OUT_TEXT_INDEX = ROOT / "02.论文" / "docx_text_index_after_ezclermont_crispr_update.txt"


def replace_exact(doc: Document, old: str, new: str) -> None:
    for p in doc.paragraphs:
        if p.text == old:
            p.text = new
            return
    raise ValueError(f"paragraph not found: {old[:120]}")


def main() -> None:
    if not BACKUP.exists():
        shutil.copy2(DOCX, BACKUP)

    doc = Document(DOCX)
    replace_exact(
        doc,
        "Phylogroup assignment for all 1,000 representative genomes was performed using the ClermonTyping method [6], which assigns E. coli isolates to phylogroups A, B1, B2, C, D, E, F, G, or cryptic/unknown categories based on a combination of quadruplex PCR markers and in-silico whole-genome typing. After correcting the original unbalanced design that severely underrepresented B2, the final re-run cohort comprised: A (n = 298), B2 (n = 244), B1 (n = 200), D (n = 103), E (n = 73), C (n = 36), F (n = 31), G (n = 11), and U (n = 4). This revised design enables group-level evaluation of B2-associated patterns and prevents D-phylogroup signals from being interpreted as a proxy for all ExPEC-associated lineages.",
        "Phylogroup assignment was performed upstream on the full set of 10,951 downloaded complete E. coli assemblies using EzClermont/ClermonTyping-compatible in-silico phylogroup classification [6]. The final 1,000-genome analytical cohort was then selected from this fully typed source set to obtain a revised, explicitly balanced phylogroup distribution: A (n = 298), B2 (n = 244), B1 (n = 200), D (n = 103), E (n = 73), C (n = 36), F (n = 31), G (n = 11), and U (n = 4). This revised design enables group-level evaluation of B2-associated patterns and prevents D-phylogroup signals from being interpreted as a proxy for all ExPEC-associated lineages.",
    )
    replace_exact(
        doc,
        "Several limitations remain. First, pathway names are represented by KEGG identifiers only; pathway-name annotation would improve biological interpretability of the 15 corrected signals. Second, phylogroup assignment currently relies on ClermonTyping output, and an independent concordance analysis using EzClermont, MLST, Mash/ANI placement against typed references, or a core-gene tree should be added before making strong taxonomic claims. Third, the residual phylogroup analysis uses ANI PCoA covariates and label-permutation sensitivity tests rather than a full core-gene tree with formal phylogenetic comparative modeling; it therefore distinguishes continuous genomic relatedness from categorical labels but does not replace tree-based evolutionary inference. Fourth, although TMHMM coverage was improved from 846 to 999 genomes after reprocessing, one B1 genome still lacks detailed TMHMM output and was treated as missing in membrane-protein analyses; DeepTMHMM was evaluated as an alternative but was not adopted for the full cohort because its runtime was incompatible with the current high-throughput design. Fifth, CRISPR spacer sharing is sequence-level but not yet matched to phage/protospacer databases, so spacer targets and ecological exposures remain unresolved. Finally, all GEM predictions remain in silico and should be experimentally validated in representative strains.",
        "Several limitations remain. First, pathway names are represented by KEGG identifiers only; pathway-name annotation would improve biological interpretability of the 15 corrected signals. Second, phylogroup assignment was performed with EzClermont/ClermonTyping-compatible classification before cohort selection; although this removes the original B2 sampling defect, a future MLST, Mash/ANI reference-placement, or core-gene tree concordance analysis would provide an additional independent phylogenomic check. Third, the residual phylogroup analysis uses ANI PCoA covariates and label-permutation sensitivity tests rather than a full core-gene tree with formal phylogenetic comparative modeling; it therefore distinguishes continuous genomic relatedness from categorical labels but does not replace tree-based evolutionary inference. Fourth, although TMHMM coverage was improved from 846 to 999 genomes after reprocessing, one B1 genome still lacks detailed TMHMM output and was treated as missing in membrane-protein analyses; DeepTMHMM was evaluated as an alternative but was not adopted for the full cohort because its runtime was incompatible with the current high-throughput design. Fifth, the available CRISPR spacer table supports sequence-level spacer sharing analysis and target-matching input generation, but no viral, prophage, plasmid, or protospacer target FASTA is currently included; therefore spacer target identities and ecological exposures remain unresolved. Finally, all GEM predictions remain in silico and should be experimentally validated in representative strains.",
    )
    doc.save(DOCX)

    if CHECKLIST.exists():
        df = pd.read_csv(CHECKLIST)
        df.loc[df["Reviewer_issue"].eq("3.2 ClermonTyping validation"), "Status"] = "clarified; optional independent check"
        df.loc[df["Reviewer_issue"].eq("3.2 ClermonTyping validation"), "Current_resolution"] = "The 1,000-genome cohort was selected after EzClermont/ClermonTyping-compatible classification of all 10,951 complete genomes."
        df.loc[df["Reviewer_issue"].eq("3.2 ClermonTyping validation"), "Data_or_action_still_needed"] = "Optional MLST, Mash/ANI reference-placement, or core-gene tree concordance if the target journal asks for an independent phylogenomic check."
        df.loc[df["Reviewer_issue"].eq("2.7 CRISPR count underpowering"), "Current_resolution"] = "Spacer-level sequence analysis is available, and unique spacer FASTA/metadata were generated for target matching."
        df.loc[df["Reviewer_issue"].eq("2.7 CRISPR count underpowering"), "Data_or_action_still_needed"] = "A viral/prophage/plasmid/protospacer target FASTA database is required for true target matching."
        df.to_csv(CHECKLIST, index=False, encoding="utf-8-sig")

    with OUT_TEXT_INDEX.open("w", encoding="utf-8") as f:
        for i, p in enumerate(Document(DOCX).paragraphs):
            f.write(f"{i:04d}: {p.text}\n")

    print(f"Updated {DOCX}")
    print(f"Backup {BACKUP}")
    print(f"Updated {CHECKLIST}")


if __name__ == "__main__":
    main()
