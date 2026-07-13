from __future__ import annotations

import shutil
from pathlib import Path

from docx import Document


ROOT = Path(r"D:\pipeline\analysis_v5")
DOCX = ROOT / "02.论文" / "Evolutionary Divergence of Genomic Defense.docx"
BACKUP = ROOT / "02.论文" / "Evolutionary Divergence of Genomic Defense.before_user_clarifications.docx"
OUT_TEXT_INDEX = ROOT / "02.论文" / "docx_text_index_after_user_clarifications.txt"


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
        "The quality-filtered source dataset used for sampling comprised 10,951 Escherichia coli genome assemblies from the NCBI Assembly database. Genome quality was assessed using CheckM2 v1.0.2 with the neural network-specific completeness model [12], and only genomes passing strict quality thresholds (completeness >= 95%, contamination <= 5%, coding density >= 0.85) were retained for downstream sampling. This source set had a median genome size of 5.13 Mb, median GC content of 50.7%, and a median of 4,863 coding sequences (Table 2). Representative genome quality metrics are available in `quality_report.tsv`; the pre-filter retrieval count should be reported if it differs from this quality-filtered source set.",
        "A total of 10,951 complete Escherichia coli genome assemblies were retrieved from the NCBI Assembly database after applying the download-time completeness filter. Genome quality was then assessed using CheckM2 v1.0.2 with the neural network-specific completeness model [12], and genomes passing strict quality thresholds (completeness >= 95%, contamination <= 5%, coding density >= 0.85) were retained for downstream sampling. This source set had a median genome size of 5.13 Mb, median GC content of 50.7%, and a median of 4,863 coding sequences (Table 2). Representative genome quality metrics are available in `quality_report.tsv`.",
    )

    replace_exact(
        doc,
        "Several limitations remain. First, pathway names are represented by KEGG identifiers only; pathway-name annotation would improve biological interpretability of the 15 corrected signals. Second, phylogroup assignment currently relies on ClermonTyping output, and an independent concordance analysis using EzClermont, MLST, Mash/ANI placement against typed references, or a core-gene tree should be added before making strong taxonomic claims. Third, the residual phylogroup analysis uses ANI PCoA covariates and label-permutation sensitivity tests rather than a full core-gene tree with formal phylogenetic comparative modeling; it therefore distinguishes continuous genomic relatedness from categorical labels but does not replace tree-based evolutionary inference. Fourth, although TMHMM coverage was improved from 846 to 999 genomes after reprocessing, one B1 genome still lacks detailed TMHMM output and was treated as missing in membrane-protein analyses; a DeepTMHMM or equivalent modern topology reanalysis would further strengthen membrane-protein inference. Fifth, CRISPR spacer sharing is sequence-level but not yet matched to phage/protospacer databases, so spacer targets and ecological exposures remain unresolved. Finally, all GEM predictions remain in silico and should be experimentally validated in representative strains.",
        "Several limitations remain. First, pathway names are represented by KEGG identifiers only; pathway-name annotation would improve biological interpretability of the 15 corrected signals. Second, phylogroup assignment currently relies on ClermonTyping output, and an independent concordance analysis using EzClermont, MLST, Mash/ANI placement against typed references, or a core-gene tree should be added before making strong taxonomic claims. Third, the residual phylogroup analysis uses ANI PCoA covariates and label-permutation sensitivity tests rather than a full core-gene tree with formal phylogenetic comparative modeling; it therefore distinguishes continuous genomic relatedness from categorical labels but does not replace tree-based evolutionary inference. Fourth, although TMHMM coverage was improved from 846 to 999 genomes after reprocessing, one B1 genome still lacks detailed TMHMM output and was treated as missing in membrane-protein analyses; DeepTMHMM was evaluated as an alternative but was not adopted for the full cohort because its runtime was incompatible with the current high-throughput design. Fifth, CRISPR spacer sharing is sequence-level but not yet matched to phage/protospacer databases, so spacer targets and ecological exposures remain unresolved. Finally, all GEM predictions remain in silico and should be experimentally validated in representative strains.",
    )

    doc.save(DOCX)

    with OUT_TEXT_INDEX.open("w", encoding="utf-8") as f:
        for i, p in enumerate(Document(DOCX).paragraphs):
            f.write(f"{i:04d}: {p.text}\n")

    print(f"Updated {DOCX}")
    print(f"Backup {BACKUP}")
    print(f"Wrote {OUT_TEXT_INDEX}")


if __name__ == "__main__":
    main()
