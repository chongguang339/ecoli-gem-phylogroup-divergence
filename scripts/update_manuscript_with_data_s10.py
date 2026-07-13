from pathlib import Path
import shutil

from docx import Document


DOCX = Path(r"D:\pipeline\analysis_v5\02.论文\Evolutionary Divergence of Genomic Defense.docx")
BACKUP = DOCX.with_name(DOCX.stem + ".before_data_s10_update.docx")


def add_table1_row(doc: Document) -> None:
    table = doc.tables[0]
    existing = [" | ".join(cell.text for cell in row.cells) for row in table.rows]
    if any("Optimized GEM model archive" in row for row in existing):
        return
    row = table.add_row()
    values = [
        "Optimized GEM model archive",
        "1,000 optimized GEMs",
        "SBML/XML models; archive SHA256 = 953f92a5d53dfe2f999d8d7ee99304d81c1e71c30e9dcb2dffff39e2bf84d7fc",
        "Data S10; exported_models.tar.gz; ecoli_1000_optimized_gem_models_manifest.csv",
    ]
    for cell, value in zip(row.cells, values):
        cell.text = value


def replace_paragraphs(doc: Document) -> None:
    replacements = {
        "Availability of data and materials: All genome accession numbers, phylogroup assignments, annotation outputs, GEM topology metrics, COBRApy/MEMOTE QC outputs, benchmark metadata, reaction-set matrices, SD01 carbon-source exchange mappings, SD01 FBA prediction/model-consistency matrices, and BacDive MG1655 phenotype-check FBA outputs are organized in the project supplementary data directory. Repository URL: pending insertion before submission.": (
            "Availability of data and materials: All genome accession numbers, phylogroup assignments, annotation outputs, GEM topology metrics, COBRApy/MEMOTE QC outputs, benchmark metadata, reaction-set matrices, SD01 carbon-source exchange mappings, SD01 FBA prediction/model-consistency matrices, BacDive MG1655 phenotype-check FBA outputs, and the 1,000 optimized GEM SBML/XML models are organized in the project supplementary data directory. The optimized GEM archive is provided as Data S10 (`exported_models.tar.gz`) with per-model and archive-level SHA256 checksums. Repository URL/DOI: to be inserted after GitHub/Zenodo or Figshare deposition before submission."
        ),
        "Supplementary Data files include: Data S1, selected 1,000-genome accession and phylogroup table; Data S2, genome assembly, COG, KEGG, GEM topology, CRISPR, SignalP, and TMHMM matrices; Data S3, statistical-test outputs with effect sizes and multiple-testing correction; Data S4, benchmark runtime and reaction-set comparison matrices; Data S5, SD01 carbon-source exchange-ID mapping tables; Data S6, SD01 770-condition FBA prediction matrices, growth-rate matrices, model inclusion and simulation status tables; Data S7, CRISPR spacer target-matching outputs; Data S8, BacDive MG1655 phenotype truth table, model-load status, FBA predictions, and accuracy metrics; Data S9, COBRApy cohort-level GEM QC tables and representative MEMOTE HTML reports.": (
            "Supplementary Data files include: Data S1, selected 1,000-genome accession and phylogroup table; Data S2, genome assembly, COG, KEGG, GEM topology, CRISPR, SignalP, and TMHMM matrices; Data S3, statistical-test outputs with effect sizes and multiple-testing correction; Data S4, benchmark runtime and reaction-set comparison matrices; Data S5, SD01 carbon-source exchange-ID mapping tables; Data S6, SD01 770-condition FBA prediction matrices, growth-rate matrices, model inclusion and simulation status tables; Data S7, CRISPR spacer target-matching outputs; Data S8, BacDive MG1655 phenotype truth table, model-load status, FBA predictions, and accuracy metrics; Data S9, COBRApy cohort-level GEM QC tables and representative MEMOTE HTML reports; Data S10, the 1,000 optimized GEM SBML/XML model archive, per-model manifest, and SHA256 checksum tables."
        ),
    }
    missing = []
    for old, new in replacements.items():
        hit = False
        for para in doc.paragraphs:
            if para.text.strip() == old:
                para.text = new
                hit = True
                break
        if not hit:
            missing.append(old[:120])
    if missing:
        raise RuntimeError("Missing target paragraphs:\n" + "\n".join(missing))


def main() -> None:
    if not BACKUP.exists():
        shutil.copy2(DOCX, BACKUP)
    doc = Document(DOCX)
    add_table1_row(doc)
    replace_paragraphs(doc)
    doc.save(DOCX)
    print(f"Updated {DOCX}")
    print(f"Backup {BACKUP}")


if __name__ == "__main__":
    main()
