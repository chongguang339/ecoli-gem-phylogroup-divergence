from __future__ import annotations

from pathlib import Path

from docx import Document
from docx.enum.table import WD_TABLE_ALIGNMENT, WD_CELL_VERTICAL_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Pt


DOCX = Path(r"D:\pipeline\analysis_v5\02.论文\Evolutionary Divergence of Genomic Defense.docx")
TEXT_INDEX = DOCX.with_name("docx_text_index_after_sd01_770_fba_update.txt")


TABLE1_ROWS = [
    ["Layer", "Parameter or operation", "Default or baseline", "Optimized setting", "Purpose"],
    ["Performance", "Alignment backend", "BLAST-based search", "DIAMOND v2.1.8", "Reduce protein-alignment runtime while retaining high-scoring annotation hits"],
    ["Performance", "Thread allocation", "Concurrent subprocesses could each request the full thread count", "Threads split between gapseq_find.sh and transporter.sh", "Avoid CPU oversubscription during reconstruction"],
    ["Performance", "Contig reuse", "Iterative grep-based scans", "Sorted set-difference using comm -23", "Reduce redundant I/O for fragmented assemblies"],
    ["Performance", "Protein hashing", "Repeated MD5 calculation and disk reads", "Unique-sequence MD5 deduplication and cache reuse", "Avoid repeated processing of identical proteins"],
    ["Workflow", "Query sequence clustering", "No mandatory identical-protein clustering", "-S identical", "Collapse 100% identical proteins before alignment"],
    ["Workflow", "Primary reference database", "Broad/global database search", "Enterobacterales fast-pass database", "Prioritize clade-relevant metabolic hits"],
    ["Model behavior", "Gap-filling medium", "Unconstrained or broad gap-filling", "Defined M9 glucose minimal medium", "Limit systematic over-prediction of substrate utilization"],
    ["Model behavior", "Taxonomic range", "auto", "all", "Maintain consistent pathway detection scope across assemblies"],
]

TABLE2_ROWS = [
    ["Dataset or analysis layer", "N or scope", "Key values", "Output used in manuscript"],
    ["NCBI source assemblies", "10,951 complete E. coli assemblies", "Download-time complete-assembly filter; CheckM2 completeness >= 95%, contamination <= 5%, coding density >= 0.85", "quality_report.tsv"],
    ["Representative cohort", "1,000 genomes", "A = 298, B2 = 244, B1 = 200, D = 103, E = 73, C = 36, F = 31, G = 11, U = 4", "selected_1000_genomes.tsv; recalculated_annotations.csv"],
    ["N50 stress-test benchmark", "20 genomes", "Actual_N50 = 22,535-675,909 bp; CDS_Count = 4,252-5,584; median speedup = 15.81-fold", "benchmark_metrics.csv; all_20_genomes_metadata.csv"],
    ["Curated-reference structural benchmark", "4 available matched reference genomes", "iB21_1397, iEC1364_W, iECDH10B_1368, iJN1463; median speedup = 12.76-fold", "modle_compare; reaction_set_compare"],
    ["SD01 carbon-source mapping", "770 aerobic/anaerobic carbon-source conditions", "738 mapped; 32 unmapped; 364 unique BiGG exchanges; 314 unique gapseq/ModelSEED exchanges", "sd01_770_carbon_conditions_exchange_id_mapping.csv"],
    ["SD01 curated-model FBA", "5 curated BiGG models", "iJO1366, iML1515, ic_1306, iUTI89_1310, iECs1327/iECSF_1327 proxy", "selected_models_770_carbon_fba_predictions_long.csv"],
    ["SD01 paired gapseq FBA", "Available default/optimized paired reconstructions", "Mean balanced accuracy default = 0.562; optimized = 0.592", "available_benchmark_770_fba_accuracy_metrics.csv"],
]

SUPP_TEXT = (
    "Supplementary Data files include: Data S1, selected 1,000-genome accession and phylogroup table; Data S2, genome "
    "assembly, COG, KEGG, GEM topology, CRISPR, SignalP, and TMHMM matrices; Data S3, statistical-test outputs with "
    "effect sizes and multiple-testing correction; Data S4, benchmark runtime and reaction-set comparison matrices; "
    "Data S5, SD01 carbon-source exchange-ID mapping tables; Data S6, SD01 770-condition FBA prediction matrices, "
    "growth-rate matrices, model-load status, and accuracy metrics; Data S7, CRISPR spacer target-matching outputs."
)


def find_para(doc: Document, prefix: str):
    for para in doc.paragraphs:
        if para.text.strip().startswith(prefix):
            return para
    raise RuntimeError(f"Could not find paragraph starting with: {prefix}")


def insert_paragraph_after(paragraph, text: str, style=None):
    new_p = paragraph._p.__class__()
    paragraph._p.addnext(new_p)
    inserted = paragraph._parent.paragraphs[paragraph._parent.paragraphs.index(paragraph) + 1]
    inserted.text = text
    if style:
        inserted.style = style
    return inserted


def move_table_after(table, paragraph):
    tbl = table._tbl
    tbl.getparent().remove(tbl)
    paragraph._p.addnext(tbl)


def add_table_after(doc: Document, paragraph, caption: str, rows: list[list[str]]) -> None:
    caption_para = doc.add_paragraph(caption)
    caption_para.style = doc.styles["Normal"]
    caption_para.runs[0].bold = True
    caption_para.paragraph_format.space_before = Pt(6)
    caption_para.paragraph_format.space_after = Pt(3)
    move_para = caption_para._p
    move_para.getparent().remove(move_para)
    paragraph._p.addnext(move_para)

    table = doc.add_table(rows=len(rows), cols=len(rows[0]))
    if "Table Grid" in [style.name for style in doc.styles]:
        table.style = "Table Grid"
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.autofit = True
    for r_idx, row in enumerate(rows):
        cells = table.rows[r_idx].cells
        for c_idx, value in enumerate(row):
            cell = cells[c_idx]
            cell.text = value
            cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
            for p in cell.paragraphs:
                p.alignment = WD_ALIGN_PARAGRAPH.CENTER if r_idx == 0 else WD_ALIGN_PARAGRAPH.LEFT
                for run in p.runs:
                    run.font.size = Pt(8)
                    if r_idx == 0:
                        run.bold = True
    move_table_after(table, caption_para)


def append_supplementary_after_availability(doc: Document) -> None:
    if any(p.text.strip().startswith("Supplementary Data files include:") for p in doc.paragraphs):
        return
    availability = find_para(doc, "Availability of data and materials:")
    new_para = doc.add_paragraph(SUPP_TEXT)
    new_para.style = availability.style
    p = new_para._p
    p.getparent().remove(p)
    availability._p.addnext(p)


def main() -> None:
    doc = Document(DOCX)
    if not any(p.text.strip().startswith("Table 1. Optimized gapseq") for p in doc.paragraphs):
        anchor1 = find_para(doc, "A complete parameter specification is provided in Table 1.")
        add_table_after(doc, anchor1, "Table 1. Optimized gapseq reconstruction parameters and rationale.", TABLE1_ROWS)
    if not any(p.text.strip().startswith("Table 2. Cohort") for p in doc.paragraphs):
        anchor2 = find_para(doc, "From this candidate pool, 1,000 genomes were selected")
        add_table_after(doc, anchor2, "Table 2. Cohort composition, benchmark datasets, and main output files.", TABLE2_ROWS)
    append_supplementary_after_availability(doc)
    doc.save(DOCX)

    refreshed = Document(DOCX)
    with TEXT_INDEX.open("w", encoding="utf-8") as handle:
        for i, para in enumerate(refreshed.paragraphs):
            if para.text.strip():
                handle.write(f"[{i}] {para.text.strip()}\n")
    print(f"Updated {DOCX}")


if __name__ == "__main__":
    main()
