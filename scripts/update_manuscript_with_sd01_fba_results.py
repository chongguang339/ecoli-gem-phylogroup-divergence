from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import shutil

from docx import Document
from docx.text.paragraph import Paragraph


DOCX = Path(r"D:\pipeline\analysis_v5\02.论文\Evolutionary Divergence of Genomic Defense.docx")
BACKUP = DOCX.with_name("Evolutionary Divergence of Genomic Defense.before_sd01_770_fba_update.docx")
TEXT_INDEX = DOCX.with_name("docx_text_index_after_sd01_770_fba_update.txt")


REPLACEMENTS = {
    "Results: We developed and deployed an optimized high-throughput annotation-and-reconstruction pipeline on 10,951": (
        "Results: We developed and deployed an optimized high-throughput annotation-and-reconstruction pipeline on "
        "10,951 high-quality E. coli genomes passing CheckM2 quality thresholds, from which 1,000 representative genomes "
        "were selected for multi-omic functional profiling and GEM reconstruction. The corrected phylogroup design "
        "substantially increased B2 representation (A = 298, B2 = 244, B1 = 200, D = 103, E = 73, C = 36, F = 31, "
        "G = 11, U = 4). In the updated benchmark, an N50-stratified 20-genome stress-test subset showed a median "
        "15.81-fold acceleration of optimized gapseq relative to the default workflow (IQR 15.40-16.35; range "
        "14.78-17.77), reducing median wall-clock time from 19,821.5 s to 1,260.5 s. Automated mapping of Monk et al. "
        "SD01 Table 5 carbon-source conditions resolved 738/770 conditions to model-readable exchange reactions "
        "(364 unique BiGG exchanges and 314 unique gapseq/ModelSEED exchanges). Across five curated BiGG models, "
        "770-condition FBA reproduced the SD01 in silico growth matrix with high agreement for iJO1366, iML1515, "
        "iECs1327/iECSF_1327, and iUTI89_1310 (mean balanced accuracy 0.964-0.984), whereas ic_1306 showed a strong "
        "aerobic false-positive burden despite high anaerobic agreement. Exploratory paired gapseq benchmarking across "
        "available default/optimized reconstructions showed modest mean improvement in balanced accuracy for optimized "
        "models (0.592 versus 0.562), but with strain- and oxygen-condition-specific tradeoffs; these SD01 data are "
        "therefore treated as an in silico consistency benchmark rather than independent experimental Biolog validation. "
        "Core-genome FastANI distances revealed that pathway composition is strongly structured by continuous genomic "
        "relatedness (KEGG Mantel Spearman r = 0.760, p = 0.001), whereas discrete phylogroup labels retained only a "
        "small KEGG residual effect after ANI correction (unique R2 = 0.017, p = 0.001). Larger residual phylogroup "
        "effects persisted for COG composition (unique R2 = 0.062) and SignalP/TMHMM surface architecture "
        "(unique R2 = 0.081)."
    ),
    "Structural reference-model benchmarking was performed for the four curated-reference genomes available": (
        "Structural reference-model benchmarking was performed for the four curated-reference genomes available in the "
        "current `modle_compare` directory: iB21_1397, iEC1364_W, iECDH10B_1368, and iJN1463. For each genome, default "
        "and optimized gapseq reconstructions were compared with the corresponding BiGG reference model by reaction "
        "counts and by the fraction of reactions removed by optimization that were nevertheless present in the curated "
        "reference. We further implemented a reproducible SD01 carbon-source FBA benchmark. SD01 Table 5 condition names "
        "were parsed from 770 aerobic/anaerobic carbon-source environments and translated into exchange reactions using "
        "the SD01 Table 3 metabolite dictionary. The resulting map contained both BiGG-style exchange identifiers "
        "(e.g., EX_glc__D_e) for curated BiGG models and gapseq/ModelSEED-style exchange identifiers "
        "(e.g., EX_cpd00027_e0) for gapseq reconstructions. FBA was run with COBRApy under an M9-like base medium, "
        "excluding glucose from the base medium and adding one carbon source at a time, with oxygen uptake enabled for "
        "C_AERO conditions and disabled for C_ANAEROBIC conditions. Prediction matrices, growth-rate matrices, and "
        "confusion-matrix metrics were exported for all successfully loaded models. Because SD01 Table 5 is explicitly "
        "an in silico growth-screen matrix from Monk et al. rather than an independent experimental Biolog table, the "
        "analysis was interpreted as model-to-literature consistency benchmarking, not as experimental phenotype "
        "validation."
    ),
    "Phenotype benchmarking remains the most important unresolved validation layer.": (
        "Phenotype benchmarking remains the most important unresolved validation layer, but the revised analysis now "
        "separates three evidence levels. First, the original MG1655 35-carbon-source M9 panel remains a single-strain "
        "proof-of-concept. Second, the new SD01 770-condition benchmark tests whether curated BiGG models and paired "
        "gapseq reconstructions reproduce the Monk et al. in silico growth-screen matrix after automated exchange-ID "
        "mapping. Third, true experimental Biolog or carbon-source data are still required before claiming independent "
        "phenotypic accuracy across strains."
    ),
    "For MG1655, the M9-constrained optimized model achieved accuracy = 82.86%": (
        "For MG1655, the M9-constrained optimized model achieved accuracy = 82.86% and sensitivity = 96.55% "
        "(TP = 28, TN = 1, FP = 5, FN = 1; 35-carbon-source panel). The SD01 Table 5 carbon-source analysis expanded "
        "the computable benchmark to 770 aerobic/anaerobic environments. Automated mapping resolved 738/770 conditions "
        "to exchange reactions, leaving 32 unmapped conditions that were retained as missing-exchange predictions. "
        "Across five curated BiGG models, carbon-source FBA showed high agreement with the SD01 in silico matrix for "
        "iECs1327/iECSF_1327 (balanced accuracy 0.989 aerobic and 0.979 anaerobic), iJO1366 (0.979 and 0.963), iML1515 "
        "(0.979 and 0.948), and iUTI89_1310 (0.989 and 0.979). The ic_1306 CFT073 model retained high anaerobic "
        "agreement (balanced accuracy 0.969) but showed poor aerobic specificity (balanced accuracy 0.557; FP = 183), "
        "identifying a model- or medium-specific false-positive burden that requires targeted curation rather than "
        "being hidden in an aggregate score."
    ),
    "Figure 1. Stratified runtime and reference-model benchmarking.": (
        "Figure 1. Stratified runtime, reference-model, and carbon-source consistency benchmarking. (A) Speedup "
        "distribution across 20 assemblies spanning Actual_N50 = 22,535-675,909 bp and CDS_Count = 4,252-5,584; the "
        "subset includes A (n = 4), B1 (n = 7), B2 (n = 4), C (n = 1), E (n = 2), cryptic/U (n = 2). (B) Default "
        "versus optimized wall-clock runtime distributions. (C) Reaction-count compaction and reaction retention in "
        "the 20-genome stress-test subset; paired XML reaction matrices further identify optimized-only and "
        "default-only reactions for every genome. (D) Four available curated-reference genomes show that most reactions "
        "removed by the optimized workflow were not present in the matched BiGG reference model. (E) Automated SD01 "
        "Table 5 carbon-source mapping and 770-condition FBA provide an in silico model-consistency benchmark, not an "
        "experimental Biolog validation dataset."
    ),
    "Of the 1,000 attempted GEM reconstructions, all 1,000 models achieved viable growth": (
        "Of the 1,000 attempted GEM reconstructions, all 1,000 models achieved viable growth (positive biomass flux under "
        "rich medium conditions) following M9-constrained gap-filling. This outcome is an expected consequence of "
        "gap-filling by design—which adds reactions until biomass production is feasible—and is therefore not interpreted "
        "as an independent quality metric; the relevant measures of model fidelity are the single-strain M9 carbon-source "
        "proof-of-concept and the SD01 770-condition in silico consistency benchmarks described in Section 3.1.3."
    ),
    "Biomass flux values among viable models range from 0.526": (
        "Biomass flux values among viable models range from 0.526 to 0.584 mmol gDW-1 h-1 under standardized rich medium "
        "conditions, with median flux = 0.5826 mmol gDW-1 h-1, closely matching the MG1655 reference value "
        "(0.5826 mmol gDW-1 h-1). This consistency indicates that the optimized workflow produced solvable biomass "
        "equations under the standardized rich-medium test condition, but it should not be treated as proof of broad "
        "phenotype accuracy; substrate-specific behavior is better assessed by the carbon-source benchmarks described "
        "above."
    ),
    "Several limitations remain. First, pathway names are represented by KEGG identifiers only;": (
        "Several limitations remain. First, pathway names are represented by KEGG identifiers only; pathway-name annotation "
        "would improve biological interpretability of the 15 corrected signals. Second, phylogroup assignment was "
        "performed with EzClermont/ClermonTyping-compatible classification before cohort selection; although this removes "
        "the original B2 sampling defect, a future MLST, Mash/ANI reference-placement, or core-gene tree concordance "
        "analysis would provide an additional independent phylogenomic check. Third, the residual phylogroup analysis uses "
        "ANI PCoA covariates and label-permutation sensitivity tests rather than a full core-gene tree with formal "
        "phylogenetic comparative modeling; it therefore distinguishes continuous genomic relatedness from categorical "
        "labels but does not replace tree-based evolutionary inference. Fourth, although TMHMM coverage was improved from "
        "846 to 999 genomes after reprocessing, one B1 genome still lacks detailed TMHMM output and was treated as missing "
        "(NA), not as a biological zero, in membrane-protein analyses; DeepTMHMM was evaluated as an alternative but was "
        "not adopted for the full cohort because its runtime was incompatible with the current high-throughput design. "
        "Fifth, CRISPR target matching depends on the completeness and taxonomic composition of the phage and plasmid "
        "databases used for spacer searches, and phage/plasmid labels can overlap when homologous mobile-element sequences "
        "are represented in both databases; therefore target absence should be interpreted as no detectable match in the "
        "queried databases rather than true absence of prior exposure. Sixth, the SD01 carbon-source benchmark uses the "
        "Monk et al. Table 5 in silico growth-screen matrix, not independent experimental Biolog data; 32/770 conditions "
        "remained unmapped automatically, iECs1327 required the available iECSF_1327/SE15 proxy, and several gapseq "
        "reconstructions showed oxygen-condition-specific sensitivity losses. Finally, all GEM predictions remain in "
        "silico and should be experimentally validated in representative strains."
    ),
    "In conclusion, this revised analysis presents an optimized and more carefully validated": (
        "In conclusion, this revised analysis presents an optimized and more carefully validated E. coli GEM reconstruction "
        "framework applied to a corrected 1,000-genome phylogroup cohort. Stratified benchmarking across 20 assemblies "
        "showed a median 15.81-fold runtime acceleration, paired XML reaction matrices identified 111.0 optimized-only "
        "reactions per genome at the median, and four available curated-reference genomes showed a median 12.76-fold "
        "acceleration with mostly non-reference-supported reactions removed by optimization. Automated SD01 carbon-source "
        "mapping translated 738/770 aerobic/anaerobic environments into model-readable exchanges and produced reproducible "
        "FBA prediction matrices for curated BiGG and available gapseq models, providing a transparent in silico "
        "consistency benchmark while preserving the need for independent experimental validation. The rebalanced cohort "
        "eliminated the earlier B2 sampling defect and substantially weakened the claim of D-specific virulence enrichment. "
        "FastANI residual analysis further showed that KEGG pathway composition is strongly tied to continuous core-genome "
        "distance, whereas discrete phylogroup labels retain smaller residual effects for KEGG, COG, GEM topology, CRISPR "
        "spacer composition, CRISPR target burden, and SignalP/TMHMM surface architecture. The main remaining validation "
        "priority is true multi-strain experimental Biolog or carbon-source data for functional accuracy benchmarking."
    ),
    "Availability of data and materials:": (
        "Availability of data and materials: All genome accession numbers, phylogroup assignments, annotation outputs, GEM "
        "topology metrics, benchmark metadata, reaction-set matrices, SD01 carbon-source exchange mappings, and FBA "
        "prediction/accuracy matrices are available as Supplementary Data Files. Genome sequences were obtained from NCBI "
        "Assembly (https://www.ncbi.nlm.nih.gov/assembly/). Modified gapseq scripts and pipeline configuration files are "
        "available at [repository URL]."
    ),
}


INSERT_AFTER_PREFIX = {
    "Structural reference-model benchmarking was performed for the four curated-reference genomes available": (
        "For the SD01 770-condition BiGG-model benchmark, five curated models were evaluated: iJO1366, iML1515, ic_1306, "
        "iUTI89_1310, and the available 1327-gene proxy iECSF_1327 for the requested iECs1327 comparison. Public BiGG "
        "SBML files were used for these curated models. For available paired gapseq reconstructions, default and optimized "
        "models were evaluated for iB21_1397, iEC1364_W, iECDH10B_1368, and MG1655/iML1515. Several local reference-model "
        "SBML files in `modle_compare` were not COBRApy-readable and were therefore recorded in the load-status table "
        "rather than forced into the analysis. All model-load status, mapped-exchange coverage, simulation status, and "
        "confusion-matrix outputs are reported in the supplementary benchmark tables."
    ),
    "For MG1655, the M9-constrained optimized model achieved accuracy = 82.86%": (
        "For the available paired gapseq reconstructions, the SD01 consistency benchmark showed modest improvement on "
        "average but not uniform superiority of the optimized workflow. Mean balanced accuracy increased from 0.562 for "
        "default gapseq reconstructions to 0.592 for optimized reconstructions, and mean MCC increased from 0.145 to "
        "0.250. The largest gains occurred for iECDH10B_1368 under both aerobic and anaerobic conditions, whereas the "
        "MG1655/iML1515 optimized reconstruction lost anaerobic sensitivity relative to the default reconstruction. These "
        "results support the optimized pipeline as a faster and more compact reconstruction workflow, but they argue "
        "against claiming universal phenotype-accuracy improvement without strain-specific curation and experimental "
        "truth data."
    ),
}


def replace_paragraph(paragraph, text: str) -> None:
    for run in paragraph.runs:
        run.text = ""
    if paragraph.runs:
        paragraph.runs[0].text = text
    else:
        paragraph.add_run(text)


def insert_paragraph_after(paragraph, text: str):
    new_p = deepcopy(paragraph._p)
    paragraph._p.addnext(new_p)
    inserted = Paragraph(new_p, paragraph._parent)
    inserted.clear()
    inserted.style = paragraph.style
    inserted.add_run(text)
    return inserted


def main() -> None:
    if not BACKUP.exists():
        shutil.copy2(DOCX, BACKUP)

    doc = Document(DOCX)

    replaced = {}
    for para in doc.paragraphs:
        text = para.text.strip()
        for prefix, new_text in REPLACEMENTS.items():
            if text.startswith(prefix):
                replace_paragraph(para, new_text)
                replaced[prefix] = True

    for prefix, insert_text in INSERT_AFTER_PREFIX.items():
        for para in doc.paragraphs:
            if para.text.strip().startswith(prefix):
                insert_paragraph_after(para, insert_text)
                break

    missing = [k for k in REPLACEMENTS if k not in replaced]
    if missing:
        raise RuntimeError("Missing replacement anchors: " + "; ".join(missing))

    doc.save(DOCX)

    refreshed = Document(DOCX)
    with TEXT_INDEX.open("w", encoding="utf-8") as handle:
        for i, para in enumerate(refreshed.paragraphs):
            if para.text.strip():
                handle.write(f"[{i}] {para.text.strip()}\n")

    print(f"Updated {DOCX}")
    print(f"Backup: {BACKUP}")
    print(f"Text index: {TEXT_INDEX}")


if __name__ == "__main__":
    main()
