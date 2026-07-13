from pathlib import Path
import shutil

from docx import Document


DOCX = Path(r"D:\pipeline\analysis_v5\02.论文\Evolutionary Divergence of Genomic Defense.docx")
BACKUP = DOCX.with_name(DOCX.stem + ".before_submission_rhetoric_cleanup.docx")


REPLACEMENTS = {
    "Phylogroup assignment was performed upstream on the full set of 10,951 downloaded complete E. coli assemblies using EzClermont/ClermonTyping-compatible in-silico phylogroup classification [7]. The run script `qc0_run_ezclermont.py` used the input directory `/home/u2024110181382/Project/02.E.coli/drep_result/high_quality_genomes`, called `ezclermont` with 10 parallel workers, parsed the `Clermont type` field from stdout/stderr, and wrote `Phylogroup_Results.csv` for downstream sampling. All 10,951 imported genomes retained valid phylogroup output after invalid-type filtering.": (
        "Phylogroup assignment was performed upstream on the full set of 10,951 downloaded complete E. coli assemblies using EzClermont/ClermonTyping-compatible in-silico phylogroup classification [7]. The ezclermont pipeline was executed via a customized Python wrapper script (`qc0_run_ezclermont.py`) that used standardized relative directories for genome assemblies, called `ezclermont` with 10 parallel workers, parsed the `Clermont type` field from stdout/stderr, and generated a phylogroup result table for downstream sampling. All 10,951 imported genomes retained valid phylogroup output after invalid-type filtering."
    ),
    "For each of the 1,000 phylogroup-stratified genomes, optimized M9-constrained models were assessed for topology, solvability, and COBRApy-derived quality metrics. The primary topology fields were number of genes, reactions, and metabolites; biomass flux under rich medium was used as a solvability check. In addition, a cohort-level QC scan quantified total reactions, total metabolites, mass-balanced reaction ratio, blocked-reaction ratio, and dead-end metabolite ratio for all reconstructed GEMs. Two representative models were further inspected with MEMOTE v0.13.0 and the resulting interactive HTML reports were retained as supplementary files. MEMOTE annotation failures were interpreted as database-namespace compliance findings, particularly for gapseq/ModelSEED identifiers, rather than as direct evidence of metabolic inconsistency.": (
        "For each of the 1,000 phylogroup-stratified genomes, optimized M9-constrained models were assessed using a two-tier GEM quality-control strategy. At the cohort scale, topology and solvability were summarized using the number of genes, reactions, and metabolites, biomass flux under a standardized rich-medium condition, and COBRApy-derived reaction/metabolite consistency metrics, including mass-balanced reaction ratio, blocked-reaction ratio, and dead-end metabolite ratio. At the representative-model scale, two phylogroup-sampled reconstructions were additionally inspected with MEMOTE v0.13.0, and the resulting interactive HTML reports were retained as supplementary files. This design emphasizes reproducible phylogroup-scale topology while retaining deeper stoichiometric and annotation audits for representative models; MEMOTE annotation findings were interpreted as database-namespace compliance results, particularly for gapseq/ModelSEED identifiers, rather than as direct evidence of metabolic inconsistency."
    ),
    "The main interpretation is bounded by four constraints. First, KEGG identifiers define pathway-level signals and require pathway-specific annotation before mechanistic claims. Second, PERMANOVA results are interpreted together with PERMDISP and ANI residual models, so PERMANOVA R2 is not treated as pure categorical phylogroup-explained variance. Third, the BacDive MG1655 panel is a small API 20E-derived phenotype check and does not replace larger Biolog or phenotype-microarray validation. Fourth, COBRApy/MEMOTE QC supports stoichiometric and topological consistency but does not replace manual model curation; blocked reactions and dead-end metabolites remain expected features of automated large-cohort reconstructions.": (
        "The main interpretation is bounded by four constraints. First, KEGG identifiers define pathway-level signals and require pathway-specific annotation before mechanistic claims. Second, PERMANOVA results are interpreted together with PERMDISP and ANI residual models, so PERMANOVA R2 is not treated as pure categorical phylogroup-explained variance. Third, the BacDive MG1655 panel is a small API 20E-derived phenotype check and does not replace larger Biolog or phenotype-microarray validation. The low specificity observed for substrate utilization, especially in the default reconstruction, is biologically and methodologically informative: automated gap-filling tends to encode metabolic potential from genomic content, whereas experimental API 20E or carbon-source assays report realized physiology under specific regulatory and culture conditions. In E. coli, genes supporting a nominal carbon-utilization route can be present but transcriptionally repressed, conditionally inducible, transporter-limited, or disconnected from biomass production under the assay environment. Therefore, false-positive substrate calls should be interpreted as the expected gap between genome-encoded metabolic capacity and expressed physiological activity, not simply as binary reconstruction failure. Fourth, COBRApy/MEMOTE QC supports stoichiometric and topological consistency but does not replace manual model curation; blocked reactions and dead-end metabolites remain expected features of automated large-cohort reconstructions."
    ),
}


def replace_text(doc: Document) -> None:
    misses = []
    for old, new in REPLACEMENTS.items():
        hit = False
        for para in doc.paragraphs:
            if para.text == old:
                para.text = new
                hit = True
                break
        if not hit:
            misses.append(old[:120])
    if misses:
        raise RuntimeError("Missing target paragraphs:\n" + "\n".join(misses))


def main() -> None:
    if not BACKUP.exists():
        shutil.copy2(DOCX, BACKUP)
    doc = Document(DOCX)
    replace_text(doc)
    doc.save(DOCX)
    print(f"Updated {DOCX}")
    print(f"Backup {BACKUP}")


if __name__ == "__main__":
    main()
