# E. coli phylogroup GEM study: reproducibility repository

The Round 4 mSystems revision adds the full-cohort MEMOTE audit, paired broad/focused models, matched-aligner control, draft-by-medium cross, reaction provenance, and corrected KEGG analysis. See [Round 4 release index](docs/ROUND4_RELEASE.md) for code-to-data mapping and the release limits. The authoritative 1,000 annotated final SBML models and manifest are in [Zenodo version 2](https://doi.org/10.5281/zenodo.22782392). Corrected supplementary tables and Round 4 evidence are in [Zenodo version 4](https://doi.org/10.5281/zenodo.24033312). The [all-versions DOI](https://doi.org/10.5281/zenodo.21341419) resolves to the latest version.

## Earlier Revision 3 content

This directory contains the code-facing additions for the mSystems resubmission of the E. coli phylogroup-resolved GEM study. The persistent data record is [Zenodo 10.5281/zenodo.21341419](https://doi.org/10.5281/zenodo.21341419).

## What changed

- Rebuilt Figure 1D from direct child reaction identifiers in authoritative paired SBML models.
- Restricted the BacDive MG1655 comparison to nine successfully mapped carbon sources; unresolved amygdalin is retained only in the mapping audit.
- Added KEGG feature-construction rules and the original 264-to-115 ontology-filter audit; Round 4 retains map01053 and updates the filtered matrix to 116 features.
- Added ANI PCoA sensitivity analyses for 5, 10, 20, and 40 positive-eigenvalue axes.
- Added the 48-component chemically defined complete-medium recipe, uptake bounds, and checksums.
- Added a 120-row six-configuration benchmark provenance table. Exact command, CPU, memory, version, checksum, and log fields are recorded for V2-V5. Legacy V1/V6 scheduler CPU and memory fields were not recorded and are explicitly marked unavailable.
- Regenerated all main and supplementary figures from authoritative tables.

## Layout

- `scripts/`: analysis, model-archive, numerical-audit, and final figure-generation scripts used for the revision. Internal manuscript-formatting utilities are intentionally excluded because they are not part of the scientific workflow.
- `data/`: compact reviewer-facing outputs, including ANI-axis sensitivity and benchmark provenance summaries.

Large assets are deposited in Zenodo rather than duplicated in GitHub. Version 2 contains the authoritative 1,000 annotated final SBML models and model manifest. Version 4 contains all 1,000 MEMOTE HTML reports and run evidence, 180 paired broad/focused reconstructions, 20 matched-aligner controls, 18 draft-by-medium controls, original per-protein annotations, exact 180-genome FASTAs and corrected Supplementary Data Sets S1 and S2. The earlier 18-model MEMOTE display is retained as a representative figure.

Scripts that reconstruct revision-wide products resolve the project root from the `ECOLI_GEM_ROOT` environment variable. The final figure script additionally accepts `ECOLI_FIGURE_WORKDIR`. These variables replace machine-specific paths; the expected input tables are identified in each script and the Zenodo workbook/package index.

## ANI axis sensitivity

```bash
python scripts/run_ani_axis_sensitivity.py \
  --ani-distance core_genome_ani_distance_matrix_by_genome_uid.csv.gz \
  --kegg-binary kegg_pathway_genome_matrix_bacterial_filtered_binary.csv \
  --output-dir ani_axis_sensitivity \
  --permutations 999 \
  --seed 20260912
```

The primary manuscript model uses 20 ANI PCoA axes. Sensitivity analyses show that 5 axes under-adjust genomic structure, whereas the unique phylogroup contribution remains small with 10, 20, or 40 axes.

## Provenance boundary

Reaction counts in the revised benchmark are direct-child SBML reaction counts. Historical JSON fields based on XML string matching are not used as reaction-count evidence. Missing legacy scheduler metadata are never imputed.
