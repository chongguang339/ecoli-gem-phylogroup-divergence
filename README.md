# Revision 3 reproducibility update

This directory contains the code-facing additions for the mSystems resubmission of the E. coli phylogroup-resolved GEM study. The persistent data record is [Zenodo 10.5281/zenodo.21341419](https://doi.org/10.5281/zenodo.21341419).

## What changed

- Rebuilt Figure 1D from direct child reaction identifiers in authoritative paired SBML models.
- Restricted the BacDive MG1655 comparison to nine successfully mapped carbon sources; unresolved amygdalin is retained only in the mapping audit.
- Added KEGG feature-construction rules and the 264-to-115 ontology-filter audit.
- Added ANI PCoA sensitivity analyses for 5, 10, 20, and 40 positive-eigenvalue axes.
- Added the 48-component chemically defined complete-medium recipe, uptake bounds, and checksums.
- Added a 120-row six-configuration benchmark provenance table. Exact command, CPU, memory, version, checksum, and log fields are recorded for V2-V5. Legacy V1/V6 scheduler CPU and memory fields were not recorded and are explicitly marked unavailable.
- Regenerated all main and supplementary figures from authoritative tables.

## Layout

- `scripts/`: analysis, model-archive, numerical-audit, and final figure-generation scripts used for the revision. Internal manuscript-formatting utilities are intentionally excluded because they are not part of the scientific workflow.
- `data/`: compact reviewer-facing outputs, including ANI-axis sensitivity and benchmark provenance summaries.
- `figures/`: final vector PDF figures.
- `KEGG_264_feature_construction_algorithm.csv`: explicit raw-feature construction rules.
- `chemically_defined_complete_medium_recipe_and_bounds.csv`: medium composition and maximum uptake bounds.

Large assets, including the 1,000 gap-fill-annotated SBML models, raw benchmark logs and JSON provenance, external-validation matrices, and 18 MEMOTE HTML reports, are deposited in the Zenodo record rather than duplicated in GitHub.

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
