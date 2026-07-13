# E. coli Phylogroup-Resolved GEM Workflow

This repository candidate contains the code and workflow-provenance files used for the manuscript:

**High-Throughput GEM Workflow Optimization for Phylogroup-Resolved Functional Divergence Analysis in *Escherichia coli***

## Contents

- `scripts/`: Python wrappers, statistical analyses, FBA benchmarking, figure generation, and manuscript-update utilities.
- `scripts/run_benchmark_modle.py`: benchmark wrapper comparing default and optimized gapseq workflows.
- `gapseq_assets/`: custom gapseq resources used by the optimized workflow.
  - `M9_glucose_medium.csv`: defined medium used for M9 glucose-constrained gap-filling.
  - `enterobacterales_gapseq_fastdb.faa`: compact Enterobacterales-focused gapseq-compatible protein database.
  - `build_gapseq_compatible_ecoli_fastdb.sh`: database construction script using a documented taxonomic regex.
- `gapseq_src_modified/`: source-script snapshot of the modified gapseq workflow used during production runs.
- `WORKFLOW_MODIFICATIONS.md`: detailed description of the gapseq parameter, database, medium, and wrapper-level modifications.
- `environment.yml`: Conda environment for downstream Python/COBRApy/statistical analyses and available command-line dependencies.
- `requirements.txt`: pip-oriented Python dependency list.
- `LICENSE`: code license.

## Optimized gapseq Command

The optimized reconstruction workflow used the following command pattern:

```bash
gapseq doall \
  -A diamond \
  -S identical \
  -E <CUSTOM_FASTDB_FASTA> \
  -m <M9_GLUCOSE_MEDIUM_CSV> \
  -f <OUTPUT_DIRECTORY> \
  -K 4 \
  <GENOME_FASTA>
```

This command is implemented in `scripts/run_benchmark_modle.py` for the reference-model and benchmark comparisons. Production cohort runs used the same database, medium, and optimized gapseq source-script logic.

## Data and Models

Processed analysis tables are available in the Supplementary Data package. The complete set of 1,000 optimized GEM SBML/XML files is provided as **Data S10**:

- `exported_models.tar.gz`
- `ecoli_1000_optimized_gem_models_manifest.csv`
- `ecoli_1000_optimized_gem_models_archive_checksum.csv`

The model archive SHA256 checksum is:

```text
953f92a5d53dfe2f999d8d7ee99304d81c1e71c30e9dcb2dffff39e2bf84d7fc
```

## External Dependencies

Some production tools are not fully captured by `environment.yml` because they require independent installation, license acceptance, large databases, or local HPC configuration. These include:

- gapseq 2.0.0 and its R dependencies
- CRISPRCasTyper or equivalent CRISPR annotation workflow
- SignalP environment present on the production server; exact executable path/version unresolved during login-node audit
- TMHMM 2.0c
- eggNOG-mapper and its database
- CheckM2 1.1.0 and its database
- FastANI 1.34 inputs and NCBI datasets

The curated production-software audit is available in `software_versions_curated.csv` / `software_versions_curated.md` in the open-science package.

## Reproducibility Note

Production runs were executed on an HPC server after compute resources were allocated through the local scheduler. Cluster-specific absolute paths are intentionally excluded from this public package. Users should replace local input/output paths with project-relative paths while preserving the documented input table names, gapseq parameters, and output structure.

## Suggested Citation

Please cite the manuscript and the archived Zenodo/Figshare/BioModels record DOI after deposition.
