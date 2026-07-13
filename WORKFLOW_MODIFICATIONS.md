# gapseq Workflow Modifications and Provenance

This study used gapseq as the base genome-scale metabolic reconstruction framework, but the production workflow was modified at the wrapper, database, medium, and source-script levels to support cohort-scale reconstruction and standardized comparisons.

## Base Reconstruction Framework

- Base tool: gapseq 2.0.0
- Optimized executable used in production: a locally installed and modified gapseq checkout.
- Public reproducibility package: this repository contains the wrapper script, modified source-script snapshot, custom medium file, and database-construction script needed to document the optimized workflow.

## Optimized gapseq Command

For each genome assembly, the optimized workflow used the following command pattern:

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

The corresponding benchmark wrapper is provided as:

```text
scripts/run_benchmark_modle.py
```

## Parameter Rationale

- `-A diamond`: uses DIAMOND 2.1.14 instead of the default BLAST-based search to reduce protein-alignment runtime.
- `-S identical`: uses identical sequence clustering/search mode for the custom reference database.
- `-E <CUSTOM_FASTDB_FASTA>`: supplies a compact Enterobacterales-focused gapseq-compatible protein database.
- `-m <M9_GLUCOSE_MEDIUM_CSV>`: constrains gap-filling with a defined M9 glucose medium.
- `-f <OUTPUT_DIRECTORY>`: writes all model outputs into sample-specific directories.
- `-K 4`: fixes per-run thread usage for benchmark comparability and to reduce CPU oversubscription on shared compute nodes.

## Custom Enterobacterales gapseq-compatible FASTA Database

The optimized workflow used a compact Enterobacterales-focused database:

```text
gapseq_assets/enterobacterales_gapseq_fastdb.faa
```

The database was generated from the gapseq bacterial sequence database using the following taxonomic filter:

```text
Escherichia|Shigella|Salmonella|Klebsiella|Enterobacter|Citrobacter
```

The public reconstruction script is provided as:

```text
gapseq_assets/build_gapseq_compatible_ecoli_fastdb.sh
```

Use it with project-relative paths:

```bash
bash gapseq_assets/build_gapseq_compatible_ecoli_fastdb.sh \
  -g <GAPSEQ_BACTERIA_SEQUENCE_DIR> \
  -o gapseq_assets/enterobacterales_gapseq_fastdb.faa \
  -r 'Escherichia|Shigella|Salmonella|Klebsiella|Enterobacter|Citrobacter'
```

## M9 Glucose Medium

The gap-filling medium is provided as:

```text
gapseq_assets/M9_glucose_medium.csv
```

The medium contains water, major ions, trace metals, ammonium, oxygen, and D-glucose as the defined carbon source.

## Modified gapseq Source-Script Snapshot

The source-script snapshot used for the optimized workflow is provided in:

```text
gapseq_src_modified/
```

Key modified or workflow-relevant files include:

- `doall.sh`
- `gapseq_find.sh`
- `transporter.sh`
- `gf.suite.R`
- `generate_GSdraft.R`

Backup files retained from development are included only to document provenance. For reproduction, use the active files without `.bak` suffixes.

## HPC Execution

Production runs were executed on an HPC server after obtaining compute resources through the local scheduler. A typical usage pattern is:

```bash
salloc <scheduler-options>
# or
sbatch <batch-script>
# or
srun <resource-options> python scripts/run_benchmark_modle.py
```

Cluster-specific absolute paths are intentionally excluded from the public workflow. Users should set project-local variables for input FASTA files, the gapseq executable, the custom FASTA database, the M9 medium, and output directories.

## Reproducibility Boundary

The deposited processed tables and SBML/XML models are the definitive outputs used in the manuscript. Re-running the complete workflow can reproduce the model-generation logic, but minor runtime differences are expected across CPU architecture, scheduler configuration, filesystem performance, and DIAMOND/gapseq minor versions.
