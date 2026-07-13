# Curated Software Versions

Versions were collected from the remote HPC server on 2026-07-07. Commands that triggered login-node resource warnings or produced no reliable version string were not forced further on the login node.

| Tool | Version/status | Path or note |
|---|---|---|
| Host | workstation | remote HPC login node |
| Kernel | Linux 4.18.0-80.7.1.el8_0.x86_64 | uname output |
| Conda | 25.11.1 | `/home/u2024110181382/miniconda3/bin/conda` |
| Python | 3.12.5 | `/home/u2024110181382/miniconda3/bin/python` |
| gapseq | 2.0.0 | `/home/u2024110181382/Project/02.E.coli/biosoft/gapseq/gapseq -v` |
| Rscript | R 4.4.3 | gapseq conda environment |
| DIAMOND | 2.1.14 | `/home/u2024110181382/.conda/envs/carveme_env_py310/bin/diamond` |
| BLAST+ | 2.5.0+ | base conda `blastp` and `makeblastdb` |
| FastANI | 1.34 | `fastani_env` conda environment |
| Prokka | 1.12 | `prokka_env` conda environment |
| CheckM2 | 1.1.0 | `checkm2` conda environment |
| EzClermont | 0.7.0 | `ez` conda environment |
| cctyper | 1.5.0 | `cctyper` conda environment |
| minced | version unresolved on login node | Java VM heap allocation failed during `--version` check |
| SignalP | environment present; executable path unresolved | `/home/u2024110181382/miniconda3/envs/signalp` |
| TMHMM | 2.0c | `/home/u2024110181382/tmhmm-2.0c/bin/tmhmm` |
| memote | 0.17.0 | `memote-conda` environment |
| COBRApy | 0.30.0 | `memote-conda` environment |
| python-libsbml | 5.20.5 | `memote-conda` environment |
| Slurm | 22.05.6 | `sbatch`, `srun`, `salloc` |
| Custom Enterobacterales fastdb | bytes=10224943; sha256=777b2e293a435ccacb0930a8ce8dbab61835c06290251889d73c60cf286a8e8c | custom gapseq FASTA database |
| M9 glucose medium | bytes=361; sha256=c44887c9c46ccec0bfaec0774e44bbaa832800845e748f901e0d73baf613bfac | defined gap-filling medium |
