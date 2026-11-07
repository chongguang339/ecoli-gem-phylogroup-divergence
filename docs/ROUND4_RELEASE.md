# Round 4 mSystems revision: reproducibility index

The revised analysis uses the final 1,000-model cohort and five additional run collections. The associated large SBML models, MEMOTE reports, evidence tables and logs belong in a **new version of the existing Zenodo record**, not in Git history. The version-specific Zenodo DOI must be inserted here after publication. The all-versions DOI cited in the manuscript is `10.5281/zenodo.21341419`; verify the live record and its files before using it as a download claim.

| Question | Code here | Zenodo upload package |
|---|---|---|
| MEMOTE of final 1,000 models | `scripts/round4/memote/` | `round4_memote_1000_reports_evidence.tar.gz`; authoritative annotated SBML archive retained from the prior record version |
| Broad versus focused library, 180 paired genomes | `scripts/round4/lineage/` | `round4_broad_focused_180_models_evidence.tar.gz` and `round4_sample180_fastas.tar.gz` |
| Matched `tblastn` versus DIAMOND, 20 genomes | `scripts/round4/benchmark/` | `round4_matched_aligner_20_evidence.tar.gz` and `round4_matched_aligner_input_fastas.zip` |
| V4/V5 draft by medium cross, 18 genomes | `scripts/round4/gapfill/` | `round4_v4v5_cross_18_evidence.tar.gz` |
| Mechanism pilot and reaction traces | `scripts/round4/gapfill/` | `round4_mechanism_controls.tar.gz`, `round4_source_tables_and_reaction_provenance.zip` |
| Corrected KEGG feature filtering | `scripts/round4/kegg/` | `round4_source_tables_and_reaction_provenance.zip` |
| Cohort and reaction provenance | `scripts/round4/provenance/` | `data/round4_manifests/` and source-tables ZIP |

The scripts under `scripts/round4/` are the original analysis scripts, preserving the absolute paths and cluster environment used for the reported runs. To rerun elsewhere, edit the paths to the downloaded data, gapseq installation and reference libraries. The run packages preserve command records and input hashes where they were generated. Historical V1 program/library bytes and V5 fast/rescue intermediate alignments were not recovered; these packages do not claim to reconstruct those unavailable historical states.

The 20 accession IDs in the historical matched-aligner benchmark do not overlap the final 1,000-model manifest. The 180-genome broad/focused sample is a subset of the final 1,000. Treat these as separate controls.

The `gapseq_assets/` focused FASTA and medium file in the earlier repository release do not have the same SHA-256 hashes as the inputs recorded by the new 180-genome run. Checksum-identical local copies are staged as `gapseq_assets/enterobacterales_gapseq_fastdb_round4.faa` (SHA-256 `777b2e293a435ccacb0930a8ce8dbab61835c06290251889d73c60cf286a8e8c`) and `gapseq_assets/M9_glucose_medium_round4.csv` (SHA-256 `c44887c9c46ccec0bfaec0774e44bbaa832800845e748f901e0d73baf613bfac`). Use these for the Round 4 configuration; keep the earlier files for historical runs. The exact nucleotide inputs for both control collections are staged in their named Zenodo packages. The installed gapseq entry script and broad-reference archive were also recovered and checksum-matched locally (SHA-256 `500b01106bb76b91f41e73b2257c773e9613222f24137f350cc6f79872233936` and `6de7305c937b4a0bbe3a5435a62c95eba2828c5ac56fc4cf1f96c3551644f81c`), but their upstream source/version and redistribution terms need review before public deposition.

The frozen checksum audit, 1,000-row raw COG category scan and its source script are staged for this release. The original 1,000 eggNOG `.emapper.annotations` have now been SHA-256-verified and staged as `round4_eggnog_annotations_1000.tar.gz` for the new Zenodo version. The exact 180 reconstruction FASTAs are likewise staged as `round4_sample180_fastas.tar.gz`. These files are **not publicly accessible until the Zenodo version is published**.

The pre-existing `run_benchmark_modle.py` wrapper is historical. Use the direct SBML reaction-count tables and matched-run summaries for the current manuscript's benchmark values.
