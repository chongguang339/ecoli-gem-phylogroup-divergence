# COG source audit and sanitized-input validation — 2026-11-07

## Confirmed from downloaded files

- `frozen_1000_file_checksum_audit.csv` has 10,000 records, including exactly 1,000 `EggNOG_Annotations` records. `eggnog_1000_raw_cog_audit.tsv` has exactly 1,000 rows with the same genome IDs, paths and recorded byte sizes.
- Each audit row marks its source file as present, has a positive annotation-row count and an empty error field. The audit totals 4,680,982 annotated protein rows (median 4,688 per genome); all reported COG_R, COG_X and COG_Y counts are zero. The accompanying JSON summary also reports 1,000 existing, size-matched source files and zero R/X counts.
- `revision4_verify_sources_20261029.py` reads the frozen manifest, finds the original `COG_category` column and writes the audit and summary in the format downloaded. It checks existence and byte size, **not current SHA-256**.
- A separate local streaming scan of the downloaded raw archive now independently re-read all 1,000 source files and all 4,680,982 protein rows. Its per-file SHA-256, annotation-row counts and R/X/Y counts matched the frozen manifest and prior audit for **all 1,000 genomes**; R, X and Y were each zero. See `INDEPENDENT_RAW_COG_RESCAN.json` and `verify_raw_cog_from_archive.py`.
- The two downloaded sanitized FASTAs match the SHA-256 in the 20-genome matched-run status table exactly: `GCA_055753585.1` = `c2c7efec59ae8f90d9f1fbe03c618cc667ac86f6a68c872b62aa46c3d91c725b`; `GCA_055766285.1` = `255437795fa5537b40125bb173174f0c98626859e696a173381be7a855ffc766`.
- `round4_source_tables_and_reaction_provenance.zip` now has 25 valid members, including the corrected submission versions of S1/S2, the COG audit files and independent rescan summary. `round4_matched_aligner_input_fastas.zip` has 23 valid members: 20 original FASTAs, two sanitized FASTAs and an audit table. All 20 actual matched-run FASTA hashes are represented.

## Remaining limitation

The original 1,000 `.emapper.annotations` and 180 exact FASTAs have now been downloaded, checksum-verified and staged for Zenodo. They are not externally accessible until the new Zenodo version is published. The upstream reason for absent R/X categories remains unknown; zero in the original annotation output does not imply biological absence.
