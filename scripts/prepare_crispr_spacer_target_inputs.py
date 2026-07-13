from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(r"D:\pipeline\analysis_v5")
SPACERS = ROOT / "00.data" / "all_genomes_crispr_spacers_detailed.csv"
META = ROOT / "00.data" / "genome_assembly_and_cog_metrics.csv"
OUT = ROOT / "03.analysis" / "crispr_target_matching_inputs"


def wrap(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    spacers = pd.read_csv(SPACERS)
    meta = pd.read_csv(META, usecols=["Genome_ID", "Phylogroup"])

    spacers["Sequence"] = (
        spacers["Sequence"]
        .astype(str)
        .str.upper()
        .str.replace(r"[^ACGTN]", "", regex=True)
    )
    spacers = spacers[spacers["Sequence"].str.len() > 0].copy()
    spacers = spacers.merge(meta, on="Genome_ID", how="left")

    per_occurrence = spacers.copy()
    per_occurrence["Occurrence_ID"] = [
        f"spacer_occ_{i:06d}" for i in range(1, len(per_occurrence) + 1)
    ]
    per_occurrence = per_occurrence[
        ["Occurrence_ID", "Genome_ID", "Phylogroup", "Array_File", "Spacer_ID", "Sequence"]
    ]
    per_occurrence.to_csv(OUT / "crispr_spacer_occurrences_for_target_matching.csv", index=False)

    unique = (
        spacers.groupby("Sequence")
        .agg(
            n_occurrences=("Sequence", "size"),
            n_genomes=("Genome_ID", "nunique"),
            n_phylogroups=("Phylogroup", "nunique"),
            genomes=("Genome_ID", lambda x: ";".join(sorted(set(map(str, x))))),
            phylogroups=("Phylogroup", lambda x: ";".join(sorted(set(map(str, x.dropna()))))),
            example_spacer_id=("Spacer_ID", "first"),
        )
        .reset_index()
        .sort_values(["n_genomes", "n_occurrences", "Sequence"], ascending=[False, False, True])
    )
    unique["Unique_Spacer_ID"] = [f"uniq_spacer_{i:06d}" for i in range(1, len(unique) + 1)]
    unique = unique[
        [
            "Unique_Spacer_ID",
            "Sequence",
            "n_occurrences",
            "n_genomes",
            "n_phylogroups",
            "phylogroups",
            "genomes",
            "example_spacer_id",
        ]
    ]
    unique.to_csv(OUT / "unique_crispr_spacers_for_target_matching.csv", index=False)

    with (OUT / "unique_crispr_spacers.fasta").open("w", encoding="ascii") as f:
        for row in unique.itertuples(index=False):
            header = (
                f">{row.Unique_Spacer_ID}|n_occurrences={row.n_occurrences}|"
                f"n_genomes={row.n_genomes}|phylogroups={row.phylogroups}"
            )
            f.write(header + "\n")
            f.write(wrap(row.Sequence) + "\n")

    with (OUT / "crispr_target_matching_README.txt").open("w", encoding="utf-8") as f:
        f.write(
            "CRISPR spacer target-matching input package\n"
            "===========================================\n\n"
            "Generated from 00.data/all_genomes_crispr_spacers_detailed.csv.\n\n"
            f"Spacer occurrences: {len(per_occurrence)}\n"
            f"Unique spacer sequences: {len(unique)}\n"
            f"Genomes represented in spacer table: {spacers['Genome_ID'].nunique()}\n\n"
            "Files:\n"
            "- unique_crispr_spacers.fasta: unique spacer query FASTA for BLAST/minimap2.\n"
            "- unique_crispr_spacers_for_target_matching.csv: sequence-level metadata.\n"
            "- crispr_spacer_occurrences_for_target_matching.csv: occurrence-level metadata.\n\n"
            "True target matching still requires an external/local target database, such as:\n"
            "- NCBI RefSeq Viral nucleotide FASTA\n"
            "- IMG/VR viral genomes\n"
            "- PHASTER/prophage sequence database\n"
            "- A curated E. coli phage/plasmid/prophage FASTA used consistently across the study\n\n"
            "Recommended command once target_db.fna is available:\n"
            "makeblastdb -in target_db.fna -dbtype nucl -parse_seqids -out target_db\n"
            "blastn -task blastn-short -query unique_crispr_spacers.fasta -db target_db "
            "-out spacer_vs_target.tsv -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' "
            "-perc_identity 90 -qcov_hsp_perc 80 -max_target_seqs 10 -num_threads 8\n\n"
            "Interpretation note: without a viral/prophage/plasmid target FASTA, spacer-level sharing can be analyzed, "
            "but phage/protospacer target identity cannot be inferred scientifically.\n"
        )

    print(f"Wrote {OUT}")
    print(f"Spacer occurrences: {len(per_occurrence)}")
    print(f"Unique spacers: {len(unique)}")
    print(f"Genomes represented: {spacers['Genome_ID'].nunique()}")
    print(f"Phylogroup missing rows: {spacers['Phylogroup'].isna().sum()}")


if __name__ == "__main__":
    main()
