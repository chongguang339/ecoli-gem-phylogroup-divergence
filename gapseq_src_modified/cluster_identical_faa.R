suppressMessages(library(Biostrings))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: cluster_identical_faa.R <input_faa> <output_faa> <map_tsv>")
}

input_faa <- args[1]
output_faa <- args[2]
map_tsv <- args[3]

if (!file.exists(input_faa)) {
  stop("Input fasta not found: ", input_faa)
}

seqs <- readAAStringSet(input_faa)
ids <- names(seqs)
seq_strings <- as.character(seqs)

# BLAST/Diamond output truncates FASTA headers at the first whitespace (space or tab).
# Use only the first word of each header as the ID so the cluster map matches
# the qseqid column in alignments.tsv during the expansion step.
ids_short <- sub("[ \t].*", "", ids)

# Identify representatives by identical sequence string (first occurrence wins)
rep_idx <- !duplicated(seq_strings)
rep_ids_short <- ids_short[rep_idx]
rep_seqs <- seqs[rep_idx]

# Map each member to its representative (using short IDs to match BLAST/Diamond output)
first_idx <- match(seq_strings, seq_strings[rep_idx])
member_rep_ids <- rep_ids_short[first_idx]

# Write mapping (rep_id = short ID matching BLAST qseqid, member_id = short ID)
map_df <- data.frame(rep_id = member_rep_ids, member_id = ids_short, stringsAsFactors = FALSE)
write.table(map_df, file = map_tsv, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Write clustered fasta
writeXStringSet(rep_seqs, filepath = output_faa, format = "fasta")
