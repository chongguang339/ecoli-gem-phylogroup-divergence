suppressMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: expand_alignments_by_map.R <alignments_in> <map_tsv> <alignments_out>")
}

align_in <- args[1]
map_tsv <- args[2]
align_out <- args[3]

if (!file.exists(align_in)) {
  stop("Alignments file not found: ", align_in)
}
if (!file.exists(map_tsv)) {
  file.copy(align_in, align_out, overwrite = TRUE)
  quit(save = "no")
}

if (file.size(align_in) == 0) {
  file.copy(align_in, align_out, overwrite = TRUE)
  quit(save = "no")
}

alignments <- fread(align_in, header = FALSE, sep = "\t", quote = "")
if (nrow(alignments) == 0) {
  file.copy(align_in, align_out, overwrite = TRUE)
  quit(save = "no")
}

map_dt <- fread(map_tsv, header = FALSE, sep = "\t")
setnames(map_dt, c("rep_id", "member_id"))

# Join on qseqid (V1) and expand for all members.
# Note: in data.table's i-based join X[Y], the join key from Y (V1) is dropped from
# the result and X's columns (rep_id, member_id) are prepended. We reconstruct V1
# from member_id and reorder columns so V1 is back in position 1.
other_cols <- setdiff(names(alignments), "V1")  # V2, V3, ... in original order
expanded <- map_dt[alignments, on = .(rep_id = V1), allow.cartesian = TRUE]
if (nrow(expanded) == 0) {
  file.copy(align_in, align_out, overwrite = TRUE)
  quit(save = "no")
}

expanded[, V1 := member_id]
expanded[, c("rep_id", "member_id") := NULL]
# Restore original column order (V1 first, then V2, V3, ...)
setcolorder(expanded, c("V1", other_cols))

fwrite(expanded, file = align_out, sep = "\t", quote = FALSE, col.names = FALSE)
