"""Join the frozen 10,951 EzClermont results, initial 1,000, and final 1,000."""
import csv
from collections import Counter
from pathlib import Path

root = Path(r"D:\pipeline\analysis_v5")
src = root / "09.mSystems_revision_data_freeze" / "00_source_snapshots"
out = root / "19.mSystems_revision4_20261028" / "01_analysis" / "initial_to_final_sampling_provenance_1000.csv"

def read(name):
    with (src / name).open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))

def key(name):
    return name.removesuffix(".fna")

pool = read("Phylogroup_Results.csv")
initial = read("Sampled_Phylogroup_Results.csv")
final = read("Selected_1000_Genomes.csv")
pool_map = {key(row["Genome"]): row for row in pool}
initial_map = {key(row["Genome"]): row for row in initial}
assert len(pool_map) == len(pool) == 10951
assert len(initial_map) == len(initial) == 1000
assert len(final) == 1000
assert len({key(row["Genome"]) for row in final}) == 1000

rows = []
for row in final:
    name = key(row["Genome"])
    original = pool_map[name]
    first = initial_map.get(name)
    rows.append({"Genome_ID": name, "Original_10951_EzClermont_Label": original["Phylogroup"],
                 "In_Initial_1000_Sample": bool(first),
                 "Initial_1000_Label": first["Phylogroup"] if first else "",
                 "Final_1000_Label": row["Phylogroup"],
                 "Label_Changed_From_Original": original["Phylogroup"] != row["Phylogroup"],
                 "Final_Selection_Status": "retained_from_initial_sample" if first else "added_after_initial_sample"})

assert sum(r["In_Initial_1000_Sample"] for r in rows) == 411
assert all(not r["Label_Changed_From_Original"] for r in rows)
with out.open("w", newline="", encoding="utf-8-sig") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
print(out)
print(Counter(r["Final_Selection_Status"] for r in rows))
