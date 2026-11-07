"""Audit previously unresolved V4/V5 reaction losses without rerunning models."""

import collections
import csv
import json
import pathlib
import xml.etree.ElementTree as ET


ROOT = pathlib.Path(r"D:\pipeline\analysis_v5")
PKG = ROOT / "19.mSystems_revision4_20261028"
RUNS = ROOT / "11.benchmark_20x6_provenance" / "remote_downloaded_20x6_runs"
ORIGINAL = RUNS / "runs"
SANITIZED = RUNS / "runs_corrected_sanitized"
OUT = PKG / "01_analysis" / "reaction_loss_resolution_20261101"
OUT.mkdir(parents=True, exist_ok=True)


def read_csv(path, delimiter=","):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


stages = read_csv(PKG / "01_analysis" / "V1_V2_and_V4_V5_reaction_stage_audit.csv")
evidence = read_csv(PKG / "01_analysis" / "V4_V5_draft_loss_reaction_evidence.csv")
seed_rows = read_csv(ROOT / "gapseq" / "gapseq" / "dat" / "seed_reactions_corrected.tsv", "\t")
seed = {row["id"]: row for row in seed_rows}


def reaction_hash(row):
    parts = []
    for item in row["stoichiometry"].split(";"):
        fields = item.split(":")
        if len(fields) < 3:
            return "ERR"
        parts.append((f"{fields[1]}[{fields[2]}]", int(float(fields[0]))))
    parts.sort()
    rev = row["reversibility"]
    if rev == "=" and parts[0][1] < 0 or rev == "<":
        parts = [(met, -coef) for met, coef in parts]
    return ";".join(f"{coef}:{met}" for met, coef in parts) + "#" + rev


hash_of = {rid: reaction_hash(row) for rid, row in seed.items() if row["gapseq.status"] in {"approved", "corrected"}}
same_hash_ids = collections.defaultdict(list)
for rid, value in hash_of.items():
    same_hash_ids[value].append(rid)


def reaction_ids(path):
    ids = set()
    for _, node in ET.iterparse(path, events=("end",)):
        if node.tag.rsplit("}", 1)[-1] == "reaction":
            rid = node.attrib.get("id", "")
            if rid.startswith("R_rxn"):
                ids.add(rid[2:10])
            node.clear()
    return ids


def model_sets(folder):
    draft_paths = list(folder.glob("*draft.xml"))
    final_paths = [p for p in folder.glob("*.xml") if not p.name.endswith("-draft.xml")]
    if len(draft_paths) != 1 or len(final_paths) != 1:
        raise RuntimeError(f"Expected one draft/final XML in {folder}: {draft_paths}, {final_paths}")
    return reaction_ids(draft_paths[0]), reaction_ids(final_paths[0])


draft_v5 = {}
for accession in sorted({r["Accession"] for r in evidence if r["evidence_class"] == "V5_qualifying_row_but_not_draft"}):
    draft_v5[accession] = model_sets(ORIGINAL / accession / "V5")[0]

qualified_rows = []
for row in evidence:
    if row["evidence_class"] != "V5_qualifying_row_but_not_draft":
        continue
    rid = row["Seed_Reaction_ID"]
    h = hash_of.get(rid)
    represented = sorted((set(same_hash_ids.get(h, [])) - {rid}) & draft_v5[row["Accession"]]) if h and h != "ERR" else []
    qualified_rows.append({
        "Accession": row["Accession"], "Seed_Reaction_ID": rid,
        "Dictionary_Status": seed.get(rid, {}).get("gapseq.status", "missing"),
        "V5_Qualifying_Rows_In_Raw_Table": row["V5_candidate_rows"],
        "V5_Post_Complex_Rows_In_Raw_Table": row["V5_post_complex_rows"],
        "Hash_Valid": h not in (None, "ERR"),
        "Approved_Equivalent_Seed_IDs": ";".join(sorted(set(same_hash_ids.get(h, [])) - {rid})) if h else "",
        "Equivalent_Present_In_V5_Draft": ";".join(represented),
        "Audit_Class": "same_stoichiometry_represented_in_V5_draft" if represented else "requires_post_raw_table_trace",
    })


def write_rows(name, rows):
    with (OUT / name).open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


write_rows("eligible_921_stoichiometry_audit.csv", qualified_rows)

missing = [row for row in evidence if row["evidence_class"] == "evidence_table_unavailable_or_malformed"]
sanitized_sets = {}
for accession in sorted({r["Accession"] for r in missing}):
    sanitized_sets[accession] = {
        ver: model_sets(SANITIZED / accession / (ver + "s"))
        for ver in ("V4", "V5")
    }
sanitized_rows = []


def evidence_by_seed(folder, wanted):
    paths = list(folder.glob("*all-Reactions.tbl"))
    assert len(paths) == 1, (folder, paths)
    rows = collections.defaultdict(list)
    with paths[0].open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t")
        assert {"dbhit", "bitscore", "status", "pathway.status", "complex.status"}.issubset(reader.fieldnames), paths[0]
        for row in reader:
            for rid in row.get("dbhit", "").split():
                if rid in wanted:
                    rows[rid].append(row)
    return rows


sanitized_evidence = {}
for accession in sanitized_sets:
    wanted = {r["Seed_Reaction_ID"] for r in missing if r["Accession"] == accession}
    sanitized_evidence[accession] = {
        ver: evidence_by_seed(SANITIZED / accession / (ver + "s"), wanted)
        for ver in ("V4", "V5")
    }


def evidence_summary(rows):
    hit = [r for r in rows if r.get("bitscore") not in (None, "", "NA")]
    eligible = [r for r in rows if
                (r.get("bitscore") not in (None, "", "NA") and float(r["bitscore"]) >= 200 and r.get("status") != "bad_blast")
                or r.get("pathway.status") in {"full", "threshold", "keyenzyme"}]
    post_complex = [r for r in eligible if r.get("complex.status") != "0"]
    best = max(hit, key=lambda r: float(r["bitscore"])) if hit else None
    return {"rows": len(rows), "hits": len(hit), "eligible": len(eligible), "post_complex": len(post_complex),
            "best_bitscore": best["bitscore"] if best else "", "best_pident": best["pident"] if best else "",
            "best_query": best["qseqid"] if best else "", "best_target": best["stitle"].split()[0] if best else ""}


for row in missing:
    accession, rid = row["Accession"], row["Seed_Reaction_ID"]
    v4_d, v4_f = sanitized_sets[accession]["V4"]
    v5_d, v5_f = sanitized_sets[accession]["V5"]
    v4_e = evidence_summary(sanitized_evidence[accession]["V4"].get(rid, []))
    v5_e = evidence_summary(sanitized_evidence[accession]["V5"].get(rid, []))
    if not v5_e["hits"]:
        evidence_class = "no_V5_alignment_in_sanitized_rerun"
    elif not v5_e["eligible"]:
        evidence_class = "V5_alignment_below_first_eligibility_in_sanitized_rerun"
    elif not v5_e["post_complex"]:
        evidence_class = "V5_complex_filter_in_sanitized_rerun"
    else:
        evidence_class = "V5_first_eligible_but_not_draft_in_sanitized_rerun"
    sanitized_rows.append({
        "Accession": accession, "Seed_Reaction_ID": rid,
        "Original_Evidence_Table_Class": row["evidence_class"],
        "Sanitized_V4_Draft": rid in v4_d,
        "Sanitized_V5_Draft": rid in v5_d,
        "Sanitized_V4_Final": rid in v4_f,
        "Sanitized_V5_Final": rid in v5_f,
        "Original_Draft_Loss_Reproduced": rid in v4_d and rid not in v5_d,
        "Sanitized_Evidence_Class": evidence_class,
        **{f"Sanitized_V4_{k}": v for k, v in v4_e.items()},
        **{f"Sanitized_V5_{k}": v for k, v in v5_e.items()},
    })
write_rows("malformed_506_sanitized_rerun_crosscheck.csv", sanitized_rows)

final_stage = [r for r in stages if r["Comparison"] == "V4_vs_V5" and r["stage_class"] != "draft_prediction_or_evidence_stage_loss" and r["A_final"] == "True" and r["B_final"] == "False"]
write_rows("V4_V5_final_stage_1032.csv", final_stage)

summary = {
    "eligible_921": len(qualified_rows),
    "eligible_921_dictionary_status": dict(collections.Counter(r["Dictionary_Status"] for r in qualified_rows)),
    "eligible_921_same_stoichiometry_represented": sum(bool(r["Equivalent_Present_In_V5_Draft"]) for r in qualified_rows),
    "eligible_921_requires_post_raw_table_trace": sum(r["Audit_Class"] == "requires_post_raw_table_trace" for r in qualified_rows),
    "malformed_506": len(sanitized_rows),
    "malformed_506_accessions": dict(collections.Counter(r["Accession"] for r in sanitized_rows)),
    "malformed_506_original_loss_reproduced_in_sanitized": sum(r["Original_Draft_Loss_Reproduced"] for r in sanitized_rows),
    "malformed_506_sanitized_evidence_class": dict(collections.Counter(r["Sanitized_Evidence_Class"] for r in sanitized_rows)),
    "V4_V5_final_stage": len(final_stage),
    "V4_V5_final_stage_exchange": sum(r["is_exchange"] == "True" for r in final_stage),
    "V4_V5_final_stage_evidence_class": dict(collections.Counter(r["Evidence_Class"] for r in final_stage)),
}
(OUT / "reaction_loss_first_pass_summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
