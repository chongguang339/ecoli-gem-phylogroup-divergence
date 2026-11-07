"""Compare matched aligner evidence tables after the paired modeling job."""

import csv
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path('/home/u2024110181382/Project/04.pipeline/revision5_v1v2_matched_20_20261102')
STATUS = ROOT / 'summary/run_status.csv'
SUMMARY = ROOT / 'summary'

with STATUS.open(newline='') as handle:
    statuses = list(csv.DictReader(handle))
by_genome = defaultdict(dict)
for row in statuses:
    by_genome[row['Accession']][row['Aligner']] = row


def evidence(path):
    tables = list(Path(path).glob('*-all-Reactions.tbl'))
    if len(tables) != 1:
        raise ValueError(f'Expected exactly one reaction evidence table in {path}: {tables}')
    with tables[0].open(newline='') as handle:
        rows = csv.DictReader((line for line in handle if not line.startswith('#')), delimiter='\t')
        counts = Counter()
        by_rxn = defaultdict(Counter)
        for row in rows:
            rxn = row['rxn']
            signature = (row['pathway'], row['subunits'], row['complex'], row['dbhit'],
                         row['src'], row['qseqid'], row['status'])
            counts[signature] += 1
            by_rxn[rxn][signature] += 1
    return tables[0], counts, by_rxn


overview = []
detail = []
for accession, aligned in by_genome.items():
    if not all(x in aligned and aligned[x]['Valid_Output'] == 'True' for x in ('tblastn', 'diamond')):
        continue
    a_path, a_counts, a_rxn = evidence(aligned['tblastn']['Output_Directory'])
    b_path, b_counts, b_rxn = evidence(aligned['diamond']['Output_Directory'])
    all_rxns = a_rxn.keys() | b_rxn.keys()
    changed = [rxn for rxn in all_rxns if a_rxn[rxn] != b_rxn[rxn]]
    overview.append({
        'Accession': accession,
        'Input_Class': aligned['tblastn']['Input_Class'],
        'Tblastn_Table': str(a_path), 'Diamond_Table': str(b_path),
        'Tblastn_Rows': sum(a_counts.values()), 'Diamond_Rows': sum(b_counts.values()),
        'Tblastn_Unique_Reactions': len(a_rxn), 'Diamond_Unique_Reactions': len(b_rxn),
        'Tblastn_Only_Reactions': len(a_rxn.keys() - b_rxn.keys()),
        'Diamond_Only_Reactions': len(b_rxn.keys() - a_rxn.keys()),
        'Shared_Reactions_Changed_Evidence': sum(rxn in a_rxn and rxn in b_rxn for rxn in changed),
        'Tblastn_Only_Evidence_Rows': sum((a_counts - b_counts).values()),
        'Diamond_Only_Evidence_Rows': sum((b_counts - a_counts).values()),
    })
    for rxn in sorted(changed):
        detail.append({
            'Accession': accession, 'Input_Class': aligned['tblastn']['Input_Class'],
            'Reaction': rxn,
            'Tblastn_Evidence_Rows': sum(a_rxn[rxn].values()),
            'Diamond_Evidence_Rows': sum(b_rxn[rxn].values()),
            'Tblastn_Only_Row_Signatures': sum((a_rxn[rxn] - b_rxn[rxn]).values()),
            'Diamond_Only_Row_Signatures': sum((b_rxn[rxn] - a_rxn[rxn]).values()),
        })

for name, rows in [('paired_evidence_overview.csv', overview),
                   ('paired_evidence_by_reaction.csv', detail)]:
    with (SUMMARY / name).open('w', newline='') as handle:
        if rows:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
print(f'Evidence comparison complete: {len(overview)} valid genome pairs; {len(detail)} changed reaction records')
