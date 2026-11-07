import csv
import json
import pathlib
import xml.etree.ElementTree as ET

import numpy as np
from scipy.spatial.distance import squareform

ROOT = pathlib.Path('/home/u2024110181382/Project/04.pipeline/revision4_finalcohort_broad_focused_20261030')
OUT = ROOT / 'summary'
sample = list(csv.DictReader((ROOT / 'inputs/sample_180.csv').open()))
labels = np.array([r['Phylogroup'] for r in sample])
n = len(sample)
z = np.load(OUT / 'paired_network_distances_180.npz')
ani = squareform(z['ani'])
H = np.eye(n) - np.ones((n, n)) / n
value, vector = np.linalg.eigh(-0.5 * H @ (ani ** 2) @ H)
order = np.argsort(value)[::-1]
positive = order[value[order] > 1e-12]
pcs = vector[:, positive] * np.sqrt(value[positive])
pair = np.tril_indices(n, -1)
X = np.column_stack([np.ones(len(pair[0]))] + [np.abs(pcs[pair[0], k] - pcs[pair[1], k]) for k in range(20)])
Q, _ = np.linalg.qr(X, mode='reduced')
same = labels[pair[0]] == labels[pair[1]]

def model_ids(path):
    ids = set()
    for _, node in ET.iterparse(path, events=('end',)):
        if node.tag.rsplit('}', 1)[-1] == 'reaction':
            rid = node.attrib.get('id')
            if rid:
                ids.add(rid)
            node.clear()
    return ids

broad, focused = [], []
for row in sample:
    base = ROOT / 'models' / row['Genome_ID']
    broad.append(model_ids(base / 'broad' / (row['Genome_ID'] + '.xml')))
    focused.append(model_ids(base / 'focused' / (row['Genome_ID'] + '.xml')))

universe = sorted(set().union(*broad))
lookup = {rid: j for j, rid in enumerate(universe)}
matrix = np.zeros((n, len(universe)), dtype=np.int16)
for i, ids in enumerate(broad):
    matrix[i, [lookup[r] for r in ids]] = 1
length = matrix.sum(axis=1)
intersection = matrix @ matrix.T
union = length[:, None] + length[None, :] - intersection
inter = intersection[pair].astype(float)
uni = union[pair].astype(float)

def r2(y):
    e = y - Q @ (Q.T @ y)
    mean = e.mean()
    sst = ((e - mean) ** 2).sum()
    return float((same.sum() * (e[same].mean() - mean) ** 2 + (~same).sum() * (e[~same].mean() - mean) ** 2) / sst)

baseline = r2(1 - inter / uni)
assert abs(baseline - 0.015246516745528035) < 1e-10
loss_count = {}
for b, f in zip(broad, focused):
    for rid in b - f:
        loss_count[rid] = loss_count.get(rid, 0) + 1
records = []
for rid, losses in loss_count.items():
    if losses < 5:
        continue
    p = matrix[:, lookup[rid]].astype(bool)
    present = int(p.sum())
    if present == n:
        continue
    a, b = p[pair[0]], p[pair[1]]
    y = 1 - (inter - (a & b)) / (uni - (a | b))
    new_r2 = r2(y)
    rec = {'Reaction_ID': rid, 'Broad_Only_Genomes': losses, 'Broad_Present_Genomes': present,
           'R2_Without_Reaction': new_r2, 'R2_Reduction': baseline - new_r2}
    for g in sorted(set(labels)):
        rec['Broad_Present_' + g] = int(p[labels == g].sum())
    records.append(rec)
records.sort(key=lambda x: x['R2_Reduction'], reverse=True)
path = OUT / 'broad_only_reaction_ablation_ranked.tsv'
with path.open('w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter='\t')
    writer.writeheader()
    writer.writerows(records)
info = {'baseline_R2': baseline, 'tested_variable_broad_only_reactions': len(records),
        'screen': 'broad-only in at least five genomes and broad-model presence varies among 180',
        'top_20': records[:20]}
(OUT / 'broad_only_reaction_ablation_ranked_summary.json').write_text(json.dumps(info, indent=2))
print(json.dumps(info, indent=2))
