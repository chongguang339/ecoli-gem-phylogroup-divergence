import csv
import json
import pathlib
import xml.etree.ElementTree as ET

import numpy as np
from scipy.spatial.distance import pdist, squareform

ROOT = pathlib.Path('/home/u2024110181382/Project/04.pipeline/revision4_finalcohort_broad_focused_20261030')
OUT = ROOT / 'summary'
sample = list(csv.DictReader((ROOT / 'inputs/sample_180.csv').open()))
labels = np.array([row['Phylogroup'] for row in sample])
data = np.load(OUT / 'paired_network_distances_180.npz')
n = len(labels)
ani = squareform(data['ani'])
H = np.eye(n) - np.ones((n, n)) / n
value, vector = np.linalg.eigh(-0.5 * H @ (ani ** 2) @ H)
order = np.argsort(value)[::-1]
positive = order[value[order] > 1e-12]
pcs = vector[:, positive] * np.sqrt(value[positive])
pair = np.tril_indices(n, -1)
X = np.column_stack([np.ones(len(pair[0]))] + [np.abs(pcs[pair[0], k] - pcs[pair[1], k]) for k in range(20)])
same = labels[pair[0]] == labels[pair[1]]

def score(sets):
    universe = sorted(set().union(*sets))
    lookup = {key: j for j, key in enumerate(universe)}
    matrix = np.zeros((n, len(universe)), dtype=bool)
    for i, ids in enumerate(sets):
        matrix[i, [lookup[key] for key in ids]] = True
    y = squareform(pdist(matrix, metric='jaccard'))[pair]
    residual = y - X @ np.linalg.lstsq(X, y, rcond=None)[0]
    mean = residual.mean()
    sst = ((residual - mean) ** 2).sum()
    a, b = residual[same].mean(), residual[~same].mean()
    r2 = (same.sum() * (a - mean) ** 2 + (~same).sum() * (b - mean) ** 2) / sst
    return {'R2': float(r2), 'same_group_effect': float(a - b)}

sets = []
gene_sets = []
focused_sets = []
focused_gene_sets = []
for row in sample:
    for config, target_all, target_gene in [('broad', sets, gene_sets), ('focused', focused_sets, focused_gene_sets)]:
        path = ROOT / 'models' / row['Genome_ID'] / config / (row['Genome_ID'] + '.xml')
        all_ids, gene_ids = set(), set()
        for _, node in ET.iterparse(path, events=('end',)):
            if node.tag.rsplit('}', 1)[-1] == 'reaction':
                rid = node.attrib.get('id')
                if rid:
                    all_ids.add(rid)
                    if any(child.tag.rsplit('}', 1)[-1] == 'geneProductAssociation' for child in node.iter()):
                        gene_ids.add(rid)
                node.clear()
        target_all.append(all_ids)
        target_gene.append(gene_ids)

five = ['R_rxn00306_c0', 'R_rxn00520_c0', 'R_rxn00934_c0', 'R_rxn01732_c0', 'R_rxn01733_c0']
two = ['R_rxn00507_c0', 'R_rxn03455_c0']
scenarios = {'baseline': [], 'five_ligases': five, 'two_oxidoreductases': two, 'all_seven': five + two}
scenarios.update({rid: [rid] for rid in five + two})
result = {}
for label, excluded in scenarios.items():
    removed = set(excluded)
    result[label] = {'excluded_reactions': excluded, 'all': score([s - removed for s in sets]),
                     'gene_evidenced': score([s - removed for s in gene_sets]),
                     'presence_in_broad': {rid: sum(rid in s for s in sets) for rid in excluded}}
baseline = result['baseline']
result['paired_intersection_only'] = {
    'description': 'For each genome, retain only reaction IDs shared by its broad and focused model; this is a set intersection sensitivity analysis, not a reconstructed model.',
    'all': score([b & f for b, f in zip(sets, focused_sets)]),
    'gene_evidenced': score([b & f for b, f in zip(gene_sets, focused_gene_sets)]),
    'broad_only_count_median': float(np.median([len(b - f) for b, f in zip(sets, focused_sets)])),
}
for label, entry in result.items():
    for kind in ('all', 'gene_evidenced'):
        entry[kind]['delta_R2_from_baseline'] = entry[kind]['R2'] - baseline[kind]['R2']
        entry[kind]['share_of_broad_focused_difference_removed'] = -entry[kind]['delta_R2_from_baseline'] / (
            baseline[kind]['R2'] - (0.006272199985216288 if kind == 'all' else 0.007974426644142333))
path = OUT / 'broad_only_reaction_ablation.json'
path.write_text(json.dumps(result, indent=2))
print(path)
print(json.dumps(result, indent=2))
