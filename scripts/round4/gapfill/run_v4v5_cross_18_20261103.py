"""Gate 18 original-evidence genomes on exact V4/V4 and V5/V5 replay, then cross gap-fill."""
import concurrent.futures as cf
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time
import xml.etree.ElementTree as ET

PROJECT = Path('/home/u2024110181382/Project/04.pipeline')
BASE = PROJECT / 'revision2_20x6_benchmark'
ROOT = PROJECT / 'revision5_v4v5_cross_18_20261103'
RSCRIPT = Path('/home/u2024110181382/miniconda3/envs/gapseq/bin/Rscript')
FILL = Path('/home/u2024110181382/Project/02.E.coli/biosoft/gapseq/src/gf.suite.R')
ROOT.mkdir(parents=True, exist_ok=True)
for sub in ('inputs', 'runs', 'summary', 'logs'):
    (ROOT / sub).mkdir(exist_ok=True)

with (BASE / 'inputs/benchmark_20_genomes.csv').open(newline='') as f:
    genomes = list(csv.DictReader(f))
with (BASE / 'inputs/benchmark_failed_2_genomes_sanitized.csv').open(newline='') as f:
    sanitized = {r['Accession'] for r in csv.DictReader(f)}
original = [r['Accession'] for r in genomes if r['Accession'] not in sanitized]
assert len(genomes) == 20 and len(original) == 18 and len(sanitized) == 2
(ROOT / 'inputs/cohort.json').write_text(json.dumps({'original_evidence_18': original,
    'sanitized_evidence_2_excluded_from_cross': sorted(sanitized)}, indent=2))


def only(folder, pattern):
    paths = list(folder.glob(pattern))
    if len(paths) != 1 or not paths[0].is_file():
        raise ValueError(f'Expected one {pattern} in {folder}: {paths}')
    return paths[0]


def final_xml(folder):
    paths = [p for p in folder.glob('*.xml') if not p.name.endswith('-draft.xml')]
    if len(paths) != 1:
        raise ValueError(f'Expected one final XML in {folder}: {paths}')
    return paths[0]


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def rxns(path):
    ids = set()
    for _, node in ET.iterparse(path, events=('end',)):
        if node.tag.rsplit('}', 1)[-1] == 'reaction':
            ids.add(node.attrib['id'])
        node.clear()
    return ids


def run(pair):
    acc, draft, medium = pair
    ddir, mdir = BASE / 'runs' / acc / draft, BASE / 'runs' / acc / medium
    inputs = [only(ddir, '*-draft.RDS'), only(ddir, '*-rxnWeights.RDS'),
              only(ddir, '*-rxnXgenes.RDS'), only(mdir, '*-medium.csv')]
    archived = final_xml(BASE / 'runs' / acc / draft) if draft == medium else None
    out = ROOT / 'runs' / acc / f'draft{draft}_medium{medium}'
    out.mkdir(parents=True, exist_ok=True)
    status_file = out / 'status.json'
    if status_file.exists():
        old = json.loads(status_file.read_text())
        if old.get('Exit_Code') == 0 and old.get('Replay_XML') and Path(old['Replay_XML']).is_file():
            return old
    cmd = [str(RSCRIPT), str(FILL), '-m', str(inputs[0]), '-n', str(inputs[3]),
           '-c', str(inputs[1]), '-g', str(inputs[2]), '-b', '100', '-f', str(out)]
    metadata = {'Accession': acc, 'Draft': draft, 'Medium': medium,
                'Inputs': {str(p): sha(p) for p in inputs}, 'Command': cmd}
    (out / 'input_command.json').write_text(json.dumps(metadata, indent=2))
    env = os.environ.copy()
    env['PATH'] = str(RSCRIPT.parent) + ':' + env['PATH']
    env['OPENBLAS_NUM_THREADS'] = '4'
    start = time.monotonic()
    with (out / 'run.log').open('w') as log:
        proc = subprocess.run(['/usr/bin/time', '-v', *cmd], stdout=log, stderr=subprocess.STDOUT, env=env)
    xml = list(out.glob('*.xml'))
    result = {'Accession': acc, 'Draft': draft, 'Medium': medium,
              'Exit_Code': proc.returncode, 'Wall_Seconds': round(time.monotonic() - start, 3),
              'Replay_XML': str(xml[0]) if len(xml) == 1 else '',
              'Archived_XML': str(archived) if archived else '',
              'Baseline_Byte_Identical': sha(xml[0]) == sha(archived) if archived and len(xml) == 1 else '',
              'Baseline_Reaction_Identical': rxns(xml[0]) == rxns(archived) if archived and len(xml) == 1 else ''}
    status_file.write_text(json.dumps(result, indent=2))
    return result


def write_csv(path, rows):
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


baseline_tasks = [(acc, version, version) for acc in original for version in ('V4', 'V5')]
with cf.ThreadPoolExecutor(max_workers=12) as pool:
    baselines = list(pool.map(run, baseline_tasks))
write_csv(ROOT / 'summary/baseline_replay.csv', baselines)
failed = [x for x in baselines if x['Exit_Code'] != 0 or x['Baseline_Byte_Identical'] is not True]
if failed:
    (ROOT / 'summary/baseline_gate_failure.json').write_text(json.dumps(failed, indent=2))
    raise SystemExit(f'Baseline gate failed for {len(failed)} / {len(baselines)} runs; cross runs not started')
print('BASELINE GATE PASSED: 36/36 byte-identical', flush=True)

cross_tasks = [(acc, 'V4', 'V5') for acc in original] + [(acc, 'V5', 'V4') for acc in original]
with cf.ThreadPoolExecutor(max_workers=12) as pool:
    crosses = list(pool.map(run, cross_tasks))
write_csv(ROOT / 'summary/cross_run_status.csv', crosses)
failed = [x for x in crosses if x['Exit_Code'] != 0 or not x['Replay_XML']]
if failed:
    (ROOT / 'summary/cross_run_failure.json').write_text(json.dumps(failed, indent=2))
    raise SystemExit(f'Cross run failed for {len(failed)} / {len(crosses)} runs')

by_key = {(r['Accession'], r['Draft'], r['Medium']): r for r in baselines + crosses}
details = []
for acc in original:
    v4draft = rxns(only(BASE / 'runs' / acc / 'V4', '*-draft.xml'))
    sets = {(d, m): rxns(Path(by_key[(acc, d, m)]['Replay_XML']))
            for d in ('V4', 'V5') for m in ('V4', 'V5')}
    target = (sets['V4', 'V4'] - sets['V5', 'V5']) - v4draft
    for reaction in sorted(target):
        v4v5 = reaction in sets['V4', 'V5']
        v5v4 = reaction in sets['V5', 'V4']
        if v4v5 and not v5v4:
            label = 'draft_sensitive'
        elif v5v4 and not v4v5:
            label = 'medium_sensitive'
        else:
            label = 'interaction_or_alternative_solution'
        details.append({'Accession': acc, 'Reaction': reaction, 'V4V4': 1,
                        'V4V5': int(v4v5), 'V5V4': int(v5v4), 'V5V5': 0,
                        'Classification': label})
write_csv(ROOT / 'summary/finalization_only_reaction_cross_classification.csv', details)
print(f'CROSS COMPLETE: 18 genomes, {len(details)} V4-only finalization-stage records', flush=True)
