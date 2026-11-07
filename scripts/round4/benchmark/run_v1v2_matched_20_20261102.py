"""Run 20 matched-reference aligner pairs; two sanitized FASTAs are flagged."""

import concurrent.futures
import csv
import hashlib
import json
import os
import pathlib
import re
import subprocess
import time
import xml.etree.ElementTree as ET

PROJECT = pathlib.Path('/home/u2024110181382/Project/04.pipeline')
ROOT = PROJECT / 'revision5_v1v2_matched_20_20261102'
ORIGINAL = PROJECT / 'revision2_20x6_benchmark/inputs/benchmark_20_genomes.csv'
SANITIZED = PROJECT / 'revision2_20x6_benchmark/inputs/benchmark_failed_2_genomes_sanitized.csv'
RUNS = ROOT / 'runs'
EXE = pathlib.Path('/home/u2024110181382/Project/02.E.coli/biosoft/gapseq/gapseq')
DB = pathlib.Path('/home/u2024110181382/miniconda3/envs/gapseq/share/gapseq/dat/seq')
PILOT = PROJECT / 'revision5_mechanism_controls_20261102/01_v1v2'
ENV_BIN = '/home/u2024110181382/miniconda3/envs/gapseq/bin'
PILOT_ACCESSION = 'GCF_055805355.1'

for folder in (ROOT, RUNS, ROOT / 'inputs', ROOT / 'summary'):
    folder.mkdir(parents=True, exist_ok=True)

with ORIGINAL.open(newline='') as handle:
    genomes = list(csv.DictReader(handle))
with SANITIZED.open(newline='') as handle:
    replacements = {r['Accession']: r['FNA_Path'] for r in csv.DictReader(handle)}
assert len(genomes) == 20 and len(replacements) == 2
for row in genomes:
    row['Input_Class'] = 'sanitized_ambiguity_codes' if row['Accession'] in replacements else 'original'
    row['FNA_Path'] = replacements.get(row['Accession'], row['FNA_Path'])
    assert pathlib.Path(row['FNA_Path']).is_file()

with (ROOT / 'inputs/cohort_manifest.csv').open('w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=list(genomes[0]))
    writer.writeheader()
    writer.writerows(genomes)


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


with (ROOT / 'inputs/shared_sha256.txt').open('w') as handle:
    for path in (EXE, DB / 'Bacteria/rev/sequences.tar.gz'):
        handle.write(f'{digest(path)}  {path}\n')


def model_paths(out):
    drafts = list(out.glob('*-draft.xml'))
    finals = [p for p in out.glob('*.xml') if not p.name.endswith('-draft.xml')]
    return drafts, finals


def elapsed_from_log(path):
    match = re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*([\d:.]+)', path.read_text(errors='replace'))
    if not match:
        return None
    parts = [float(x) for x in match.group(1).split(':')]
    return sum(value * 60 ** power for power, value in enumerate(reversed(parts)))


def run_one(row):
    accession, fna = row['Accession'], pathlib.Path(row['FNA_Path'])
    records = []
    for aligner in ('tblastn', 'diamond'):
        if accession == PILOT_ACCESSION:
            out = PILOT / f'{accession}_modified_{aligner}_v1db_16core'
            status_path = out / 'matched20_status.json'
            status = {'Accession': accession, 'Input_Class': row['Input_Class'], 'FNA_SHA256': digest(fna),
                      'Aligner': aligner, 'Output_Directory': str(out), 'Reused_Completed_Pilot': True,
                      'Exit_Code': 0, 'Wall_Seconds': elapsed_from_log(out / 'run.log')}
            log_text = (out / 'run.log').read_text(errors='replace')
            assert str(fna) in log_text, f'Pilot input differs from cohort manifest: {accession}'
        else:
            out = RUNS / accession / aligner
            out.mkdir(parents=True, exist_ok=True)
            status_path = out / 'matched20_status.json'
            if status_path.exists():
                prior = json.loads(status_path.read_text())
                drafts, finals = model_paths(out)
                if prior.get('Valid_Output') and len(drafts) == 1 and len(finals) == 1:
                    records.append(prior)
                    continue
            cmd = [str(EXE), 'doall', '-A', aligner, '-D', str(DB), '-O', '-K', '16', '-f', str(out), str(fna)]
            env = os.environ.copy()
            env['PATH'] = ENV_BIN + ':' + env['PATH']
            env['OPENBLAS_NUM_THREADS'] = '4'
            env['GAPSEQ_CACHE_DIR'] = str(out / 'cache')
            with (out / 'command.json').open('w') as handle:
                json.dump({'Command': cmd, 'FNA_SHA256': digest(fna), 'Input_Class': row['Input_Class']}, handle, indent=2)
            start = time.monotonic()
            with (out / 'run.log').open('w') as handle:
                proc = subprocess.run(['/usr/bin/time', '-v', *cmd], stdout=handle, stderr=subprocess.STDOUT, env=env)
            status = {'Accession': accession, 'Input_Class': row['Input_Class'], 'FNA_SHA256': digest(fna),
                      'Aligner': aligner, 'Output_Directory': str(out), 'Reused_Completed_Pilot': False,
                      'Exit_Code': proc.returncode, 'Wall_Seconds': time.monotonic() - start}
        drafts, finals = model_paths(out)
        status['Valid_Output'] = status['Exit_Code'] == 0 and len(drafts) == 1 and len(finals) == 1
        status['Draft_XML'] = str(drafts[0]) if len(drafts) == 1 else ''
        status['Final_XML'] = str(finals[0]) if len(finals) == 1 else ''
        status_path.write_text(json.dumps(status, indent=2))
        records.append(status)
    return records


all_status = []
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
    for records in executor.map(run_one, genomes):
        all_status.extend(records)
        print('FINISHED', records[0]['Accession'], [(r['Aligner'], r['Valid_Output']) for r in records], flush=True)

with (ROOT / 'summary/run_status.csv').open('w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=list(dict.fromkeys(k for row in all_status for k in row)))
    writer.writeheader()
    writer.writerows(all_status)


def reactions(path):
    ids = set()
    for _, node in ET.iterparse(path, events=('end',)):
        if node.tag.rsplit('}', 1)[-1] == 'reaction':
            ids.add(node.attrib['id'])
        node.clear()
    return ids


pairs = []
for row in genomes:
    accession = row['Accession']
    by_align = {r['Aligner']: r for r in all_status if r['Accession'] == accession}
    if not all(r['Valid_Output'] for r in by_align.values()):
        continue
    for stage, key in (('draft', 'Draft_XML'), ('final', 'Final_XML')):
        a, b = reactions(by_align['tblastn'][key]), reactions(by_align['diamond'][key])
        pairs.append({'Accession': accession, 'Input_Class': row['Input_Class'], 'Stage': stage,
                      'Tblastn_Count': len(a), 'Diamond_Count': len(b),
                      'Tblastn_Only': len(a-b), 'Diamond_Only': len(b-a),
                      'Jaccard': len(a & b) / len(a | b)})
with (ROOT / 'summary/paired_reaction_counts.csv').open('w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=list(pairs[0]))
    writer.writeheader()
    writer.writerows(pairs)
print('VALID PAIRS', len(pairs) // 2, 'OF', len(genomes), flush=True)
