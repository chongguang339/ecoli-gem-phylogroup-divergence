import csv,random,pathlib,collections,hashlib
root=pathlib.Path('/home/u2024110181382/Project/04.pipeline')
out=root/'revision4_finalcohort_broad_focused_20261030'
freeze=root/'revision_20260905_data_freeze'
master=list(csv.DictReader((freeze/'frozen_1000_master_manifest.csv').open(encoding='utf-8-sig')))
audit=list(csv.DictReader((freeze/'frozen_1000_file_checksum_audit.csv').open(encoding='utf-8-sig')))
fasta={r['Genome_ID']:r for r in audit if r['File_Type']=='FASTA'}
ani=root/'final_cohort_1000/analysis/core_genome_ani_distance_matrix.csv'
with ani.open(newline='') as h:
 rd=csv.reader(h); header=next(rd)[1:]; dist={row[0]:dict(zip(header,map(float,row[1:]))) for row in rd}
quotas={'A':25,'B1':25,'B2':25,'C':20,'D':25,'E':25,'F':20,'G':11,'U':4}
by=collections.defaultdict(list)
for r in master: by[r['Final_Phylogroup']].append(r['Genome_ID'])
rng=random.Random(42); selected=[]
for group,quota in quotas.items():
 ids=sorted(by[group]); assert len(ids)>=quota,(group,len(ids))
 assert all(i in dist for i in ids),(group,'ANI row missing')
 medoid=min(ids,key=lambda i:sum(dist[i][j] for j in ids)/len(ids))
 if len(ids)<=quota: picks=ids
 else:
  ranked=sorted(ids,key=lambda i:(dist[medoid][i],i)); bins=[ranked[k*len(ranked)//3:(k+1)*len(ranked)//3] for k in range(3)]
  wanted=[quota//3+(k<quota%3) for k in range(3)]
  picks=sorted(x for bucket,n in zip(bins,wanted) for x in rng.sample(bucket,n))
 for i in picks:
  assert i in fasta and pathlib.Path(fasta[i]['Path']).is_file(),i
  selected.append({'Genome_ID':i,'Phylogroup':group,'FASTA':fasta[i]['Path'],'Frozen_FASTA_SHA256':fasta[i]['SHA256'],'ANI_Medoid':medoid,'Distance_to_Medoid':dist[medoid][i]})
assert len(selected)==180 and len({r['Genome_ID'] for r in selected})==180
p=out/'inputs/sample_180.csv'
with p.open('w',newline='') as h:
 w=csv.DictWriter(h,fieldnames=list(selected[0]));w.writeheader();w.writerows(selected)
with (out/'inputs/sample_180_ids.txt').open('w') as h:
 h.write(''.join(r['Genome_ID']+'\n' for r in selected))
print('sample',len(selected),dict(collections.Counter(r['Phylogroup'] for r in selected)),flush=True)
