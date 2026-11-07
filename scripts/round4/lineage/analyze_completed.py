import csv, pathlib, collections, statistics, json, xml.etree.ElementTree as ET
import numpy as np
from scipy.spatial.distance import pdist
from scipy.stats import wilcoxon
P=pathlib.Path('/home/u2024110181382/Project/04.pipeline')
R=P/'revision4_finalcohort_broad_focused_20261030'; O=R/'summary'
sample=list(csv.DictReader((R/'inputs/sample_180.csv').open()))
status=list(csv.DictReader((O/'run_status.tsv').open(),delimiter='\t'))
assert len(sample)==180 and len(status)==360
assert all(x['Status']=='SUCCESS' for x in status)
def reactions(path):
 ids=set(); gene=set(); context=ET.iterparse(path,events=('end',))
 for _,e in context:
  if e.tag.rsplit('}',1)[-1]=='reaction':
   rid=e.attrib.get('id')
   if rid:
    ids.add(rid)
    if any(z.tag.rsplit('}',1)[-1]=='geneProductAssociation' for z in e.iter()):gene.add(rid)
   e.clear()
 return ids,gene
rows=[];sets={'broad':[],'focused':[]};gsets={'broad':[],'focused':[]}
for s in sample:
 id=s['Genome_ID']; pair={}
 for c in ('broad','focused'):
  p=R/'models'/id/c/(id+'.xml'); assert p.is_file() and p.stat().st_size>1000000,p
  rx,g=reactions(p);pair[c]=(rx,g);sets[c].append(rx);gsets[c].append(g)
 b,bg=pair['broad'];f,fg=pair['focused']
 rows.append({'Genome_ID':id,'Phylogroup':s['Phylogroup'],'Broad_Reactions':len(b),'Focused_Reactions':len(f),'Broad_Only':len(b-f),'Focused_Only':len(f-b),'Shared':len(b&f),'Jaccard':len(b&f)/len(b|f),'Broad_Gene_Evidenced':len(bg),'Focused_Gene_Evidenced':len(fg),'Broad_Gene_Only':len(bg-fg),'Focused_Gene_Only':len(fg-bg)})
with (O/'paired_reaction_counts_180.tsv').open('w',newline='') as h:
 w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def summary(vals):
 a=np.asarray(vals,dtype=float);return {'median':float(np.median(a)),'mean':float(a.mean()),'q25':float(np.quantile(a,.25)),'q75':float(np.quantile(a,.75)),'min':float(a.min()),'max':float(a.max())}
by=collections.defaultdict(list)
for r in rows:by[r['Phylogroup']].append(r)
reports={g:{'n':len(rs),**{k:summary([float(r[k]) for r in rs]) for k in ('Broad_Reactions','Focused_Reactions','Broad_Only','Focused_Only','Jaccard','Broad_Gene_Only')}} for g,rs in by.items()}
reports['ALL']={'n':len(rows),**{k:summary([float(r[k]) for r in rows]) for k in ('Broad_Reactions','Focused_Reactions','Broad_Only','Focused_Only','Jaccard','Broad_Gene_Only')}}
reports['ALL']['wilcoxon_broad_vs_focused_reaction_count_p']=float(wilcoxon([r['Broad_Reactions'] for r in rows],[r['Focused_Reactions'] for r in rows]).pvalue)
(O/'paired_reaction_summary.json').write_text(json.dumps(reports,indent=2))
# Model-set Jaccard distances and a group-label PERMANOVA, using the same sampled IDs for both configurations.
def distvec(ss):
 universe=sorted(set().union(*ss));idx={x:i for i,x in enumerate(universe)}
 a=np.zeros((len(ss),len(universe)),dtype=bool)
 for j,s in enumerate(ss):a[j,[idx[x] for x in s]]=1
 return pdist(a,metric='jaccard')
def permanova(v,labels,nperm=999,seed=42):
 n=len(labels);D=np.zeros((n,n));D[np.triu_indices(n,1)]=v;D=D+D.T
 labs=np.asarray(labels);uniq=np.unique(labs);total=(D**2).sum()/(2*n)
 def calc(x):
  within=sum((D[np.ix_(np.where(x==g)[0],np.where(x==g)[0])]**2).sum()/(2*sum(x==g)) for g in uniq)
  between=total-within;return (between/(len(uniq)-1))/(within/(n-len(uniq))),between/total
 f,r2=calc(labs);rng=np.random.default_rng(seed);ge=sum(calc(rng.permutation(labs))[0]>=f for _ in range(nperm))
 return {'F':float(f),'R2':float(r2),'p_permutation':(ge+1)/(nperm+1),'permutations':nperm}
labels=[r['Phylogroup'] for r in rows]
with (P/'final_cohort_1000/analysis/core_genome_ani_distance_matrix.csv').open() as h:
 reader=csv.reader(h);head=next(reader);am={r[0]:dict(zip(head[1:],map(float,r[1:]))) for r in reader}
ani=np.array([am[a['Genome_ID']][b['Genome_ID']] for i,a in enumerate(sample) for b in sample[i+1:]])
network={}
distance_vectors={}
for c in ('broad','focused'):
 network[c]={}
 for typ,ss in [('all',sets[c]),('gene_evidenced',gsets[c])]:
  v=distvec(ss);distance_vectors[c+'_'+typ]=v;network[c][typ]={'permanova':permanova(v,labels),'ANI_distance_Pearson_r':float(np.corrcoef(v,ani)[0,1]),'median_Jaccard_distance':float(np.median(v))}
np.savez_compressed(O/'paired_network_distances_180.npz',ani=ani,**distance_vectors)
(O/'sampled_lineage_network_summary.json').write_text(json.dumps(network,indent=2))
# MEMOTE annotated final models: all 1000 and phylogroup-specific distributions.
mem=list(csv.DictReader((P/'model_10.29/memote_1000_output/memote_all_models.csv').open(encoding='utf-8-sig')))
master=list(csv.DictReader((P/'revision_20260905_data_freeze/frozen_1000_master_manifest.csv').open(encoding='utf-8-sig')))
group={r['Genome_ID']:r['Final_Phylogroup'] for r in master}
mr=[]
for r in mem:
 id=r['model_id'];score=r['total_score'];mr.append({'Genome_ID':id,'Phylogroup':group.get(id,''),'Exit_Code':r['exit_code'],'Total_Score':score,'Report_Exists':pathlib.Path(r['report_html']).is_file()})
with (O/'memote_1000_scores_by_phylogroup.tsv').open('w',newline='') as h:
 w=csv.DictWriter(h,fieldnames=list(mr[0]),delimiter='\t');w.writeheader();w.writerows(mr)
assert len(mr)==1000 and len(set(r['Genome_ID'] for r in mr))==1000
ms={}
for g in sorted(set(group.values())):
 v=[float(r['Total_Score']) for r in mr if r['Phylogroup']==g and r['Exit_Code']=='0' and r['Total_Score']]
 ms[g]={'n_expected':sum(x==g for x in group.values()),'n_success':len(v),'scores':summary(v) if v else None}
ms['ALL']={'n_expected':1000,'n_success':sum(r['Exit_Code']=='0' and bool(r['Total_Score']) for r in mr),'n_reports':sum(r['Report_Exists'] for r in mr),'scores':summary([float(r['Total_Score']) for r in mr if r['Exit_Code']=='0' and r['Total_Score']])}
(O/'memote_1000_summary.json').write_text(json.dumps(ms,indent=2))
print('PAIR_ALL',reports['ALL']);print('LINEAGE',network);print('MEMOTE_ALL',ms['ALL'],flush=True)
