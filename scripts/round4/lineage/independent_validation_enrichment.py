import csv,hashlib,pathlib,collections,json
import libsbml
from scipy.stats import fisher_exact
P=pathlib.Path('/home/u2024110181382/Project/04.pipeline/revision4_finalcohort_broad_focused_20261030');O=P/'summary'
sample=list(csv.DictReader((P/'inputs/sample_180.csv').open()));old={r['Genome_ID']:r for r in csv.DictReader((O/'paired_reaction_counts_180.tsv').open(),delimiter='\t')}
checks=[]; unique=collections.defaultdict(set);names={}
for s in sample:
 id=s['Genome_ID'];h=hashlib.sha256()
 with open(s['FASTA'],'rb') as f:
  for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
 hash_ok=h.hexdigest()==s['Frozen_FASTA_SHA256']; pair={}
 for conf in ('broad','focused'):
  path=P/'models'/id/conf/(id+'.xml');doc=libsbml.readSBMLFromFile(str(path));model=doc.getModel()
  if model is None:raise RuntimeError('no SBML model '+str(path))
  ids={r.getId() for r in model.getListOfReactions()}; bad=doc.getNumErrors(libsbml.LIBSBML_SEV_ERROR)
  if len(ids)!=model.getNumReactions():raise RuntimeError('duplicate reaction IDs '+str(path))
  for r in model.getListOfReactions():
   if r.getId() and r.getName():names.setdefault(r.getId(),r.getName())
  pair[conf]=ids
  checks.append({'Genome_ID':id,'Phylogroup':s['Phylogroup'],'Config':conf,'FASTA_SHA256_Match':hash_ok,'SBML_Error_Count':bad,'SBML_Reactions':len(ids),'Prior_Count_Match':len(ids)==int(old[id][conf.title()+'_Reactions'])})
 unique[s['Phylogroup']].update(pair['broad']-pair['focused'])
 if not hash_ok:raise RuntimeError('FASTA hash mismatch '+id)
with (O/'independent_fasta_sbml_validation.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(checks[0]),delimiter='\t');w.writeheader();w.writerows(checks)
assert len(checks)==360 and all(r['FASTA_SHA256_Match'] and r['Prior_Count_Match'] for r in checks)
# Reuse independent libSBML reaction lists for a group-vs-rest prevalence audit.
loss={};labels={s['Genome_ID']:s['Phylogroup'] for s in sample}
for s in sample:
 id=s['Genome_ID']; b=libsbml.readSBMLFromFile(str(P/'models'/id/'broad'/(id+'.xml'))).getModel(); f=libsbml.readSBMLFromFile(str(P/'models'/id/'focused'/(id+'.xml'))).getModel()
 loss[id]={r.getId() for r in b.getListOfReactions()}-{r.getId() for r in f.getListOfReactions()}
allrx=sorted(set().union(*loss.values()));groups=sorted(set(labels.values()));tests=[]
for rx in allrx:
 for g in groups:
  yes=sum(rx in loss[i] for i in loss if labels[i]==g);no=sum(rx in loss[i] for i in loss if labels[i]!=g)
  ng=sum(labels[i]==g for i in loss);nrest=len(loss)-ng
  if yes+no<5:continue
  _,p=fisher_exact([[yes,ng-yes],[no,nrest-no]])
  tests.append({'Reaction_ID':rx,'Reaction_Name':names.get(rx,''),'Phylogroup':g,'Group_Loss':yes,'Group_N':ng,'Other_Loss':no,'Other_N':nrest,'Risk_Difference':yes/ng-no/nrest,'Fisher_p':p})
order=sorted(range(len(tests)),key=lambda i:tests[i]['Fisher_p']);m=len(order);q=[1.0]*m;best=1.0
for rank in range(m-1,-1,-1):
 i=order[rank];best=min(best,tests[i]['Fisher_p']*m/(rank+1));q[i]=best
for row,adj in zip(tests,q):row['BH_q']=adj
with (O/'broad_only_reaction_group_enrichment.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(tests[0]),delimiter='\t');w.writeheader();w.writerows(tests)
summary={'fasta_hash_match':sum(r['FASTA_SHA256_Match'] for r in checks),'sbml_count_match':sum(r['Prior_Count_Match'] for r in checks),'sbml_models':len(checks),'models_with_sbml_errors':sum(r['SBML_Error_Count']>0 for r in checks),'tested_reaction_group_pairs':len(tests),'BH_q_lt_0_05':sum(r['BH_q']<.05 for r in tests)}
(O/'independent_validation_enrichment_summary.json').write_text(json.dumps(summary,indent=2));print(summary,flush=True)
