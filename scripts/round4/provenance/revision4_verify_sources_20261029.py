import csv, os, pathlib, collections, json
root=pathlib.Path('/home/u2024110181382/Project/04.pipeline')
out=root/'revision4_source_inventory_20261029'
out.mkdir(exist_ok=True)
manifest=root/'revision_20260905_data_freeze'/'frozen_1000_file_checksum_audit.csv'
rows=list(csv.DictReader(manifest.open(newline='',encoding='utf-8-sig')))
ann=[r for r in rows if r['File_Type']=='EggNOG_Annotations']
assert len(ann)==1000,len(ann)
results=[]
for r in ann:
 p=pathlib.Path(r['Path']); counts=collections.Counter(); n=0; error=''
 if not p.is_file(): error='MISSING'
 else:
  try:
   idx=None
   with p.open(encoding='utf-8',errors='replace') as h:
    for line in h:
     if line.startswith('#query'):
      cols=line.lstrip('#').rstrip('\n').split('\t'); idx=cols.index('COG_category') if 'COG_category' in cols else None; break
    if idx is None: error='NO_COG_HEADER'
    else:
     for line in h:
      if not line or line.startswith('#'): continue
      x=line.rstrip('\n').split('\t'); n+=1
      if len(x)>idx and x[idx]!='-':
       for c in set(x[idx]):
        if c.isalpha() and c.isupper(): counts[c]+=1
  except Exception as e: error=type(e).__name__+':'+str(e)
 results.append(dict(Genome_ID=r['Genome_ID'],Path=str(p),Exists=p.is_file(),Bytes=p.stat().st_size if p.is_file() else 0,Expected_Bytes=r['Bytes'],Annotated_Rows=n,COG_R=counts['R'],COG_X=counts['X'],COG_Y=counts['Y'],Error=error))
with (out/'eggnog_1000_raw_cog_audit.tsv').open('w',newline='',encoding='utf-8') as h:
 w=csv.DictWriter(h,fieldnames=list(results[0]),delimiter='\t');w.writeheader();w.writerows(results)
summary={'records':len(results),'exists':sum(x['Exists'] for x in results),'bytes_match':sum(str(x['Bytes'])==x['Expected_Bytes'] for x in results),'nonzero_R_genomes':sum(x['COG_R']>0 for x in results),'nonzero_X_genomes':sum(x['COG_X']>0 for x in results),'total_R':sum(x['COG_R'] for x in results),'total_X':sum(x['COG_X'] for x in results),'errors':collections.Counter(x['Error'] for x in results)}
with (out/'eggnog_1000_raw_cog_summary.json').open('w') as h:json.dump(summary,h,indent=2)
print(summary,flush=True)
