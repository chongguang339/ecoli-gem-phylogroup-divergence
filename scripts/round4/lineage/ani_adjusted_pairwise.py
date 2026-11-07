import csv,json,pathlib,numpy as np
from scipy.spatial.distance import squareform
R=pathlib.Path('/home/u2024110181382/Project/04.pipeline/revision4_finalcohort_broad_focused_20261030');O=R/'summary'
a=np.load(O/'paired_network_distances_180.npz'); sample=list(csv.DictReader((R/'inputs/sample_180.csv').open()));labels=np.array([x['Phylogroup'] for x in sample]);n=len(labels)
ani=squareform(a['ani']);H=np.eye(n)-np.ones((n,n))/n;B=-.5*H@(ani**2)@H;val,vec=np.linalg.eigh(B);ix=np.argsort(val)[::-1];ok=ix[val[ix]>1e-12];pcs=vec[:,ok]*np.sqrt(val[ok]);assert pcs.shape[1]>=20
pairs=np.tril_indices(n,-1);X=np.column_stack([np.ones(len(pairs[0]))]+[np.abs(pcs[pairs[0],i]-pcs[pairs[1],i]) for i in range(20)])
res={}
for name in ('broad_all','focused_all','broad_gene_evidenced','focused_gene_evidenced'):
 y=squareform(a[name])[pairs];beta=np.linalg.lstsq(X,y,rcond=None)[0];e=y-X@beta;mean=e.mean();sst=((e-mean)**2).sum()
 def est(lbl):
  same=lbl[pairs[0]]==lbl[pairs[1]];u=e[same].mean();v=e[~same].mean();return (same.sum()*(u-mean)**2+(~same).sum()*(v-mean)**2)/sst,u-v
 obs,effect=est(labels);rng=np.random.default_rng(42);ext=sum(est(rng.permutation(labels))[0]>=obs for _ in range(999))
 res[name]={'pairwise_residual_same_group_R2':float(obs),'same_group_distance_effect':float(effect),'label_permutation_p':(ext+1)/1000,'ANI_PCoA_axes':20}
(O/'ani_adjusted_pairwise_sensitivity.json').write_text(json.dumps(res,indent=2));print(json.dumps(res,indent=2),flush=True)
