"""Exact multi-bit basin-retention curves, two starting-state measures.
Synchronous published model only, not physical intervention or topology causality.
"""
import json,itertools,sys,hashlib
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).parent))
import flower_model as F
ROOT=Path(__file__).resolve().parents[1]
def main():
 n=12;states=np.array(list(itertools.product([0,1],repeat=n)),np.uint8);powers=1<<np.arange(n-1,-1,-1);succ=F.step_batch(states)@powers;assert len(states)==4096
 labels=[]
 for start in range(4096):
  seen={};cur=start
  while cur not in seen:seen[cur]=len(seen);cur=int(succ[cur])
  labels.append(tuple(sorted(s for s,i in seen.items() if i>=seen[cur])))
 unique=sorted(set(labels));lab=np.array([unique.index(a) for a in labels]);count=np.array([int(x).bit_count() for x in range(4096)]);ids=np.arange(4096);rows=[]
 for name,bits in F.PUBLISHED_FP.items():
  fixed=int(bits,2);assert int(succ[fixed])==fixed;group=lab[fixed];members=ids[lab==group];hit_fixed=lab==group;dist=count[ids^fixed]
  fp=[float(hit_fixed[dist==k].mean()) for k in range(13)]
  sums=np.zeros(13);totals=np.zeros(13)
  for start in members:
   dist=count[ids^start];sums+=np.bincount(dist,weights=hit_fixed,minlength=13);totals+=np.bincount(dist,minlength=13)
  basin=(sums/totals).tolist();assert abs((1-basin[1])-json.loads((ROOT/'results/flower_model.json').read_text())['labeled'][name]['fragility'])<.0001
  # Independent random bit flips, exact combinatorial weighted measure.
  import math
  qs=np.linspace(0,.5,501);fp_prob=np.array([sum(math.comb(n,k)*q**k*(1-q)**(n-k)*fp[k] for k in range(13)) for q in qs]);basin_prob=np.array([sum(math.comb(n,k)*q**k*(1-q)**(n-k)*basin[k] for k in range(13)) for q in qs])
  def crossing(a):
   ix=np.where(a<=.5)[0];return float(qs[ix[0]]) if len(ix) else None
  rows.append({'organ':name,'fixed_state':bits,'basin_size':len(members),'fixed_start_retention_by_exact_radius':fp,'uniform_basin_start_retention_by_exact_radius':basin,'fixed_start_first_half_retention_q':crossing(fp_prob),'basin_start_first_half_retention_q':crossing(basin_prob),'independent_flip_grid_q':[0,.01,.05,.1,.2,.5],'fixed_start_retention_on_grid':[float(sum(math.comb(n,k)*q**k*(1-q)**(n-k)*fp[k] for k in range(13))) for q in [0,.01,.05,.1,.2,.5]],'basin_start_retention_on_grid':[float(sum(math.comb(n,k)*q**k*(1-q)**(n-k)*basin[k] for k in range(13))) for q in [0,.01,.05,.1,.2,.5]]})
 out={'model':'12-node published floral threshold Boolean network, synchronous update','source':'https://pageperso.lis-lab.fr/~sylvain.sene/files/publi_pres/rgs18.pdf','n_states':4096,'n_attractors':len(unique),'measurement':'probability of same eventual attractor after exact k distinct bit flips, enumerated without sampling','rows':rows,'model_code_sha256':hashlib.sha256((ROOT/'experiments/flower_model.py').read_bytes()).hexdigest(),'limits':['Artificial uniform state-space measure, not distribution of biological embryonic states','Independent bit noise is not calibrated gene-expression noise','Perturbation occurs before relaxation, not continuously during development','Synchronous update; no causal inference that topology alone determines fragility','Half-retention crossings quantized at q=.001']}
 (ROOT/'results/flower_perturbation_radius.json').write_text(json.dumps(out,indent=2)+'\n')
 for r in rows:print(r['organ'],'k1 fixed/basin',r['fixed_start_retention_by_exact_radius'][1],r['uniform_basin_start_retention_by_exact_radius'][1],'q50',r['fixed_start_first_half_retention_q'],r['basin_start_first_half_retention_q'])
if __name__=='__main__':main()
