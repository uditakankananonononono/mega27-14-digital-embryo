"""Fixed-time, fixed-coordinate linear amplification diagnostic."""
import json,hashlib
from pathlib import Path
import numpy as np
from scipy.linalg import expm,svdvals,eigvalsh
ROOT=Path(__file__).resolve().parents[1]
def compute():
 pth=ROOT/'experiments/spatial_ode_nonnormal_plan.json';plan=json.loads(pth.read_text());hashes={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in plan['inputs']};old=json.loads((ROOT/plan['inputs'][2]).read_text());rows=[]
 for f in plan['inputs'][:2]:
  J=np.loadtxt(ROOT/f,delimiter=',');assert J.shape==(len(old['species_ids']),)*2 and np.isfinite(J).all()
  rows.append({'matrix':f,'spectral_abscissa':float(np.linalg.eigvals(J).real.max()),'euclidean_logarithmic_norm':float(eigvalsh((J+J.T)/2)[-1]),'propagator_norms':[{'time':t,'operator_2_norm':float(svdvals(expm(t*J))[0])} for t in [0,1,10,100]]})
 return {'plan_sha256':hashlib.sha256(pth.read_bytes()).hexdigest(),'input_sha256':hashes,'native_jacobian_finite':old['native_jacobian_finite'],'species_order':old['species_ids'],'linearizations':rows,'largest_fixed_grid_norm':max(r['operator_2_norm'] for c in rows for r in c['propagator_norms']),'scope':plan['scope']}
if __name__=='__main__':
 j=compute();(ROOT/'results/spatial_ode_nonnormal_audit.json').write_text(json.dumps(j,indent=2)+'\n');print({k:v for k,v in j.items() if k!='species_order'})
