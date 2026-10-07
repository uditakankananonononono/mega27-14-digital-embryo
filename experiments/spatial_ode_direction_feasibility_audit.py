"""Nonnegative endpoint geometry of a linear maximizing direction, not dynamics."""
import json,hashlib
from pathlib import Path
import numpy as np
from scipy.linalg import expm,svd
ROOT=Path(__file__).resolve().parents[1]
def compute():
 pth=ROOT/'experiments/spatial_ode_direction_feasibility_plan.json';p=json.loads(pth.read_text());hashes={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in p['inputs']};station=json.loads((ROOT/p['inputs'][2]).read_text());jac=json.loads((ROOT/p['inputs'][3]).read_text())
 cells=station['residuals'][-1]['species'];x=np.array([r['concentration'] for r in cells]);ids=jac['species_ids'];assert [r['id'] for r in cells]==ids;assert len(x)==len(ids) and np.all(x>=0);rows=[]
 for f in p['inputs'][:2]:
  J=np.loadtxt(ROOT/f,delimiter=',');u,s,vt=svd(expm(10*J));v=vt[0];signs=[]
  for sign in [1,-1]:
   z=sign*v;idx=np.where(z<0)[0];limits=x[idx]/(-z[idx]);i=idx[np.argmin(limits)];alpha=float(limits.min());signs.append({'sign':sign,'negative_components':int(len(idx)),'maximum_nonnegative_amplitude':alpha,'amplitude_over_endpoint_norm':alpha/float(np.linalg.norm(x)),'limiting_species':ids[i],'limiting_concentration':float(x[i]),'limiting_direction_component':float(z[i]),'minimum_final_concentration':float(np.min(x+alpha*z))})
  rows.append({'matrix':f,'time':10,'maximum_linear_norm':float(s[0]),'unit_right_singular_vector':v.tolist(),'signs':signs})
 return {'plan_sha256':hashlib.sha256(pth.read_bytes()).hexdigest(),'input_sha256':hashes,'endpoint_norm':float(np.linalg.norm(x)),'minimum_endpoint_concentration':float(np.min(x)),'linearizations':rows,'limits':p['limits']}
if __name__=='__main__':
 j=compute();(ROOT/'results/spatial_ode_direction_feasibility_audit.json').write_text(json.dumps(j,indent=2)+'\n');print({k:v for k,v in j.items() if k!='linearizations'});print([{k:v for k,v in r.items() if k!='unit_right_singular_vector'} for r in j['linearizations']])
