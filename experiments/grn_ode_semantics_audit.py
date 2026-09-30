"""Defined finite-horizon Hill model audit, not equivalence to Boolean dynamics."""
import hashlib,itertools,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from digitalembryo.segment_polarity import successor

def rhs(x,K,n,a):
 h=x**n/(K**n+x**n)
 wg,en,hh,ptc,ci=[h[...,i] for i in range(5)]
 return np.stack([np.minimum(1,ci+a*wg),wg,en,ci*(1-hh),np.zeros_like(wg)],axis=-1)-x

def integrate(x,K,n,a,dt,horizon=100):
 cur=x.copy()
 for _ in range(round(horizon/dt)):
  k1=rhs(cur,K,n,a);k2=rhs(cur+dt*k1/2,K,n,a);k3=rhs(cur+dt*k2/2,K,n,a);k4=rhs(cur+dt*k3,K,n,a)
  cur+=dt*(k1+2*k2+2*k3+k4)/6
 return cur

def main():
 raw=(ROOT/'experiments/grn_ode_semantics_plan.json').read_bytes();p=json.loads(raw);x=np.array(list(itertools.product(p['initial_grid'],repeat=5)),float)
 rows=[]
 for K,n,a in itertools.product(p['Ks'],p['ns'],p['strengths']):
  finals=[integrate(x,K,n,a,dt,p['horizon']) for dt in p['steps']]
  for dt,xf in zip(p['steps'],finals):
   rows.append({'K':K,'n':n,'a':a,'dt':dt,'initial_count':len(x),'terminal_wg_on_fraction':float(np.mean(xf[:,0]>.5)),'max_rhs_residual':float(np.max(np.abs(rhs(xf,K,n,a)))),'n_distinct_rounded_6':len(np.unique(xf.round(6),axis=0)),'min_state_coordinate':float(xf.min()),'max_state_coordinate':float(xf.max()),'dt_halving_max_endpoint_difference':float(np.max(np.abs(finals[0]-finals[1])))})
 out={'plan_sha256':hashlib.sha256(raw).hexdigest(),'boolean_wg_on_state':[1,1,1,0,0],'boolean_successor':list(successor((1,1,1,0,0))),'ode_rhs_at_same_state_K05_n4_a0':rhs(np.array([1,1,1,0,0.]),.5,4,0).tolist(),'rows':rows,'limits':p['limits']}
 (ROOT/'results/grn_ode_semantics_audit.json').write_text(json.dumps(out,indent=2)+'\n')
 print('Rows',len(rows),'max solver difference',max(r['dt_halving_max_endpoint_difference'] for r in rows),'max residual',max(r['max_rhs_residual'] for r in rows))
 for K in p['Ks']:
  print(K,[(n,a,next(r['terminal_wg_on_fraction'] for r in rows if r['K']==K and r['n']==n and r['a']==a and r['dt']==.05)) for n in p['ns'] for a in [.85,.9,1]])
if __name__=='__main__':main()
