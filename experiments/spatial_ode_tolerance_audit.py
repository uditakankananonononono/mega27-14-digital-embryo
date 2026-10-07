"""Frozen tighter-tolerance comparison; no retuning model or disturbances."""
import json,hashlib
from pathlib import Path
import spatial_ode_nonlinear_audit as N
ROOT=Path(__file__).resolve().parents[1]
def compute(upstream):
 p=ROOT/'experiments/spatial_ode_tolerance_plan.json';old=ROOT/'results/spatial_ode_nonlinear_audit.json';base=json.loads(old.read_text());tight=N.compute(upstream,rtol=1e-10,atol=1e-16);rows=[]
 for a,b in zip(base['perturbations'],tight['perturbations']):
  rows.append({'fractional_increase':a['fractional_increase'],'gain_max_abs_change':max(abs(x['observed_gain']-y['observed_gain']) for x,y in zip(a['times'],b['times'])),'linear_error_max_abs_change':max(abs(x['linear_relative_error']-y['linear_relative_error']) for x,y in zip(a['times'],b['times'])),'baseline_drift_max_abs_change':max(abs(x['baseline_drift_norm']-y['baseline_drift_norm']) for x,y in zip(a['times'],b['times']))})
 return {'plan_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'previous_result_sha256':hashlib.sha256(old.read_bytes()).hexdigest(),'tight_rtol':1e-10,'tight_atol':1e-16,'comparison':rows,'tight_results':tight,'limits':json.loads(p.read_text())['limits']}
if __name__=='__main__':
 import argparse
 a=argparse.ArgumentParser();a.add_argument('--upstream',required=True);v=a.parse_args();j=compute(v.upstream);(ROOT/'results/spatial_ode_tolerance_audit.json').write_text(json.dumps(j,indent=2)+'\n');print(j['comparison'])
