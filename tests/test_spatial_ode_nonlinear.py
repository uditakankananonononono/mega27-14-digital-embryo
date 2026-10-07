import json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def test_frozen_paired_nonlinear_scope():
 j=json.loads((ROOT/'results/spatial_ode_nonlinear_audit.json').read_text());assert j['plan_sha256']==hashlib.sha256((ROOT/'experiments/spatial_ode_nonlinear_plan.json').read_bytes()).hexdigest()
 assert [r['fractional_increase'] for r in j['perturbations']]==[0.001,0.01]
 for r in j['perturbations']:
  assert [z['time'] for z in r['times']]==[0,1,10,100]
  assert all(z['minimum_concentration']>=0 for z in r['times'])
  assert r['times'][-1]['observed_gain']<1 and r['times'][-1]['linear_relative_error']<.01
