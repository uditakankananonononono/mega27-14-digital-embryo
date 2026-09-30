import json,sys,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
from grn_ode_semantics_audit import rhs,integrate,successor

def test_boolean_persistence_is_not_hill_equilibrium():
 x=(1,1,1,0,0);assert successor(x)==x
 assert rhs(np.array(x,float),.5,4,0)[0]==-1
 assert rhs(np.array([0,0,0,0,1.]),.5,4,0)[4]==-1

def test_plan_saved_settings_and_numerical_checks():
 j=json.loads((ROOT/'results/grn_ode_semantics_audit.json').read_text());raw=(ROOT/'experiments/grn_ode_semantics_plan.json').read_bytes()
 assert j['plan_sha256']==hashlib.sha256(raw).hexdigest() and len(j['rows'])==120
 assert all(r['initial_count']==243 and r['dt_halving_max_endpoint_difference']<2e-8 for r in j['rows'])
 assert all(r['terminal_wg_on_fraction']==0 for r in j['rows'] if r['a']==0)
 assert all(r['terminal_wg_on_fraction']==8/9 for r in j['rows'] if r['a']==.85 and r['K']==.2)
 assert all(r['terminal_wg_on_fraction']==0 for r in j['rows'] if r['a']==1 and r['K']==.65 and r['n'] in [2,4])

def test_rk4_halving_for_explicit_grid():
 import itertools
 x=np.array(list(itertools.product([0,.5,1],repeat=5)),float)
 fine=integrate(x,.2,4,.85,.05);coarse=integrate(x,.2,4,.85,.1)
 assert np.max(np.abs(fine-coarse))<2e-8
 assert np.mean(fine[:,0]>.5)==8/9
