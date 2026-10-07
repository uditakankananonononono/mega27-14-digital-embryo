import json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def test_frozen_tolerance_identity():
 j=json.loads((ROOT/'results/spatial_ode_tolerance_audit.json').read_text());assert j['plan_sha256']==hashlib.sha256((ROOT/'experiments/spatial_ode_tolerance_plan.json').read_bytes()).hexdigest()
 assert j['tight_rtol']==1e-10 and j['tight_atol']==1e-16
 assert all(r['gain_max_abs_change']<2e-7 for r in j['comparison'])
