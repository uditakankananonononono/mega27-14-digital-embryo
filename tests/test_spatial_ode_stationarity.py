import json,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def test_stationarity_snapshot_geometry():
 j=json.loads((ROOT/'results/spatial_ode_stationarity_audit.json').read_text())
 assert j['plan_sha256']==hashlib.sha256((ROOT/'experiments/spatial_ode_stationarity_plan.json').read_bytes()).hexdigest()
 a=np.loadtxt(ROOT/'results/spatial_ode_stationarity_trajectory.csv',delimiter=',',skiprows=1)
 assert a.shape==(2201,21) and np.isfinite(a).all()
 assert np.array_equal(a[:,0],np.arange(2201)*5)
 assert j['prior_overlap_max_abs_error']==0
 assert j['residuals'][0]['max_abs_derivative']>1e-8
 assert j['residuals'][1]['max_abs_derivative']<=1e-8
 assert j['stationary_at_11000_under_frozen_threshold']
 assert np.isclose(j['ptc_cell2_minus_cell4_1100'],a[220,14]-a[220,16],atol=1e-25)
 assert np.isclose(j['ptc_cell2_minus_cell4_11000'],a[-1,14]-a[-1,16],atol=1e-25)
 assert abs(j['ptc_cell2_minus_cell4_11000'])<1e-15
