import hashlib,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def test_jacobian_failure_and_independent_difference_checks():
 j=json.loads((ROOT/'results/spatial_ode_jacobian_audit.json').read_text())
 assert j['plan_sha256']==hashlib.sha256((ROOT/'experiments/spatial_ode_jacobian_plan.json').read_bytes()).hexdigest()
 assert not j['native_jacobian_finite'] and j['native_nonfinite_entries']==1320
 assert j['spectral_abscissa'] is None and j['positive_modes_gt1e8'] is None
 assert len(j['native_nonfinite_columns'])==10
 assert j['prior_endpoint_max_abs_error']==0
 for label,c in zip(['1e5','1e6'],j['finite_difference_checks']):
  a=np.loadtxt(ROOT/f'results/spatial_ode_jacobian_fd{label}.csv',delimiter=',');assert a.shape==(132,132) and np.isfinite(a).all()
  ev=np.linalg.eigvals(a)
  assert c['negative_modes_ltminus1e8']==132 and c['positive_modes_gt1e8']==0
  assert np.isclose(ev.real.max(),c['finite_difference_spectral_abscissa'],atol=1e-9)
