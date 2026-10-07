import importlib.util,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
s=importlib.util.spec_from_file_location('non',ROOT/'experiments/spatial_ode_nonnormal_audit.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
def test_frozen_grid_and_archived_identity():
 j=m.compute();old=json.loads((ROOT/'results/spatial_ode_nonnormal_audit.json').read_text());assert j['input_sha256']==old['input_sha256'];assert not j['native_jacobian_finite']
 for a,b in zip(j['linearizations'],old['linearizations']):
  assert [v['time'] for v in a['propagator_norms']]==[0,1,10,100]
  assert np.allclose([v['operator_2_norm'] for v in a['propagator_norms']],[v['operator_2_norm'] for v in b['propagator_norms']],rtol=1e-9)
  assert a['spectral_abscissa']<0 and a['euclidean_logarithmic_norm']>0
