import importlib.util,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
s=importlib.util.spec_from_file_location('direction',ROOT/'experiments/spatial_ode_direction_feasibility_audit.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
def test_geometry_replay():
 j=m.compute();old=json.loads((ROOT/'results/spatial_ode_direction_feasibility_audit.json').read_text());assert j['input_sha256']==old['input_sha256']
 for a,b in zip(j['linearizations'],old['linearizations']):
  assert np.isclose(np.linalg.norm(a['unit_right_singular_vector']),1)
  assert sorted(x['negative_components'] for x in a['signs'])==[49,82]
  assert np.allclose(sorted(x['maximum_nonnegative_amplitude'] for x in a['signs']), sorted(x['maximum_nonnegative_amplitude'] for x in b['signs']),rtol=1e-3,atol=0)
