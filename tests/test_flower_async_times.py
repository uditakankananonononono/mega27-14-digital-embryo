import json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def test_absorption_time_boundary_and_finite_horizon():
 j=json.loads((ROOT/'results/flower_async_absorption_times.json').read_text());a=np.load(ROOT/'results/flower_async_absorption_times.npz')
 assert a['expected_updates'].min()>=-1e-9;assert a['variance_updates'].min()>=-1e-8
 for c in j['closed_classes']:assert np.isclose(a['expected_updates'][c[0]],0)
 for r in j['rows']:
  vals=list(r['recovered_original_fate_by_update_horizon'].values());assert all(x<=y+1e-10 for x,y in zip(vals,vals[1:]));assert vals[-1]<=r['async_fixed_one_bit_retention']+1e-10
 r={x['organ']:x for x in j['rows']}
 assert r['petal']['recovered_original_fate_by_update_horizon']['12']>r['carpel']['recovered_original_fate_by_update_horizon']['12']
 assert r['petal']['recovered_original_fate_by_update_horizon']['24']<r['carpel']['recovered_original_fate_by_update_horizon']['24']
