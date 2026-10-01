import json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
from grn_initial_measure_audit import compute
def test_exact_reproduction_and_frozen_measures():
 j=json.loads((ROOT/'results/grn_initial_measure_audit.json').read_text());assert compute()==j
 assert len(j['rows'])==44 and j['n_census_states']==16384
 assert all(0<=v<=1+1e-12 for r in j['rows'] for v in r['probabilities'].values())
 assert all(r['support_states']==106 for r in j['rows'] if 'hamming' in r['measure'])
def test_uniform_rows_match_release_census():
 j=json.loads((ROOT/'results/grn_initial_measure_audit.json').read_text());old=json.loads((ROOT/'results/grn_release_audit.json').read_text())
 for r in j['rows']:
  if r['measure']!='uniform' or r['node'] is None:continue
  other=next(x for x in old['rows'] if (x['node'],x['value'],x['forced_updates'])==(r['node'],r['value'],r['forced_updates']))
  assert np.isclose(r['probabilities']['exact_original_reference'],other['released_reference_fraction'])
def test_endpoint_and_near_reference_counts():
 j=json.loads((ROOT/'results/grn_initial_measure_audit.json').read_text())
 assert all(r['probabilities']['cycle']==0 and r['probabilities']['exact_original_reference']==r['probabilities']['any_wg_on_fixed_point'] for r in j['rows'])
 r=next(r for r in j['rows'] if r['node'] is None and 'hamming' in r['measure']);assert np.isclose(r['probabilities']['exact_original_reference'],29/106)
