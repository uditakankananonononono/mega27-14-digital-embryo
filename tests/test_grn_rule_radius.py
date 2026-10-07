import json,math,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
from grn_rule_radius_audit import compute
def test_all_shells_and_global_totals():
 j=json.loads((ROOT/'results/grn_rule_radius_audit.json').read_text());assert compute()==j;assert j['shell_sizes']==[math.comb(14,k) for k in range(15)]
 old=json.loads((ROOT/'results/grn_rule_basin_audit.json').read_text());assert sum(j['baseline']['return_counts'])==128
 for r,p in zip(j['rows'],old['rows']):
  assert all(r[k]==p[k] for k in ['target','regulator','replacement'])
  if r['reference_valid']:assert sum(r['return_counts'])==p['exact_original_reference_starts'] and r['return_counts'][0]==1
  else:assert r['return_counts'] is None and r['fractions'] is None
 assert j['baseline']['return_counts'][1:3]==[7,21]
 assert sum(r['reference_valid'] for r in j['rows'])==38
