import itertools,json,sys
from pathlib import Path
import pytest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
from grn_rule_boundary_audit import compute,summarize

def test_boundary_reconciles_existing_basin_counts():
 j=json.loads((ROOT/'results/grn_rule_boundary_audit.json').read_text());assert compute()==j
 old=json.loads((ROOT/'results/grn_rule_basin_audit.json').read_text())
 assert j['baseline']['basin_size']==128 and j['baseline']['exit_edges']==896 and j['baseline']['exit_fraction']==.5
 for r,p in zip(j['rows'],old['rows']):
  assert r['reference_valid']==p['original_reference_is_fixed']
  if r['reference_valid']:
   assert r['basin_size']==p['exact_original_reference_starts']
   assert r['directed_edges']==14*r['basin_size']==r['exit_edges']+r['internal_edges']
   assert sum(x['exit_edges'] for x in r['per_node'])==r['exit_edges']
  else:assert r['exit_fraction'] is None and r['per_node'] is None

def test_small_graph_boundary_oracle():
 for graph in itertools.product(range(4),repeat=4):
  for ref in range(4):
   j=summarize(graph,ref,['a','b'])
   if graph[ref]!=ref:assert not j['reference_valid'];continue
   basin=[]
   for start in range(4):
    seen=set();cur=start
    while cur!=ref and cur not in seen:seen.add(cur);cur=graph[cur]
    if cur==ref:basin.append(start)
   exits=sum((s^(1<<i)) not in basin for s in basin for i in range(2))
   assert j['exit_edges']==exits and j['exit_fraction']==exits/(2*len(basin))

def test_invalid_boundary_geometry():
 with pytest.raises(ValueError):summarize([0,1,2],0,['a','b'])
