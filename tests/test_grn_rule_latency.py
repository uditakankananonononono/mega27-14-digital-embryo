import itertools,json,sys
from pathlib import Path
import pytest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
from grn_rule_latency_audit import compute,hit_times

def test_latency_matches_shell_and_basin_records():
 j=json.loads((ROOT/'results/grn_rule_latency_audit.json').read_text());assert compute()==j
 old=json.loads((ROOT/'results/grn_rule_radius_audit.json').read_text())
 for r,p in zip([j['baseline']]+j['rows'],[old['baseline']]+old['rows']):
  assert r['reference_valid']==p['reference_valid']
  if r['reference_valid']:
   assert r['eventual_counts']==p['return_counts']
   assert r['cumulative']['32']==r['eventual_counts']
   assert r['shell_histograms'][0]=={'0':1}
   assert max(int(k) for h in r['shell_histograms'] for k in h)==r['max_latency']
  else:assert r['max_latency'] is None and r['shell_histograms'] is None
 assert j['baseline']['max_latency']==3

def test_exhaustive_small_graph_hit_time_oracle():
 for n in range(1,5):
  for graph in itertools.product(range(n),repeat=n):
   for ref in range(n):
    if graph[ref]!=ref:continue
    expected=[]
    for start in range(n):
     cur=start;seen=set();steps=0
     while cur!=ref and cur not in seen:
      seen.add(cur);cur=graph[cur];steps+=1
     expected.append(steps if cur==ref else -1)
    assert hit_times(graph,ref)==expected


def test_invalid_reference_endpoint():
 with pytest.raises(ValueError,match='not fixed'):hit_times([1,0],0)
