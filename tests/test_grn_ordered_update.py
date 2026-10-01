import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
from grn_ordered_update_audit import compute,graph_for
import fragility_index as F
def test_complete_saved_ordered_census():
 j=json.loads((ROOT/'results/grn_ordered_update_audit.json').read_text());assert compute()==j and len(j['rows'])==9
 assert [r['counts']['original_reference'] for r in j['rows']]==[1024,16384,8192,1024,16384,4096,1024,16384,8192]
 assert all(r['n_initial']==16384 and r['counts']['sweep_cycle']==0 for r in j['rows'])
def test_sweep_fixed_points_are_original_rule_fixed_points():
 j=json.loads((ROOT/'results/grn_ordered_update_audit.json').read_text());rules=F.parse_bnet(F.BNET);free=list(rules)
 for r in j['rows']:
  if r['forced_node'] is not None:continue
  graph=graph_for(rules,free,r['node_order'])
  for code,nxt in enumerate(graph):
   if code!=nxt:continue
   state={v:(code>>i)&1 for i,v in enumerate(free)};state.update(F.EXT);successor=F.successor(state,rules,F.EXT)
   assert all(successor[v]==state[v] for v in free)
def test_microstep_forcing_remains_forced():
 rules=F.parse_bnet(F.BNET);free=list(rules)
 for node,val in [('v_wg',1),('v_CIA',1)]:
  graph=graph_for(rules,free,free,(node,val));bit=1<<free.index(node)
  assert all(bool(s&bit)==bool(val) for s in graph)
