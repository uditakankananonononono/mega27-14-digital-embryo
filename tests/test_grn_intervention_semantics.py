import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
import fragility_index as F
from grn_maintained_clamp_audit import graph_labels

def test_historical_clamp_is_not_maintained():
 j=json.loads((ROOT/'results/grn_intervention_semantics_audit.json').read_text());r=F.parse_bnet(F.BNET)
 nxt=F.successor(j['initial_state'],r,F.EXT)
 assert nxt==j['first_unclamped_successor'] and nxt['v_en']==0!=j['intended_clamp']
 assert len(j['full_exit_one_nodes'])==7 and 'v_wg' in j['full_exit_one_nodes']

def test_exact_functional_graph_classification():
 assert graph_labels([0,0,3,2],1)==['wg_off_fixed','wg_off_fixed','cycle','cycle']
 assert graph_labels([1,1,1,2],1)==['wg_on_fixed']*4

def test_maintained_clamp_fractions_and_tautology():
 j=json.loads((ROOT/'results/grn_maintained_clamp_audit.json').read_text());rows={(r['node'],r['value']):r for r in j['rows']}
 assert len(rows)==28 and j['baseline_wg_on_fixed_states']==128
 assert all(sum(r['outcomes'].values())==r['initial_state_count']==8192 for r in rows.values())
 assert rows['v_CIA',1]['wg_on_fixed_fraction']==.25
 assert rows['v_wg',1]['wg_on_fixed_fraction']==1
 assert rows['v_wg',0]['wg_on_fixed_fraction']==0
 # Source rule independently checks maintained-clamp trajectories.
 rules=F.parse_bnet(F.BNET);free=j['node_order']
 for node,val in [('v_CIA',1),('v_ci',1),('v_EN_protein',0)]:
  live=[v for v in free if v!=node];wins=0
  import itertools
  for bits in itertools.product([0,1],repeat=len(live)):
   s=dict(zip(live,bits));s[node]=val;s.update(F.EXT);seen=set()
   while True:
    key=tuple(s[v] for v in free)
    if key in seen:break
    seen.add(key);n=F.successor(s,rules,F.EXT);n[node]=val
    if all(n[v]==s[v] for v in free):wins+=int(n['v_wg']==1);break
    s=n
  assert wins==rows[node,val]['outcomes']['wg_on_fixed']
