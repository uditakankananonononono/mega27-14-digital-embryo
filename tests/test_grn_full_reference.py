import json,sys,hashlib,itertools
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
from grn_full_reference_audit import terminal_ids,compute
import fragility_index as F

def test_terminal_ids_graph_examples():
 assert terminal_ids([0,0,3,2])==[0,0,-1,-1]
 assert terminal_ids([1,1,1,2])==[1,1,1,1]

def test_full_reference_saved_and_projected_false_positives():
 j=json.loads((ROOT/'results/grn_full_reference_audit.json').read_text());assert compute()==j
 assert j['baseline_exact_reference_starts']==128
 assert j['plan_sha256']==hashlib.sha256((ROOT/'experiments/grn_full_reference_plan.json').read_bytes()).hexdigest()
 rows={(r['node'],r['value']):r for r in j['rows']};assert len(rows)==28
 for r in rows.values():
  f=r['fractions'];assert f['full_reference']<=f['projected_reference']<=f['wg_on_fixed']
  if r['value']!=r['reference_value']:assert f['full_reference']==0
 assert rows['v_wg',1]['fractions']['full_reference']==1
 for key in [('v_PTC_protein',0),('v_hh',1),('v_ptc',0)]:
  assert rows[key]['fractions']['wg_on_fixed']==1/128 and rows[key]['fractions']['projected_reference']==0

def test_independent_scalar_maintained_reference_trajectories():
 j=json.loads((ROOT/'results/grn_full_reference_audit.json').read_text());rules=F.parse_bnet(F.BNET);free=j['node_order'];ref=j['original_reference_state']
 for node,val in [('v_PTC_protein',0),('v_CIA',1),('v_wg',1)]:
  live=[v for v in free if v!=node];full=projected=marker=0
  for bits in itertools.product([0,1],repeat=len(live)):
   s=dict(zip(live,bits));s[node]=val;s.update(F.EXT);seen=set()
   while True:
    key=tuple(s[v] for v in free)
    if key in seen:break
    seen.add(key);nxt=F.successor(s,rules,F.EXT);nxt[node]=val
    if all(nxt[v]==s[v] for v in free):
     full+=all(nxt[v]==ref[v] for v in free);projected+=all(nxt[v]==ref[v] for v in live);marker+=nxt['v_wg']==1;break
    s=nxt
  r=next(r for r in j['rows'] if r['node']==node and r['value']==val)
  assert [full,projected,marker]==[r['counts'][k] for k in ['full_reference','projected_reference','wg_on_fixed']]
