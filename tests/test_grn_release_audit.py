import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
from grn_release_audit import compute
import fragility_index as F

def test_exact_release_saved_counts():
 j=json.loads((ROOT/'results/grn_release_audit.json').read_text());assert compute()==j and len(j['rows'])==112
 assert j['baseline_reference_count']==128 and all(r['n_initial']==16384 for r in j['rows'])
 rows={(r['node'],r['value'],r['forced_updates']):r for r in j['rows']}
 assert [rows['v_wg',1,k]['released_reference_fraction'] for k in [1,2,4,8]]==[.0625,.125,.5,1]
 assert [rows['v_CIA',1,k]['released_reference_fraction'] for k in [1,2,4,8]]==[.03125,.0625,.25,.25]
 assert all(r['exact_reference_at_release']<=r['eventual_reference_after_release'] for r in j['rows'])

def test_independent_scalar_release_from_all_states():
 j=json.loads((ROOT/'results/grn_release_audit.json').read_text());free=j['node_order'];rules=F.parse_bnet(F.BNET);ref=j['reference_integer']
 for node,val,duration in [('v_wg',1,2),('v_CIA',1,4)]:
  count=0
  for code in range(16384):
   s={v:(code>>i)&1 for i,v in enumerate(free)};s.update(F.EXT);s[node]=val
   for _ in range(duration):s=F.successor(s,rules,F.EXT);s[node]=val
   seen=set()
   while True:
    key=tuple(s[v] for v in free)
    if key in seen:break
    seen.add(key);nxt=F.successor(s,rules,F.EXT)
    if all(nxt[v]==s[v] for v in free):
     count+=sum(nxt[v]<<i for i,v in enumerate(free))==ref;break
    s=nxt
  row=next(r for r in j['rows'] if r['node']==node and r['value']==val and r['forced_updates']==duration);assert count==row['eventual_reference_after_release']
