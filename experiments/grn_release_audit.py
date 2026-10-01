"""Exact finite transient forcing followed by restoration of original dynamics."""
import json,hashlib,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'experiments'))
import fragility_index as F
from grn_full_reference_audit import terminal_ids

def compute():
 raw=(ROOT/'experiments/grn_release_plan.json').read_bytes();rules=F.parse_bnet(F.BNET);free=[v for v in rules if v not in F.EXT];states=np.arange(1<<len(free),dtype=np.int32);graph=[]
 for code in states:
  s={v:(int(code)>>i)&1 for i,v in enumerate(free)};s.update(F.EXT);n=F.successor(s,rules,F.EXT);graph.append(sum(n[v]<<i for i,v in enumerate(free)))
 graph=np.array(graph);term=np.array(terminal_ids(graph));wg=1<<free.index('v_wg');refs=[int(s) for s in states if graph[s]==s and s&wg];assert len(refs)==1;ref=refs[0];rows=[]
 for i,node in enumerate(free):
  bit=1<<i
  for val in [0,1]:
   cur=(states&~bit)|(val<<i)
   for step in range(1,9):
    cur=(graph[cur]&~bit)|(val<<i)
    if step in [1,2,4,8]:rows.append({'node':node,'value':val,'forced_updates':step,'n_initial':len(states),'exact_reference_at_release':int(np.sum(cur==ref)),'eventual_reference_after_release':int(np.sum(term[cur]==ref)),'released_reference_fraction':float(np.mean(term[cur]==ref)),'distinct_release_states':len(np.unique(cur))})
 return {'plan_sha256':hashlib.sha256(raw).hexdigest(),'source_model_sha256':hashlib.sha256(F.BNET.encode()).hexdigest(),'node_order':free,'reference_integer':ref,'baseline_reference_count':int(np.sum(term==ref)),'rows':rows,'limits':json.loads(raw)['limits']}
def main():
 j=compute();(ROOT/'results/grn_release_audit.json').write_text(json.dumps(j,indent=2)+'\n')
 for node in ['v_CIA','v_CI_protein','v_ci','v_EN_protein','v_wg','v_PTC_protein']:
  print(node,[(r['value'],r['forced_updates'],r['released_reference_fraction']) for r in j['rows'] if r['node']==node])
if __name__=='__main__':main()
