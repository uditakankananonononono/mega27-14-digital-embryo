"""Enumerate every external binary context with endpoint validity retained."""
import hashlib,itertools,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'experiments'))
import fragility_index as F
from grn_full_reference_audit import terminal_ids
def compute():
 raw=(ROOT/'experiments/grn_external_regime_plan.json').read_bytes();p=json.loads(raw);old=json.loads((ROOT/'results/grn_full_reference_audit.json').read_text());ref=old['original_reference_integer'];rules=F.parse_bnet(F.BNET);free=[v for v in rules if v not in F.EXT];states=np.arange(1<<len(free));wg=1<<free.index('v_wg');rows=[]
 for bits in itertools.product([0,1],repeat=3):
  ext=dict(zip(p['externals'],bits));graph=[]
  for code in states:
   s={v:(int(code)>>i)&1 for i,v in enumerate(free)};s.update(ext);n=F.successor(s,rules,ext);graph.append(sum(n[v]<<i for i,v in enumerate(free)))
  graph=np.array(graph);term=np.array(terminal_ids(graph));fixed=[int(s) for s in states if graph[s]==s]
  for node,val,duration in [(None,None,0)]+[tuple(c) for c in p['conditions']]:
   cur=states.copy()
   if node is not None:
    i=free.index(node);bit=1<<i;cur=(cur&~bit)|(val<<i)
    for _ in range(duration):cur=(graph[cur]&~bit)|(val<<i)
   ends=term[cur];rows.append({'externals':ext,'node':node,'value':val,'forced_updates':duration,'n_initial':len(states),'original_internal_reference_is_fixed':bool(graph[ref]==ref),'fixed_endpoints':fixed,'counts':{'original_internal_reference':int(np.sum(ends==ref)),'wg_on_fixed':int(np.sum((ends>=0)&((ends&wg)!=0))),'cycle':int(np.sum(ends<0))}})
 return {'plan_sha256':hashlib.sha256(raw).hexdigest(),'model_sha256':hashlib.sha256(F.BNET.encode()).hexdigest(),'original_internal_reference_integer':ref,'node_order':free,'rows':rows,'limits':p['limits']}
if __name__=='__main__':
 j=compute();(ROOT/'results/grn_external_regime_audit.json').write_text(json.dumps(j,indent=2)+'\n');print(json.dumps(j['rows'],indent=2))
