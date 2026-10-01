"""Weight the same complete state census under frozen nonuniform initial measures."""
import hashlib,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'experiments'))
import fragility_index as F
from grn_full_reference_audit import terminal_ids
def compute():
 raw=(ROOT/'experiments/grn_initial_measure_plan.json').read_bytes();p=json.loads(raw);rules=F.parse_bnet(F.BNET);free=[v for v in rules if v not in F.EXT];states=np.arange(1<<len(free));graph=[]
 for code in states:
  s={v:(int(code)>>i)&1 for i,v in enumerate(free)};s.update(F.EXT);n=F.successor(s,rules,F.EXT);graph.append(sum(n[v]<<i for i,v in enumerate(free)))
 graph=np.array(graph);term=np.array(terminal_ids(graph));wg=1<<free.index('v_wg');ref=next(int(s) for s in states if graph[s]==s and s&wg)
 pop=np.array([int(s).bit_count() for s in states]);near=np.array([(int(s)^ref).bit_count()<=2 for s in states]);weights={'uniform':np.ones(len(states))/len(states),'independent_bit_p025':.25**pop*.75**(len(free)-pop),'independent_bit_p075':.75**pop*.25**(len(free)-pop),'uniform_hamming_distance_le2_from_reference':near/near.sum()}
 rows=[]
 for node,val,duration in [(None,None,0)]+[tuple(c) for c in p['conditions']]:
  cur=states.copy()
  if node is not None:
   bit=1<<free.index(node);cur=(cur&~bit)|(val<<free.index(node))
   for _ in range(duration):cur=(graph[cur]&~bit)|(val<<free.index(node))
  ends=term[cur];masks={'exact_original_reference':ends==ref,'any_wg_on_fixed_point':(ends>=0)&((ends&wg)!=0),'cycle':ends<0}
  for name,w in weights.items():
   assert np.isclose(w.sum(),1)
   rows.append({'node':node,'value':val,'forced_updates':duration,'measure':name,'support_states':int(np.sum(w>0)),'probabilities':{k:float(np.sum(w[mask])) for k,mask in masks.items()}})
 return {'plan_sha256':hashlib.sha256(raw).hexdigest(),'model_sha256':hashlib.sha256(F.BNET.encode()).hexdigest(),'reference_integer':ref,'node_order':free,'n_census_states':len(states),'rows':rows,'limits':p['limits']}
if __name__=='__main__':
 j=compute();(ROOT/'results/grn_initial_measure_audit.json').write_text(json.dumps(j,indent=2)+'\n')
 for r in j['rows']:print(r)
