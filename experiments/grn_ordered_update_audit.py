"""Exact complete-sweep graphs under three fixed sequential update orders."""
import hashlib,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'experiments'))
import fragility_index as F
from grn_full_reference_audit import terminal_ids
def graph_for(rules,free,order,forced=None):
 graph=[]
 for code in range(1<<len(free)):
  s={v:(code>>i)&1 for i,v in enumerate(free)};s.update(F.EXT)
  if forced:s[forced[0]]=forced[1]
  for node in order:
   code_fn=rules[node][0];s[node]=int(bool(eval(code_fn,{'__builtins__':{}},s)))
   if forced:s[forced[0]]=forced[1]
  graph.append(sum(s[v]<<i for i,v in enumerate(free)))
 return np.array(graph)
def compute():
 raw=(ROOT/'experiments/grn_ordered_update_plan.json').read_bytes();p=json.loads(raw);rules=F.parse_bnet(F.BNET);free=[v for v in rules if v not in F.EXT];ref=json.loads((ROOT/'results/grn_full_reference_audit.json').read_text())['original_reference_integer'];states=np.arange(1<<len(free));wg=1<<free.index('v_wg');orders={'source_order':free,'reverse_source_order':free[::-1],'fixed_permutation_seed2183':np.random.default_rng(2183).permutation(free).tolist()};rows=[]
 for name,order in orders.items():
  graph=graph_for(rules,free,order);term=np.array(terminal_ids(graph));assert graph[ref]==ref
  for node,val,duration in [(None,None,0)]+[tuple(c) for c in p['conditions']]:
   cur=states.copy()
   if node is not None:
    i=free.index(node);cur=(cur&~(1<<i))|(val<<i);forcedgraph=graph_for(rules,free,order,(node,val))
    for _ in range(duration):cur=forcedgraph[cur]
   ends=term[cur];rows.append({'schedule':name,'node_order':order,'forced_node':node,'forced_value':val,'forced_sweeps':duration,'n_initial':len(states),'graph_sha256':hashlib.sha256(graph.astype('<u2').tobytes()).hexdigest(),'counts':{'original_reference':int(np.sum(ends==ref)),'wg_on_fixed_sweep_endpoint':int(np.sum((ends>=0)&((ends&wg)!=0))),'sweep_cycle':int(np.sum(ends<0))}})
 return {'plan_sha256':hashlib.sha256(raw).hexdigest(),'model_sha256':hashlib.sha256(F.BNET.encode()).hexdigest(),'reference_integer':ref,'rows':rows,'limits':p['limits']}
if __name__=='__main__':
 j=compute();(ROOT/'results/grn_ordered_update_audit.json').write_text(json.dumps(j,indent=2)+'\n')
 for r in j['rows']:print(r['schedule'],r['forced_node'],r['counts'])
