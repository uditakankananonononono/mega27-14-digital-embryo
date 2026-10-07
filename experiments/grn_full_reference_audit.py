"""Enumerate exact attractor identity, avoiding a clamped phenotype marker tautology."""
import hashlib,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'experiments'))
import fragility_index as F

def terminal_ids(graph):
 from numbers import Integral
 if len(graph)==0 or any(isinstance(x,bool) or not isinstance(x,Integral) or x<0 or x>=len(graph) for x in graph):raise ValueError('invalid functional graph')
 labels=[None]*len(graph)
 for start in range(len(graph)):
  if labels[start] is not None:continue
  path=[];pos={};cur=start
  while labels[cur] is None and cur not in pos:
   pos[cur]=len(path);path.append(cur);cur=graph[cur]
  if labels[cur] is not None:label=labels[cur]
  else:label=cur if len(path)-pos[cur]==1 else -1
  for x in path:labels[x]=label
 return labels

def compute():
 rules=F.parse_bnet(F.BNET);free=[v for v in rules if v not in F.EXT];states=range(1<<len(free));graph=[]
 for code in states:
  state={v:(code>>i)&1 for i,v in enumerate(free)};state.update(F.EXT);nxt=F.successor(state,rules,F.EXT)
  graph.append(sum(nxt[v]<<i for i,v in enumerate(free)))
 wg=1<<free.index('v_wg');refs=[x for x in states if graph[x]==x and x&wg];assert len(refs)==1;ref=refs[0];baseline=terminal_ids(graph);rows=[]
 for i,node in enumerate(free):
  bit=1<<i
  for val in [0,1]:
   altered=[(graph[(s&~bit)|(val<<i)]&~bit)|(val<<i) for s in states];labels=terminal_ids(altered);initial=[s for s in states if (s>>i)&1==val];ends=[labels[s] for s in initial]
   count={'full_reference':sum(x==ref for x in ends),'projected_reference':sum(x>=0 and (x&~bit)==(ref&~bit) for x in ends),'wg_on_fixed':sum(x>=0 and x&wg!=0 for x in ends),'cycle':sum(x<0 for x in ends)}
   rows.append({'node':node,'value':val,'reference_value':(ref>>i)&1,'n_initial':len(initial),'counts':count,'fractions':{k:v/len(initial) for k,v in count.items()},'distinct_fixed_endpoints':sorted(set(x for x in ends if x>=0)),'graph_sha256':hashlib.sha256(b''.join(x.to_bytes(2,'big') for x in altered)).hexdigest()})
 raw=(ROOT/'experiments/grn_full_reference_plan.json').read_bytes()
 return {'plan_sha256':hashlib.sha256(raw).hexdigest(),'model_sha256':hashlib.sha256(F.BNET.encode()).hexdigest(),'node_order':free,'original_reference_integer':ref,'original_reference_state':{v:(ref>>i)&1 for i,v in enumerate(free)},'baseline_exact_reference_starts':baseline.count(ref),'rows':rows,'limits':json.loads(raw)['limits']}
def main():
 out=compute();(ROOT/'results/grn_full_reference_audit.json').write_text(json.dumps(out,indent=2)+'\n')
 print('reference',out['original_reference_state'])
 for r in out['rows']:print(r['node'],r['value'],r['fractions'])
if __name__=='__main__':main()
