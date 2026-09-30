"""Exact finite functional-graph census of maintained binary interventions.
Original historical initial-restriction results remain untouched.
"""
import itertools,json,sys,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'experiments'))
import fragility_index as F

def graph_labels(next_map,wg_bit):
 labels=[None]*len(next_map)
 for start in range(len(next_map)):
  if labels[start] is not None:continue
  path=[];pos={};cur=start
  while labels[cur] is None and cur not in pos:
   pos[cur]=len(path);path.append(cur);cur=next_map[cur]
  if labels[cur] is not None:label=labels[cur]
  else:
   cycle=path[pos[cur]:];label=('wg_on_fixed' if cur&wg_bit else 'wg_off_fixed') if len(cycle)==1 else 'cycle'
  for state in path:labels[state]=label
 return labels

def main():
 rules=F.parse_bnet(F.BNET);free=[v for v in rules if v not in F.EXT];n=len(free);states=range(1<<n)
 encode=lambda s:sum(int(s[v])<<i for i,v in enumerate(free))
 successor=[]
 for state in states:
  s={v:(state>>i)&1 for i,v in enumerate(free)};s.update(F.EXT)
  successor.append(encode(F.successor(s,rules,F.EXT)))
 wg_bit=1<<free.index('v_wg');baseline=graph_labels(successor,wg_bit)
 assert baseline.count('wg_on_fixed')==128
 old=json.loads((ROOT/'results/grn_fragility.json').read_text());rows=[]
 for i,node in enumerate(free):
  bit=1<<i
  for val in [0,1]:
   # Persist clamp before every rule evaluation and after every transition.
   graph=[(successor[(state&~bit)|(val<<i)]&~bit)|(val<<i) for state in states]
   assert all(((nxt>>i)&1)==val for nxt in graph)
   labels=graph_labels(graph,wg_bit);initial=[s for s in states if ((s>>i)&1)==val]
   counts={k:sum(labels[s]==k for s in initial) for k in ['wg_on_fixed','wg_off_fixed','cycle']}
   assert sum(counts.values())==8192
   key=node+('=ON' if val else '=OFF');historic=old['knockout_predictions'][key]['basin_fraction']
   row={'node':node,'value':val,'initial_state_count':len(initial),'outcomes':counts,'wg_on_fixed_fraction':counts['wg_on_fixed']/len(initial),'historical_initial_restriction_fraction':historic,'graph_sha256':hashlib.sha256(bytes().join(x.to_bytes(2,'big') for x in graph)).hexdigest()}
   rows.append(row);print(key,row['wg_on_fixed_fraction'],historic,flush=True)
 out={'source_bnet_sha256':hashlib.sha256(F.BNET.encode()).hexdigest(),'externals':F.EXT,'node_order':free,'baseline_wg_on_fixed_states':128,'baseline_fraction':128/16384,'definition':'Fraction of uniformly weighted states in the 13-node nonclamped initial space reaching any wg=1 fixed point, under transition with chosen node forced before evaluation and after every step. Not exact biological WT state or a measured intervention.','rows':rows,'limits':['Binary maintained clamp and literal Boolean rules; no physical validation','Uniform reduced-space initial measure differs from baseline full state space','Fixed-point wg marker is not a complete wild-type phenotype','Historical field named knockout_predictions was initial restriction, remains unchanged']}
 (ROOT/'results/grn_maintained_clamp_audit.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
