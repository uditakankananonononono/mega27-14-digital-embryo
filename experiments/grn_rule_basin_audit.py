"""Exact counterfactual basin enumeration with inadmissible reference endpoints retained."""
import hashlib,json,re,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'experiments'))
import fragility_index as F
from grn_full_reference_audit import terminal_ids
def graph_for(lines,nodes):
 codes=np.arange(1<<len(nodes));env={n:(codes>>i)&1 for i,n in enumerate(nodes)};env.update(F.EXT);graph=np.zeros(len(codes),dtype=np.int64)
 for i,node in enumerate(nodes):
  expr=lines[node]
  # NumPy boolean arrays preserve the original parenthesized &, |, ! grammar.
  local={n:np.asarray(v,dtype=bool) for n,v in env.items()};expr=re.sub(r'\b([01])\b',lambda m:'CONST'+m[1],expr).replace('!','~');local.update(CONST0=np.bool_(False),CONST1=np.bool_(True));val=np.asarray(eval(compile(expr,'<vector-bnet>','eval'),{'__builtins__':{}},local),dtype=np.int64);graph|=val<<i
 return graph.tolist()
def compute():
 old=json.loads((ROOT/'results/grn_full_reference_audit.json').read_text());nodes=old['node_order'];ref=old['original_reference_integer'];lines={}
 for line in F.BNET.splitlines()[1:]:
  if line.strip():target,expr=line.split(',',1);lines[target.strip()]=expr.strip()
 basegraph=graph_for(lines,nodes);base=terminal_ids(basegraph).count(ref);assert base==old['baseline_exact_reference_starts'];rows=[]
 for r in json.loads((ROOT/'results/grn_rule_validity_audit.json').read_text())['rows']:
  altered=dict(lines);altered[r['target']]=re.sub(r'\b'+re.escape(r['regulator'])+r'\b',str(r['replacement']),lines[r['target']]);graph=graph_for(altered,nodes);labels=terminal_ids(graph);valid=graph[ref]==ref;assert valid==r['original_reference_is_fixed'];rows.append({k:r[k] for k in ['target','regulator','replacement','original_reference_is_fixed']}|{'n_initial':len(graph),'exact_original_reference_starts':labels.count(ref) if valid else None,'cycle_starts':labels.count(-1),'fixed_starts':sum(x>=0 for x in labels),'distinct_fixed_endpoints':sorted(set(x for x in labels if x>=0)),'graph_sha256':hashlib.sha256(b''.join(x.to_bytes(2,'big') for x in graph)).hexdigest()})
 raw=(ROOT/'experiments/grn_rule_basin_plan.json').read_bytes();return {'plan_sha256':hashlib.sha256(raw).hexdigest(),'baseline_original_reference_starts':base,'rows':rows,'limits':json.loads(raw)['limits']}
if __name__=='__main__':
 j=compute();(ROOT/'results/grn_rule_basin_audit.json').write_text(json.dumps(j,indent=2)+'\n');valid=[r['exact_original_reference_starts'] for r in j['rows'] if r['original_reference_is_fixed']];print('baseline',j['baseline_original_reference_starts'],'valid',len(valid),'range',min(valid),max(valid),'changed',sum(v!=j['baseline_original_reference_starts'] for v in valid),'cycle range',min(r['cycle_starts'] for r in j['rows']),max(r['cycle_starts'] for r in j['rows']))
