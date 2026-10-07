"""Exact basin-directed flip boundaries, not biological perturbation frequencies."""
import hashlib,json,re,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
from grn_rule_basin_audit import F,graph_for,terminal_ids

def summarize(graph,ref,nodes):
 if len(graph)!=(1<<len(nodes)):raise ValueError('state geometry')
 if graph[ref]!=ref:return {'reference_valid':False,'basin_size':None,'directed_edges':None,'exit_edges':None,'internal_edges':None,'exit_fraction':None,'per_node':None,'reference_flip_exits':None}
 labels=np.asarray(terminal_ids(graph));inside=labels==ref;members=np.flatnonzero(inside);per=[]
 for i,node in enumerate(nodes):
  exits=int(np.sum(~inside[members^(1<<i)]));per.append({'node':node,'exit_edges':exits,'exit_fraction':exits/len(members)})
 edges=len(nodes)*len(members);exits=sum(r['exit_edges'] for r in per)
 return {'reference_valid':True,'basin_size':len(members),'directed_edges':edges,'exit_edges':exits,'internal_edges':edges-exits,'exit_fraction':exits/edges,'per_node':per,'reference_flip_exits':[node for i,node in enumerate(nodes) if not inside[ref^(1<<i)]]}

def compute():
 raw=(ROOT/'experiments/grn_rule_boundary_plan.json').read_bytes();p=json.loads(raw);old=json.loads((ROOT/'results/grn_full_reference_audit.json').read_text());ref=old['original_reference_integer'];nodes=old['node_order'];lines={}
 for line in F.BNET.splitlines()[1:]:
  if line.strip():target,expr=line.split(',',1);lines[target.strip()]=expr.strip()
 baseline=summarize(graph_for(lines,nodes),ref,nodes);rows=[]
 for r in json.loads((ROOT/'results/grn_rule_validity_audit.json').read_text())['rows']:
  altered=dict(lines);altered[r['target']]=re.sub(r'\b'+re.escape(r['regulator'])+r'\b',str(r['replacement']),lines[r['target']]);rows.append({k:r[k] for k in ['target','regulator','replacement']}|summarize(graph_for(altered,nodes),ref,nodes))
 return {'plan_sha256':hashlib.sha256(raw).hexdigest(),'baseline':baseline,'rows':rows,'limits':p['limits']}
if __name__=='__main__':
 j=compute();(ROOT/'results/grn_rule_boundary_audit.json').write_text(json.dumps(j,indent=2)+'\n');print('baseline',j['baseline']);v=[r for r in j['rows'] if r['reference_valid']];print('range',min(r['exit_fraction'] for r in v),max(r['exit_fraction'] for r in v),'different',sum(r['exit_fraction']!=j['baseline']['exit_fraction'] for r in v))
 for r in v:
  if r['basin_size']!=128:print(r['target'],r['regulator'],r['replacement'],r['basin_size'],r['exit_edges'],r['exit_fraction'])
