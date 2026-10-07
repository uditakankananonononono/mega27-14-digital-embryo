"""Exact reference first-hit times, conditional on endpoint admissibility."""
import hashlib,json,re,sys
from collections import Counter
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
from grn_rule_basin_audit import F,graph_for,terminal_ids

def hit_times(graph,ref):
 if graph[ref]!=ref:raise ValueError('reference not fixed')
 labels=np.asarray(terminal_ids(graph));returns=labels==ref
 cur=np.arange(len(graph));times=np.full(len(graph),-1,dtype=int);times[ref]=0
 graph=np.asarray(graph);step=0
 while np.any(returns&(times<0)):
  step+=1
  if step>len(graph):raise AssertionError('unresolved reference basin')
  cur=graph[cur];times[returns&(times<0)&(cur==ref)]=step
 return times.tolist()

def compute():
 raw=(ROOT/'experiments/grn_rule_latency_plan.json').read_bytes();p=json.loads(raw)
 old=json.loads((ROOT/'results/grn_full_reference_audit.json').read_text());ref=old['original_reference_integer'];nodes=old['node_order'];shells=[(s^ref).bit_count() for s in range(1<<len(nodes))];sizes=[shells.count(k) for k in range(len(nodes)+1)];lines={}
 for line in F.BNET.splitlines()[1:]:
  if line.strip():target,expr=line.split(',',1);lines[target.strip()]=expr.strip()
 def summarize(graph):
  if graph[ref]!=ref:return {'reference_valid':False,'shell_histograms':None,'cumulative':None,'eventual_counts':None,'max_latency':None}
  times=hit_times(graph,ref);hist=[dict(sorted(Counter(times[s] for s,k in enumerate(shells) if k==radius and times[s]>=0).items())) for radius in range(len(nodes)+1)]
  return {'reference_valid':True,'shell_histograms':[{str(k):v for k,v in h.items()} for h in hist],'cumulative':{str(h):[sum(v for k,v in bins.items() if k<=h) for bins in hist] for h in p['horizons']},'eventual_counts':[sum(h.values()) for h in hist],'max_latency':max(times)}
 baseline=summarize(graph_for(lines,nodes));rows=[]
 for r in json.loads((ROOT/'results/grn_rule_validity_audit.json').read_text())['rows']:
  altered=dict(lines);altered[r['target']]=re.sub(r'\b'+re.escape(r['regulator'])+r'\b',str(r['replacement']),lines[r['target']]);rows.append({k:r[k] for k in ['target','regulator','replacement']}|summarize(graph_for(altered,nodes)))
 return {'plan_sha256':hashlib.sha256(raw).hexdigest(),'shell_sizes':sizes,'baseline':baseline,'rows':rows,'limits':p['limits']}
if __name__=='__main__':
 j=compute();(ROOT/'results/grn_rule_latency_audit.json').write_text(json.dumps(j,indent=2)+'\n');print('baseline',j['baseline']);valid=[r for r in j['rows'] if r['reference_valid']];print('valid',len(valid),'max_latency_range',min(r['max_latency'] for r in valid),max(r['max_latency'] for r in valid))
