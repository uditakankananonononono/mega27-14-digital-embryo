"""Exact shell-conditioned basin counts under saved rule counterfactuals."""
import hashlib,json,re,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
from grn_rule_basin_audit import F,graph_for,terminal_ids
def compute():
 raw=(ROOT/'experiments/grn_rule_radius_plan.json').read_bytes();old=json.loads((ROOT/'results/grn_full_reference_audit.json').read_text());ref=old['original_reference_integer'];nodes=old['node_order'];lines={}
 for line in F.BNET.splitlines()[1:]:
  if line.strip():t,e=line.split(',',1);lines[t.strip()]=e.strip()
 shells=[(s^ref).bit_count() for s in range(1<<len(nodes))];sizes=[shells.count(k) for k in range(len(nodes)+1)]
 def summarize(graph):
  labels=terminal_ids(graph);valid=graph[ref]==ref;counts=[sum(labels[s]==ref for s,r in enumerate(shells) if r==k) for k in range(len(nodes)+1)] if valid else None
  return {'reference_valid':valid,'return_counts':counts,'fractions':[c/n for c,n in zip(counts,sizes)] if valid else None}
 baseline=summarize(graph_for(lines,nodes));rows=[]
 for r in json.loads((ROOT/'results/grn_rule_validity_audit.json').read_text())['rows']:
  altered=dict(lines);altered[r['target']]=re.sub(r'\b'+re.escape(r['regulator'])+r'\b',str(r['replacement']),lines[r['target']]);rows.append({k:r[k] for k in ['target','regulator','replacement']}|summarize(graph_for(altered,nodes)))
 return {'plan_sha256':hashlib.sha256(raw).hexdigest(),'shell_sizes':sizes,'baseline':baseline,'rows':rows,'limits':json.loads(raw)['limits']}
if __name__=='__main__':
 j=compute();(ROOT/'results/grn_rule_radius_audit.json').write_text(json.dumps(j,indent=2)+'\n');print('baseline',j['baseline']);valid=[r for r in j['rows'] if r['reference_valid']]
 for k in [1,2]:print('radius',k,'minmax',min(r['fractions'][k] for r in valid),max(r['fractions'][k] for r in valid),'different',sum(r['return_counts'][k]!=j['baseline']['return_counts'][k] for r in valid))
