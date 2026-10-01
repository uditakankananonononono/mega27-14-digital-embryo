"""Test reference validity before any counterfactual basin claim."""
import hashlib,json,re,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'experiments'))
import fragility_index as F
def compute():
 raw=(ROOT/'experiments/grn_rule_validity_plan.json').read_bytes();old=json.loads((ROOT/'results/grn_full_reference_audit.json').read_text());state=dict(old['original_reference_state']);state.update(F.EXT);lines={}
 for line in F.BNET.splitlines()[1:]:
  if line.strip():target,expr=line.split(',',1);lines[target.strip()]=expr.strip()
 rows=[]
 for target,expr in lines.items():
  for regulator in sorted(set(re.findall(r'v_[A-Za-z0-9_]+',expr))):
   for val in [0,1]:
    altered=dict(lines);altered[target]=re.sub(r'\b'+re.escape(regulator)+r'\b',str(val),expr);text='targets,factors\n'+'\n'.join(t+', '+e for t,e in altered.items());rules=F.parse_bnet(text);nxt=F.successor(state,rules,F.EXT);changed=[v for v in old['node_order'] if nxt[v]!=state[v]]
    rows.append({'target':target,'regulator':regulator,'replacement':val,'original_reference_is_fixed':not changed,'changed_internal_coordinates':changed,'modified_model_sha256':hashlib.sha256(text.encode()).hexdigest()})
 return {'plan_sha256':hashlib.sha256(raw).hexdigest(),'original_model_sha256':hashlib.sha256(F.BNET.encode()).hexdigest(),'rows':rows,'limits':json.loads(raw)['limits']}
if __name__=='__main__':
 j=compute();(ROOT/'results/grn_rule_validity_audit.json').write_text(json.dumps(j,indent=2)+'\n');print('variants',len(j['rows']),'invalid',sum(not r['original_reference_is_fixed'] for r in j['rows']))
