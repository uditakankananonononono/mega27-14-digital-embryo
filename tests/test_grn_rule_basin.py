import json,re,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
from grn_rule_basin_audit import compute,graph_for,F
def test_vector_successors_equal_independent_scalar_parser_all_variants():
 old=json.loads((ROOT/'results/grn_full_reference_audit.json').read_text());nodes=old['node_order'];lines={}
 for line in F.BNET.splitlines()[1:]:
  if line.strip():target,expr=line.split(',',1);lines[target.strip()]=expr.strip()
 for r in json.loads((ROOT/'results/grn_rule_validity_audit.json').read_text())['rows']:
  changed=dict(lines);changed[r['target']]=re.sub(r'\b'+re.escape(r['regulator'])+r'\b',str(r['replacement']),lines[r['target']]);g=graph_for(changed,nodes);rules=F.parse_bnet('targets,factors\n'+'\n'.join(n+', '+e for n,e in changed.items()))
  for code in range(len(g)):
   state={v:(code>>i)&1 for i,v in enumerate(nodes)};state.update(F.EXT);nxt=F.successor(state,rules,F.EXT);assert g[code]==sum(nxt[v]<<i for i,v in enumerate(nodes))
def test_rule_basin_complete_and_invalid_not_zeroed():
 j=json.loads((ROOT/'results/grn_rule_basin_audit.json').read_text());assert compute()==j and len(j['rows'])==58 and j['baseline_original_reference_starts']==128
 invalid=[r for r in j['rows'] if not r['original_reference_is_fixed']];valid=[r for r in j['rows'] if r['original_reference_is_fixed']];assert len(invalid)==20 and all(r['exact_original_reference_starts'] is None for r in invalid)
 assert len(valid)==38 and sum(r['exact_original_reference_starts']!=128 for r in valid)==8
 assert all(r['fixed_starts']+r['cycle_starts']==16384 for r in j['rows'])
