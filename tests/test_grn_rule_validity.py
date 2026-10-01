import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
from grn_rule_validity_audit import compute
def test_full_counterfactual_reference_validity():
 j=json.loads((ROOT/'results/grn_rule_validity_audit.json').read_text());assert compute()==j
 assert len(j['rows'])==58 and sum(not r['original_reference_is_fixed'] for r in j['rows'])==20
 assert all(r['changed_internal_coordinates']==([] if r['original_reference_is_fixed'] else [r['target']]) for r in j['rows'])
 assert len({(r['target'],r['regulator'],r['replacement']) for r in j['rows']})==58
