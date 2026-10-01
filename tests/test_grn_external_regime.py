import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
from grn_external_regime_audit import compute
def test_all_regimes_and_validity():
 j=json.loads((ROOT/'results/grn_external_regime_audit.json').read_text());assert compute()==j and len(j['rows'])==24
 base=[r for r in j['rows'] if r['node'] is None];assert len(base)==8 and sum(r['original_internal_reference_is_fixed'] for r in base)==3
 assert all(r['counts']['original_internal_reference']==0 for r in j['rows'] if not r['original_internal_reference_is_fixed'])
 assert all(r['counts']['cycle']==0 and r['counts']['original_internal_reference']==r['counts']['wg_on_fixed'] for r in j['rows'])
def test_original_regime_reproduces_release():
 j=json.loads((ROOT/'results/grn_external_regime_audit.json').read_text());rows=[r for r in j['rows'] if tuple(r['externals'].values())==(1,0,0)]
 assert [r['counts']['original_internal_reference'] for r in rows]==[128,16384,4096]
 assert all(r['n_initial']==16384 for r in j['rows'])
