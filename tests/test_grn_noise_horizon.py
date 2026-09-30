import json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
import grn_noise_horizon_audit as N
import fragility_index as F

def test_vector_lookup_matches_scalar_every_state():
 free,table=N.transition_table();rules=F.parse_bnet(F.BNET)
 for code in range(len(table)):
  state={v:(code>>i)&1 for i,v in enumerate(free)};state.update(F.EXT)
  nxt=F.successor(state,rules,F.EXT)
  assert table[code]==sum(nxt[v]<<i for i,v in enumerate(free))

def test_frozen_plan_and_saved_endpoint_separation():
 import hashlib
 planraw=(ROOT/'experiments/grn_noise_horizon_plan.json').read_bytes();p=json.loads(planraw);j=json.loads((ROOT/'results/grn_noise_horizon_audit.json').read_text())
 assert j['plan_sha256']==hashlib.sha256(planraw).hexdigest()
 assert len(j['rows'])==4*5*3
 free,table=N.transition_table();wg=1<<free.index('v_wg')
 ref=j['original_wg_on_fixed_state_integers'];assert len(ref)==1
 assert all(table[x]==x and x&wg for x in ref)
 for r in j['rows']:
  assert r['terminal_reference_fixed_fraction']<=r['terminal_wg_fraction']
  assert r['early_stop_reference_fixed_fraction']<=r['early_stop_wg_fraction']
  if r['q']==0:assert r['terminal_wg_fraction']==r['early_stop_wg_fraction']==r['terminal_reference_fixed_fraction']
 assert any(r['terminal_wg_fraction']>10*r['terminal_reference_fixed_fraction'] for r in j['rows'] if r['q']==.1)

def test_noiseless_uniform_sample_matches_exact_basin():
 free,table=N.transition_table();states=np.arange(len(table),dtype=np.uint16);cur=states.copy()
 for _ in range(100):cur=table[cur]
 wg=1<<free.index('v_wg');assert int(np.sum((cur&wg)>0))==128
 j=json.loads((ROOT/'results/grn_noise_horizon_audit.json').read_text())
 for rep in range(4):
  rng=np.random.default_rng(930514+rep);cur=rng.integers(0,len(table),4096,dtype=np.uint16)
  for _ in range(100):cur=table[cur]
  row=next(r for r in j['rows'] if r['q']==0 and r['horizon']==100 and r['replicate']==rep)
  assert np.mean((cur&wg)>0)==row['terminal_wg_fraction']
