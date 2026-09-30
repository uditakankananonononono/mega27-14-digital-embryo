"""Archive consistency and probability integrity for enumerated Markov absorption."""
import json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def test_six_closed_classes_and_normalized_absorption():
 j=json.loads((ROOT/'results/flower_async_absorption.json').read_text());a=np.load(ROOT/'results/flower_async_absorption_probabilities.npz');H=a['probabilities']
 assert H.shape==(4096,6);assert j['n_closed_recurrent_classes']==6
 assert all(len(c)==1 for c in j['closed_classes'])
 assert np.max(abs(H.sum(axis=1)-1))<1e-9;assert H.min()>-1e-10
 assert j['linear_residual_max']<1e-10
 for col,c in enumerate(j['closed_classes']):assert np.isclose(H[c[0],col],1)
def test_fixed_start_schedule_ranking_reversal():
 j=json.loads((ROOT/'results/flower_async_absorption.json').read_text());r={x['organ']:x for x in j['rows']}
 assert r['carpel']['synchronous_fixed_one_bit_retention']<r['petal']['synchronous_fixed_one_bit_retention']
 assert r['carpel']['async_fixed_one_bit_retention']>r['petal']['async_fixed_one_bit_retention']
 assert np.isclose(sum(x['async_uniform_initial_absorption_mass'] for x in r.values()),1)
