import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def test_comparator_and_mechanism_scope():
 text=(ROOT/'paper/sections/conclusion.tex').read_text();info=(ROOT/'paper/sections/info_threshold.tex').read_text()
 assert 'not pattern images' in text and 'does not identify whether a living embryo' in text
 assert '2.8566' in info and 'not a minimum biological information threshold' in info
 j=json.loads((ROOT/'results/cf_decoder.json').read_text())
 assert j['models']['liu_threshold_lnA']['rmse_pctEL']<4.29<j['models']['liu_threshold_trueD']['rmse_pctEL']
 ledger=json.loads((ROOT/'results/manuscript_claim_consistency.json').read_text())
 assert len(ledger)>=10
