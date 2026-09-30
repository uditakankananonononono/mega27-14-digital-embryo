import json,sys,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
from turing_metric_audit import run

def test_metric_audit_reproduces_and_exposes_intercept():
 before=(ROOT/'results/turing_metric_audit.json').read_bytes();run();assert (ROOT/'results/turing_metric_audit.json').read_bytes()==before
 j=json.loads(before);assert j['n']==171
 assert j['plan_sha256']==hashlib.sha256((ROOT/'experiments/turing_metric_plan.json').read_bytes()).hexdigest()
 for seed in [0,1]:
  r={x['model']:x for x in j['rows'] if x['seed']==seed}
  assert r['ridge_with_intercept']['global']['ordinary_rmse']<r['linear']['global']['ordinary_rmse']
  assert r['ridge_with_intercept']['global']['relative_rmse']>r['linear']['global']['relative_rmse']
  assert r['ridge_no_intercept']['equal_fold_mean']['relative_rmse']>1.23
  assert r['train_mean']['global']['ordinary_rmse']<r['linear']['global']['ordinary_rmse']
