"""Verify fold/global metric weighting and missing ridge intercept on saved cache."""
import json,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def feats(r):
 c=np.array(r['curve']);n=len(c);return [np.argmax(c)/(n-1),np.mean(c>=.5*c.max()),c.max(),np.sum(c*np.arange(n))/(c.sum()+1e-12)/(n-1)]
def metric(p,y):return {'relative_rmse':float(np.sqrt(np.mean(((p-y)/(np.abs(y)+1e-9))**2))),'ordinary_rmse':float(np.sqrt(np.mean((p-y)**2)))}
def run():
 raw=(ROOT/'results/turing_sims_cache.json').read_bytes();rows=json.loads(raw);planraw=(ROOT/'experiments/turing_metric_plan.json').read_bytes();plan=json.loads(planraw)
 y=np.array([(r['k2_sim']-r['k2_min'])/(r['k2_max']-r['k2_min']) for r in rows]);linear=np.array([(r['k2_lin']-r['k2_min'])/(r['k2_max']-r['k2_min']) for r in rows]);F=np.array([feats(r) for r in rows]);out=[]
 for seed in plan['seeds']:
  fold=np.array_split(np.random.default_rng(seed).permutation(len(y)),6);preds={k:np.zeros(len(y)) for k in ['linear','train_mean','ridge_no_intercept','ridge_with_intercept']};folds=[]
  for i,te in enumerate(fold):
   tr=np.setdiff1d(np.arange(len(y)),te);mu,sd=F[tr].mean(0),F[tr].std(0)+1e-9;Z=(F[tr]-mu)/sd;T=(F[te]-mu)/sd
   w0=np.linalg.solve(Z.T@Z+np.eye(4),Z.T@y[tr]);w1=np.linalg.solve(Z.T@Z+np.eye(4),Z.T@(y[tr]-y[tr].mean()))
   for k,p in [('linear',linear[te]),('train_mean',np.full(len(te),y[tr].mean())),('ridge_no_intercept',T@w0),('ridge_with_intercept',T@w1+y[tr].mean())]:
    preds[k][te]=p;folds.append({'fold':i,'model':k,'n':len(te),**metric(p,y[te])})
  for k,p in preds.items():out.append({'seed':seed,'model':k,'global':metric(p,y),'equal_fold_mean':{m:float(np.mean([r[m] for r in folds if r['model']==k])) for m in ['relative_rmse','ordinary_rmse']}})
 result={'cache_sha256':hashlib.sha256(raw).hexdigest(),'plan_sha256':hashlib.sha256(planraw).hexdigest(),'n':len(y),'target_min':float(y.min()),'target_max':float(y.max()),'target_mean':float(y.mean()),'rows':out,'limits':plan['limits']};(ROOT/'results/turing_metric_audit.json').write_text(json.dumps(result,indent=2)+'\n')
 for r in out:print(r)
if __name__=='__main__':run()
