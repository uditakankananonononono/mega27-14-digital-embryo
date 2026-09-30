"""Paired first-unchanged stopping vs continuous-noise finite-horizon outcomes."""
import hashlib,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'experiments'))
import fragility_index as F

def transition_table():
 rules=F.parse_bnet(F.BNET);free=[v for v in rules if v not in F.EXT];n=len(free);table=[]
 for state in range(1<<n):
  s={v:(state>>i)&1 for i,v in enumerate(free)};s.update(F.EXT);nxt=F.successor(s,rules,F.EXT)
  table.append(sum(nxt[v]<<i for i,v in enumerate(free)))
 return free,np.array(table,dtype=np.uint16)

def main():
 planraw=(ROOT/'experiments/grn_noise_horizon_plan.json').read_bytes();plan=json.loads(planraw)
 free,table=transition_table();n=len(free);wg_bit=1<<free.index('v_wg');states=np.arange(len(table),dtype=np.uint16)
 reference=states[(table==states)&((states&wg_bit)>0)];assert len(reference)>0
 rows=[]
 for rep in range(plan['replicates']):
  initial=np.random.default_rng(plan['seed']+rep).integers(0,len(table),4096,dtype=np.uint16)
  for q in plan['q']:
   rng=np.random.default_rng(plan['seed']+100+rep);cur=initial.copy();early=initial.copy();stopped=np.zeros(len(cur),dtype=bool);first=np.zeros(len(cur),dtype=np.int32)
   for step in range(1,max(plan['horizons'])+1):
    # Identical uniforms across q within a replicate; no early RNG termination.
    flip=(rng.random((len(cur),n))<q);mask=np.sum(flip*(1<<np.arange(n)),axis=1,dtype=np.uint16)
    nxt=table[cur]^mask;new=(nxt==cur)&~stopped;early[new]=nxt[new];first[new]=step;stopped|=new;cur=nxt
    if step in plan['horizons']:
     early_or_terminal=np.where(stopped,early,cur)
     row={'replicate':rep,'q':q,'horizon':step,'n_paths':len(cur),'terminal_wg_fraction':float(np.mean((cur&wg_bit)>0)),'terminal_reference_fixed_fraction':float(np.mean(np.isin(cur,reference))),'early_stop_wg_fraction':float(np.mean((early_or_terminal&wg_bit)>0)),'early_stop_reference_fixed_fraction':float(np.mean(np.isin(early_or_terminal,reference))),'first_unchanged_seen_fraction':float(np.mean(stopped)),'median_first_unchanged_if_seen':float(np.median(first[stopped])) if stopped.any() else None};rows.append(row)
 keys=['terminal_wg_fraction','terminal_reference_fixed_fraction','early_stop_wg_fraction','early_stop_reference_fixed_fraction']
 summary=[{'q':q,'horizon':h,'replicate_means':{k:float(np.mean([r[k] for r in rows if r['q']==q and r['horizon']==h])) for k in keys},'replicate_population_sds':{k:float(np.std([r[k] for r in rows if r['q']==q and r['horizon']==h])) for k in keys}} for q in plan['q'] for h in plan['horizons']]
 out={'summary':summary,'plan_sha256':hashlib.sha256(planraw).hexdigest(),'source_bnet_sha256':hashlib.sha256(F.BNET.encode()).hexdigest(),'node_order':free,'original_wg_on_fixed_state_integers':reference.tolist(),'rows':rows,'limits':['Four Monte Carlo repeats of 4096 iid starts; repeated horizons and q are dependent','Reference fixed-point identity is a model endpoint, not verified biological WT','Continuing noise q>0 never makes an unchanged step permanently absorbing','Historical grn_update_rules.json preserved; new random initial measure differs from its systematic quarter sample']}
 (ROOT/'results/grn_noise_horizon_audit.json').write_text(json.dumps(out,indent=2)+'\n')
 for q in plan['q']:
  for h in plan['horizons']:
   x=[r for r in rows if r['q']==q and r['horizon']==h]
   print(q,h,{k:round(float(np.mean([r[k] for r in x])),6) for k in ['terminal_wg_fraction','terminal_reference_fixed_fraction','early_stop_wg_fraction','early_stop_reference_fixed_fraction','first_unchanged_seen_fraction']},flush=True)
if __name__=='__main__':main()
