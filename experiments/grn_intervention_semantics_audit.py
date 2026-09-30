"""Show historical clamp function changes only initial measure; no result overwrite."""
import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'experiments'))
import fragility_index as F

def main():
 rules=F.parse_bnet(F.BNET);free=[v for v in rules if v not in F.EXT]
 s={v:0 for v in free};s.update(F.EXT);s['v_en']=1
 nxt=F.successor(s,rules,F.EXT)
 assert nxt['v_en']==0
 old=json.loads((ROOT/'results/grn_fragility.json').read_text())
 ones=[k for k,v in old['per_node_exit_fraction'].items() if v==1];zeros=[k for k,v in old['per_node_exit_fraction'].items() if v==0]
 assert len(ones)==len(zeros)==7
 j={'historical_source':'experiments/fragility_index.py: wt_basin initializes clamp but attractor_of updates without it','node':'v_en','intended_clamp':1,'initial_state':s,'first_unclamped_successor':nxt,'maintained_clamp_first_successor':dict(nxt,v_en=1),'full_exit_one_nodes':ones,'full_exit_zero_nodes':zeros,'interpretation':'Historical knockout_predictions are conditional initial-state outcomes, not maintained-clamp experiments. Old numerical results preserved.','noise_scope':'update_rule_sensitivity.py stops on an unchanged step or at 300 sweeps, then labels terminal wg. No eventual stochastic absorption guarantee.'}
 (ROOT/'results/grn_intervention_semantics_audit.json').write_text(json.dumps(j,indent=2)+'\n');print(json.dumps(j,indent=2))
if __name__=='__main__':main()
