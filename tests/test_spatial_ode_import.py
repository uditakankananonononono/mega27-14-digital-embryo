import hashlib,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
from spatial_ode_import_audit import table

def test_saved_spatial_replay_frozen_checks():
 j=json.loads((ROOT/'results/spatial_ode_import_audit.json').read_text())
 assert hashlib.sha256((ROOT/'experiments/spatial_ode_import_plan.json').read_bytes()).hexdigest()==j['plan_sha256']
 assert len(j['comparisons'])==40 and all(r['within_frozen_threshold'] for r in j['comparisons'])
 assert j['published_argmax_matches']=={'wg':True,'ptc':False,'en':True,'hh':True}
 t=np.loadtxt(ROOT/'results/spatial_ode_import_trajectory.csv',delimiter=',',skiprows=1)
 assert t.shape==(221,21) and np.array_equal(t[:,0],np.arange(221)*5) and np.all(np.isfinite(t))
 assert j['endpoints']['ptc']['argmax_cell_zero_index']==3
 assert abs(j['endpoints']['ptc']['values'][1]-j['endpoints']['ptc']['values'][3])<1e-9
 for i,g in enumerate(['en','ci','hh','ptc','wg']):assert np.array_equal(t[-1,1+i*4:1+(i+1)*4],j['endpoints'][g]['values'])

def test_reference_column_parser(tmp_path):
 p=tmp_path/'ref.tsv';p.write_text('time\ten_0_0\n0\t1\n5\t.5\n')
 t=table(p);assert np.array_equal(t['time'],[0,5]) and np.array_equal(t['en_0_0'],[1,.5])
