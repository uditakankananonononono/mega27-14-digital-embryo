"""Pinned independent spatial SBML trajectory replay, not biological validation."""
import argparse,hashlib,json,subprocess,csv
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def table(path):
 with open(path) as f:
  rows=list(csv.reader(f,delimiter='\t'));names=rows[0];data=np.array([[float(x) for x in r] for r in rows[1:] if r],float)
 if len(set(names))!=len(names) or not np.all(np.isfinite(data)):raise ValueError('reference table')
 return {n:data[:,i] for i,n in enumerate(names)}

def run(upstream):
 import roadrunner
 raw=(ROOT/'experiments/spatial_ode_import_plan.json').read_bytes();p=json.loads(raw);upstream=Path(upstream)
 if subprocess.check_output(['git','-C',str(upstream),'rev-parse','HEAD']).decode().strip()!=p['source_commit']:raise ValueError('source commit')
 base=upstream/'vonDassow2000';sbml=base/p['sbml'];assert hashlib.sha256(sbml.read_bytes()).hexdigest()==p['sbml_sha256']
 names=[f'{g}_0_{c}' for g in ['en','ci','hh','ptc','wg'] for c in range(4)];r=roadrunner.RoadRunner(str(sbml));r.integrator.absolute_tolerance=1e-12;r.integrator.relative_tolerance=1e-6;r.integrator.stiff=False;r.timeCourseSelections=['time']+names;result=np.asarray(r.simulate(0,1100,221));assert result.shape==(221,21) and np.all(np.isfinite(result))
 comp=[];hashes={};reference_ptc_endpoints={};times=np.arange(221)*5
 for source in ['tellurium','copasi']:
  path=base/f'timecourse1/timecourse1.{source}.tsv';ref=table(path);hashes[source]=hashlib.sha256(path.read_bytes()).hexdigest();assert np.array_equal(ref['time' if source=='tellurium' else 'Time'],times)
  reference_ptc_endpoints[source]=[float(ref[f'ptc_0_{c}' if source=='tellurium' else f'[ptc_0,{c}]'][-1]) for c in range(4)]
  for i,name in enumerate(names,1):
   key=name if source=='tellurium' else '['+name.replace('_0_', '_0,')+']'
   expected=ref[key];error=float(np.max(np.abs(result[:,i]-expected)));scale=max(1e-12,float(np.max(np.abs(expected))));comp.append({'source':source,'species':name,'max_abs_error':error,'scale_normalized_max_error':error/scale,'within_frozen_threshold':error/scale<=1e-3})
 endpoints={g:{'values':result[-1,1+k*4:1+(k+1)*4].tolist(),'argmax_cell_zero_index':int(np.argmax(result[-1,1+k*4:1+(k+1)*4]))} for k,g in enumerate(['en','ci','hh','ptc','wg'])}
 expected={'wg':1,'ptc':1,'en':2,'hh':2};out={'plan_sha256':hashlib.sha256(raw).hexdigest(),'source_url':p['source_url'],'source_commit':p['source_commit'],'sbml_sha256':p['sbml_sha256'],'reference_sha256':hashes,'reference_ptc_endpoints':reference_ptc_endpoints,'roadrunner_version':roadrunner.__version__,'samples':221,'species':names,'minimum_concentration':float(result[:,1:].min()),'comparisons':comp,'endpoints':endpoints,'published_argmax_matches':{g:endpoints[g]['argmax_cell_zero_index']==c for g,c in expected.items()},'limits':p['limits']}
 (ROOT/'results/spatial_ode_import_audit.json').write_text(json.dumps(out,indent=2)+'\n');np.savetxt(ROOT/'results/spatial_ode_import_trajectory.csv',result,delimiter=',',header=','.join(['time']+names),comments='');print('comparisons_max',max(c['scale_normalized_max_error'] for c in comp),'fails',sum(not c['within_frozen_threshold'] for c in comp));print('endpoints',endpoints);return out
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--upstream',required=True);x=a.parse_args();run(x.upstream)
