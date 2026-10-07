"""Unchanged selected spatial ODE setting: finite residual, not asymptotic proof."""
import argparse,hashlib,json,subprocess
from pathlib import Path
import numpy as np
import roadrunner
ROOT=Path(__file__).resolve().parents[1]
def compute(upstream):
 pth=ROOT/'experiments/spatial_ode_stationarity_plan.json';p=json.loads(pth.read_bytes());files=list(Path(upstream).rglob(p['sbml']))
 if subprocess.check_output(['git','-C',str(upstream),'rev-parse','HEAD']).decode().strip()!=p['source_commit']:raise ValueError('source commit')
 if len(files)!=1:raise ValueError('source path ambiguity')
 f=files[0]
 if hashlib.sha256(f.read_bytes()).hexdigest()!=p['sbml_sha256']:raise ValueError('source digest')
 names=[f'{g}_0_{c}' for g in ['en','ci','hh','ptc','wg'] for c in range(4)];r=roadrunner.RoadRunner(str(f));r.integrator.absolute_tolerance=1e-12;r.integrator.relative_tolerance=1e-6;r.integrator.stiff=False;r.timeCourseSelections=['time']+names
 a=np.asarray(r.simulate(0,1100,221));snap=[]
 def residual(t):
  ids=list(r.model.getFloatingSpeciesIds());x=np.asarray(r.model.getFloatingSpeciesConcentrations());dx=np.asarray(r.getRatesOfChange());return {'time':t,'floating_species':len(ids),'max_abs_derivative':float(np.max(np.abs(dx))),'max_concentration_scaled_derivative':float(np.max(np.abs(dx)/np.maximum(1e-12,np.abs(x)))),'species':[{'id':i,'concentration':float(v),'derivative':float(d)} for i,v,d in zip(ids,x,dx)]}
 snap.append(residual(1100));b=np.asarray(r.simulate(1100,11000,1981));snap.append(residual(11000));trajectory=np.vstack([a,b[1:]])
 if trajectory.shape!=(2201,21) or not np.all(np.isfinite(trajectory)):raise ValueError('trajectory geometry/nonfinite')
 prior=np.loadtxt(ROOT/'results/spatial_ode_import_trajectory.csv',delimiter=',',skiprows=1)
 endpoints={g:{'values_1100':a[-1,1+k*4:1+(k+1)*4].tolist(),'values_11000':b[-1,1+k*4:1+(k+1)*4].tolist(),'argmax_1100':int(np.argmax(a[-1,1+k*4:1+(k+1)*4])),'argmax_11000':int(np.argmax(b[-1,1+k*4:1+(k+1)*4]))} for k,g in enumerate(['en','ci','hh','ptc','wg'])}
 j={'plan_sha256':hashlib.sha256(pth.read_bytes()).hexdigest(),'sbml_sha256':p['sbml_sha256'],'source_commit':p['source_commit'],'prior_overlap_max_abs_error':float(np.max(np.abs(a-prior))),'residuals':snap,'endpoint_mRNA_max_abs_drift':float(np.max(np.abs(a[-1,1:]-b[-1,1:]))),'endpoints':endpoints,'ptc_cell2_minus_cell4_1100':float(a[-1,14]-a[-1,16]),'ptc_cell2_minus_cell4_11000':float(b[-1,14]-b[-1,16]),'stationary_at_11000_under_frozen_threshold':snap[-1]['max_abs_derivative']<=1e-8,'roadrunner_version':roadrunner.__version__,'limits':p['limits']}
 (ROOT/'results/spatial_ode_stationarity_audit.json').write_text(json.dumps(j,indent=2)+'\n');np.savetxt(ROOT/'results/spatial_ode_stationarity_trajectory.csv',trajectory,delimiter=',',header=','.join(['time']+names),comments='');print({k:v for k,v in j.items() if k!='residuals'});print([(s['time'],s['max_abs_derivative'],s['max_concentration_scaled_derivative']) for s in snap]);return j
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--upstream',required=True);x=a.parse_args();compute(x.upstream)
