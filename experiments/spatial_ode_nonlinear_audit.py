"""Paired nonlinear versus linear deviations on two frozen feasible directions."""
import json,hashlib,subprocess
from pathlib import Path
import numpy as np
import roadrunner
from scipy.linalg import expm
ROOT=Path(__file__).resolve().parents[1]
def compute(upstream):
 pth=ROOT/'experiments/spatial_ode_nonlinear_plan.json';p=json.loads(pth.read_text());u=Path(upstream);f=u/'vonDassow2000/vonDassow2000_1x4.timecourse1.xml';assert subprocess.check_output(['git','-C',str(u),'rev-parse','HEAD']).decode().strip()==p['source_commit'];assert hashlib.sha256(f.read_bytes()).hexdigest()==p['sbml_sha256']
 old=ROOT/'results/spatial_ode_stationarity_audit.json';jac=ROOT/'results/spatial_ode_jacobian_fd1e6.csv';cells=json.loads(old.read_text())['residuals'][-1]['species'];x=np.array([s['concentration'] for s in cells]);J=np.loadtxt(jac,delimiter=',');times=[0,1,10,100]
 def trajectory(z):
  r=roadrunner.RoadRunner(str(f));assert list(r.model.getFloatingSpeciesIds())==[s['id'] for s in cells];r.integrator.relative_tolerance=1e-8;r.integrator.absolute_tolerance=1e-14;r.integrator.stiff=False;r.model.setFloatingSpeciesConcentrations(z);r.model.setTime(11000);states=[z.copy()]
  for a,b in zip(times[:-1],times[1:]):
   r.simulate(11000+a,11000+b,2);states.append(np.asarray(r.model.getFloatingSpeciesConcentrations()).copy())
  return states
 base=trajectory(x);rows=[]
 for eps in [0.001,0.01]:
  delta=eps*x;pert=trajectory(x+delta);details=[]
  for t,b,z in zip(times,base,pert):
   observed=z-b;linear=expm(t*J)@delta;details.append({'time':t,'observed_gain':float(np.linalg.norm(observed)/np.linalg.norm(delta)),'linear_gain':float(np.linalg.norm(linear)/np.linalg.norm(delta)),'linear_relative_error':float(np.linalg.norm(observed-linear)/np.linalg.norm(linear)),'minimum_concentration':float(z.min()),'baseline_drift_norm':float(np.linalg.norm(b-x))})
  rows.append({'fractional_increase':eps,'initial_deviation_norm':float(np.linalg.norm(delta)),'times':details})
 return {'plan_sha256':hashlib.sha256(pth.read_bytes()).hexdigest(),'source_commit':p['source_commit'],'sbml_sha256':p['sbml_sha256'],'endpoint_sha256':hashlib.sha256(old.read_bytes()).hexdigest(),'jacobian_sha256':hashlib.sha256(jac.read_bytes()).hexdigest(),'roadrunner_version':roadrunner.__version__,'perturbations':rows,'limits':p['limits']}
if __name__=='__main__':
 import argparse
 a=argparse.ArgumentParser();a.add_argument('--upstream',required=True);v=a.parse_args();j=compute(v.upstream);(ROOT/'results/spatial_ode_nonlinear_audit.json').write_text(json.dumps(j,indent=2)+'\n');print(j)
