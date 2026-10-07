"""Local finite-state Jacobian evidence, not nonlinear/global stability proof."""
import argparse,hashlib,json,subprocess
from pathlib import Path
import numpy as np
import roadrunner
ROOT=Path(__file__).resolve().parents[1]
def compute(upstream):
 pth=ROOT/'experiments/spatial_ode_jacobian_plan.json';p=json.loads(pth.read_bytes());u=Path(upstream)
 if subprocess.check_output(['git','-C',str(u),'rev-parse','HEAD']).decode().strip()!=p['source_commit']:raise ValueError('source pin')
 f=u/'vonDassow2000'/p['sbml']
 if hashlib.sha256(f.read_bytes()).hexdigest()!=p['sbml_sha256']:raise ValueError('SBML digest')
 r=roadrunner.RoadRunner(str(f));r.integrator.absolute_tolerance=1e-12;r.integrator.relative_tolerance=1e-6;r.integrator.stiff=False;r.simulate(0,1100,221);r.simulate(1100,11000,1981)
 ids=list(r.model.getFloatingSpeciesIds());x=np.asarray(r.model.getFloatingSpeciesConcentrations()).copy();dx=np.asarray(r.getRatesOfChange()).copy();prior=json.loads((ROOT/'results/spatial_ode_stationarity_audit.json').read_text());px=np.array([s['concentration'] for s in prior['residuals'][-1]['species']]);J=np.asarray(r.getFullJacobian());n=len(ids)
 if J.shape!=(n,n):raise ValueError('Jacobian geometry')
 native_finite=bool(np.all(np.isfinite(J)));ev=np.linalg.eigvals(J) if native_finite else np.array([],complex);checks=[];matrices=[]
 r.model.setFloatingSpeciesConcentrations(x);dx=np.asarray(r.getRatesOfChange()).copy()
 for rel in [1e-5,1e-6]:
  F=np.zeros_like(J);forward=0
  for i in range(n):
   h=rel*max(1e-8,abs(x[i]));v=x.copy();v[i]+=h;r.model.setFloatingSpeciesConcentrations(v);a=np.asarray(r.getRatesOfChange()).copy()
   if x[i]>=h:
    v=x.copy();v[i]-=h;r.model.setFloatingSpeciesConcentrations(v);b=np.asarray(r.getRatesOfChange()).copy();F[:,i]=(a-b)/(2*h)
   else:F[:,i]=(a-dx)/h;forward+=1
  r.model.setFloatingSpeciesConcentrations(x)
  if not np.all(np.isfinite(F)):raise ValueError('finite difference nonfinite')
  matrices.append(F);e=float(np.max(np.abs(F-J))) if native_finite else None;ee=np.linalg.eigvals(F);checks.append({'eigenvalues':[{'real':float(z.real),'imag':float(z.imag)} for z in sorted(ee,key=lambda z:z.real,reverse=True)],'positive_modes_gt1e8':int(np.sum(ee.real>1e-8)),'negative_modes_ltminus1e8':int(np.sum(ee.real<-1e-8)),'near_zero_real_modes':int(np.sum(np.abs(ee.real)<=1e-8)),'relative_step':rel,'step_floor_scale':1e-8,'forward_columns_near_zero':forward,'max_abs_jacobian_error':e,'relative_max_error':e/max(1e-12,float(np.max(np.abs(J)))) if native_finite else None,'finite_difference_spectral_abscissa':float(ee.real.max())})
 out={'plan_sha256':hashlib.sha256(pth.read_bytes()).hexdigest(),'source_commit':p['source_commit'],'sbml_sha256':p['sbml_sha256'],'floating_species':n,'species_ids':ids,'prior_endpoint_max_abs_error':float(np.max(np.abs(x-px))),'endpoint_max_abs_derivative':float(np.max(np.abs(dx))),'native_jacobian_finite':native_finite,'native_nonfinite_entries':int(np.sum(~np.isfinite(J))),'native_nonfinite_columns':[ids[i] for i in range(n) if np.any(~np.isfinite(J[:,i]))],'finite_difference_step_max_abs_difference':float(np.max(np.abs(matrices[0]-matrices[1]))),'jacobian_max_abs_entry':float(np.max(np.abs(J))) if native_finite else None,'spectral_abscissa':float(ev.real.max()) if native_finite else None,'positive_modes_gt1e8':int(np.sum(ev.real>1e-8)) if native_finite else None,'negative_modes_ltminus1e8':int(np.sum(ev.real<-1e-8)) if native_finite else None,'near_zero_real_modes':int(np.sum(np.abs(ev.real)<=1e-8)) if native_finite else None,'eigenvalues':[{'real':float(z.real),'imag':float(z.imag)} for z in sorted(ev,key=lambda z:z.real,reverse=True)],'finite_difference_checks':checks,'roadrunner_version':roadrunner.__version__,'limits':p['limits']}
 (ROOT/'results/spatial_ode_jacobian_audit.json').write_text(json.dumps(out,indent=2)+'\n');np.savetxt(ROOT/'results/spatial_ode_jacobian.csv',J,delimiter=',');np.savetxt(ROOT/'results/spatial_ode_jacobian_fd1e5.csv',matrices[0],delimiter=',');np.savetxt(ROOT/'results/spatial_ode_jacobian_fd1e6.csv',matrices[1],delimiter=',');print({k:v for k,v in out.items() if k not in ['species_ids','eigenvalues','finite_difference_checks']});print([{k:v for k,v in c.items() if k!='eigenvalues'} for c in checks]);print('top eigenvalues',out['eigenvalues'][:8]);return out
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--upstream',required=True);x=a.parse_args();compute(x.upstream)
