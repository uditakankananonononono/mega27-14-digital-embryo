"""Exact finite-state random-single-node asynchronous absorption probabilities.
Sparse linear solves, not trajectory sampling. Update rule is not biology-calibrated.
"""
import sys,itertools,json,hashlib
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix,eye
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import spsolve
sys.path.insert(0,str(Path(__file__).parent));import flower_model as F
ROOT=Path(__file__).resolve().parents[1]
def main():
 n=12;N=1<<n;states=np.array(list(itertools.product([0,1],repeat=n)),np.uint8);powers=1<<np.arange(n-1,-1,-1);updated=F.step_batch(states);ids=np.arange(N);dest=np.empty((N,n),int)
 for node in range(n):dest[:,node]=ids+(updated[:,node].astype(int)-states[:,node])*powers[node]
 P=coo_matrix((np.full(N*n,1/n),(np.repeat(ids,n),dest.ravel())),shape=(N,N)).tocsr();assert np.allclose(P.sum(axis=1),1)
 count,labels=connected_components(P,directed=True,connection='strong');closed=np.ones(count,bool)
 for source in range(N):
  for target in dest[source]:
   if labels[source]!=labels[target]:closed[labels[source]]=False
 classes=[ids[labels==label].tolist() for label in np.where(closed)[0]];recurrent=np.array(sorted(x for c in classes for x in c));transient=np.setdiff1d(ids,recurrent)
 fixed={name:int(bits,2) for name,bits in F.PUBLISHED_FP.items()};assert all([v] in classes for v in fixed.values())
 Q=P[transient][:,transient];B=np.zeros((N,len(classes)));rhs=np.zeros((len(transient),len(classes)))
 for col,c in enumerate(classes):B[c,col]=1;rhs[:,col]=np.asarray(P[transient][:,c].sum(axis=1)).ravel()
 B[transient]=spsolve(eye(len(transient),format='csr')-Q,rhs)
 residual=float(np.max(np.abs(B-P@B)));assert residual<1e-10 and np.allclose(B.sum(axis=1),1,atol=1e-9);assert B.min()>-1e-10
 sync=updated@powers;sl=[]
 for start in ids:
  seen={};cur=int(start)
  while cur not in seen:seen[cur]=len(seen);cur=int(sync[cur])
  sl.append(tuple(sorted(s for s,i in seen.items() if i>=seen[cur])))
 rows=[]
 for name,f in fixed.items():
  col=classes.index([f]);neighbors=f^powers;members=np.array([i for i,a in enumerate(sl) if a==(f,)])
  stay=B[neighbors,col];retained=np.mean([B[x^powers,col].mean() for x in members])
  rows.append({'organ':name,'fixed_state':F.PUBLISHED_FP[name],'async_fixed_one_bit_retention':float(stay.mean()),'async_fixed_one_bit_nodewise':dict(zip(F.NODES,map(float,stay))),'synchronous_fixed_one_bit_retention':float(np.mean([sl[x]==(f,) for x in neighbors])),'synchronous_basin_size':len(members),'async_probability_from_uniform_synchronous_basin':float(B[members,col].mean()),'async_one_bit_from_uniform_synchronous_basin':float(retained),'async_uniform_initial_absorption_mass':float(B[:,col].mean())})
 out={'source':'https://pageperso.lis-lab.fr/~sylvain.sene/files/publi_pres/rgs18.pdf','model_code_sha256':hashlib.sha256((ROOT/'experiments/flower_model.py').read_bytes()).hexdigest(),'update':'At each step select one of12nodes uniformly, apply published threshold to that node only; self loops allowed','n_states':N,'n_closed_recurrent_classes':len(classes),'closed_classes':classes,'linear_residual_max':residual,'probability_row_sum_max_error':float(np.max(abs(B.sum(axis=1)-1))),'rows':rows,'limits':['Random single-node update is a modeling choice, not empirically calibrated timing','Uniform synchronous basin starts are held identical to compare kernels, not an asynchronous basin definition','Stochastic schedule absorption probability differs from deterministic basin size','No causal topology or physical gene perturbation claim']}
 (ROOT/'results/flower_async_absorption.json').write_text(json.dumps(out,indent=2)+'\n');np.savez_compressed(ROOT/'results/flower_async_absorption_probabilities.npz',probabilities=B,states=states)
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
