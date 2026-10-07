import itertools,sys
from pathlib import Path
import pytest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'experiments'))
from grn_full_reference_audit import terminal_ids
def oracle(graph,start):
 seen=set();cur=start
 while cur not in seen:
  seen.add(cur)
  if graph[cur]==cur:return cur
  cur=graph[cur]
 return -1
def test_all_small_functional_graphs_match_independent_trajectory_oracle():
 total=0
 for n in range(1,6):
  for graph in itertools.product(range(n),repeat=n):
   assert terminal_ids(graph)==[oracle(graph,s) for s in range(n)];total+=1
 assert total==3413
@pytest.mark.parametrize('graph',[[],[-1],[1],[0,2],[0,False],[0,.5]])
def test_invalid_graph_rejected_not_mislabelled(graph):
 with pytest.raises(ValueError):terminal_ids(graph)
def test_numpy_graph_supported():
 import numpy as np
 assert terminal_ids(np.array([1,1,3,2]))==[1,1,-1,-1]
