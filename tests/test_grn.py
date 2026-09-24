import numpy as np
from digitalembryo.grn import (count_fixed_points, derrida_sensitivity,
                               boolean_network_derrida, simulate_toggle)


def test_toggle_bistable_two_basins():
    fates = count_fixed_points(a1=2.5, a2=2.5, n=3.0, grid=20)
    # symmetric bistable switch: both fates reachable, neither dominates
    assert fates["high_u"] > 50 and fates["high_v"] > 50
    assert abs(fates["high_u"] - fates["high_v"]) / sum(fates.values()) < 0.3


def test_toggle_trajectory_converges():
    _, U, V = simulate_toggle(0.1, 2.0)
    assert abs(U[-1] - U[-50]) < 1e-3 and abs(V[-1] - V[-50]) < 1e-3
    assert V[-1] > U[-1]  # started high-v, must stay in high-v basin


def test_derrida_measured_vs_theory():
    theory = derrida_sensitivity(K=3, p=0.5)
    measured = boolean_network_derrida(N=12, K=3, p=0.5, n_nets=4, n_pairs=16, seed=7)
    assert abs(measured - theory) < 0.35  # annealed approx is loose but close


def test_ordered_vs_chaotic_regimes():
    assert derrida_sensitivity(K=1) < 1.0   # ordered
    assert derrida_sensitivity(K=4) > 1.0   # chaotic (2*4*0.25 = 2)
