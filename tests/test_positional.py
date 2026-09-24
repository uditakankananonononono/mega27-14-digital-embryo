import numpy as np
from digitalembryo.positional import (exponential_gradient,
                                      french_flag_boundaries,
                                      positional_error_vs_noise,
                                      mutual_information_positions)


def test_gradient_shape():
    x = np.linspace(0, 1, 50)
    M = exponential_gradient(x, M0=1.0, lam=0.2)
    assert np.all(np.diff(M) < 0)
    assert abs(M[0] - 1.0) < 1e-9
    assert abs(np.log(M[10] / M[0]) + x[10] / 0.2) < 1e-9


def test_flag_boundaries_monotone():
    b = french_flag_boundaries(thresholds=(0.5, 0.2))
    assert b[0] < b[1]  # higher threshold closer to source
    assert abs(b[0] - 0.2 * np.log(2)) < 0.02


def test_error_grows_with_noise_and_distance():
    rows = positional_error_vs_noise(noises=(0.05, 0.2), trials=200, seed=2)
    r = {(d["noise"], d["threshold"]): d for d in rows}
    assert r[(0.2, 0.5)]["x_std"] > r[(0.05, 0.5)]["x_std"]
    assert r[(0.05, 0.2)]["x_std"] > r[(0.05, 0.5)]["x_std"]


def test_mutual_information_positive_bounded():
    mi = mutual_information_positions(trials=400, seed=4)
    assert 0.5 < mi < 3.0  # bits; bounded by log2(8 readout bins)
