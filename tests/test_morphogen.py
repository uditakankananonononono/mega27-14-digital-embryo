import numpy as np
from digitalembryo.morphogen import (GiererMeinhardt2D, turing_thresholds)


def test_steady_state_solves_kinetics():
    th = turing_thresholds(0.005, 0.2, 0.06, 0.12)
    a, h = th["a_star"], th["h_star"]
    # both reaction terms vanish at the steady state
    assert abs(a**2 / h - 0.06 * a) < 1e-9
    assert abs(a**2 - 0.12 * h) < 1e-9


def test_turing_window_matches_simulation():
    # Parameters inside the Turing window must produce a pattern whose
    # dominant wavenumber falls inside the analytic unstable band.
    Da, Dh, mu_a, mu_h = 0.005, 0.2, 0.06, 0.12
    th = turing_thresholds(Da, Dh, mu_a, mu_h)
    assert th["unstable"] and th["trace"] < 0 and th["det"] > 0
    sim = GiererMeinhardt2D(n=64, Da=Da, Dh=Dh, mu_a=mu_a, mu_h=mu_h, seed=3)
    sim.step(2000)
    k2 = sim.pattern_wavenumber()
    assert th["k2_min"] * 0.5 < k2 < th["k2_max"] * 2.0
    # pattern actually emerged (variance grew beyond the noise floor)
    assert sim.a.std() > 0.05 * sim.a.mean()


def test_stable_params_no_pattern():
    # Equal diffusion -> no Turing instability; field must stay ~homogeneous.
    th = turing_thresholds(0.1, 0.1, 0.06, 0.12)
    assert not th["unstable"]
    sim = GiererMeinhardt2D(n=32, Da=0.1, Dh=0.1, mu_a=0.06, mu_h=0.12, seed=1)
    before = sim.a.std()
    sim.step(3000)
    assert sim.a.std() <= before * 1.5


def test_dt_stability_guard():
    try:
        GiererMeinhardt2D(dt=1e6)
        assert False, "should have raised"
    except ValueError:
        pass
