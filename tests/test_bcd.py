"""Hermetic tests for the real-data Bicoid module (fixture = 3 real embryos)."""
import numpy as np
import pytest

from digitalembryo import bcd

FIX = np.load(__file__.replace("test_bcd.py", "fixtures/bcd_gradients.npz"))


def _grad(i):
    x, y = FIX[f"g{i}_x"], FIX[f"g{i}_y"]
    m = x < 0.6
    return x[m], y[m]


def test_single_exponential_recovers_published_length_scale():
    # Gregor 2007 / Liu 2013: lambda ~ 0.15-0.25 egg lengths
    for i in range(3):
        x, y = _grad(i)
        s = bcd.fit_single(x, y)
        assert 0.05 < s["lam"] < 0.40, s
        assert s["A"] > 0


def test_fit_is_deterministic():
    x, y = _grad(0)
    a = bcd.fit_single(x, y)
    b = bcd.fit_single(x, y)
    assert a == b


def test_double_fit_nests_single():
    x, y = _grad(1)
    s = bcd.fit_single(x, y)
    d = bcd.fit_double(x, y)
    assert d is not None
    assert d["rss"] <= s["rss"] * 1.0001  # nested model must not fit worse
    assert d["l1"] <= d["l2"]


def test_f_test_bounds():
    x, y = _grad(2)
    s = bcd.fit_single(x, y)
    d = bcd.fit_double(x, y)
    p = bcd.f_test(s["rss"], d["rss"], 2, 4, len(x))
    assert 0.0 <= p <= 1.0


def test_bootstrap_ci_contains_mle():
    x, y = _grad(0)
    s = bcd.fit_single(x, y)
    lo, hi, mean = bcd.bootstrap_lambda(x, y, n_boot=50, rng=0)
    assert lo < s["lam"] < hi
    assert lo > 0
