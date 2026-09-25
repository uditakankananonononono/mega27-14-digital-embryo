import numpy as np

from digitalembryo.gapgene import AXIS, boundary_position, load_if_sessions


def test_sessions_have_real_embryos():
    ss = load_if_sessions()
    assert len(ss) == 11
    assert sum(len(s["embryos"]) for s in ss) == 2180


def test_boundary_recovery_synthetic():
    p = 1.0 / (1 + np.exp((AXIS - 0.45) * 60))  # sharp step at 0.45
    assert abs(boundary_position(p) - 0.45) < 0.01


def test_boundary_rejects_flat():
    assert not np.isfinite(boundary_position(np.ones_like(AXIS)))


def test_hb_boundary_increases_with_dose():
    from digitalembryo.gapgene import hb_boundaries
    meds = {}
    for s in load_if_sessions():
        if "Hb" in s["genes"] and s["line"] in ("1XA", "2XA2IIIA"):
            b = hb_boundaries(s)
            meds.setdefault(s["line"], []).append(np.nanmedian(b))
    assert np.mean(meds["2XA2IIIA"]) > np.mean(meds["1XA"])
