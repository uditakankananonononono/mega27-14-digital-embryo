"""Multi-cell id-021 port: wild-type stripe-pair fixed point under the right SLP prepattern."""
import sys, pathlib
import pytest
mpbn = pytest.importorskip("mpbn")
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "experiments"))
from segment_polarity_mc2 import build_multicell, pattern_of

def _patterns(slp):
    m = mpbn.MPBooleanNetwork(build_multicell(4, slp))
    return {pattern_of(fp, 4) for fp in m.fixedpoints()}

def test_wildtype_fixed_point_exists():
    assert "WE.." in _patterns((1, 0, 0, 1))

def test_wildtype_needs_slp_off_after_en():
    # SLP on in the cell posterior to en kills the wild-type pair
    assert "WE.." not in _patterns((1, 0, 1, 1))

def test_all_off_is_always_fixed():
    assert "...." in _patterns((0, 0, 0, 0))
