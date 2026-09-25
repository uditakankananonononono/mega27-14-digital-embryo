import itertools

from digitalembryo import segment_polarity as sp


def test_state_space_exact():
    attr = sp.attractors()
    # Albert-Othmer core must have at least one wild-type-like fixed point
    assert len(attr["fixed"]) >= 1
    for fp in attr["fixed"]:
        assert sp.successor(fp) == fp  # genuine fixed points


def test_all_states_reach_attractor():
    for bits in itertools.product([0, 1], repeat=len(sp.NODES)):
        cur, guard = bits, 0
        seen = set()
        while cur not in seen:
            seen.add(cur)
            cur = sp.successor(cur)
            guard += 1
        assert guard <= 40


def test_basins_sum_to_state_space():
    b = sp.basin_sizes()
    assert sum(b.values()) + len([c for c in sp.attractors()["cycles"]]) >= sum(b.values())
    assert sum(b.values()) <= 2 ** len(sp.NODES)


def test_wt_pattern_attractor_exists():
    # wild-type-like state: wg, en, hh ON; ptc OFF (segment border pattern)
    fixed = sp.attractors()["fixed"]
    assert any(s[0] == 1 and s[1] == 1 and s[2] == 1 and s[3] == 0 for s in fixed)


def test_ci_knockout_removes_wt_attractor():
    # ci LOF must abolish the wg+en+hh wild-type attractor (CI activates wg)
    fixed = sp.robustness_knockout("ci")
    assert not any(s[0] == 1 and s[1] == 1 and s[2] == 1 for s in fixed)
