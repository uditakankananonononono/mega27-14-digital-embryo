import pathlib

from digitalembryo.bnet import attractor_analysis, parse_bnet, successor

BNET = (pathlib.Path(__file__).parent.parent / "data/raw/grn_models/id021_bodysegmentation.bnet").read_text()
RULES = parse_bnet(BNET)
EXT = {"v_hh_external": 0, "v_WG_external": 0, "v_SLP": 0}


def test_parse_all_vars():
    assert len(RULES) == 14
    assert "v_wg" in RULES and "v_CIA" in RULES


def test_fixed_points_are_fixed():
    a = attractor_analysis(RULES, EXT)
    assert len(a["fixed"]) >= 1
    for fp in a["fixed"]:
        nxt = successor(fp, RULES, EXT)
        free = [v for v in RULES if v not in EXT]
        assert all(nxt[v] == fp[v] for v in free)


def test_bistability_condition():
    # published: SLP=0 with WG signal gives 2 stable states (bistability)
    a = attractor_analysis(RULES, {"v_hh_external": 0, "v_WG_external": 1, "v_SLP": 0})
    assert len(a["fixed"]) == 2


def test_no_cycles_sync():
    for hh in (0, 1):
        for wg in (0, 1):
            for slp in (0, 1):
                a = attractor_analysis(RULES, {"v_hh_external": hh, "v_WG_external": wg, "v_SLP": slp})
                assert a["cycles"] == []
