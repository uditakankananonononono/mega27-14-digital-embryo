"""Second-method cross-check of Study B attractors with AEON (biodivine) and mpbn.
Fixed points are update-semantics invariant, so any correct engine must return
the same fixed-point count per condition as our exact synchronous enumeration."""
import json, pathlib, sys
ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from digitalembryo.bnet import parse_bnet, attractor_analysis

bnet_text = (ROOT / "data/raw/grn_models/id021_bodysegmentation.bnet").read_text()
rules = parse_bnet(bnet_text)
committed = json.load(open(ROOT / "results/grn_real.json"))

CONDS = {}
for c in committed["conditions"]:
    hh, wg, slp = (p[-1] for p in c.split("_"))
    CONDS[c] = {"v_hh_external": int(hh), "v_WG_external": int(wg), "v_SLP": int(slp)}

ours = {}
for cond, ext in CONDS.items():
    a = attractor_analysis(rules, ext)
    ours[cond] = len(a["fixed"])
out = {"ours_counts": ours, "committed_counts": {c: committed["conditions"][c]["n_fixed"] for c in CONDS}}
out["ours_match_committed"] = all(out["committed_counts"][c] == ours[c] for c in ours)

try:
    from biodivine_aeon import BooleanNetwork, AsynchronousGraph, Attractors
    aeon_counts = {}
    for cond, ext in CONDS.items():
        rg = BooleanNetwork.from_bnet(bnet_text)
        for k, v in ext.items():
            rg = rg.fix_variable(k, v)
        n_fp = 0
        for a in Attractors.attractors(AsynchronousGraph(rg)):
            if a.is_singleton():
                n_fp += 1
        aeon_counts[cond] = n_fp
    out["aeon_fixed_counts"] = aeon_counts
    out["aeon_agrees"] = all(aeon_counts[c] == ours[c] for c in ours)
except Exception as e:
    out["aeon_error"] = f"{type(e).__name__}: {e}"[:300]

try:
    from mpbn import MPBNet
    mpbn_counts = {}
    for cond, ext in CONDS.items():
        m = MPBNet(bnet_text)
        for k, v in ext.items():
            m[k] = v
        mpbn_counts[cond] = len(list(m.fixedpoints()))
    out["mpbn_fixed_counts"] = mpbn_counts
    out["mpbn_agrees"] = all(mpbn_counts[c] == ours[c] for c in ours)
except Exception as e:
    out["mpbn_error"] = f"{type(e).__name__}: {e}"[:300]

json.dump(out, open(ROOT / "results/grn_crosscheck.json", "w"), indent=1)
print(json.dumps(out, indent=1)[:1500])
