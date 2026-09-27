"""embryosim - CLI tool for digital-embryo analysis on real data.

Subcommands:
  fit-bcd      Fit SDD exponential to real Bcd gradients (Liu 2013 dataset)
  hb-boundary  Hb boundary positions vs Bcd dosage
  turing       Run a Gierer-Meinhardt pattern vs its analytic Turing band
  attractors   Exact attractor enumeration of a .bnet GRN model
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]


def cmd_fit_bcd(args):
    from .bcd import fit_single, load_live_gradients

    recs = load_live_gradients()
    if args.line:
        recs = [r for r in recs if r["line"] == args.line]
    lams = []
    for r in recs:
        m = r["x"] < 0.6
        if m.sum() >= 15:
            lams.append(fit_single(r["x"][m], r["intensity"][m])["lam"])
    lams = np.array(lams)
    out = {"n_gradients": len(lams), "lambda_median_EL": float(np.median(lams)),
           "lambda_iqr": [float(np.percentile(lams, 25)), float(np.percentile(lams, 75))],
           "dataset": "Liu et al. PNAS 2013 (Zenodo 4942019)"}
    print(json.dumps(out, indent=2))


def cmd_hb_boundary(args):
    from .gapgene import hb_boundaries, load_if_sessions

    rows = []
    for s in load_if_sessions():
        if "Hb" not in s["genes"]:
            continue
        b = hb_boundaries(s)
        b = b[np.isfinite(b)]
        if len(b) >= 5:
            rows.append({"line": s["line"], "n": len(b),
                         "hb_boundary_median_EL": float(np.median(b))})
    print(json.dumps(rows, indent=2))


def cmd_turing(args):
    from .morphogen import GiererMeinhardt2D, turing_thresholds

    th = turing_thresholds(args.da, args.dh, args.mu_a, args.mu_h)
    out = {"turing_unstable": th["unstable"]}
    if th["unstable"]:
        sim = GiererMeinhardt2D(n=args.n, Da=args.da, Dh=args.dh,
                                mu_a=args.mu_a, mu_h=args.mu_h, seed=args.seed)
        sim.step(args.steps)
        k2 = float(sim.pattern_wavenumber())
        out.update({"analytic_band_k2": [th["k2_min"], th["k2_max"]],
                    "simulated_dominant_k2": k2,
                    "inside_band": bool(th["k2_min"] <= k2 <= th["k2_max"])})
    print(json.dumps(out, indent=2))


def cmd_attractors(args):
    from .bnet import attractor_analysis, parse_bnet

    text = pathlib.Path(args.model).read_text()
    rules = parse_bnet(text)
    ext = {}
    for kv in args.fix or []:
        k, v = kv.split("=")
        ext[k] = int(v)
    a = attractor_analysis(rules, ext)
    print(json.dumps({"n_states": a["n_states"], "n_fixed": len(a["fixed"]),
                      "n_cycles": len(a["cycles"]), "fixed_points": a["fixed"]},
                     indent=2))


def cmd_fragility(args):
    """Developmental fragility index on a .bnet GRN model."""
    from .bnet import parse_bnet
    import itertools, json as _json
    import numpy as _np
    rules = parse_bnet(pathlib.Path(args.model).read_text())
    ext = {}
    for pin in (args.fix or []):
        k, v = pin.split("="); ext[k] = int(v)
    free = [n for n in rules if n not in ext]
    from .bnet import successor as _succ
    def step(state):
        ns = _succ(state, rules, ext)
        return ns
    # enumerate
    n = len(free)
    succ = {}
    for bits in itertools.product([0, 1], repeat=n):
        st = dict(ext); st.update(dict(zip(free, bits)))
        ns = step(st)
        succ[tuple(bits)] = tuple(ns[f] for f in free)
    def root_of(s):
        seen = {}
        cur = s
        while cur not in seen:
            seen[cur] = len(seen)
            cur = succ[cur]
        return cur
    # WT = wg-ON fixed point
    wt = None
    for s in succ:
        if succ[s] == s:
            st = dict(zip(free, s))
            if st.get("v_wg", 0) == 1:
                wt = s
    if wt is None:
        print(_json.dumps({"error": "no wg-ON fixed point"})); return
    basin = [s for s in succ if root_of(s) == wt]
    exits = 0; total = 0
    for s in basin:
        for b in range(n):
            t = list(s); t[b] ^= 1; t = tuple(t)
            total += 1
            if root_of(t) != wt:
                exits += 1
    print(_json.dumps({"model": args.model, "n_states": 2 ** n,
                       "wt_basin_states": len(basin),
                       "wt_basin_fraction": len(basin) / 2 ** n,
                       "fragility_index": exits / total if total else None},
                      indent=2))


def cmd_info_threshold(args):
    """Minimum-information threshold: decode CF position from degraded Bcd."""
    from .decoder import embryo_table, grouped_folds
    import numpy as _np
    import json as _json
    rows = embryo_table()
    X = _np.array([r["profile"] for r in rows])
    y = _np.array([r["cf"] for r in rows])
    folds = grouped_folds(rows, 5)
    rmses = []
    for tr, te in folds:
        Xm = _np.column_stack([_np.ones(len(X)), X])
        w = _np.linalg.solve(Xm[tr].T @ Xm[tr] + 1.0 * _np.eye(Xm.shape[1]), Xm[tr].T @ y[tr])
        pred = Xm[te] @ w
        rmses.append(float(_np.sqrt(_np.mean((pred - y[te]) ** 2))))
    print(_json.dumps({"full_gradient_rmse_pct_EL_mean": 100 * float(_np.mean(rmses)),
                       "liu_2013_benchmark_pct_EL": 4.94}, indent=2))


def main(argv=None):
    p = argparse.ArgumentParser(prog="embryosim", description=__doc__)
    sub = p.add_subparsers(required=True)
    f = sub.add_parser("fit-bcd", help="fit SDD exponential to real Bcd gradients")
    f.add_argument("--line", default=None, help="fly line name, e.g. 2XA")
    f.set_defaults(fn=cmd_fit_bcd)
    h = sub.add_parser("hb-boundary", help="Hb boundary vs dosage")
    h.set_defaults(fn=cmd_hb_boundary)
    t = sub.add_parser("turing", help="Turing pattern vs analytic band")
    t.add_argument("--da", type=float, default=0.005)
    t.add_argument("--dh", type=float, default=0.2)
    t.add_argument("--mu-a", type=float, default=0.06)
    t.add_argument("--mu-h", type=float, default=0.12)
    t.add_argument("--n", type=int, default=48)
    t.add_argument("--steps", type=int, default=2500)
    t.add_argument("--seed", type=int, default=0)
    t.set_defaults(fn=cmd_turing)
    a = sub.add_parser("attractors", help="exact attractors of a .bnet model")
    a.add_argument("--model", default=str(ROOT / "data/raw/grn_models/id021_bodysegmentation.bnet"))
    a.add_argument("--fix", nargs="*", help="external pins like v_hh_external=1")
    a.set_defaults(fn=cmd_attractors)
    fr = sub.add_parser("fragility", help="developmental fragility index of a .bnet model")
    fr.add_argument("--model", default=str(ROOT / "data/raw/grn_models/id021_bodysegmentation.bnet"))
    fr.add_argument("--fix", nargs="*", default=["v_hh_external=1", "v_WG_external=0", "v_SLP=0"])
    fr.set_defaults(fn=cmd_fragility)
    it = sub.add_parser("info-threshold", help="minimum-information decode of CF position")
    it.set_defaults(fn=cmd_info_threshold)
    args = p.parse_args(argv)
    args.fn(args)


if __name__ == "__main__":
    main()
