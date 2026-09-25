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
    args = p.parse_args(argv)
    args.fn(args)


if __name__ == "__main__":
    main()
