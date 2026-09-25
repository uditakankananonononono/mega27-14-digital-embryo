#!/usr/bin/env python3
"""Sweep mu_a across the Turing band: measured vs linear wavenumber."""
import json
import multiprocessing as mp
import pathlib

import numpy as np

from digitalembryo.turing_ext import run_point

ROOT = pathlib.Path(__file__).resolve().parent.parent


def _w(mu_a):
    return run_point(0.005, 0.2, float(mu_a), 0.12, n=48, steps=2500)


def main():
    mu_as = np.linspace(0.015, 0.085, 15)
    with mp.Pool(2) as pool:
        rows = [r for r in pool.map(_w, mu_as) if r]
    eps = np.array([r["epsilon"] for r in rows])
    dev = np.array([r["dev"] for r in rows])
    # correction law: dev ~ c * eps^gamma (log-log fit on positive deviations)
    m = dev > 0
    law = None
    if m.sum() >= 4:
        g, c = np.polyfit(np.log(eps[m]), np.log(dev[m]), 1)
        law = {"gamma": float(g), "c": float(np.exp(c)), "n_pos": int(m.sum())}
    out = {"model": "Gierer-Meinhardt 2D, Da=0.005 Dh=0.2 mu_h=0.12 rho=1",
           "n_points": len(rows), "law": law, "rows": rows}
    (ROOT / "results" / "turing_sweep.json").write_text(json.dumps(out, indent=1))
    for r in rows:
        print(f"mu_a={r['mu_a']:.3f} eps={r['epsilon']:+.3f} k2_lin={r['k2_lin']:7.3f} "
              f"k2_sim={r['k2_sim']:7.3f} dev={r['dev']:+.3f}")
    print("law:", law)


if __name__ == "__main__":
    main()
