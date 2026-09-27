#!/usr/bin/env python3
"""Boolean vs continuous ODE (verdict #4): do the attractor conclusions survive?

Hill-function ODE counterpart of the Albert-Othmer 5-node core
(activations AND -> product, OR -> sum capped; degradation linear).
Scan initial conditions on a grid; classify the wg-ON attractor basin
fraction and compare with the Boolean basin under the same core.
"""
import json
import pathlib

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
NODES = ["wg", "en", "hh", "ptc", "ci"]
IDX = {n: i for i, n in enumerate(NODES)}
# edges (source, target, sign) from digitalembryo.segment_polarity
EDGES = [("ci", "wg", +1), ("wg", "en", +1), ("en", "hh", +1),
         ("hh", "ptc", -1), ("ptc", "ci", -1), ("ci", "ptc", +1)]

K = 0.5   # Hill threshold
NH = 4.0  # Hill coefficient


def hill(x, sign):
    h = (x ** NH) / (K ** NH + x ** NH)
    return h if sign > 0 else 1.0 - h


def rhs(x):
    dx = np.zeros(5)
    for i, t in enumerate(NODES):
        acts = [hill(x[IDX[s]], sg) for s, tt, sg in EDGES if tt == t and sg > 0]
        reps = [hill(x[IDX[s]], sg) for s, tt, sg in EDGES if tt == t and sg < 0]
        # production: OR over activators; AND with (1 - repression) per repressor
        a = min(1.0, sum(acts)) if acts else 0.0
        r = 1.0
        for rr in reps:
            r *= rr  # hill(x, -1) already returns 1-h
        dx[i] = a * r - x[i]
    return dx


def integrate(x0, t_end=60.0, dt=0.05):
    x = np.array(x0, dtype=float)
    for _ in range(int(t_end / dt)):
        k1 = rhs(x)
        x = x + dt * k1
        x = np.clip(x, 0, 1.2)
    return x


def main():
    grid = np.linspace(0, 1, 5)
    on, tot = 0, 0
    finals = []
    import itertools
    for x0 in itertools.product(grid, repeat=5):
        xf = integrate(x0)
        finals.append(xf)
        if xf[IDX["wg"]] > 0.5:
            on += 1
        tot += 1
    finals = np.array(finals)
    # count distinct attractors (rounded)
    uniq = np.unique(np.round(finals, 1), axis=0)
    # Boolean reference: enumerate the 5-node core
    import sys
    sys.path.insert(0, str(ROOT / "src"))
    from digitalembryo.segment_polarity import basin_sizes, attractors
    at = attractors()
    bs = basin_sizes()
    out = {
        "design": "Hill ODE counterpart of the 5-node Albert-Othmer core; 5^5=3125 grid of initial conditions",
        "hill_K": K, "hill_n": NH,
        "ode_wg_on_fraction": on / tot,
        "ode_n_distinct_attractors_approx": int(len(uniq)),
        "ode_wg_on_attractor_examples": [list(map(float, u)) for u in uniq if u[IDX["wg"]] > 0.5][:4],
        "boolean_fixed_points": [list(map(int, a)) for a in at["fixed"]],
        "boolean_n_cycles": len(at["cycles"]),
        "boolean_basin_sizes": {str(k): v for k, v in bs.items()} if isinstance(bs, dict) else str(bs),
        "reading": "if the ODE preserves a wg-ON attractor with nonzero basin, the bistability conclusion survives continuous dynamics",
    }
    (ROOT / "results" / "ode_vs_boolean.json").write_text(json.dumps(out, indent=1))
    print(json.dumps({k: out[k] for k in ("ode_wg_on_fraction", "ode_n_distinct_attractors_approx", "boolean_fixed_points")}, indent=1))


if __name__ == "__main__":
    main()
