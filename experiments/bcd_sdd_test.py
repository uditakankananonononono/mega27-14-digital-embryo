#!/usr/bin/env python3
"""Full-dataset SDD test over all quality-filtered Liu 2013 Bcd gradients.

Outputs results/bcd_sdd.json with per-gradient fits and population stats.
"""
import json
import multiprocessing as mp
import pathlib

import numpy as np

from digitalembryo.bcd import bootstrap_lambda, f_test, fit_double, fit_single, load_live_gradients

ROOT = pathlib.Path(__file__).resolve().parent.parent


def _work(rec):
    x, y = rec["x"], rec["intensity"]
    m = x < 0.6
    if m.sum() < 15:
        return None
    x, y = x[m], y[m]
    s = fit_single(x, y)
    d = fit_double(x, y)
    out = {
        "line": rec["line"], "embryo": rec["embryo"], "side": rec["side"],
        "n": int(len(x)), "A": s["A"], "lam": s["lam"], "rss1": s["rss"],
        "egg_length_um": rec["egg_length_um"],
    }
    if d is not None:
        out.update({"rss2": d["rss"], "l1": d["l1"], "l2": d["l2"],
                    "p_f": f_test(s["rss"], d["rss"], 2, 4, len(x))})
    return out


DOSAGE = {"1IIA": 1, "1IIC": 1, "1XA": 1, "2IIA": 1, "2IIB": 1,
          "1XA1IIA": 2, "1IIA1IIIA": 2, "2IIIA": 2, "2XA": 1, "2IIIB": 2,
          "2IIC": 2, "1XA2IIA": 2, "2XA1IIA": 2, "2XA2IIA": 3, "2IIA2IIIA": 3,
          "2XA2IIIA": 3, "2XA2IIIB": 3, "2XA1IIA2IIIA": 4, "2XA2IIA2IIIA": 4,
          "2XA2IIC2IIIA": 4}


def main():
    recs = load_live_gradients()
    with mp.Pool(2) as pool:
        rows = [r for r in pool.map(_work, recs) if r]
    lams = np.array([r["lam"] for r in rows])
    pvals = np.array([r.get("p_f", 1.0) for r in rows])
    eggs = np.array([r["egg_length_um"] for r in rows])
    amps = np.array([r["A"] for r in rows])
    doses = np.array([DOSAGE.get(r["line"], np.nan) for r in rows])

    per_dose = {}
    for d in sorted(set(doses[~np.isnan(doses)].astype(int))):
        m = doses == d
        per_dose[str(d)] = {
            "n": int(m.sum()),
            "lam_median": float(np.median(lams[m])),
            "lam_iqr": [float(np.percentile(lams[m], 25)), float(np.percentile(lams[m], 75))],
            "A_median": float(np.median(amps[m])),
            "egg_median_um": float(np.median(eggs[m])),
        }
    valid = ~np.isnan(doses)
    # Spearman: does lambda scale with egg length within matched dosage?
    from scipy.stats import spearmanr
    rho_egg, p_egg = spearmanr(lams[valid], eggs[valid])
    rho_dose, p_dose = spearmanr(lams[valid], doses[valid])
    rho_A_dose, p_A_dose = spearmanr(amps[valid], doses[valid])

    # bootstrap CI on median lambda of the 2XA flagship line
    rec2xa = [r for r in load_live_gradients() if r["line"] == "2XA"]
    x, y = rec2xa[0]["x"], rec2xa[0]["intensity"]
    m = x < 0.6
    lo, hi, mean = bootstrap_lambda(x[m], y[m], n_boot=200, rng=0)

    out = {
        "dataset": "Liu et al. PNAS 2013, Zenodo record 4942019, LiveImaging.mat",
        "n_gradients": len(rows),
        "n_lines": len({r["line"] for r in rows}),
        "lam_population": {"median": float(np.median(lams)),
                           "iqr": [float(np.percentile(lams, 25)), float(np.percentile(lams, 75))],
                           "cv": float(np.std(lams) / np.mean(lams))},
        "f_test_reject_fraction": {"p05": float((pvals < 0.05).mean()),
                                   "p05_bonferroni": float((pvals < 0.05 / len(pvals)).mean())},
        "per_dosage": per_dose,
        "spearman": {"lam_vs_egg": [float(rho_egg), float(p_egg)],
                     "lam_vs_dosage": [float(rho_dose), float(p_dose)],
                     "A_vs_dosage": [float(rho_A_dose), float(p_A_dose)]},
        "flagship_2xa_embryo0_lambda_boot_ci95": [lo, hi],
        "rows": rows,
    }
    (ROOT / "results").mkdir(exist_ok=True)
    (ROOT / "results" / "bcd_sdd.json").write_text(json.dumps(out, indent=1))
    print(json.dumps({k: v for k, v in out.items() if k != "rows"}, indent=1))


if __name__ == "__main__":
    main()
