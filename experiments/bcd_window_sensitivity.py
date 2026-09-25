#!/usr/bin/env python3
"""Window sensitivity of the apparent Bcd length scale + DCLS robustness."""
import json
import pathlib

import numpy as np
from scipy import stats

from digitalembryo.bcd import load_live_gradients, loglinear_lambda

ROOT = pathlib.Path(__file__).resolve().parent.parent
DOS = {"1IIA": 1, "1IIC": 1, "1XA": 1, "2IIA": 1, "2IIB": 1, "1XA1IIA": 2,
       "1IIA1IIIA": 2, "2IIIA": 2, "2XA": 1, "2IIIB": 2, "2IIC": 2,
       "1XA2IIA": 2, "2XA1IIA": 2, "2XA2IIA": 3, "2IIA2IIIA": 3, "2XA2IIIA": 3,
       "2XA2IIIB": 3, "2XA1IIA2IIIA": 4, "2XA2IIA2IIIA": 4, "2XA2IIC2IIIA": 4}
WINDOWS = [(0.0, 0.4), (0.0, 0.5), (0.0, 0.6), (0.05, 0.45), (0.1, 0.6), (0.2, 0.8)]

def main():
    recs = [r for r in load_live_gradients() if r["line"] in DOS]
    out = {"dataset": "Liu 2013 Zenodo 4942019 (dosage lines only)", "windows": {}}
    for lo, hi in WINDOWS:
        lams, doses = [], []
        for r in recs:
            lam = loglinear_lambda(r["x"], r["intensity"], lo, hi)
            if np.isfinite(lam) and 0 < lam < 1:
                lams.append(lam)
                doses.append(DOS[r["line"]])
        lams, doses = np.array(lams), np.array(doses)
        rho, p = stats.spearmanr(doses, lams)
        out["windows"][f"{lo}-{hi}"] = {
            "n": int(len(lams)), "lam_median": float(np.median(lams)),
            "lam_iqr": [float(np.percentile(lams, 25)), float(np.percentile(lams, 75))],
            "dcls_rho": float(rho), "dcls_p": float(p)}
        print(f"[{lo:.2f},{hi:.2f}] n={len(lams)} lam={np.median(lams):.4f} DCLS rho={rho:+.3f} p={p:.2e}")
    (ROOT / "results" / "bcd_windows.json").write_text(json.dumps(out, indent=1))

if __name__ == "__main__":
    main()
