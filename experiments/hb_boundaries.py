#!/usr/bin/env python3
"""Hb boundary positions vs Bcd dosage (Liu 2013 reproduction + per-embryo stats)."""
import json
import pathlib

import numpy as np
from scipy import stats

from digitalembryo.gapgene import hb_boundaries, load_if_sessions

ROOT = pathlib.Path(__file__).resolve().parent.parent
DOS = {"1XA": 1, "2XA": 1, "2XA2IIIA": 3, "1XA BNT": 1, "2XA BNT": 2}  # bcd-GFP dose (Table S1)


def main():
    ss = load_if_sessions()
    per_session, all_b, all_d = [], [], []
    for s in ss:
        if "Hb" not in s["genes"]:
            continue
        b = hb_boundaries(s)
        b = b[np.isfinite(b)]
        if len(b) < 5:
            continue
        d = DOS.get(s["line"])
        row = {"line": s["line"], "genes": s["genes"], "n": int(len(b)),
               "median_EL": float(np.median(b)),
               "iqr": [float(np.percentile(b, 25)), float(np.percentile(b, 75))],
               "bcd_dose": d}
        per_session.append(row)
        if d:
            all_b.extend(b.tolist())
            all_d.extend([d] * len(b))
    rho, p = stats.spearmanr(all_d, all_b)
    out = {"sessions": per_session, "n_embryos": len(all_b),
           "boundary_vs_dosage_spearman": [float(rho), float(p)],
           "note": "Hb boundary shifts posterior with Bcd dosage (Liu 2013 reproduction)"}
    (ROOT / "results" / "hb_boundaries.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
