#!/usr/bin/env python3
"""CF-position decoding benchmark: Liu 2013 threshold model vs ridge vs CNN."""
import json
import pathlib

import numpy as np

from digitalembryo.bcd import loglinear_lambda
from digitalembryo.decoder import cnn_decoder, embryo_table, grouped_folds, liu_baseline, ridge_features

DOS = {"1IIA": 1, "1IIC": 1, "1XA": 1, "2IIA": 1, "2IIB": 1, "1XA1IIA": 2,
       "1IIA1IIIA": 2, "2IIIA": 2, "2XA": 1, "2IIIB": 2, "2IIC": 2,
       "1XA2IIA": 2, "2XA1IIA": 2, "2XA2IIA": 3, "2IIA2IIIA": 3, "2XA2IIIA": 3,
       "2XA2IIIB": 3, "2XA1IIA2IIIA": 4, "2XA2IIA2IIIA": 4, "2XA2IIC2IIIA": 4}
# truncated lines carry no clean integer dose; use their bcd copy estimate from
# Liu Table S1 naming (BNT/BT/BN series ~2 copies of the transgene + endogenous)
for extra in ["1IIA BNT", "1XA BNT", "1XA BT", "1XA BN", "2XA BNT", "2XA BT",
              "2XA BN", "2XA1IIA BT", "2XA1IIA BN"]:
    DOS.setdefault(extra, 2)


def liu_true_dose(rows, train, test):
    """Liu 2013 exactly: xCF = Sx*ln(D)+c with the true genetic dose D."""
    D = np.array([np.log(DOS[r["line"]]) for r in rows])
    y = np.array([r["cf"] for r in rows])
    a, b = np.polyfit(D[train], y[train], 1)
    return a * D[test] + b


def cnn_seeded(rows, train, test):
    preds = []
    for seed in (0, 1, 2):
        preds.append(cnn_decoder(rows, train, test, seed=seed))
    return np.mean(preds, axis=0)

ROOT = pathlib.Path(__file__).resolve().parent.parent


def rmse(a, b):
    return float(np.sqrt(np.mean((np.asarray(a) - np.asarray(b)) ** 2)))


def main():
    rows = embryo_table()
    folds = grouped_folds(rows, k=5, seed=0)
    y = np.array([r["cf"] for r in rows])
    res = {"n_embryos": len(rows), "n_lines": len({r["line"] for r in rows}),
           "cv": "5-fold, disjoint fly lines", "models": {}}
    for name, fn in [("liu_threshold_lnA", liu_baseline),
                     ("liu_threshold_trueD", liu_true_dose),
                     ("ridge_lnA_lambda", ridge_features),
                     ("cnn1d_profile_3seed", cnn_seeded)]:
        preds = np.full(len(rows), np.nan)
        for train, test in folds:
            preds[test] = fn(rows, train, test)
        r = rmse(preds, y)
        ss = 1 - np.sum((preds - y) ** 2) / np.sum((y - y.mean()) ** 2)
        res["models"][name] = {"rmse_EL": r, "rmse_pctEL": 100 * r, "r2": float(ss)}
        print(f"{name:20s} RMSE={100*r:.2f}%EL  R2={ss:.3f}")
    # paired bootstrap: CNN vs Liu-trueD RMSE difference
    rng = np.random.default_rng(0)
    pa = np.full(len(rows), np.nan); pb = np.full(len(rows), np.nan)
    for train, test in folds:
        pa[test] = liu_true_dose(rows, train, test)
        pb[test] = cnn_seeded(rows, train, test)
    idx = np.arange(len(rows)); diffs = []
    for _ in range(2000):
        s_ = rng.choice(idx, size=len(idx), replace=True)
        diffs.append(rmse(pa[s_], y[s_]) - rmse(pb[s_], y[s_]))
    lo, hi = np.percentile(diffs, [2.5, 97.5])
    res["paired_bootstrap_cnn_minus_liu"] = {
        "rmse_gain_mean_pctEL": float(100 * np.mean(diffs)),
        "ci95_pctEL": [float(100 * lo), float(100 * hi)],
        "p_gain_le_0": float(np.mean(np.array(diffs) <= 0))}
    print(f"CNN gain over Liu-trueD: {100*np.mean(diffs):.3f}%EL CI95 [{100*lo:.3f},{100*hi:.3f}]")
    res["published_bar"] = {
        "liu2013_Sx_pctEL": 10.5,
        "liu2013_relC_at_CF_sd_pct": "16+-5",
        "note": "Liu et al. PNAS 2013 threshold model; bar = its CV RMSE on identical data"}
    (ROOT / "results" / "cf_decoder.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
