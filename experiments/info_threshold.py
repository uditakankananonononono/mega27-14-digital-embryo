#!/usr/bin/env python3
"""Information-threshold sweep: minimum Bicoid-gradient information needed to
decode cephalic-furrow position.

Verdict reframe 2 (decoder limits -> minimum-information threshold): degrade the
input gradient three ways (additive noise, anterior-window truncation, sparse
sampling) and measure held-out decode RMSE under the same disjoint-line 5-fold
protocol as cf_decoder.py. The critical threshold is where RMSE crosses the
Liu 2013 true-dose benchmark (4.94% EL).
"""
import json
import pathlib

import numpy as np

from digitalembryo.decoder import GRID, embryo_table, grouped_folds

ROOT = pathlib.Path(__file__).resolve().parent.parent
LIU_TRUED_PCT = 4.937638103679565  # archived benchmark, results/cf_decoder.json
MEAN_PCT = None  # filled from data (predict-train-mean baseline)


def ridge_decode(X, y, train, test, lam=1.0):
    Xm = np.column_stack([np.ones(len(X)), X])
    w = np.linalg.solve(Xm[train].T @ Xm[train] + lam * np.eye(Xm.shape[1]),
                        Xm[train].T @ y[train])
    return Xm[test] @ w


def sweep(rows, degrade, levels, seed=0):
    """Run 5-fold disjoint-line CV at each degradation level."""
    y = np.array([r["cf"] for r in rows])
    P = np.array([r["profile"] for r in rows])  # log-intensity on GRID
    folds = grouped_folds(rows, k=5, seed=seed)
    out = []
    for lev in levels:
        errs = []
        for train, test in folds:
            Xt = degrade(P, lev, train, test)
            pred = ridge_decode(Xt, y, train, test)
            errs.append(np.sqrt(np.mean((pred - y[test]) ** 2)))
        out.append({"level": lev, "rmse_pctEL": float(np.mean(errs)) * 100.0})
    return out


def main():
    rows = embryo_table()
    y = np.array([r["cf"] for r in rows])
    folds = grouped_folds(rows, k=5, seed=0)
    mean_err = np.mean([np.sqrt(np.mean((y[train].mean() - y[test]) ** 2)) for train, test in folds]) * 100.0

    rng = np.random.default_rng(7)

    def noise_arm(sig):
        def f(P, lev, train, test):
            scale = P[train].std()
            return P + rng.normal(0, lev * scale, P.shape)
        return f

    def trunc_arm(P, w, train, test):
        X = P.copy()
        X[:, GRID > w] = 0.0
        return X

    def sparse_arm(P, k, train, test):
        k = int(k)
        if k >= P.shape[1]:
            return P.copy()
        idx = np.linspace(0, P.shape[1] - 1, k).round().astype(int)
        X = np.zeros_like(P)
        X[:, idx] = P[:, idx]
        return X

    res = {
        "design": "5-fold disjoint-fly-line CV; ridge decoder on 64-point log-Bcd profile; three degradation arms",
        "n_embryos": len(rows), "n_lines": len({r['line'] for r in rows}),
        "benchmarks": {"liu_truedose_pctEL": LIU_TRUED_PCT, "predict_mean_pctEL": float(mean_err)},
        "noise": sweep(rows, noise_arm(None), [0.0, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0]),
        "truncation": sweep(rows, trunc_arm, [0.6, 0.55, 0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1]),
        "sparse": sweep(rows, sparse_arm, [64, 32, 16, 12, 8, 6, 4, 2, 1]),
    }
    # critical thresholds: first level crossing the Liu benchmark
    for arm in ("noise", "truncation", "sparse"):
        rows_arm = res[arm]
        cross = next((r["level"] for r in rows_arm if r["rmse_pctEL"] > LIU_TRUED_PCT), None)
        res[f"{arm}_critical_level_vs_liu"] = cross
    out = ROOT / "results" / "info_threshold.json"
    out.write_text(json.dumps(res, indent=2))
    print(json.dumps({k: res[k] for k in ("noise_critical_level_vs_liu", "truncation_critical_level_vs_liu", "sparse_critical_level_vs_liu")}, indent=2))
    print("full:", res["truncation"])
    print("sparse:", res["sparse"])
    print("noise:", res["noise"])


if __name__ == "__main__":
    main()
