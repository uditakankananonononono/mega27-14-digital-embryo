"""Real Bicoid gradient analysis - Liu, Morrison & Gregor PNAS 2013 data.

SDD model: B(x) = A * exp(-x/lambda) (synthesis-diffusion-degradation steady
state). Discovery test: per-embryo model comparison of the pure exponential
against a two-component exponential, with an F-test and bootstrap CIs on
the length scale lambda.
"""
from __future__ import annotations

import pathlib

import numpy as np
import scipy.io
import scipy.optimize
import scipy.stats

DATA = pathlib.Path(__file__).resolve().parents[2] / "data" / "raw" / "liu2013"


def load_live_gradients(path=None, quality_filter=True):
    """Return list of per-embryo Bcd gradient records from LiveImaging.mat."""
    path = path or DATA / "LiveImaging.mat"
    m = scipy.io.loadmat(str(path), squeeze_me=True, struct_as_record=False)
    records = []
    for line in np.atleast_1d(m["FlyLines"]):
        name = str(getattr(line, "FlyLineName"))
        embryos = np.atleast_1d(getattr(line, "Embryos"))
        for idx, e in enumerate(embryos):
            orient = getattr(e, "Orientation", None)
            lr = float(getattr(orient, "LR", np.nan)) if orient is not None else np.nan
            ap = float(getattr(orient, "AP", np.nan)) if orient is not None else np.nan
            if quality_filter and not (lr == 0 and ap == 0):
                continue  # Liu et al. selection: LR=0, AP=0 only
            g = getattr(e, "Gradient")
            for side in ("left", "right"):
                xy = np.asarray(getattr(g, side), dtype=float)
                if xy.ndim != 2 or xy.shape[1] != 2 or xy.shape[0] < 10:
                    continue
                records.append(
                    {
                        "line": name,
                        "embryo": idx,
                        "side": side,
                        "x": xy[:, 0],
                        "intensity": xy[:, 1],
                        "egg_length_um": float(getattr(e, "EggLength", np.nan)),
                        "cf": float(getattr(e, "CF", np.nan)),
                    }
                )
    return records


def _exp1(x, a, lam):
    return a * np.exp(-x / lam)


def _exp2(x, a1, l1, a2, l2):
    return a1 * np.exp(-x / l1) + a2 * np.exp(-x / l2)


def fit_single(x, y):
    """MLE (least squares, Gaussian error) fit of the SDD single exponential."""
    y = np.clip(y, 1e-3, None)
    p0 = [y.max(), 0.15]
    p, _ = scipy.optimize.curve_fit(_exp1, x, y, p0=p0, maxfev=20000)
    rss = float(np.sum((y - _exp1(x, *p)) ** 2))
    return {"A": p[0], "lam": p[1], "rss": rss, "k": 2}


def fit_double(x, y):
    """Two-component exponential fit (transport + retention modes)."""
    s = fit_single(x, y)
    p0 = [s["A"] * 0.8, s["lam"] * 0.6, s["A"] * 0.2, s["lam"] * 2.5]
    bounds = ([0, 1e-3, 0, 1e-3], [np.inf, 2.0, np.inf, 2.0])
    try:
        p, _ = scipy.optimize.curve_fit(_exp2, x, y, p0=p0, bounds=bounds, maxfev=40000)
    except RuntimeError:
        return None
    rss = float(np.sum((y - _exp2(x, *p)) ** 2))
    # order components by length scale: l1 short, l2 long
    order = np.argsort([p[1], p[3]])
    p = [p[0], p[1], p[2], p[3]] if order[0] == 0 else [p[2], p[3], p[0], p[1]]
    return {"A1": p[0], "l1": p[1], "A2": p[2], "l2": p[3], "rss": rss, "k": 4}


def f_test(rss1, rss2, k1, k2, n):
    """Nested-model F test: does the k2-parameter model beat the k1 one?"""
    if rss2 <= 0 or rss1 <= rss2:
        return 1.0
    f = ((rss1 - rss2) / (k2 - k1)) / (rss2 / (n - k2))
    return float(1 - scipy.stats.f.cdf(f, k2 - k1, n - k2))


def bootstrap_lambda(x, y, n_boot=200, rng=None):
    """Residual bootstrap CI for the SDD length scale lambda."""
    rng = np.random.default_rng(rng)
    s = fit_single(x, y)
    resid = y - _exp1(x, s["A"], s["lam"])
    lams = []
    for _ in range(n_boot):
        yb = _exp1(x, s["A"], s["lam"]) + rng.choice(resid, size=len(resid), replace=True)
        try:
            lams.append(fit_single(x, np.clip(yb, 1e-3, None))["lam"])
        except RuntimeError:
            continue
    lo, hi = np.percentile(lams, [2.5, 97.5])
    return float(lo), float(hi), float(np.mean(lams))


def loglinear_lambda(x, y, lo=0.2, hi=0.8):
    """Liu et al. estimator: linear fit of ln C(x) over [lo, hi] EL; lambda=-1/slope."""
    x, y = np.asarray(x), np.asarray(y)
    m = (x >= lo) & (x <= hi) & (y > 0)
    if m.sum() < 8:
        return np.nan
    slope = np.polyfit(x[m], np.log(y[m]), 1)[0]
    return float(-1.0 / slope) if slope < 0 else np.nan
