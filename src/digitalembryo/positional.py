"""Project C: positional information and axis formation.

C1: Wolpert French-flag: morphogen gradient interpreted by thresholds.
    Quantifies positional error vs gradient noise (information-theoretic bound).
C2: Bicoid-style exponential gradient + gap-gene readout (Hunchback-like),
    benchmarked against measured Drosophila embryo precision (~1% egg length,
    Gregor et al. 2007 - cited in paper, not retrained).
"""
from __future__ import annotations
import numpy as np


def exponential_gradient(x, M0=1.0, lam=0.2):
    """Bicoid-like exponential decay from anterior source; x in [0,1] egg length."""
    return M0 * np.exp(-x / lam)


def french_flag_boundaries(M0=1.0, lam=0.2, thresholds=(0.5, 0.2), n=200):
    """Positions where the gradient crosses each threshold (blue/white/red)."""
    x = np.linspace(0, 1, n)
    M = exponential_gradient(x, M0, lam)
    out = []
    for t in thresholds:
        idx = int(np.argmin(np.abs(M - t)))
        out.append(x[idx])
    return np.array(out)


def positional_error_vs_noise(lam=0.2, thresholds=(0.5, 0.2),
                              noises=(0.01, 0.05, 0.1, 0.2), trials=400, seed=0):
    """Positional error (std of boundary position, in egg-length units) as a
    function of multiplicative morphogen readout noise. Theory: for an
    exponential gradient, dx/x ~ noise * lam / x  -> error grows with distance.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for sigma in noises:
        for t in thresholds:
            xb_theory = lam * np.log(1.0 / t)  # boundary at M0*exp(-x/lam)=t
            xs = []
            for _ in range(trials):
                noisy_t = t * np.exp(sigma * rng.standard_normal())
                xs.append(lam * np.log(1.0 / noisy_t))
            xs = np.array(xs)
            rows.append(dict(noise=sigma, threshold=t, x_theory=xb_theory,
                             x_mean=float(xs.mean()), x_std=float(xs.std()),
                             theory_std=float(xb_theory * sigma * lam / max(xb_theory, 1e-9))))
    return rows


def mutual_information_positions(lam=0.2, M0=1.0, noise=0.1, n_pos=64, trials=600, seed=1):
    """Mutual information (bits) between true position and threshold-readout
    of a noisy exponential gradient - small discrete-MI estimator."""
    rng = np.random.default_rng(seed)
    x_true = rng.uniform(0, 1, trials)
    M = exponential_gradient(x_true, M0, lam) * np.exp(noise * rng.standard_normal(trials))
    bins = np.quantile(M, np.linspace(0, 1, 9))
    m_binned = np.clip(np.digitize(M, bins[1:-1]), 0, 7)
    x_binned = np.clip((x_true * n_pos).astype(int), 0, n_pos - 1)
    joint = np.zeros((n_pos, 8))
    for xi, mi in zip(x_binned, m_binned):
        joint[xi, mi] += 1
    joint /= joint.sum()
    px = joint.sum(1, keepdims=True); pm = joint.sum(0, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(joint > 0, joint / (px @ pm), 1.0)
        mi = float(np.nansum(joint * np.log2(ratio)))
    return mi
