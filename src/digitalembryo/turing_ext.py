"""Beyond-onset wavenumber renormalization for Gierer-Meinhardt patterns.

Linear theory predicts the fastest-growing mode k2_lin (max of the exact
dispersion relation). The nonlinear steady-state pattern wavenumber deviates;
we quantify the deviation as a function of distance from onset.
"""
from __future__ import annotations

import numpy as np

from .morphogen import GiererMeinhardt2D, turing_thresholds


def dispersion(Da, Dh, mu_a, mu_h, rho, k2):
    """Exact max growth rate sigma(k^2) from the 2x2 linearised RD system."""
    th = turing_thresholds(Da, Dh, mu_a, mu_h, rho)
    a_star, h_star = th["a_star"], th["h_star"]
    fa = 2 * rho * a_star / h_star - mu_a
    fh = -rho * a_star**2 / h_star**2
    ga = 2 * rho * a_star
    gh = -mu_h
    tr_k = fa + gh - (Da + Dh) * k2
    p_k = Da * Dh * k2**2 - (Dh * fa + Da * gh) * k2 + (fa * gh - fh * ga)
    disc = tr_k**2 - 4 * p_k
    disc = np.maximum(disc, 0.0)
    return (tr_k + np.sqrt(disc)) / 2


def linear_fastest_k2(Da, Dh, mu_a, mu_h, rho=1.0):
    """k^2 maximising the exact dispersion relation inside the unstable band."""
    th = turing_thresholds(Da, Dh, mu_a, mu_h, rho)
    if not th["unstable"]:
        return np.nan
    grid = np.linspace(th["k2_min"], th["k2_max"], 4000)
    sig = dispersion(Da, Dh, mu_a, mu_h, rho, grid)
    return float(grid[int(np.argmax(sig))])


def onset_distance(Da, Dh, mu_a, mu_h, rho=1.0):
    """epsilon: normalised distance from the Turing-onset boundary in mu_a.

    epsilon = 1 - mu_a/mu_a_crit where mu_a_crit is the largest mu_a still
    unstable (bisection). epsilon=0 is onset, larger = deeper in the band.
    """
    if not turing_thresholds(Da, Dh, mu_a, mu_h, rho)["unstable"]:
        return np.nan
    lo, hi = mu_a, mu_a * 64.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if turing_thresholds(Da, Dh, mid, mu_h, rho)["unstable"]:
            lo = mid
        else:
            hi = mid
    mu_crit = lo
    return float(1.0 - mu_a / mu_crit), float(mu_crit)


def run_point(Da, Dh, mu_a, mu_h, rho=1.0, n=48, steps=2500, seed=0):
    """One sweep point: simulate to steady state, return measured k2 + theory."""
    th = turing_thresholds(Da, Dh, mu_a, mu_h, rho)
    if not th["unstable"]:
        return None
    k2_lin = linear_fastest_k2(Da, Dh, mu_a, mu_h, rho)
    eps, mu_crit = onset_distance(Da, Dh, mu_a, mu_h, rho)
    sim = GiererMeinhardt2D(n=n, Da=Da, Dh=Dh, mu_a=mu_a, mu_h=mu_h, rho=rho, seed=seed)
    sim.step(steps)
    k2_sim = float(sim.pattern_wavenumber())
    return {"Da": Da, "Dh": Dh, "mu_a": mu_a, "mu_h": mu_h,
            "k2_lin": k2_lin, "k2_sim": k2_sim, "epsilon": eps,
            "mu_crit": mu_crit, "dev": (k2_sim - k2_lin) / k2_lin}
