"""Project B: gene-regulatory-network differentiation.

B1: Gardner-Collins mutual-inhibition toggle switch (bistable cell-fate circuit).
    du/dt = a1/(1+v^n) - u
    dv/dt = a2/(1+u^n) - v
B2: Boolean ensemble - random GRNs classified as ordered vs chaotic via the
    Derrida annealed approximation, validated by direct simulation.
"""
from __future__ import annotations
import numpy as np


def toggle_rhs(state, a1, a2, n):
    u, v = state
    return np.array([a1 / (1.0 + v**n) - u, a2 / (1.0 + u**n) - v])


def simulate_toggle(u0, v0, a1=2.5, a2=2.5, n=3.0, dt=0.05, t_end=200.0):
    """Integrate the toggle switch; returns (t, U, V) trajectories."""
    steps = int(t_end / dt)
    U = np.empty(steps); V = np.empty(steps); T = np.arange(steps) * dt
    s = np.array([u0, v0], dtype=float)
    for i in range(steps):
        s = s + dt * toggle_rhs(s, a1, a2, n)
        s = np.clip(s, 0, None)
        U[i], V[i] = s
    return T, U, V


def count_fixed_points(a1=2.5, a2=2.5, n=3.0, grid=200):
    """Count stable fixed points by basin sampling on a grid of ICs."""
    fates = []
    for u0 in np.linspace(0.05, max(a1, a2) * 1.2, grid):
        for v0 in np.linspace(0.05, max(a1, a2) * 1.2, grid):
            s = np.array([u0, v0])
            for _ in range(4000):
                s = s + 0.05 * toggle_rhs(s, a1, a2, n)
                s = np.clip(s, 0, None)
            fates.append("high_u" if s[0] > s[1] else "high_v")
    fates = np.array(fates)
    return {"high_u": int((fates == "high_u").sum()),
            "high_v": int((fates == "high_v").sum())}


def derrida_sensitivity(K, p=0.5):
    """Annealed approximation: average sensitivity of a Boolean network.
    lambda = 2 K p (1-p). lambda > 1 -> chaotic, < 1 -> ordered.
    """
    return 2 * K * p * (1 - p)


def boolean_network_derrida(N=12, K=3, p=0.5, n_nets=8, n_pairs=24, seed=0):
    """Measured Derrida coefficient: mean Hamming divergence after one step
    from 1-bit-flipped initial states, averaged over random networks.
    Returns measured lambda to compare against derrida_sensitivity(K, p).
    """
    rng = np.random.default_rng(seed)
    lambdas = []
    for _ in range(n_nets):
        funcs = rng.random((N, 2**K)) < p
        inputs = rng.integers(0, N, size=(N, K))
        for _ in range(n_pairs):
            x = rng.integers(0, 2, size=N)
            y = x.copy(); y[rng.integers(N)] ^= 1
            def step(z):
                idx = (z[inputs] * (1 << np.arange(K - 1, -1, -1))).sum(axis=1)
                return funcs[np.arange(N), idx].astype(int)
            d = int((step(x) ^ step(y)).sum())
            lambdas.append(d / 1.0)  # divergence per unit initial distance
    return float(np.mean(lambdas))
