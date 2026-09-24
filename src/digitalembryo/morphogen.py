"""Project A: Turing pattern formation in a 2D embryonic field.

Gierer-Meinhardt activator-inhibitor kinetics:
    da/dt = rho * a^2 / h - mu_a * a + Da * lap(a) + sigma_a
    dh/dt = rho * a^2     - mu_h * h + Dh * lap(h)
Linearised Turing (diffusion-driven) instability analysis is provided so the
simulation can be validated against the analytic threshold - math rigor gate.
"""
from __future__ import annotations
import numpy as np


def turing_thresholds(Da, Dh, mu_a, mu_h, rho=1.0):
    """Analytic Turing-instability conditions at the homogeneous steady state.

    Returns dict with steady state, Jacobian trace/det, critical wavenumber
    range (k2_min, k2_max) of unstable modes, and whether instability exists.
    Derivation: linearise f,g at (a*,h*); require tr(J)<0, det(J)>0 (stable
    without diffusion) and that the dispersion relation det(J - k^2 D) has a
    positive root band (unstable with diffusion).
    """
    # steady state: h-eq gives h = rho a^2/mu_h; a-eq gives h = rho a/mu_a
    # -> a* = mu_h/mu_a, h* = rho mu_h / mu_a^2
    a_star = mu_h / mu_a
    h_star = rho * a_star / mu_a
    # Jacobian of reaction terms at (a*, h*)
    fa = 2 * rho * a_star / h_star - mu_a
    fh = -rho * a_star**2 / h_star**2
    ga = 2 * rho * a_star
    gh = -mu_h
    tr = fa + gh
    det = fa * gh - fh * ga
    # Dispersion: p(k^2) = det(J - k^2 D) = Da Dh k^4 - (Dh fa + Da gh) k^2 + det
    B = Dh * fa + Da * gh
    disc = B * B - 4 * Da * Dh * det
    unstable = tr < 0 and det > 0 and B > 0 and disc > 0
    if unstable:
        k2_max = (B + np.sqrt(disc)) / (2 * Da * Dh)
        k2_min = (B - np.sqrt(disc)) / (2 * Da * Dh)
    else:
        k2_min = k2_max = float("nan")
    return dict(a_star=a_star, h_star=h_star, trace=tr, det=det,
                B=B, unstable=bool(unstable), k2_min=k2_min, k2_max=k2_max)


class GiererMeinhardt2D:
    """Explicit-Euler finite-difference solver, periodic boundaries."""

    def __init__(self, n=64, Da=0.005, Dh=0.2, mu_a=0.06, mu_h=0.12,
                 rho=1.0, dt=None, seed=0):
        self.n = n
        self.Da, self.Dh, self.mu_a, self.mu_h, self.rho = Da, Dh, mu_a, mu_h, rho
        th = turing_thresholds(Da, Dh, mu_a, mu_h, rho)
        rng = np.random.default_rng(seed)
        self.a = th["a_star"] * (1.0 + 0.01 * rng.standard_normal((n, n)))
        self.h = th["h_star"] * (1.0 + 0.01 * rng.standard_normal((n, n)))
        # CFL-style stability bound for explicit diffusion: dt <= dx^2/(4D)
        max_d = max(Da, Dh)
        cfl = 1.0 / (4.0 * max_d + 1e-12)
        self.dt = dt if dt is not None else 0.8 * cfl
        if self.dt > cfl:
            raise ValueError("dt violates explicit-diffusion stability bound")

    @staticmethod
    def _lap(u):
        return (-4 * u + np.roll(u, 1, 0) + np.roll(u, -1, 0)
                + np.roll(u, 1, 1) + np.roll(u, -1, 1))

    def step(self, steps=1):
        for _ in range(steps):
            a, h = self.a, self.h
            act = self.rho * a * a / np.maximum(h, 1e-9)
            da = act - self.mu_a * a + self.Da * self._lap(a)
            dh = self.rho * a * a - self.mu_h * h + self.Dh * self._lap(h)
            self.a = np.clip(a + self.dt * da, 0.0, None)
            self.h = np.clip(h + self.dt * dh, 1e-9, None)
        return self.a, self.h

    def pattern_wavenumber(self):
        """Dominant k^2 of the activator field via radial FFT power."""
        f = np.fft.fft2(self.a - self.a.mean())
        p = np.abs(f) ** 2
        kx = np.fft.fftfreq(self.n) * 2 * np.pi
        ky = np.fft.fftfreq(self.n) * 2 * np.pi
        KX, KY = np.meshgrid(kx, ky, indexing="ij")
        k2 = KX**2 + KY**2
        p = p.ravel(); k2 = k2.ravel()
        p[k2 == 0] = 0
        return k2[np.argmax(p)]
