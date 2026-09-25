import json, pathlib, sys
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))
from digitalembryo.morphogen import turing_thresholds
fig, ax = plt.subplots(figsize=(7.5, 4.2))
for mu_a in [0.025, 0.04, 0.06, 0.08]:
    Da, Dh, mu_h = 0.005, 0.2, 0.12
    tr = mu_a - mu_h; det = mu_a * mu_h
    B = Dh * mu_a - Da * mu_h
    k2 = np.linspace(0.01, 6, 400)
    # max eigenvalue of M = J - k2 D, computed directly
    Ja, Jh = mu_a, -mu_h
    Jfh = -mu_a**2 / mu_h; Jga = 2 * mu_h
    tau = tr - (Da + Dh) * k2
    Delta = det - B * k2 + Da * Dh * k2**2
    disc = np.sqrt(np.maximum(tau**2 - 4 * Delta, 0))
    re = np.where(tau**2 - 4 * Delta >= 0, 0.5 * (tau + disc), 0.5 * tau)
    ax.plot(k2, re, label=f"$\\mu_a$={mu_a}")
ax.axhline(0, color='k', lw=0.6)
ax.set_xlabel("$k^2$"); ax.set_ylabel("Re $\\sigma_+(k^2)$")
ax.legend(); fig.tight_layout()
fig.savefig('paper/figures/fig_dispersion.pdf')
print('saved')
