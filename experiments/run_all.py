"""Run all three digital-embryo studies; writes results JSON + figures."""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from digitalembryo.morphogen import GiererMeinhardt2D, turing_thresholds
from digitalembryo.grn import (count_fixed_points, simulate_toggle,
                               derrida_sensitivity, boolean_network_derrida)
from digitalembryo.positional import (positional_error_vs_noise,
                                      mutual_information_positions,
                                      french_flag_boundaries)

OUT = os.path.join(os.path.dirname(__file__), "..", "results")
FIG = os.path.join(OUT, "figures")
os.makedirs(FIG, exist_ok=True)
res = {}

# --- Project A: Turing patterns, simulation vs analytic threshold ---
th = turing_thresholds(0.005, 0.2, 0.06, 0.12)
sim = GiererMeinhardt2D(n=64, seed=3)
sim.step(6000)
k2 = sim.pattern_wavenumber()
res["A_turing"] = dict(analytic={k: (v if isinstance(v, float) else v)
                                 for k, v in th.items()},
                       sim_dominant_k2=float(k2),
                       pattern_std=float(sim.a.std()),
                       inside_band=bool(th["k2_min"]*0.5 < k2 < th["k2_max"]*2.0))
fig, ax = plt.subplots(1, 2, figsize=(8, 3.5))
ax[0].imshow(sim.a, cmap="viridis"); ax[0].set_title("Activator field (Turing pattern)")
kk = np.linspace(max(th["k2_min"]*0.2, 1e-6), th["k2_max"]*3, 200)
B, Da, Dh, det = th["B"], 0.005, 0.2, th["det"]
ax[1].plot(kk, Da*Dh*kk**2 - B*kk + det)
ax[1].axhline(0, color="k", lw=0.5); ax[1].axvline(k2, color="r", ls="--", label="sim dominant $k^2$")
ax[1].set_xlabel("$k^2$"); ax[1].set_ylabel("dispersion det(J-$k^2$D)"); ax[1].legend()
fig.tight_layout(); fig.savefig(os.path.join(FIG, "A_turing.png"), dpi=150); plt.close(fig)

# --- Project B: toggle switch fates + Derrida ---
fates = count_fixed_points(grid=24)
T, U, V = simulate_toggle(0.1, 2.0)
derrida = {str(K): dict(theory=derrida_sensitivity(K),
                        measured=boolean_network_derrida(K=K, n_nets=4, n_pairs=12, seed=K))
           for K in (1, 2, 3, 4)}
res["B_grn"] = dict(fates=fates, derrida=derrida)
fig, ax = plt.subplots(1, 2, figsize=(8, 3.5))
ax[0].plot(T, U, label="u (fate A)"); ax[0].plot(T, V, label="v (fate B)")
ax[0].set_xlabel("time"); ax[0].legend(); ax[0].set_title("Toggle-switch differentiation")
Ks = [1, 2, 3, 4]
ax[1].plot(Ks, [derrida[str(k)]["theory"] for k in Ks], "o-", label="annealed theory")
ax[1].plot(Ks, [derrida[str(k)]["measured"] for k in Ks], "s--", label="measured")
ax[1].axhline(1, color="k", lw=0.5); ax[1].set_xlabel("K inputs"); ax[1].set_ylabel("Derrida $\\lambda$")
ax[1].legend(); ax[1].set_title("Ordered vs chaotic GRNs")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "B_grn.png"), dpi=150); plt.close(fig)

# --- Project C: positional information ---
err = positional_error_vs_noise(trials=300, seed=2)
mis = [mutual_information_positions(noise=nz, trials=500, seed=5) for nz in (0.05, 0.1, 0.2, 0.4)]
res["C_positional"] = dict(error_rows=err, mi_bits=list(map(float, mis)),
                           flag_boundaries=list(map(float, french_flag_boundaries())))
fig, ax = plt.subplots(1, 2, figsize=(8, 3.5))
for t in (0.5, 0.2):
    xs = [d["noise"] for d in err if d["threshold"] == t]
    ys = [d["x_std"] for d in err if d["threshold"] == t]
    ax[0].plot(xs, ys, "o-", label=f"threshold {t}")
ax[0].set_xlabel("morphogen noise"); ax[0].set_ylabel("boundary std (egg length)")
ax[0].legend(); ax[0].set_title("Positional error vs noise")
ax[1].plot([0.05, 0.1, 0.2, 0.4], mis, "s-")
ax[1].set_xlabel("noise"); ax[1].set_ylabel("positional information (bits)")
ax[1].set_title("Information carried by gradient")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "C_positional.png"), dpi=150); plt.close(fig)

with open(os.path.join(OUT, "results.json"), "w") as f:
    json.dump(res, f, indent=2)
print(json.dumps({k: (v if not isinstance(v, dict) else list(v)) for k, v in res.items()}, indent=1))
print("A inside_band:", res["A_turing"]["inside_band"], "| B fates:", fates)
