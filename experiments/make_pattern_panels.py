"""2x2 panel of steady-state Gierer-Meinhardt activator fields across the mu_a sweep."""
import pathlib, sys
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))
from digitalembryo.morphogen import GiererMeinhardt2D
fig, axes = plt.subplots(2, 2, figsize=(8, 8))
for ax, mu_a in zip(axes.ravel(), [0.025, 0.04, 0.055, 0.075]):
    sim = GiererMeinhardt2D(n=48, Da=0.005, Dh=0.2, mu_a=mu_a, mu_h=0.12, seed=3)
    sim.step(2500)
    ax.imshow(sim.a, cmap='viridis')
    ax.set_title(f"$\\mu_a$={mu_a}  ($k^2_{{sim}}$={sim.pattern_wavenumber():.2f})")
    ax.set_xticks([]); ax.set_yticks([])
fig.tight_layout()
fig.savefig('paper/figures/fig_turing_panels.pdf')
print('saved')
