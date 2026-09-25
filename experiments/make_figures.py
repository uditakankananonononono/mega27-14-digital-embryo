#!/usr/bin/env python3
"""Generate all paper figures from committed results/*.json and raw data."""
import json
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
FIG = ROOT / "paper" / "figures"
FIG.mkdir(parents=True, exist_ok=True)
RES = ROOT / "results"

import sys
sys.path.insert(0, str(ROOT / "src"))


def fig_bcd_gradients():
    from digitalembryo.bcd import fit_single, load_live_gradients
    recs = [r for r in load_live_gradients() if r["line"] == "2XA"][:6]
    fig, axes = plt.subplots(2, 3, figsize=(9, 5.2), sharex=True)
    for ax, r in zip(axes.ravel(), recs):
        m = r["x"] < 0.6
        x, y = r["x"][m], r["intensity"][m]
        s = fit_single(x, y)
        ax.semilogy(x, y, ".", ms=2, color="0.4")
        xs = np.linspace(0, 0.6, 200)
        ax.semilogy(xs, s["A"] * np.exp(-xs / s["lam"]), "r-", lw=1.2)
        ax.set_title(f"embryo {r['embryo']} {r['side']}  $\\lambda$={s['lam']:.3f} EL", fontsize=8)
    fig.suptitle("Real Bcd-GFP gradients (Liu 2013, line 2XA) with SDD exponential fits")
    fig.tight_layout()
    fig.savefig(FIG / "fig_bcd_gradients.pdf")


def fig_lambda_dosage():
    d = json.loads((RES / "bcd_sdd.json").read_text())
    rows = [r for r in d["rows"]]
    DOS = {"1IIA": 1, "1IIC": 1, "1XA": 1, "2IIA": 1, "2IIB": 1, "1XA1IIA": 2,
           "1IIA1IIIA": 2, "2IIIA": 2, "2XA": 1, "2IIIB": 2, "2IIC": 2,
           "1XA2IIA": 2, "2XA1IIA": 2, "2XA2IIA": 3, "2IIA2IIIA": 3, "2XA2IIIA": 3,
           "2XA2IIIB": 3, "2XA1IIA2IIIA": 4, "2XA2IIA2IIIA": 4, "2XA2IIC2IIIA": 4}
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    for dose, c in [(1, "C0"), (2, "C1"), (3, "C2"), (4, "C3")]:
        v = [r["lam"] for r in rows if DOS.get(r["line"]) == dose]
        jit = np.random.default_rng(0).normal(0, 0.03, len(v))
        ax.scatter(np.full(len(v), dose) + jit, v, s=3, alpha=0.25, color=c)
        ax.errorbar(dose, np.median(v), yerr=np.std(v) / np.sqrt(len(v)),
                    fmt="ks", ms=6, capsize=3, zorder=5)
    ax.set_xlabel("bcd-GFP dosage (functional copies)")
    ax.set_ylabel(r"gradient length scale $\lambda$ (EL)")
    ax.set_title("Dosage-coupled length scale (DCLS): per-embryo $\\lambda$ vs dosage")
    fig.tight_layout()
    fig.savefig(FIG / "fig_lambda_dosage.pdf")


def fig_windows():
    d = json.loads((RES / "bcd_windows.json").read_text())
    ws = d["windows"]
    keys = sorted(ws, key=lambda k: float(k.split("-")[0]) + float(k.split("-")[1]) / 100)
    mid = [(float(k.split("-")[0]) + float(k.split("-")[1])) / 2 for k in keys]
    lam = [ws[k]["lam_median"] for k in keys]
    rho = [ws[k]["dcls_rho"] for k in keys]
    fig, ax1 = plt.subplots(figsize=(5.4, 3.6))
    ax1.plot(mid, lam, "o-", label=r"median $\lambda$ (EL)")
    ax1.axhline(0.165, color="k", ls="--", lw=1, label="Liu 2013 corrected $\\lambda$=0.165")
    ax1.set_xlabel("fit-window midpoint (EL)")
    ax1.set_ylabel(r"apparent $\lambda$ (EL)")
    ax2 = ax1.twinx()
    ax2.plot(mid, rho, "s--", color="C3", label=r"DCLS $\rho$")
    ax2.set_ylabel(r"DCLS Spearman $\rho$", color="C3")
    ax2.axhline(0, color="C3", lw=0.5, alpha=0.4)
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, fontsize=7, loc="upper left")
    ax1.set_title("Window-sensitivity law and anterior weighting of DCLS")
    fig.tight_layout()
    fig.savefig(FIG / "fig_windows.pdf")


def fig_cf_decoder():
    d = json.loads((RES / "cf_decoder.json").read_text())
    m = d["models"]
    names = ["liu_threshold_trueD", "liu_threshold_lnA", "ridge_lnA_lambda", "cnn1d_profile_3seed"]
    labels = ["Liu 2013\nthreshold (true D)", "Liu 2013\nthreshold (ln A)", "ridge\n(ln A, $\\lambda$)", "1D-CNN\n(full profile)"]
    vals = [100 * m[n]["rmse_EL"] for n in names]
    r2s = [m[n]["r2"] for n in names]
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    bars = ax.bar(labels, vals, color=["0.6", "0.6", "C2", "C1"])
    for b, r2 in zip(bars, r2s):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.05,
                f"$R^2$={r2:.2f}", ha="center", fontsize=8)
    ax.set_ylabel("CF-position RMSE (%EL), held-out lines")
    ax.set_title("Benchmark: CNN decoder vs Liu 2013 threshold model")
    fig.tight_layout()
    fig.savefig(FIG / "fig_cf_decoder.pdf")


def fig_gap_decode():
    d = json.loads((RES / "gap_decode.json").read_text())
    keys = [k for k in d if ":" in k]
    keys.sort()
    vals = [d[k]["median_abs_err_pctEL"] for k in keys]
    fig, ax = plt.subplots(figsize=(6.2, 3.4))
    ax.barh(keys, vals, color="C4")
    ax.axvline(1.0, color="k", ls="--", lw=1, label="Petkova 2019 4-gene $\\sim$1% EL")
    ax.set_xlabel("median |decode error| (%EL)")
    ax.legend(fontsize=7)
    ax.set_title("Positional decoding from gap-gene pairs (held-out embryos)")
    fig.tight_layout()
    fig.savefig(FIG / "fig_gap_decode.pdf")


def fig_grn_learning():
    d = json.loads((RES / "grn_learning_curve.json").read_text())
    n = [c["n_train"] for c in d["curve"]]
    g = [np.mean(c["acc_gnn"]) for c in d["curve"]]
    m = [np.mean(c["acc_mlp"]) for c in d["curve"]]
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    ax.plot(n, g, "o-", label="GNN (regulatory graph)")
    ax.plot(n, m, "s--", label="MLP (flat vector)")
    ax.set_xlabel("training states (class-balanced)")
    ax.set_ylabel("held-out accuracy")
    ax.set_ylim(0.75, 1.01)
    ax.legend(fontsize=8)
    ax.set_title("Basin-boundary learning curve, Cell Collective id-021")
    fig.tight_layout()
    fig.savefig(FIG / "fig_grn_learning.pdf")


def fig_turing():
    d = json.loads((RES / "turing_sweep.json").read_text())
    eps = [r["epsilon"] for r in d["rows"]]
    dev = [r["dev"] for r in d["rows"]]
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    ax.plot(eps, dev, "o-")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel(r"distance from onset $\varepsilon = 1 - \mu_a/\mu_{a,crit}$")
    ax.set_ylabel(r"$(k^2_{sim} - k^2_{lin})/k^2_{lin}$")
    ax.set_title("Beyond-onset wavenumber upshift (Gierer-Meinhardt)")
    fig.tight_layout()
    fig.savefig(FIG / "fig_turing_dev.pdf")


if __name__ == "__main__":
    fig_bcd_gradients(); print("1")
    fig_lambda_dosage(); print("2")
    fig_windows(); print("3")
    fig_cf_decoder(); print("4")
    fig_gap_decode(); print("5")
    fig_grn_learning(); print("6")
    fig_turing(); print("7")
