#!/usr/bin/env python3
"""CNN surrogate benchmark: predict steady-state pattern wavenumber from the
linear dispersion curve; head-to-head vs the linear fastest-mode predictor.
"""
import json
import multiprocessing as mp
import pathlib

import numpy as np

from digitalembryo.morphogen import GiererMeinhardt2D, turing_thresholds
from digitalembryo.turing_ext import dispersion, linear_fastest_k2

ROOT = pathlib.Path(__file__).resolve().parent.parent
GRID = 64


def _point(i):
    rng = np.random.default_rng(1000 + i)
    Da = 10 ** rng.uniform(-2.6, -2.0)   # 0.0025 - 0.01
    Dh = 10 ** rng.uniform(-0.9, -0.3)   # 0.126 - 0.5
    mu_a = 10 ** rng.uniform(-2.0, -1.0)  # 0.01 - 0.1
    mu_h = 10 ** rng.uniform(-1.2, -0.6)  # 0.063 - 0.25
    th = turing_thresholds(Da, Dh, mu_a, mu_h)
    if not th["unstable"]:
        return None
    band = np.linspace(th["k2_min"], th["k2_max"], GRID)
    sig = dispersion(Da, Dh, mu_a, mu_h, 1.0, band)
    k2_lin = linear_fastest_k2(Da, Dh, mu_a, mu_h)
    sim = GiererMeinhardt2D(n=48, Da=Da, Dh=Dh, mu_a=mu_a, mu_h=mu_h, seed=i)
    sim.step(2500)
    k2_sim = float(sim.pattern_wavenumber())
    return {"curve": (sig / (sig.max() + 1e-12)).tolist(), "k2_min": th["k2_min"],
            "k2_max": th["k2_max"], "k2_lin": k2_lin, "k2_sim": k2_sim,
            "params": [Da, Dh, mu_a, mu_h]}


def cnn_fit_predict(train, test, epochs=400, seed=0):
    import torch

    torch.manual_seed(seed)
    X = torch.tensor([r["curve"] for r in ROWS], dtype=torch.float32)
    # predict relative position of k2 within the band (bounded target)
    y = torch.tensor([(r["k2_sim"] - r["k2_min"]) / (r["k2_max"] - r["k2_min"])
                      for r in ROWS], dtype=torch.float32)
    net = torch.nn.Sequential(
        torch.nn.Conv1d(1, 16, 5, padding=2), torch.nn.ReLU(), torch.nn.MaxPool1d(2),
        torch.nn.Conv1d(16, 32, 5, padding=2), torch.nn.ReLU(), torch.nn.MaxPool1d(2),
        torch.nn.Flatten(), torch.nn.Linear(32 * 16, 64), torch.nn.ReLU(),
        torch.nn.Linear(64, 1))
    opt = torch.optim.Adam(net.parameters(), lr=3e-3)
    lossf = torch.nn.MSELoss()
    for _ in range(epochs):
        opt.zero_grad()
        loss = lossf(net(X[train].unsqueeze(1)).squeeze(-1), y[train])
        loss.backward()
        opt.step()
    with torch.no_grad():
        frac = net(X[test].unsqueeze(1)).squeeze(-1).numpy()
    return np.array([ROWS[i]["k2_min"] + f * (ROWS[i]["k2_max"] - ROWS[i]["k2_min"])
                     for i, f in zip(test, frac)])


if __name__ == "__main__":
    global ROWS
    with mp.Pool(2) as pool:
        ROWS = [r for r in pool.map(_point, range(72)) if r]
    rng = np.random.default_rng(0)
    idx = rng.permutation(len(ROWS))
    k = 6
    folds = [(np.setdiff1d(idx, idx[i::k]), idx[i::k]) for i in range(k)]
    y = np.array([r["k2_sim"] for r in ROWS])
    ylin = np.array([r["k2_lin"] for r in ROWS])
    pred_cnn = np.full(len(ROWS), np.nan)
    for train, test in folds:
        pred_cnn[test] = np.mean([cnn_fit_predict(train, test, seed=s) for s in (0, 1)], axis=0)
    rel = lambda p: float(np.sqrt(np.mean(((p - y) / y) ** 2)))
    out = {
        "n_sims": len(ROWS),
        "cv": "6-fold over parameter sets",
        "linear_fastest_mode_relRMSE": rel(ylin),
        "cnn_dispersion_relRMSE": rel(pred_cnn),
        "improvement": rel(ylin) / rel(pred_cnn),
    }
    (ROOT / "results" / "turing_cnn.json").write_text(json.dumps(out, indent=1))
    np.save(ROOT / "results" / "turing_cnn_preds.npy",
            np.column_stack([y, ylin, pred_cnn]))
    print(json.dumps(out, indent=1))
