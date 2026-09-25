"""Positional decoding: predict cephalic-furrow position from Bcd gradients.

Head-to-head against Liu et al. PNAS 2013 threshold model xCF = Sx*ln(D)+c.
"""
from __future__ import annotations

import numpy as np

from .bcd import fit_single, load_live_gradients

GRID = np.linspace(0.0, 0.6, 64)


def embryo_table():
    """One row per embryo: mean gradient on a fixed grid, fit params, CF."""
    recs = load_live_gradients()
    by_emb = {}
    for r in recs:
        by_emb.setdefault((r["line"], r["embryo"]), []).append(r)
    rows = []
    for (line, idx), sides in by_emb.items():
        cf = sides[0]["cf"]
        if not np.isfinite(cf):
            continue
        profs, fits = [], []
        for s in sides:
            x, y = s["x"], s["intensity"]
            m = (x >= 0) & (x <= 0.6) & (y > 0)
            if m.sum() < 15:
                continue
            profs.append(np.interp(GRID, x[m], np.log(y[m])))
            fits.append(fit_single(x[m], y[m]))
        if not profs:
            continue
        A = float(np.mean([f["A"] for f in fits]))
        lam = float(np.mean([f["lam"] for f in fits]))
        rows.append({"line": line, "embryo": idx, "cf": float(cf) / 100.0,  # stored in %EL in the .mat
                     "A": A, "lam": lam, "profile": np.mean(profs, axis=0)})
    return rows


def liu_baseline(rows, train, test):
    """Liu 2013: xCF = Sx * ln(D_proxy) + c, D proxied by gradient amplitude A."""
    D = np.array([np.log(r["A"]) for r in rows])
    y = np.array([r["cf"] for r in rows])
    a, b = np.polyfit(D[train], y[train], 1)
    return a * D[test] + b


def ridge_features(rows, train, test):
    X = np.array([[np.log(r["A"]), r["lam"]] for r in rows])
    y = np.array([r["cf"] for r in rows])
    Xm = np.column_stack([np.ones(len(rows)), X])
    lam_r = 1e-6
    w = np.linalg.solve(Xm[train].T @ Xm[train] + lam_r * np.eye(3), Xm[train].T @ y[train])
    return Xm[test] @ w


def cnn_decoder(rows, train, test, epochs=300, seed=0, lr=5e-3):
    """Small 1D-CNN on the log-intensity profile (fixed grid)."""
    import torch

    torch.manual_seed(seed)
    X = torch.tensor(np.stack([r["profile"] for r in rows]), dtype=torch.float32)
    X = (X - X[train].mean()) / (X[train].std() + 1e-6)
    y = torch.tensor([r["cf"] for r in rows], dtype=torch.float32)

    net = torch.nn.Sequential(
        torch.nn.Conv1d(1, 16, 5, padding=2), torch.nn.ReLU(), torch.nn.MaxPool1d(2),
        torch.nn.Conv1d(16, 32, 5, padding=2), torch.nn.ReLU(), torch.nn.MaxPool1d(2),
        torch.nn.Flatten(), torch.nn.Linear(32 * 16, 64), torch.nn.ReLU(),
        torch.nn.Linear(64, 1))
    opt = torch.optim.Adam(net.parameters(), lr=lr)
    lossf = torch.nn.MSELoss()
    Xt, yt = X[train].unsqueeze(1), y[train]
    net.train()
    for _ in range(epochs):
        opt.zero_grad()
        loss = lossf(net(Xt).squeeze(-1), yt)
        loss.backward()
        opt.step()
    net.eval()
    with torch.no_grad():
        return net(X[test].unsqueeze(1)).squeeze(-1).numpy()


def grouped_folds(rows, k=5, seed=0):
    """k folds with disjoint fly lines (held-out-line generalization)."""
    rng = np.random.default_rng(seed)
    lines = np.array(sorted({r["line"] for r in rows}))
    rng.shuffle(lines)
    rlines = np.array([r["line"] for r in rows])
    folds = []
    for i in range(k):
        test_lines = set(lines[i::k])
        test = np.array([j for j, ln in enumerate(rlines) if ln in test_lines])
        train = np.array([j for j, ln in enumerate(rlines) if ln not in test_lines])
        folds.append((train, test))
    return folds
