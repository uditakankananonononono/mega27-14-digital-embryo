#!/usr/bin/env python3
"""4-model comparison on cached Turing sims (verdict #19): predict the relative
position of the simulated wavenumber within the linearly unstable band.
Models: (1) linear theory alone (fastest linear mode), (2) CNN on dispersion
curve (the retracted n=58 winner), (3) MLP on the 4 physical parameters,
(4) ridge on hand features of the curve (peak position, width, height).
6-fold CV x 2 seeds, relative RMSE. All honest, same protocol.
"""
import json
import numpy as np

ROWS = json.load(open("results/turing_sims_cache.json"))

def targets(rows):
    return np.array([(r["k2_sim"] - r["k2_min"]) / (r["k2_max"] - r["k2_min"]) for r in rows])

def relrmse(pred, true):
    return float(np.sqrt(np.mean(((pred - true) / (np.abs(true) + 1e-9)) ** 2)))

def curve_feats(r):
    c = np.array(r["curve"]); n = len(c)
    pk = float(np.argmax(c)) / (n - 1)
    half = c >= 0.5 * c.max()
    width = float(half.sum()) / n
    centroid = float((c * np.arange(n)).sum() / (c.sum() + 1e-12)) / n
    return [pk, width, float(c.max()), centroid]

def kfolds(n, k=6, seed=0):
    rng = np.random.default_rng(seed)
    idx = rng.permutation(n)
    return np.array_split(idx, k)

def main():
    y = targets(ROWS)
    # model 1: linear theory prediction of relative position
    y_lin = np.array([(r["k2_lin"] - r["k2_min"]) / (r["k2_max"] - r["k2_min"]) for r in ROWS])
    res = {"linear_theory_relRMSE": relrmse(y_lin, y)}
    # models 2-4 with CV
    import torch
    F = np.array([curve_feats(r) for r in ROWS])
    P = np.array([np.log10(r["params"]) for r in ROWS])
    C = np.array([r["curve"] for r in ROWS], dtype=np.float32)
    cnn_scores, mlp_scores, ridge_scores = [], [], []
    for seed in (0, 1):
        for te in kfolds(len(ROWS), 6, seed):
            tr = np.setdiff1d(np.arange(len(ROWS)), te)
            torch.manual_seed(seed)
            # CNN
            net = torch.nn.Sequential(
                torch.nn.Conv1d(1, 16, 5, padding=2), torch.nn.ReLU(), torch.nn.MaxPool1d(2),
                torch.nn.Conv1d(16, 32, 5, padding=2), torch.nn.ReLU(), torch.nn.MaxPool1d(2),
                torch.nn.Flatten(), torch.nn.Linear(32 * 16, 64), torch.nn.ReLU(),
                torch.nn.Linear(64, 1))
            opt = torch.optim.Adam(net.parameters(), lr=1e-3)
            Xtr = torch.tensor(C[tr]).unsqueeze(1); ytr = torch.tensor(y[tr], dtype=torch.float32)
            for ep in range(300):
                opt.zero_grad(); loss = ((net(Xtr).squeeze(-1) - ytr) ** 2).mean(); loss.backward(); opt.step()
            with torch.no_grad():
                pred = net(torch.tensor(C[te]).unsqueeze(1)).squeeze(-1).numpy()
            cnn_scores.append(relrmse(pred, y[te]))
            # MLP on params
            mlp = torch.nn.Sequential(torch.nn.Linear(4, 32), torch.nn.ReLU(), torch.nn.Linear(32, 32), torch.nn.ReLU(), torch.nn.Linear(32, 1))
            opt = torch.optim.Adam(mlp.parameters(), lr=1e-3)
            Ptr = torch.tensor(P[tr], dtype=torch.float32)
            for ep in range(300):
                opt.zero_grad(); loss = ((mlp(Ptr).squeeze(-1) - ytr) ** 2).mean(); loss.backward(); opt.step()
            with torch.no_grad():
                pred = mlp(torch.tensor(P[te], dtype=torch.float32)).squeeze(-1).numpy()
            mlp_scores.append(relrmse(pred, y[te]))
            # ridge on curve features
            mu, sd = F[tr].mean(0), F[tr].std(0) + 1e-9
            Ztr = (F[tr] - mu) / sd; Zte = (F[te] - mu) / sd
            A = Ztr.T @ Ztr + 1.0 * np.eye(Ztr.shape[1])
            w = np.linalg.solve(A, Ztr.T @ y[tr])
            ridge_scores.append(relrmse(Zte @ w, y[te]))
    res.update({
        "cnn_curve_relRMSE": round(float(np.mean(cnn_scores)), 4),
        "mlp_params_relRMSE": round(float(np.mean(mlp_scores)), 4),
        "ridge_curve_features_relRMSE": round(float(np.mean(ridge_scores)), 4),
        "linear_theory_relRMSE": round(res["linear_theory_relRMSE"], 4),
        "protocol": "6-fold CV x 2 seeds; target = relative position of k2_sim within [k2_min, k2_max]; n=171 cached sims",
        "reading": "if no model beats linear theory, the beyond-onset shift is not learnable from these inputs - a boundary result, not a modeling failure"})
    json.dump(res, open("results/turing_4model.json", "w"), indent=1)
    print(json.dumps(res, indent=1))

if __name__ == "__main__":
    main()
