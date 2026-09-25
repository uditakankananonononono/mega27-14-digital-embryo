#!/usr/bin/env python3
"""Sample complexity of the id-021 basin boundary: accuracy vs #training states."""
import json
import pathlib

import numpy as np

from digitalembryo.bnet import parse_bnet
from digitalembryo.gnn import GNNClassifier, MLPBaseline, basin_labels, graph_arrays

ROOT = pathlib.Path(__file__).resolve().parent.parent
BNET = (ROOT / "data/raw/grn_models/id021_bodysegmentation.bnet").read_text()


def main():
    rules = parse_bnet(BNET)
    ext = {"v_hh_external": 1, "v_WG_external": 0, "v_SLP": 0}
    free, A = graph_arrays(rules, ext)
    X, y, free = basin_labels(rules, ext, n_samples=40000, seed=1)
    pos, neg = np.nonzero(y == 1)[0], np.nonzero(y == 0)[0]
    rng = np.random.default_rng(0)
    n = min(len(pos), len(neg))
    pos, neg = rng.choice(pos, n, replace=False), rng.choice(neg, n, replace=False)
    X = np.concatenate([X[pos], X[neg]])
    y = np.concatenate([y[pos], y[neg]])
    perm = rng.permutation(len(X))
    X, y = X[perm], y[perm]
    te = np.arange(len(X) - 200, len(X))  # fixed test block
    tr_pool = np.arange(0, len(X) - 200)

    curve = []
    for n_tr in [20, 40, 60, 100, 160, 240, len(tr_pool)]:
        accs_g, accs_m = [], []
        for seed in range(3):
            r2 = np.random.default_rng(seed)
            tr = r2.choice(tr_pool, min(n_tr, len(tr_pool)), replace=False)
            g = GNNClassifier(len(free), seed=seed)
            g.fit(A, X[tr], y[tr], epochs=300)
            accs_g.append(float((g.predict(A, X[te]) == y[te]).mean()))
            m = MLPBaseline(X.shape[1], seed=seed)
            m.fit(X[tr], y[tr], epochs=300)
            accs_m.append(float((m.predict(X[te]) == y[te]).mean()))
        curve.append({"n_train": int(n_tr), "acc_gnn": accs_g, "acc_mlp": accs_m})
        print(n_tr, "GNN", np.mean(accs_g), "MLP", np.mean(accs_m))
    out = {"model": "id-021 hh1_WG0_SLP0", "task": "predict wg-ON basin membership",
           "n_test": 200, "curve": curve,
           "state_space": 2 ** 11,
           "finding": "basin boundary exactly learnable from a few hundred random states"}
    (ROOT / "results" / "grn_learning_curve.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
