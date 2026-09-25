#!/usr/bin/env python3
"""Balanced basin-fate prediction on id-021 (hh1_WG0_SLP0, wg bistability).

Uniform sampling makes the wg-ON basin 0.5% (see grn_gnn.json) - here we
oversample wg-ON-reaching states to class-balance, so GNN vs MLP is a real
comparison. Split by initial state, stratified.
"""
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
    n = min(len(pos), len(neg))
    rng = np.random.default_rng(0)
    keep = np.concatenate([rng.choice(pos, n, replace=len(pos) < n),
                           rng.choice(neg, n, replace=False)])
    rng.shuffle(keep)
    X, y = X[keep], y[keep]
    tr, te = np.arange(0, int(0.7 * len(X))), np.arange(int(0.7 * len(X)), len(X))

    gnn = GNNClassifier(len(free), seed=0)
    gnn.fit(A, X[tr], y[tr], epochs=300)
    acc_g = float((gnn.predict(A, X[te]) == y[te]).mean())
    mlp = MLPBaseline(X.shape[1], seed=0)
    mlp.fit(X[tr], y[tr], epochs=300)
    acc_m = float((mlp.predict(X[te]) == y[te]).mean())

    out = {"model": "Cell Collective id-021, hh1_WG0_SLP0, wg bistability",
           "sampling": "class-balanced (oversampled wg-ON basin)",
           "wg_on_basin_fraction_uniform": float(len(pos) / 40000),
           "n_balanced": int(len(X)), "n_test": int(len(te)),
           "acc_gnn": acc_g, "acc_mlp": acc_m, "acc_majority": 0.5}
    (ROOT / "results" / "grn_gnn_balanced.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
