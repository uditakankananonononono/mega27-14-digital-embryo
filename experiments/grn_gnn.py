#!/usr/bin/env python3
"""GNN vs MLP: predict basin fate (WG-protein ON at attractor) from initial
state on the real Cell Collective id-021 network (bistable condition)."""
import json
import pathlib

import numpy as np

from digitalembryo.bnet import parse_bnet
from digitalembryo.gnn import GNNClassifier, MLPBaseline, basin_labels, graph_arrays

ROOT = pathlib.Path(__file__).resolve().parent.parent
BNET = (ROOT / "data/raw/grn_models/id021_bodysegmentation.bnet").read_text()


def main():
    rules = parse_bnet(BNET)
    ext = {"v_hh_external": 1, "v_WG_external": 0, "v_SLP": 0}  # bistable in wg
    free, A = graph_arrays(rules, ext)
    X, y, free = basin_labels(rules, ext, n_samples=2000, seed=0)
    # replicate initial-state vector across nodes as per-node feature plane
    rng = np.random.default_rng(0)
    idx = rng.permutation(len(X))
    tr, te = idx[:1400], idx[1400:]

    gnn = GNNClassifier(len(free), seed=0)
    gnn.fit(A, X[tr], y[tr], epochs=300)
    acc_gnn = float((gnn.predict(A, X[te]) == y[te]).mean())

    mlp = MLPBaseline(X.shape[1], seed=0)
    mlp.fit(X[tr], y[tr], epochs=300)
    acc_mlp = float((mlp.predict(X[te]) == y[te]).mean())

    maj = float(max(y[te].mean(), 1 - y[te].mean()))
    out = {"model": "Cell Collective id-021, condition hh1_WG0_SLP0 (bistable in wg)",
           "n_samples": len(X), "n_test": len(te),
           "positive_rate": float(y.mean()),
           "acc_majority": maj, "acc_mlp": acc_mlp, "acc_gnn": acc_gnn}
    (ROOT / "results" / "grn_gnn.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
