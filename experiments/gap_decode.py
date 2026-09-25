#!/usr/bin/env python3
"""Positional decoding from gap-gene pairs (Liu 2013 IF data).

Per embryo, per position: expression vector g(x) in R^d (d=2 genes).
Decoder: kNN regression on held-out embryos. Published comparison points:
Dubuis 2013 / Petkova 2019 report ~1%% EL for 4 gap genes jointly, several %%
for single genes. Here: pairs.
"""
import json
import pathlib

import numpy as np

from digitalembryo.gapgene import AXIS, load_if_sessions

ROOT = pathlib.Path(__file__).resolve().parent.parent
STEP = 10  # use every 10th of the 1000 axis points -> 100 positions


def session_matrix(session):
    """embryos x positions x genes (dorsal+ventral mean), positions subsampled."""
    gi = range(len(session["genes"]))
    mats = []
    for e in session["embryos"]:
        d, v = e["dorsal"], e["ventral"]
        if d.ndim != 2 or d.shape[1] != len(session["genes"]):
            continue
        prof = np.nanmean(np.stack([d, v]), axis=0)  # 1000 x genes
        mats.append(prof[::STEP])
    return np.array(mats)  # embryos x 100 x genes


def knn_decode(train_X, train_pos, test_X, k=5):
    """train_X: (Ntr, genes) with positions train_pos; predict pos for test."""
    d = ((test_X[:, None, :] - train_X[None, :, :]) ** 2).sum(-1)
    nn = np.argsort(d, axis=1)[:, :k]
    return train_pos[nn].mean(1)


def main():
    out = {}
    for s in load_if_sessions():
        if len(s["genes"]) != 2:
            continue
        M = session_matrix(s)
        if len(M) < 30:
            continue
        emb, npos, ng = M.shape
        pos = AXIS[::STEP]
        X = M.reshape(emb * npos, ng)
        P = np.tile(pos, emb)
        E = np.repeat(np.arange(emb), npos)
        ok = np.isfinite(X).all(1)
        X, P, E = X[ok], P[ok], E[ok]
        rng = np.random.default_rng(0)
        emb_ids = rng.permutation(emb)
        te_e = set(emb_ids[: emb // 5].tolist())
        te = np.array([e in te_e for e in E])
        pred = knn_decode(X[~te], P[~te], X[te])
        err = np.abs(pred - P[te])
        key = f"{s['line']}:{','.join(s['genes'])}"
        out[key] = {"n_embryos": int(emb), "n_test_pts": int(te.sum()),
                    "median_abs_err_pctEL": float(100 * np.median(err)),
                    "mean_abs_err_pctEL": float(100 * err.mean())}
        print(key, f"median|err|={100*np.median(err):.2f}%EL  mean={100*err.mean():.2f}%EL  (n_emb={emb})")
    out["published_context"] = {
        "petkova2019_4gapgenes": "~1% EL joint decoding (Cell 176:844)",
        "dubuis2013_single_gene": "single gap genes carry ~2-3 bits; %EL errors several x larger than joint"}
    (ROOT / "results" / "gap_decode.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
