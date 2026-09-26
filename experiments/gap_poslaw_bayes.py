#!/usr/bin/env python3
"""POSLAW-2: near-optimal (Gaussian MAP) decoding replication of POSLAW-1.

kNN is a lower bound on decoding precision (POSLAW-1 caveat). Here: per-position
Gaussian model - for each position j, expression vector ~ N(mu_j, Sigma_j) fit on
train embryos; decode test points by max posterior with flat prior over positions.
Diagonal-regularized covariance. Same sessions, subsets, folds, nulls as POSLAW-1.
Question: do the d^-2 scaling law and the dosage-shallowing discovery replicate
under near-optimal decoding, or were they kNN artifacts?
"""
import itertools, json, pathlib
import numpy as np
from digitalembryo.gapgene import AXIS, load_if_sessions

STEP = 20
SEEDS = 2

def session_matrix(session):
    mats = []
    for e in session["embryos"]:
        d, v = e["dorsal"], e["ventral"]
        if d.ndim != 2 or d.shape[1] != len(session["genes"]):
            continue
        mats.append(np.nanmean(np.stack([d, v]), axis=0)[::STEP])
    return np.array(mats)

def bayes_decode(trM, teM):
    """trM: (n_tr, npos, g); teM: (n_te, npos, g). Returns predicted pos index per test point."""
    n_tr, npos, g = trM.shape
    mu = np.nanmean(trM, axis=0)                      # npos x g
    # pooled covariance per position with shrinkage to diagonal
    cen = trM - mu[None]
    cov = np.array([np.cov(cen[:, j, :].T) for j in range(npos)])  # npos x g x g
    if g == 1:
        cov = cov.reshape(npos, 1, 1)
    reg = 1e-3 * np.trace(cov, axis1=1, axis2=2).mean()
    cov = cov + reg * np.eye(g)
    prec = np.linalg.inv(cov)                          # npos x g x g
    logdet = np.linalg.slogdet(cov)[1]
    te = teM.reshape(-1, g)                            # points x g
    diff = te[:, None, :] - mu[None, :, :]             # points x npos x g
    maha = np.einsum('pjg,jgk,pjk->pj', diff, prec, diff)
    ll = -0.5 * (maha + logdet[None, :])
    return np.nanargmax(ll, axis=1)

def decode_err(M, genes_idx, seed, permute=False):
    emb, npos, _ = M.shape
    X = M[:, :, genes_idx].copy()
    # impute missing values with the per-(position, gene) mean over embryos that have them
    fill = np.nanmean(X, axis=0)
    fill = np.where(np.isfinite(fill), fill, np.nanmean(X))
    inds = np.where(~np.isfinite(X))
    X[inds] = np.take(fill, inds[1] * X.shape[2] + inds[2]) if False else fill[inds[1], inds[2]]
    rng = np.random.default_rng(seed)
    ids = rng.permutation(emb)
    nte = max(1, emb // 5)
    tr, te = X[ids[nte:]], X[ids[:nte]]
    pred_idx = bayes_decode(tr, te)
    true_idx = np.tile(np.arange(npos), len(te))
    if permute:
        rng2 = np.random.default_rng(1000 + seed)
        true_idx = np.roll(true_idx, rng2.integers(1, npos))
    return float(np.median(np.abs(pred_idx - true_idx) * (100.0 / npos)))

def main():
    out = {"decoder": "Gaussian MAP (near-optimal), shrinkage-reg covariance",
           "sessions": {}}
    for s in load_if_sessions():
        if len(s["genes"]) < 3:
            continue
        M = session_matrix(s)
        if len(M) < 30:
            continue
        name = f"{s.get('line','?')}:{'+'.join(s['genes'])}"
        rec = {"n_embryos": int(len(M)), "by_size": {}}
        d = len(s["genes"])
        for sz in range(1, d + 1):
            errs, nulls = [], []
            for sub in itertools.combinations(range(d), sz):
                for seed in range(SEEDS):
                    e = decode_err(M, list(sub), seed)
                    if np.isfinite(e):
                        errs.append(e)
                        nulls.append(decode_err(M, list(sub), seed, permute=True))
            rec["by_size"][str(sz)] = {"median_abs_err_pctEL": float(np.median(errs)),
                                       "null_median_pctEL": float(np.median(nulls)),
                                       "n_runs": len(errs)}
        sizes = sorted(int(k) for k in rec["by_size"])
        xs = np.array(sizes, float)
        ys = np.array([rec["by_size"][str(k)]["median_abs_err_pctEL"] for k in sizes])
        if len(sizes) >= 2 and (ys > 0).all():
            a, b = np.polyfit(np.log(xs), np.log(ys), 1)
            rec["scaling_exponent_alpha"] = float(-a)
            rec["extrapolated_err_at_4genes_pctEL"] = float(np.exp(np.polyval([a, b], np.log(4.0))))
        out["sessions"][name] = rec
        print(name, {k: round(v["median_abs_err_pctEL"], 2) for k, v in rec["by_size"].items()},
              "nulls", {k: round(v["null_median_pctEL"], 1) for k, v in rec["by_size"].items()},
              "alpha=", round(rec.get("scaling_exponent_alpha", float("nan")), 3),
              "extrap@4=", round(rec.get("extrapolated_err_at_4genes_pctEL", float("nan")), 2), flush=True)
    pathlib.Path("results").mkdir(exist_ok=True)
    pathlib.Path("results/gap_poslaw_bayes.json").write_text(json.dumps(out, indent=1))
    print("saved results/gap_poslaw_bayes.json")

if __name__ == "__main__":
    main()
