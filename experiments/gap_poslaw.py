#!/usr/bin/env python3
"""POSLAW-1: positional-decoding error scaling with gap-gene count.

Within each imaging session, decode position (kNN, k=5, embryo-held-out folds)
using every gene subset of size 1..d; median absolute error (%EL) vs subset size.
Fit log-log slope (error ~ d^-alpha). Published reference points:
Dubuis 2013 (single gap gene: several %EL), Petkova 2019 (4 gap genes: ~1%EL).
Controls: circular position permutation per embryo (structure-preserved null).
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

def knn_decode(trX, trP, teX, k=5, chunk=400):
    preds = []
    for i in range(0, len(teX), chunk):
        d2 = ((teX[i:i+chunk, None, :] - trX[None, :, :]) ** 2).sum(-1)
        preds.append(trP[np.argsort(d2, axis=1)[:, :k]].mean(1))
    return np.concatenate(preds)

def decode_err(M, genes_idx, seed, permute=False):
    emb, npos, _ = M.shape
    pos = AXIS[::STEP]
    X = M[:, :, genes_idx].reshape(emb * npos, len(genes_idx))
    P = np.tile(pos, emb); E = np.repeat(np.arange(emb), npos)
    ok = np.isfinite(X).all(1); X, P, E = X[ok], P[ok], E[ok]
    if permute:
        rng = np.random.default_rng(1000 + seed)
        for e in np.unique(E):
            m = E == e
            P[m] = np.roll(P[m], rng.integers(1, npos))
    rng = np.random.default_rng(seed)
    ids = rng.permutation(emb)
    te_e = set(ids[: emb // 5].tolist())
    te = np.array([e in te_e for e in E])
    if te.sum() == 0 or (~te).sum() == 0:
        return np.nan
    err = np.abs(knn_decode(X[~te], P[~te], X[te]) - P[te]) * 100.0  # fraction EL -> %EL
    return float(np.median(err))

def main():
    out = {"reference": {"dubuis_2013_single_gene": "several %EL",
                          "petkova_2019_four_genes_pctEL": 1.0},
           "sessions": {}}
    for s in load_if_sessions():
        if len(s["genes"]) < 3:
            continue
        M = session_matrix(s)
        if len(M) < 30:
            continue
        name = f"{s.get('line','?')}:{'+'.join(s['genes'])}"
        rec = {"n_embryos": int(len(M)), "genes": s["genes"], "by_size": {}}
        d = len(s["genes"])
        for sz in range(1, d + 1):
            errs, nulls = [], []
            for sub in itertools.combinations(range(d), sz):
                for seed in range(SEEDS):
                    e = decode_err(M, list(sub), seed)
                    if np.isfinite(e):
                        errs.append(e)
                        nulls.append(decode_err(M, list(sub), 0, permute=True)) if seed == 0 else None
            rec["by_size"][str(sz)] = {
                "median_abs_err_pctEL": float(np.median(errs)),
                "iqr": [float(np.percentile(errs, 25)), float(np.percentile(errs, 75))],
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
              "alpha=", round(rec.get("scaling_exponent_alpha", float("nan")), 3),
              "extrap@4=", round(rec.get("extrapolated_err_at_4genes_pctEL", float("nan")), 2))
    pathlib.Path("results").mkdir(exist_ok=True)
    pathlib.Path("results/gap_poslaw.json").write_text(json.dumps(out, indent=1))
    print("saved results/gap_poslaw.json")

if __name__ == "__main__":
    main()
