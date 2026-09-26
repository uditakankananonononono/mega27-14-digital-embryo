#!/usr/bin/env python3
"""POSLAW-4: is the pole rescue (POSLAW-3) decoder-dependent?

POSLAW-2 showed the dosage effect on ABSOLUTE precision was kNN-specific while
the exponent-shallowing replicated under MAP. POSLAW-3's pole rescue (3-gene
flattens the 1-gene U-shaped error profile) was measured with kNN only. Here
the identical position-resolved band analysis is rerun with the Gaussian MAP
decoder from POSLAW-2. If the pole rescue survives MAP, it is a data property;
if it vanishes, it was a kNN artifact (like absolute precision at 2x).
Same sessions, STEP=20, embryo-held-out 80/20, 5 bands.
"""
import json, pathlib
import numpy as np
from digitalembryo.gapgene import AXIS, load_if_sessions
from gap_poslaw import session_matrix
from gap_poslaw_bayes import bayes_decode

STEP = 20
NBANDS = 5

def band_errs_bayes(M, genes_idx, seed=0):
    emb, npos, _ = M.shape
    X = M[:, :, genes_idx].copy()
    fill = np.nanmean(X, axis=0)
    fill = np.where(np.isfinite(fill), fill, np.nanmean(X))
    inds = np.where(~np.isfinite(X))
    X[inds] = fill[inds[1], inds[2]]
    rng = np.random.default_rng(seed)
    ids = rng.permutation(emb)
    nte = max(1, emb // 5)
    tr, te = X[ids[nte:]], X[ids[:nte]]
    pred_idx = bayes_decode(tr, te)
    true_idx = np.tile(np.arange(npos), len(te))
    err = np.abs(pred_idx - true_idx) * (100.0 / npos)
    bands = np.minimum((true_idx * NBANDS / npos).astype(int), NBANDS - 1)
    return err, bands

def main():
    out = {"decoder": "Gaussian MAP", "bands": NBANDS,
           "band_edges_pctEL": [round(100 * i / NBANDS, 1) for i in range(NBANDS + 1)],
           "sessions": {}}
    for s in load_if_sessions():
        if len(s["genes"]) < 3:
            continue
        M = session_matrix(s)
        if len(M) < 30:
            continue
        name = f"{s.get('line','?')}:{'+'.join(s['genes'])}"
        rec = {"n_embryos": int(len(M)), "by_band": []}
        singles = []
        for g in range(len(s["genes"])):
            singles.append(band_errs_bayes(M, [g]))
        r3 = band_errs_bayes(M, list(range(len(s["genes"]))))
        for b in range(NBANDS):
            meds1 = [float(np.median(e[bnd == b])) for e, bnd in singles if (bnd == b).sum() >= 3]
            med3 = float(np.median(r3[0][r3[1] == b])) if (r3[1] == b).sum() >= 3 else float("nan")
            rec["by_band"].append({"band": b,
                "single_gene_median_pctEL": float(np.median(meds1)) if meds1 else None,
                "three_gene_median_pctEL": med3,
                "gain_ratio": (float(np.median(meds1)) / med3) if meds1 and med3 and med3 > 0 else None})
        out["sessions"][name] = rec
        print(name, flush=True)
        for row in rec["by_band"]:
            print("  band", row["band"], "1g:", None if row["single_gene_median_pctEL"] is None else round(row["single_gene_median_pctEL"], 2),
                  "3g:", None if row["three_gene_median_pctEL"] is None else round(row["three_gene_median_pctEL"], 2),
                  "gain:", None if row["gain_ratio"] is None else round(row["gain_ratio"], 2), flush=True)
    pathlib.Path("results/gap_poslaw_posres_bayes.json").write_text(json.dumps(out, indent=1))
    print("saved results/gap_poslaw_posres_bayes.json")

if __name__ == "__main__":
    main()
