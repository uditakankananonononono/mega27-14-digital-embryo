#!/usr/bin/env python3
"""POSLAW-3: position-resolved decoding error - WHERE does the joint-decoding
gain live along the AP axis?

Same sessions/data/harness as POSLAW-1 (gap_poslaw.py). For 1-gene vs 3-gene
decoding, per-sample absolute errors are binned by TRUE position into 5 bands
across the axis. Hypothesis: the joint-decoding gain concentrates mid-axis,
where gap-gene expression boundaries carry position information; ends of the
axis are decodable from single genes already.
"""
import json, pathlib
import numpy as np
from digitalembryo.gapgene import AXIS, load_if_sessions
from gap_poslaw import session_matrix, knn_decode

STEP = 20
NBANDS = 5

def band_errs(M, genes_idx, seed=0):
    emb, npos, _ = M.shape
    pos = AXIS[::STEP]
    X = M[:, :, genes_idx].reshape(emb * npos, len(genes_idx))
    P = np.tile(pos, emb); E = np.repeat(np.arange(emb), npos)
    ok = np.isfinite(X).all(1); X, P, E = X[ok], P[ok], E[ok]
    rng = np.random.default_rng(seed)
    ids = rng.permutation(emb)
    te_e = set(ids[: emb // 5].tolist())
    te = np.array([e in te_e for e in E])
    if te.sum() == 0 or (~te).sum() == 0:
        return None
    pred = knn_decode(X[~te], P[~te], X[te])
    err = np.abs(pred - P[te]) * 100.0
    bands = np.minimum((P[te] * NBANDS).astype(int), NBANDS - 1)
    return err, bands

def main():
    out = {"bands": NBANDS, "band_edges_pctEL": [round(100 * i / NBANDS, 1) for i in range(NBANDS + 1)],
           "sessions": {}}
    for s in load_if_sessions():
        if len(s["genes"]) < 3:
            continue
        M = session_matrix(s)
        if len(M) < 30:
            continue
        name = f"{s.get('line','?')}:{'+'.join(s['genes'])}"
        rec = {"n_embryos": int(len(M)), "by_band": []}
        e1_all = np.zeros((NBANDS, 0)); e3_all = np.zeros((NBANDS, 0))
        # best single gene = Hb (canonical) plus all singles; use median over singles per band
        singles, triple = [], list(range(len(s["genes"])))
        for g in range(len(s["genes"])):
            r = band_errs(M, [g])
            if r: singles.append(r)
        r3 = band_errs(M, triple)
        for b in range(NBANDS):
            meds1 = [float(np.median(e[bnd == b])) for e, bnd in singles if (bnd == b).sum() >= 5]
            med3 = float(np.median(r3[0][r3[1] == b])) if (r3[1] == b).sum() >= 5 else float("nan")
            rec["by_band"].append({"band": b,
                "single_gene_median_pctEL": float(np.median(meds1)) if meds1 else None,
                "three_gene_median_pctEL": med3,
                "gain_ratio": (float(np.median(meds1)) / med3) if meds1 and med3 and med3 > 0 else None})
        out["sessions"][name] = rec
        print(name)
        for row in rec["by_band"]:
            print("  band", row["band"], "1g:", None if row["single_gene_median_pctEL"] is None else round(row["single_gene_median_pctEL"], 2),
                  "3g:", None if row["three_gene_median_pctEL"] is None else round(row["three_gene_median_pctEL"], 2),
                  "gain:", None if row["gain_ratio"] is None else round(row["gain_ratio"], 2))
    pathlib.Path("results/gap_poslaw_posres.json").write_text(json.dumps(out, indent=1))
    print("saved results/gap_poslaw_posres.json")

if __name__ == "__main__":
    main()
