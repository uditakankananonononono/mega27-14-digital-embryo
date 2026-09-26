# POSLAW-1: positional-decoding error scaling with gap-gene count (2026-09-26)

Experiment: `experiments/gap_poslaw.py` (kNN k=5, embryo-held-out 20% folds, 2 seeds,
50 positions per embryo, every subset of size 1..3, circular-permutation nulls).
Data: Liu 2013 IF sessions with >=3 gap genes (Hb, Eve, Kr) at 1x and 2x Bcd dosage.

## Results (median abs decoding error, % egg length)

| session | 1 gene | 2 genes | 3 genes | exponent alpha | extrap @4 genes |
|---|---|---|---|---|---|
| 1XA (1x Bcd, n=244) | 21.6 | 6.2 | 2.4 | 1.98 | 1.4 |
| 2XA (2x Bcd, n=142) | 21.0 | 8.0 | 4.2 | 1.46 | 2.8 |

Permutation nulls (positions circularly shifted per embryo): 24-26 %EL at every
subset size. Real decoding beats the null 1.1x (1 gene), 3-4x (2 genes), ~10x
(3 genes). HONEST CAVEAT: single-gap-gene decoding on this data is barely better
than a structure-preserved null (21.6 vs 24.4 %EL) - the multi-gene synergy, not
single-gene precision, is where the positional information lives here.

## Claims and caveats

- Error falls as ~ d^-2 at 1x dosage; the d=4 extrapolation lands at 1.4 %EL,
  consistent with Petkova 2019's ~1 %EL for joint 4-gap-gene decoding (Cell 176:844)
  on different data with Bayes-optimal decoding. Ours is a kNN lower bound on
  precision, so the agreement is conservative. NOT a beat claim: different dataset,
  simpler decoder, extrapolation by one step only.
- DISCOVERY CANDIDATE: the 2x-dosage line has a SHALLOWER scaling exponent
  (1.46 vs 1.98) and worse 3-gene precision (4.2 vs 2.4 %EL). Doubling Bcd does not
  add positional information through the gap layer - it degrades integration.
  Connects to the repo's DCLS result (Bcd length scale is dosage-dependent).
- Caveats: only two 3-gene sessions exist; single-gene 21 %EL is worse than
  published single-gene estimates (several %EL) because this IF data is coarser;
  exponent fit is on 3 points.

## POSLAW-2 (2026-09-26): near-optimal decoder replication
experiments/gap_poslaw_bayes.py - Gaussian MAP decoder (per-position mean/cov with
shrinkage; train-mean imputation), same sessions/folds/nulls.

| session | 1 gene | 2 genes | 3 genes | alpha | extrap @4 |
|---|---|---|---|---|---|
| 1XA (kNN) | 21.6 | 6.2 | 2.4 | 1.98 | 1.4 |
| 1XA (Bayes) | 18.0 | 2.0 | 2.0 | 2.13 | 0.78 |
| 2XA (kNN) | 21.0 | 8.0 | 4.2 | 1.46 | 2.8 |
| 2XA (Bayes) | 13.0 | 4.0 | 2.0 | 1.70 | 1.23 |

ROBUST: the dosage-shallowing of the scaling exponent replicates across decoders
(1.98->1.46 kNN; 2.13->1.70 Bayes). FRAGILE: the worse 3-gene precision at 2x
dosage (2.4->4.2 under kNN) disappears under Bayes (2.0 both) - decoder-dependent,
downgraded from finding to kNN-specific observation. CAVEAT: 50 positions quantize
errors to 2% EL steps; the Bayes 2/3-gene values sit at the quantization floor.
