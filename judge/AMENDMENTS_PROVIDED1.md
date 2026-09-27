# Lane-14 AMENDMENT QUEUE - LOCKED BEFORE EXECUTION
Verdict: PROVIDED round 1 (WhatsApp 10:42:00 IST, wamid...OTVCOTJCRgA=, verbatim in judge/round_provided1_verdict_whatsapp.txt).
Locked: 2026-09-27 10:44 IST, against 53pp build 2374c15. No execution began before this lock.

## Standing directive (hers, verbatim prefix)
"NO NEGATIVE RESULTS SHOULD BE THERE" - every negative becomes a positive-framed finding (mechanism/limit/threshold), never hidden.

## Weakness items (her numbering 2-19) mapped vs current paper
- #2 biological novelty unclear: PARTIAL (Study C is real Bicoid/gap-gene biology) -> reframe + tie each result to fate stability/transitions/robustness
- #3 real developmental system for Study A: NEW (compare Turing surrogate vs limb/pigmentation/embryonic patterning data)
- #4 Boolean vs continuous ODE: NEW (attractor survival under ODE)
- #5 ML question explicit: PARTIAL (sample-complexity section exists) -> sharpen to "observations needed to infer hidden landscapes"
- #6 second developmental network: NEW (generalization test)
- #7 benchmark vs random/simple/standard tools: PARTIAL (GNN vs MLP) -> add random sampling + standard network analysis
- #8 biological interpretation per result: EDITORIAL pass
- #9 uncertainty analysis (edges/rules/thresholds): NEW
- #10 sync vs async/stochastic updates: NEW
- #11 basin-size sensitivity: NEW
- #12 Study C grounding: PARTIAL (error in %EL) -> tie to embryo-to-embryo variation/developmental precision
- #13 headline discovery: EDITORIAL (pick flagship)
- #14 ML vs biological baselines: PARTIAL (kNN baseline exists) -> add marker/threshold baselines
- #15 generalization other organism/system: NEW
- #16 perturbation prediction: PARTIAL (bcd dosage) -> add gene-change -> landscape-change prediction
- #17 package as tool (CLI/docs/examples): NEW
- #18 rediscovery framing: EDITORIAL (hidden-property quantification, not "reproduced bistability")
- #19 CNN failure interpretation: PARTIAL (surrogate negative disclosed) -> 4-model comparison (CNN/PINN/GNN/direct-parameter)

## Five negative-result reframes (her NEGATIVE RESULT 1-3 + named reframes)
1. CNN wavenumber failure -> "local spatial representations insufficient for emergent patterns" + information-need question
2. Decoder limits -> minimum-information threshold (information-bottleneck sweep: full/noisy/partial/sparse gradient)
3. Boolean simplification -> robustness-to-model-uncertainty result
4. AI-vs-mechanistic comparison -> explicit
5. Small wild-type basin (0.77%) -> developmental fragility index

## Redesign spine (her strongest-pivot)
Title direction: "Computational Limits and Robustness of Predicting Developmental Cell-Fate Landscapes"
One story: A (why image AI fails emergent patterns) -> B (topology creates robust fate landscapes) -> C (minimum information to decode states).
RULE: title pivot lands only after the supporting analyses exist (fragility index + information threshold first) - no claiming a landscape before computing it.

## NOT-to-do (hers, locked)
- No "CNN failed, but we tried another model" engineering framing
- No "our model predicts development" overreach
- No adding ML models without biological interpretation

## First deliverables (cheap, existing data)
1. Developmental fragility index from existing basin enumeration (reframe 5 / #11 sensitivity seed)
2. Information-threshold sweep on Study C decoder (reframe 2): full -> noisy -> partial -> sparse Bicoid gradient, find critical threshold
