# mega27-14-digital-embryo
Three computational embryogenesis studies with analytic validation:
- **A**: Turing reaction-diffusion (Gierer-Meinhardt) - simulated pattern wavenumber validated against the analytic instability band.
- **B**: GRN differentiation - Gardner toggle switch bistability + Derrida ordered/chaotic transition (theory vs measurement).
- **C**: Positional information - Wolpert French flag, Bicoid-like gradients, error scaling and mutual information.
Run: `pip install -e . && pytest && python experiments/run_all.py`

## CLI tool (embryosim)
`pip install -e .` exposes the `embryosim` command:
- `embryosim fit-bcd [--line 2XA]` - fit SDD exponentials to real Bcd gradients (Liu 2013)
- `embryosim hb-boundary` - Hb boundary positions vs Bcd dosage across fly lines
- `embryosim turing [--da 0.005 --dh 0.2 ...]` - Gierer-Meinhardt pattern vs analytic Turing band
- `embryosim attractors [--model FILE.bnet --fix v_SLP=0 ...]` - exact attractor enumeration
- `embryosim fragility` - developmental fragility index (reproduces the paper's 0.781% WT basin, F=0.500)
- `embryosim info-threshold` - restricted-input CF decode (separate full-profile CNN benchmark: 2.38% EL RMSE vs the 4.94% EL true-dose Liu threshold comparator, but the measured-amplitude threshold comparator in the same benchmark is 2.8566% EL and ridge is 2.55%, so the margin depends on the comparator chosen; see results/cf_decoder.json. The single-point (4.29% EL) and anterior-window (2.99% EL) restricted-input decoders do not beat the measured-amplitude comparator. Not a biological minimum-information claim)
Every command prints JSON with the dataset/provenance fields named inline.

### Fixed-coordinate transient linear amplification
The full finite-difference linearizations have negative spectral abscissa, but raw-concentration Euclidean propagator norms at frozen times1/10/100 are9.641/36.778/6.206. Positive symmetric-part logarithmic norm5.283 supplies a distinct non-normal transient-growth diagnostic. No nonlinear/positivity-feasible disturbance or biological robustness is certified; nativeJacobian failure remains. See `results/spatial_ode_nonnormal_audit.json`.
The time10 leading singular direction's initial nonnegative amplitude ceiling is1.497e-6 in its less restrictive sign (4.861e-7 of endpointnorm);opposite sign ceiling is~4.7e-79. Exact limiting species change across finite-difference steps,not stable biological bottlenecks. This is geometry only,no nonlinearfeasible transientgrowth measured. See `results/spatial_ode_direction_feasibility_audit.json`.
