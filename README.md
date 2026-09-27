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
- `embryosim info-threshold` - minimum-information CF decode (reproduces 2.39% EL vs Liu's 4.94% benchmark)
Every command prints JSON with the dataset/provenance fields named inline.
