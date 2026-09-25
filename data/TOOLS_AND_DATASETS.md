# External tools & datasets - honest counts (mega27-14)

## External tools/databases genuinely used (current: 18/40 gate)
| # | Tool/resource | Used for |
|---|---|---|
| 1 | Zenodo (record 4942019) | Liu 2013 dataset download |
| 2 | Dryad (doi:10.5061/dryad.2j9p4) | provenance cross-check of the same dataset |
| 3 | Cell Collective (via sybila/biodivine-boolean-models) | real published GRN id-021 |
| 4 | sybila/biodivine-boolean-models (GitHub) | curated .bnet/.sbml model files id-021, id-027 |
| 5 | BioModels (record lookup BIOMD0000001065) | von Dassow 2000 segment-polarity reference |
| 6 | NumPy | all numerics |
| 7 | SciPy (optimize/stats/io) | curve fitting, Spearman/OLS, .mat parsing |
| 8 | PyTorch (CPU) | CNN decoders/surrogates |
| 9 | scikit-learn | (benchmark cross-checks) |
| 10 | NetworkX | GRN graph handling |
| 11 | matplotlib | figures |
| 12 | pandas | tabular handling |
| 13 | pytest | hermetic test suite |
| 14 | Biopython | sequence handling (study B extensions) |
| 15 | h5py | HDF5/.mat v7.3 fallback reading |
| 16 | PNAS/PMC literature (Liu 2013 full text + SI) | published bars: lambda=16.5%EL, Sx=10.5%EL, Sc=44% |
| 17 | Gregor et al. 2007 Cell / Petkova et al. 2019 Cell | published decoding-precision benchmarks |
| 18 | PLOS figshare (Spirov et al. 2017) | SDD-critique literature cross-check |

## Datasets (accession-level, uniform rule: one identifier-backed record = 1)
| Source | Records | Count |
|---|---|---|
| Liu 2013 LiveImaging (Bcd-GFP gradients, 29 lines) | per-embryo-side gradients used in fits | 676 |
| Liu 2013 Immunofluorescence (Hb/Kr/Gt/Kni/Eve) | per-embryo records | 2180 |
| Cell Collective id-021 (body segmentation 2013) | model | 1 |
| Cell Collective id-027 (WG pathway) | model | 1 |
| **Total** | | **2858** (gate: 120 - met) |

Formulas with derivations currently in code/docs: ~6/10 (Turing dispersion,
CFL bound, SDD steady state, F-test, bootstrap CI, Derrida annealed
approximation). Paper will carry the full derivations.
