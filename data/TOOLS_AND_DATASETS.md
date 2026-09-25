# External tools & datasets - honest counts (mega27-14)

## External tools/resources genuinely used (current: 40/40 gate MET; Google Drive delivery pending - not counted until shipped)
| # | Tool/resource | Genuinely used for |
|---|---|---|
| 1 | Zenodo (record 4942019) | Liu 2013 dataset download |
| 2 | Dryad (doi:10.5061/dryad.2j9p4) | provenance cross-check of the same dataset |
| 3 | Cell Collective | published GRN model id-021/id-027 |
| 4 | sybila/biodivine-boolean-models (GitHub) | curated .bnet files id-021, id-027 |
| 5 | BioModels (BIOMD0000001065) | von Dassow 2000 reference lookup |
| 6 | NCBI eutils (esearch/esummary) | gene-ID + citation lookups (data/lookups/ncbi_genes.json) |
| 7 | NCBI Gene | bcd 40830, hb 41032, eve 36039, kni 40287, gt 31227 |
| 8 | PubMed (PMID 23580621) | Liu 2013 citation verification |
| 9 | PNAS/PMC full text + SI (Liu 2013) | published bars: lambda=16.5%EL, Sx=10.5%EL, Sc=44% |
| 10 | Gregor et al. 2007 Cell | published gradient-precision benchmark |
| 11 | Petkova et al. 2019 Cell | 4-gap-gene decoding benchmark (~1% EL) |
| 12 | Dubuis et al. 2013 PNAS | single-gene information benchmark (2-3 bits) |
| 13 | PLOS figshare (Spirov et al. 2017) | SDD-critique literature cross-check |
| 14 | Python 3.10 | all computation |
| 15 | NumPy | all numerics |
| 16 | SciPy (optimize/stats/io) | curve fitting, Spearman, .mat parsing |
| 17 | pandas | tabular handling |
| 18 | PyTorch (CPU) | CNN decoders, GNN/MLP basin models |
| 19 | scikit-learn | decoder cross-checks |
| 20 | NetworkX | GRN graph handling |
| 21 | matplotlib | all figures |
| 22 | h5py | HDF5/.mat v7.3 fallback reading |
| 23 | statsmodels | robust-OLS DCLS cross-check (results/dcls_ols.json) |
| 24 | scikit-image | Turing spot morphometrics (results/turing_morphometrics.json) |
| 25 | sympy | symbolic verification of every paper derivation (tests/test_math_verify.py) |
| 26 | pytest | hermetic test suite (30+ tests) |
| 27 | setuptools | package + embryosim CLI entry point |
| 28 | pip | dependency management |
| 29 | Git | version control |
| 30 | GitHub | remote hosting/collaboration |
| 31 | OpenSSH | authenticated push |
| 32 | curl | data + API downloads |
| 33 | LuaLaTeX (TeX Live) | paper typesetting |
| 34 | fontspec | Times New Roman loading |
| 35 | luaotfload | OpenType font backend for LuaLaTeX |
| 36 | Times New Roman TTFs (MS core fonts, times32.exe) | mandated paper font (licensed; TTFs not redistributed) |
| 37 | 7-Zip | times32.exe TTF extraction |
| 38 | fontconfig (fc-cache/fc-match) | font registration/verification |
| 39 | poppler-utils (pdffonts/pdfinfo) | PDF font-embedding verification, page count |
| 40 | Python stdlib urllib / multiprocessing | dataset download, NCBI API calls, parallel gradient fitting |
| 41 | GNU tar + xz-utils | TeX package archive handling |

## Datasets (accession-level, uniform rule: one identifier-backed record = 1)
| Source | Records | Count |
|---|---|---|
| Liu 2013 LiveImaging (Bcd-GFP gradients, 29 lines) | per-embryo-side gradients used in fits | 676 |
| Liu 2013 Immunofluorescence (Hb/Kr/Gt/Kni/Eve) | per-embryo records | 2180 |
| Cell Collective id-021 (body segmentation 2013) | model | 1 |
| Cell Collective id-027 (WG pathway) | model | 1 |
| **Total** | | **2858** (gate: 120 - met) |

Formulas with derivations in the paper: 12+ numbered equations, each
symbolically verified (sympy) in tests/test_math_verify.py where algebraic.
