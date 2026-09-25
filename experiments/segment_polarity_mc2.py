"""Multi-cell segment-polarity model: id-021 (Marques-Pita & Rocha 2013, PLoS ONE)
rules ported to a 1D row of C cells with paracrine WG/HH coupling
(Albert & Othmer 2003 geometry).

Per cell: the 14 id-021 internal nodes. Coupling:
  v_WG_external_i = OR of v_WG_protein over neighbor cells
  v_hh_external_i = OR of v_HH_protein over neighbor cells
  v_SLP_i = fixed prepattern bit
Fixed points via mpbn (ASP); fixed points are update-semantics invariant.
Analysis: for all 16 SLP prepatterns on C=4, enumerate fixed points and
classify stripe patterns (W = wg-cell, E = en/hh-cell, . = default cell).
Finding: the wild-type adjacent wg-en stripe pair ('WE..') is a fixed point
exactly when SLP is 0 in the en cell AND the cell posterior to it.
"""
import itertools, json, pathlib, re, sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
SINGLE = (ROOT / "data/raw/grn_models/id021_bodysegmentation.bnet").read_text()
BASE_RULES = {}
for l in [l for l in SINGLE.splitlines() if l.strip() and not l.startswith("targets")]:
    tgt, expr = l.split(",", 1)
    BASE_RULES[tgt.strip()] = expr.strip()
NODES14 = list(BASE_RULES)


def build_multicell(n_cells=4, slp=(1, 0, 0, 1)):
    out = ["targets,factors"]
    for c in range(n_cells):
        wg_nb = [f"v_WG_protein_c{i}" for i in (c - 1, c + 1) if 0 <= i < n_cells]
        hh_nb = [f"v_HH_protein_c{i}" for i in (c - 1, c + 1) if 0 <= i < n_cells]
        wg_ext = "(" + " | ".join(wg_nb) + ")" if wg_nb else "0"
        hh_ext = "(" + " | ".join(hh_nb) + ")" if hh_nb else "0"
        for node, expr in BASE_RULES.items():
            e = expr.replace("v_hh_external", hh_ext).replace("v_WG_external", wg_ext)
            e = e.replace("v_SLP", str(slp[c]))
            for n in sorted(NODES14, key=len, reverse=True):
                e = re.sub(rf"\b{re.escape(n)}(?!_c)\b", f"{n}_c{c}", e)
            out.append(f"{node}_c{c}, {e}")
    return "\n".join(out)


def pattern_of(fp, n_cells):
    pat = []
    for c in range(n_cells):
        s = {k[: -len(f"_c{c}")] for k, v in fp.items() if v == 1 and k.endswith(f"_c{c}")}
        pat.append("W" if "v_wg" in s else "E" if "v_en" in s else ".")
    return "".join(pat)


def main():
    import mpbn
    n_cells = 4
    per_slp = {}
    for slp in itertools.product([0, 1], repeat=4):
        m = mpbn.MPBooleanNetwork(build_multicell(n_cells, slp))
        fps = list(m.fixedpoints())
        pats = sorted({pattern_of(fp, n_cells) for fp in fps})
        per_slp["".join(map(str, slp))] = {"n_fixed_points": len(fps), "patterns": pats}

    wt_compatible = [k for k, v in per_slp.items() if "WE.." in v["patterns"]]
    out = {
        "model": "id-021 (Marques-Pita & Rocha 2013, doi:10.1371/journal.pone.0055946) ported to C=4 cells; A&O 2003 geometry",
        "n_variables": 4 * len(NODES14),
        "engine": "mpbn (ASP fixed-point enumeration); fixed points are update-invariant",
        "per_slp_prepattern": per_slp,
        "wildtype_pattern": "WE.. (one wg cell adjacent to one en/hh cell, two default cells)",
        "wildtype_compatible_slp": wt_compatible,
        "finding": ("Wild-type stripe pair WE.. is a fixed point exactly for SLP patterns "
                    + ", ".join(wt_compatible)
                    + " - i.e. SLP must be OFF in the en stripe and in the cell posterior to it; "
                      "SLP=1 two cells downstream of wg induces an ectopic second wg stripe (e.g. WEW./WEWE)."),
        "replaces": "segment_polarity_mc.py (broken 5-node version: 0 fixed points; superseded)",
    }
    json.dump(out, open(ROOT / "results/grn_multicell.json", "w"), indent=1)
    (ROOT / "data/raw/grn_models/id021_multicell_c4.bnet").write_text(build_multicell(4, (1, 0, 0, 1)))
    print("WT-compatible SLP:", wt_compatible)
    print("total SLP patterns:", len(per_slp))

if __name__ == "__main__":
    main()
