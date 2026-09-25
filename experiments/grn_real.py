#!/usr/bin/env python3
"""Study B with the real published network: Cell Collective id-021 (body
segmentation in Drosophila, 2013; Albert-Othmer lineage), attractors under
all 8 external-input conditions (hh_ext, WG_ext, SLP)."""
import json
import pathlib

from digitalembryo.bnet import attractor_analysis, parse_bnet

ROOT = pathlib.Path(__file__).resolve().parent.parent
BNET = (ROOT / "data/raw/grn_models/id021_bodysegmentation.bnet").read_text()


def main():
    rules = parse_bnet(BNET)
    out = {"model": "Cell Collective id-021 BODY-SEGMENTATION-IN-DROSOPHILA-2013",
           "doi": "10.1371/journal.pone.0055946", "conditions": {}}
    for hh in (0, 1):
        for wg in (0, 1):
            for slp in (0, 1):
                ext = {"v_hh_external": hh, "v_WG_external": wg, "v_SLP": slp}
                a = attractor_analysis(rules, ext)
                key = f"hh{hh}_WG{wg}_SLP{slp}"
                out["conditions"][key] = {
                    "n_fixed": len(a["fixed"]),
                    "n_cycles": len(a["cycles"]),
                    "fixed": a["fixed"],
                    "cycle_lengths": [len(c) for c in a["cycles"]]}
                print(key, "fixed:", len(a["fixed"]), "cycles:", [len(c) for c in a["cycles"]])
                for fp in a["fixed"]:
                    print("   ", {k: v for k, v in sorted(fp.items())})
    (ROOT / "results" / "grn_real.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
