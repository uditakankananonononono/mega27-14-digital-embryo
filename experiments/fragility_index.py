#!/usr/bin/env python3
"""Developmental fragility index on the published id-021 segment-polarity network.

Verdict reframe 5 + weakness #11/#16: instead of reporting the small wild-type
basin as a bare negative, quantify HOW fragile it is:
  1. exact WT-basin enumeration (synchronous, condition hh1_WG0_SLP0),
  2. single-bit-flip fragility: fraction of one-bit perturbations of basin
     states whose trajectory exits the wg-ON basin (overall and per node),
  3. gene-knockout prediction: clamp each free node to 0 / 1, re-enumerate,
     report basin fold-change (intervention predictions),
  4. edge-removal sensitivity: delete each regulatory edge, re-enumerate.
"""
import itertools
import re
import json
import pathlib

from digitalembryo.bnet import parse_bnet, successor

ROOT = pathlib.Path(__file__).resolve().parent.parent
BNET = (ROOT / "data/raw/grn_models/id021_bodysegmentation.bnet").read_text()
EXT = {"v_hh_external": 1, "v_WG_external": 0, "v_SLP": 0}


def attractor_of(state, rules, ext, free, max_steps=200):
    cur, seen = state, set()
    for _ in range(max_steps):
        key = tuple(cur[v] for v in free)
        if key in seen:
            return cur, "cycle"
        seen.add(key)
        nxt = successor(cur, rules, ext)
        if all(nxt[v] == cur[v] for v in free):
            return nxt, "fixed"
        cur = nxt
    return cur, "unknown"


def wt_basin(rules, ext, free, clamp=None):
    """Set of initial states (as tuples over free-nonclamped vars) reaching a wg-ON fixed point."""
    clamp = clamp or {}
    live = [v for v in free if v not in clamp]
    basin = set()
    for bits in itertools.product([0, 1], repeat=len(live)):
        s = dict(zip(live, bits))
        s.update(clamp)
        s.update(ext)
        att, kind = attractor_of(s, rules, ext, free)
        if kind == "fixed" and att.get("v_wg", 0) == 1:
            basin.add(bits)
    return basin, len(live)


def main():
    rules = parse_bnet(BNET)
    free = [v for v in rules if v not in EXT]
    basin, n_live = wt_basin(rules, EXT, free)
    total = 2 ** n_live
    frac = len(basin) / total

    # single-bit-flip fragility
    flips_exit = {v: 0 for v in free}
    flips_tot = 0
    exit_tot = 0
    idx = {v: i for i, v in enumerate(free)}
    for bits in basin:
        for v in free:
            b2 = list(bits)
            b2[idx[v]] ^= 1
            s = dict(zip(free, b2))
            s.update(EXT)
            att, kind = attractor_of(s, rules, EXT, free)
            inside = kind == "fixed" and att.get("v_wg", 0) == 1
            flips_tot += 1
            if not inside:
                exit_tot += 1
                flips_exit[v] += 1
    fragility = exit_tot / flips_tot
    per_node = {v: flips_exit[v] / len(basin) for v in free}

    # knockout / constitutive-on predictions
    knock = {}
    base = frac
    for v in free:
        for val in (0, 1):
            b, nl = wt_basin(rules, EXT, free, clamp={v: val})
            knock[f"{v}={'OFF' if val == 0 else 'ON'}"] = {
                "basin_fraction": len(b) / (2 ** nl),
                "fold_vs_wt": (len(b) / (2 ** nl)) / base if base else None}

    # edge-removal sensitivity (patch raw bnet text, re-parse)
    raw_lines = {}
    for l in BNET.splitlines():
        if l.strip() and not l.startswith("targets"):
            tgt, expr = l.split(",", 1)
            raw_lines[tgt.strip()] = expr.strip()
    edge_sens = {}
    for tgt, expr in raw_lines.items():
        for src_v in list(raw_lines) + list(EXT):
            if re.search(rf"\b{re.escape(src_v)}\b", expr):
                mod_lines = dict(raw_lines)
                mod_lines[tgt] = re.sub(rf"\b{re.escape(src_v)}\b", "0", expr)
                text = "targets,factors\n" + "\n".join(f"{t}, {e}" for t, e in mod_lines.items())
                try:
                    r2 = parse_bnet(text)
                    b, nl = wt_basin(r2, EXT, free)
                    edge_sens[f"{src_v}->{tgt}"] = {"basin_fraction": len(b) / (2 ** nl),
                                          "fold_vs_wt": (len(b) / (2 ** nl)) / base if base else None}
                except Exception as e:
                    edge_sens[f"{src_v}->{tgt}"] = {"error": str(e)[:80]}

    out = {
        "model": "Cell Collective id-021, condition hh1_WG0_SLP0",
        "n_free_nodes": len(free), "state_space": total,
        "wt_basin_states": len(basin), "wt_basin_fraction": frac,
        "fragility_index": fragility,
        "fragility_definition": "fraction of single-bit flips of WT-basin states whose trajectory exits the wg-ON basin",
        "per_node_exit_fraction": per_node,
        "knockout_predictions": knock,
        "edge_removal_sensitivity": edge_sens,
    }
    (ROOT / "results" / "grn_fragility.json").write_text(json.dumps(out, indent=1))
    print(json.dumps({k: out[k] for k in ("wt_basin_states", "wt_basin_fraction", "fragility_index")}, indent=1))
    print("worst nodes:", sorted(per_node.items(), key=lambda kv: -kv[1])[:4])
    kills = {k: v for k, v in knock.items() if v.get("fold_vs_wt") == 0}
    print("knockouts abolishing WT basin:", list(kills))


if __name__ == "__main__":
    main()
