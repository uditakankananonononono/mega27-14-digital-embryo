#!/usr/bin/env python3
"""Verdict #6 generalization test: developmental fragility index on a SECOND
developmental network - the gap-gene cross-repression module (hb-kr-gt-kni
style: anterior morphogen activates hb; hb/kr/gt/kni cross-repress in a chain;
each gene represses its posterior neighbor and the next-but-one).
4 nodes -> 16 states, exact enumeration, same fragility definition as id-021:
fraction of single-bit perturbations of WT-basin states exiting the basin.
WT pattern: single anterior-domain state hb ON, rest OFF (anterior gap domain).
"""
import json
import numpy as np

NODES = ["hb", "kr", "gt", "kni"]
# edges (source, target, sign): classic gap cross-repression chain
# bcd->hb activation is a fixed prepattern input (anterior), modeled as const 1
EDGES = [("hb", "kr", -1), ("kr", "gt", -1), ("gt", "kni", -1),
         ("hb", "gt", -1), ("kr", "kni", -1),
         ("kni", "gt", -1), ("gt", "kr", -1)]  # mutual repression of adjacent
BCD = 1  # anterior prepattern present

def step(state):
    out = []
    for i, t in enumerate(NODES):
        s = 0
        for src, tgt, sg in EDGES:
            if tgt == t:
                v = (state >> NODES.index(src)) & 1
                s += sg * v
        if t == "hb":
            s += BCD  # morphogen activation
        out.append(1 if s > 0 else 0)
    return sum(b << i for i, b in enumerate(out))

def main():
    n = 16
    succ = np.array([step(s) for s in range(n)])
    # attractors (fixed points + cycles) via pointer following
    belong = -np.ones(n, dtype=int)
    for s0 in range(n):
        if belong[s0] != -1: continue
        path, s = [], s0
        while belong[s] == -1 and s not in path:
            path.append(s); s = int(succ[s])
        root = s
        for p in path: belong[p] = root
    # WT state: hb ON only = 0b0001 = 1
    WT = 1
    wt_root = belong[WT]
    wt_basin = np.nonzero(belong == wt_root)[0]
    exits, total = 0, 0
    per_bit_exit = np.zeros(4); per_bit_tot = np.zeros(4)
    for s in wt_basin:
        for b in range(4):
            t = int(s) ^ (1 << b)
            total += 1; per_bit_tot[b] += 1
            if belong[t] != wt_root:
                exits += 1; per_bit_exit[b] += 1
    out = {
        "network": "gap-gene cross-repression module (hb/kr/gt/kni), anterior Bcd prepattern",
        "n_states": n,
        "wt_state": "hb-ON only (anterior gap domain)",
        "wt_basin_states": [int(x) for x in wt_basin],
        "wt_basin_fraction": round(len(wt_basin) / n, 4),
        "fragility_index": round(exits / total, 4) if total else None,
        "per_gene_exit_fraction": {NODES[b]: round(float(per_bit_exit[b] / per_bit_tot[b]), 4) for b in range(4)},
        "comparison": {"id021_single_cell": {"wt_basin_fraction": 0.0078, "fragility": 0.500}},
        "reading": "generalization: does the fragility framework's structure (rare WT basin, gene-specific exit asymmetry) recur in an unrelated developmental module",
    }
    json.dump(out, open("results/grn_fragility_gapgene.json", "w"), indent=1)
    print(json.dumps(out, indent=1))

if __name__ == "__main__":
    main()
