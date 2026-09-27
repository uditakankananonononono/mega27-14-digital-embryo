"""Cross-organism generalization test (verdict #15): Mendoza & Alvarez-Buylla
threshold Boolean network of Arabidopsis floral organ identity.

Model parameters transcribed from the published weight matrix and threshold
vector (original 12-node network; nodes EMF1, TFL1, LFY, AP1, CAL, LUG, UFO,
BFU, AG, AP3, PI, SUP), update rule x_i(t+1) = H(sum_j w_ij x_j - theta_i),
H(z)=1 iff z>0. Source: Mendoza & Alvarez-Buylla (1998) as reprinted in full
in Ruz & Goles (2018), Fig. 1 (pageperso.lis-lab.fr/~sylvain.sene/files/publi_pres/rgs18.pdf).

Validation targets from the same source: synchronous update must yield exactly
13 attractors (6 fixed points + 7 limit cycles of length 2), and the six fixed
points must match the published organ-identity states.
"""
import itertools, json
import numpy as np

NODES = ["EMF1","TFL1","LFY","AP1","CAL","LUG","UFO","BFU","AG","AP3","PI","SUP"]
W = np.array([
 [ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],   # EMF1
 [ 1, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0],   # TFL1
 [-2,-1, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0],   # LFY
 [-1, 0, 5, 0, 0, 0, 0, 0,-1, 0, 0, 0],   # AP1
 [ 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0],   # CAL
 [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],   # LUG
 [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],   # UFO
 [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0],   # BFU
 [ 0,-2, 1,-2, 0,-1, 0, 0, 0, 0, 0, 0],   # AG
 [ 0, 0, 3, 0, 0, 0, 2, 1, 0, 0, 0,-2],   # AP3
 [ 0, 0, 4, 0, 0, 0, 1, 1, 0, 0, 0,-1],   # PI
 [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],   # SUP
], dtype=int)
THETA = np.array([0, 0, 3, -1, 1, 0, 0, 1, -1, 0, 0, 0])

PUBLISHED_FP = {
 "sepal":         "000100000000",
 "petal":         "000100010110",
 "carpel":        "000000001000",
 "stamen":        "000000011110",
 "inflorescence": "110000000000",
 "mutant":        "110000010110",
}

def step_batch(X):
    return (X @ W.T - THETA > 0).astype(np.uint8)

def main():
    n = len(NODES)
    states = np.array(list(itertools.product([0,1], repeat=n)), dtype=np.uint8)
    powers = (1 << np.arange(n-1, -1, -1)).astype(np.int64)
    succ = step_batch(states) @ powers  # succ[i] = integer index of successor state
    def attractor(s):
        seen = {}; cur = s
        while cur not in seen:
            seen[cur] = len(seen); cur = succ[cur]
        cyc = [c for c, k in seen.items() if k >= seen[cur]]
        return tuple(sorted(cyc))
    att_of = {}
    attractors = {}
    for s in range(2**n):
        a = attractor(s)
        att_of[s] = a
        attractors.setdefault(a, []).append(s)
    n_fixed = sum(1 for a in attractors if len(a) == 1)
    n_cyc2 = sum(1 for a in attractors if len(a) == 2)
    # validate published fixed points
    fp_found = {NODES[0]: None}
    fp_match = {}
    for name, bits in PUBLISHED_FP.items():
        s = int(bits, 2)
        ok = int(succ[s]) == s
        fp_match[name] = bool(ok)
    # basin sizes and per-fixed-point fragility
    frag = {}
    basin_sizes = {}
    for a, members in attractors.items():
        basin_sizes[a] = len(members)
        if len(a) == 1:
            exits = 0; total = 0
            for s in members:
                for b in range(n):
                    t = s ^ (1 << (n-1-b))
                    total += 1
                    if att_of[t] != a:
                        exits += 1
            frag[a] = exits / total
    # label fixed points
    labeled = {}
    for name, bits in PUBLISHED_FP.items():
        s = int(bits, 2)
        a = attractor(s)
        labeled[name] = {"state": bits, "is_fixed": fp_match[name],
                         "basin": len(attractors.get(a, [])),
                         "fragility": round(frag.get(a, float("nan")), 4)}
    out = {
     "n_states": 2**n, "n_attractors": len(attractors),
     "n_fixed_points": n_fixed, "n_limit_cycles_len2": n_cyc2,
     "published_fixed_points_match": fp_match,
     "labeled": labeled,
     "comparison_segment_polarity": {"wt_basin_frac": 128/16384, "fragility": 0.500},
     "gap_gene_cascade": {"fragility": 0.0, "structure": "single global attractor"},
    }
    json.dump(out, open("results/flower_model.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != "labeled"}, indent=1))
    print(json.dumps(labeled, indent=1))

if __name__ == "__main__":
    main()
