"""Study B core: the Drosophila segment-polarity gene network.

Topology taken from Albert & Othmer, J. Theor. Biol. 223:1-18 (2003),
the canonical published Boolean model of the segment-polarity GRN
(single-cell, 4 genes/proteins: wg/WG, en/EN, hh/HH, ptc/PTC, ci/CI, PH,
plus SMO). Simplified single-cell Boolean core used here: 5 nodes with
published interaction signs; synchronous update; attractor enumeration by
exhaustive state-space search (2^5 = 32 states, exact, no sampling).
"""
from __future__ import annotations

import itertools

import numpy as np

# Albert-Othmer-inspired single-cell Boolean core.
# Edges: (source, target, sign) sign=+1 activation, -1 repression.
# Sources: Table 1 / Fig. 1 of Albert & Othmer 2003 (wg->en, en->hh,
# hh->ptc(repress), ptc->|ci, ci->wg(+), ptc self via HH binding).
NODES = ["wg", "en", "hh", "ptc", "ci"]
EDGES = [
    ("ci", "wg", +1),    # CI activates wg
    ("wg", "en", +1),    # WG (paracrine) maintains en
    ("en", "hh", +1),    # EN activates hh
    ("hh", "ptc", -1),   # HH binds PTC, removing repression (ptc activity down)
    ("ptc", "ci", -1),   # PTC represses ci (via SMO sequestration)
    ("ci", "ptc", +1),   # CI (activator form) induces ptc transcription
]


def successor(state, edges=EDGES, nodes=NODES):
    """Synchronous Boolean successor: node ON iff signed input sum > 0,
    OFF iff < 0; ties keep current state (hysteresis convention - models
    protein half-life persistence in the reduced single-cell core)."""
    idx = {n: i for i, n in enumerate(nodes)}
    nxt = list(state)
    for t_i, t in enumerate(nodes):
        s = 0
        for src, tgt, sign in edges:
            if tgt == t:
                s += sign * state[idx[src]]
        if s > 0:
            nxt[t_i] = 1
        elif s < 0:
            nxt[t_i] = 0
    return tuple(nxt)


def attractors(nodes=NODES, edges=EDGES):
    """Exact attractor enumeration over the full state space."""
    n = len(nodes)
    attr = {"fixed": [], "cycles": []}
    seen_global = set()
    for bits in itertools.product([0, 1], repeat=n):
        if bits in seen_global:
            continue
        path, cur, seen = [], bits, {}
        while cur not in seen:
            seen[cur] = len(path)
            path.append(cur)
            cur = successor(cur, nodes=nodes, edges=edges)
        seen_global.update(path)
        if seen[cur] == len(path) - 1:
            attr["fixed"].append(cur)
        else:
            attr["cycles"].append(path[seen[cur]:])
    return attr


def basin_sizes(nodes=NODES, edges=EDGES):
    """Basin size (number of states flowing to each attractor)."""
    n = len(nodes)
    attr_states = {a for a in attractors(nodes, edges)["fixed"]}
    basins = {a: 0 for a in attr_states}
    for bits in itertools.product([0, 1], repeat=n):
        cur, guard = bits, 0
        while cur not in attr_states and guard < 2**n + 2:
            cur = successor(cur, nodes=nodes, edges=edges)
            guard += 1
        if cur in basins:
            basins[cur] += 1
    return basins


def robustness_knockout(node, nodes=NODES, edges=EDGES):
    """Attractor set after fixing `node` OFF (loss-of-function)."""
    idx = {n: i for i, n in enumerate(nodes)}

    def succ_lof(state):
        s = successor(state, nodes=nodes, edges=edges)
        s = list(s)
        s[idx[node]] = 0
        return tuple(s)

    fixed = set()
    for bits in itertools.product([0, 1], repeat=len(nodes)):
        if bits[idx[node]] != 0:
            continue
        cur, guard = bits, 0
        seen = set()
        while cur not in seen and guard < 40:
            seen.add(cur)
            cur = succ_lof(cur)
            guard += 1
        fixed.add(cur)
    return sorted(fixed)
