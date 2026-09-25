"""SUPERSEDED (kept as record): first multi-cell attempt; 0 fixed points, broken.
Working model: experiments/segment_polarity_mc2.py (see results/grn_multicell.json).

Study B: multi-cell segment-polarity Boolean model (Albert & Othmer 2003).

1D row of C cells x 5 nodes (wg, en, hh, ptc, ci). Paracrine coupling:
WG and HH signal to immediate neighbours. Exact attractor enumeration over
all 2^(5C) states via vectorized successor + pointer jumping.
Convention: default-off threshold (node ON iff signed input sum > 0).
"""
from __future__ import annotations

import numpy as np

NODES = ["wg", "en", "hh", "ptc", "ci"]


def successor_all(n_cells=4):
    """Vectorized successor map over all 2^(5*n_cells) states -> uint32 array."""
    n_bits = 5 * n_cells
    states = np.arange(2**n_bits, dtype=np.uint32)

    def bit(node, cell):
        b = (states >> np.uint32(cell * 5 + NODES.index(node))) & np.uint32(1)
        return b.astype(np.int16)

    def neigh(node, cell):
        left = bit(node, cell - 1) if cell > 0 else 0
        right = bit(node, cell + 1) if cell < n_cells - 1 else 0
        return bit(node, cell) + left + right

    out = np.zeros(2**n_bits, dtype=np.uint32)
    for c in range(n_cells):
        hh_sig = (neigh("hh", c) > 0).astype(np.int16)        # HH signal received
        wg_n = bit("ci", c)                                   # CI-A -> wg
        en_n = (neigh("wg", c) > 0).astype(np.int16)          # paracrine WG -> en
        hh_n = en_n                                           # en -> hh
        ptc_n = (bit("ci", c) & (hh_sig == 0)).astype(np.int16)  # CI-A -> ptc when no HH
        ci_n = ((bit("en", c) == 0) & ((hh_sig > 0) | (bit("ptc", c) == 0))).astype(np.int16)  # en -| ci; PTC unbound -| ci
        for node, val in [("wg", wg_n), ("en", en_n), ("hh", hh_n),
                          ("ptc", ptc_n), ("ci", ci_n)]:
            out |= (val.astype(np.uint32) << np.uint32(c * 5 + NODES.index(node)))
    return out


def fixed_points(succ):
    """All states mapping to themselves."""
    return np.nonzero(succ == np.arange(len(succ), dtype=np.uint32))[0]


def basin_map(succ):
    """Pointer-jump: attractor (fixed point) reached by every state."""
    n = len(succ)
    root = succ.copy()
    for _ in range(40):
        new = root[root]
        if np.array_equal(new, root):
            break
        root = new
    return root


def decode(state, n_cells=4):
    return [[(state >> (c * 5 + i)) & 1 for i in range(5)] for c in range(n_cells)]


def stripe_summary(state, n_cells=4):
    d = decode(state, n_cells)
    return {n: [row[i] for row in d] for i, n in enumerate(NODES)}
