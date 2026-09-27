#!/usr/bin/env python3
"""Update-rule sensitivity (verdict weakness #10): does the wild-type basin
result survive asynchronous and stochastic updates?

Synchronous: exact enumeration (128/2048, results/grn_fragility.json).
Asynchronous: random-order single-node updates, full enumeration of initial
  states with replicated random schedules.
Stochastic: synchronous updates with bit-flip noise probability q per step.
"""
import itertools
import json
import pathlib

import numpy as np

from digitalembryo.bnet import parse_bnet, successor

ROOT = pathlib.Path(__file__).resolve().parent.parent
BNET = (ROOT / "data/raw/grn_models/id021_bodysegmentation.bnet").read_text()
EXT = {"v_hh_external": 1, "v_WG_external": 0, "v_SLP": 0}


def async_trajectory(s, rules, free, rng, max_steps=300):
    cur = dict(s)
    for _ in range(max_steps):
        order = rng.permutation(free)
        moved = False
        for tgt in order:
            env = dict(cur)
            env.update(EXT)
            val = int(bool(eval(rules[tgt][0], {"__builtins__": {}}, env)))
            if val != cur[tgt]:
                cur[tgt] = val
                moved = True
        if not moved:
            break
    return cur


def noisy_trajectory(s, rules, free, q, rng, max_steps=300):
    cur = dict(s)
    for _ in range(max_steps):
        nxt = successor(cur, rules, EXT)
        for v in free:
            if rng.random() < q:
                nxt[v] ^= 1
        if all(nxt[v] == cur[v] for v in free):
            cur = nxt
            break
        cur = nxt
    return cur


def main():
    rules = parse_bnet(BNET)
    free = [v for v in rules if v not in EXT]
    states = [dict(zip(free, bits)) for bits in itertools.product([0, 1], repeat=len(free))]

    # synchronous exact
    sync_wt = 0
    for s in states:
        cur, seen = s, set()
        for _ in range(200):
            key = tuple(cur[v] for v in free)
            if key in seen:
                break
            seen.add(key)
            nxt = successor(cur, rules, EXT)
            if all(nxt[v] == cur[v] for v in free):
                cur = nxt
                break
            cur = nxt
        if cur.get("v_wg", 0) == 1:
            sync_wt += 1

    rng = np.random.default_rng(0)
    async_counts = []
    for rep in range(8):
        c = 0
        rng_r = np.random.default_rng(rep)
        for s in states:
            if async_trajectory(s, rules, free, rng_r).get("v_wg", 0) == 1:
                c += 1
        async_counts.append(c)

    noise = {}
    for q in (0.001, 0.01, 0.05, 0.1):
        counts = []
        for rep in range(4):
            c = 0
            rng_r = np.random.default_rng(100 + rep)
            for s in states[::4]:  # quarter sample for speed
                if noisy_trajectory(s, rules, free, q, rng_r).get("v_wg", 0) == 1:
                    c += 1
            counts.append(c * 4)
        noise[str(q)] = {"mean_states": float(np.mean(counts)), "sd": float(np.std(counts))}

    out = {
        "model": "Cell Collective id-021, condition hh1_WG0_SLP0",
        "synchronous_exact_wt_states": sync_wt,
        "synchronous_exact_fraction": sync_wt / len(states),
        "async_wt_states_over_8_reps": async_counts,
        "async_fraction_mean": float(np.mean(async_counts) / len(states)),
        "async_fraction_sd": float(np.std(async_counts) / len(states)),
        "stochastic_noise_sweep": noise,
        "reading": "basin conclusion is update-rule dependent if async fraction diverges materially from synchronous 0.781%",
    }
    (ROOT / "results" / "grn_update_rules.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
