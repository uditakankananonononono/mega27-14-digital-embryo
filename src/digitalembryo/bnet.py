"""Minimal .bnet (BoolNet format) parser + exact synchronous dynamics."""
from __future__ import annotations

import itertools
import re


def parse_bnet(text):
    """Return {var: (expr_fn, [inputs])} from 'targets,factors' bnet text."""
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    if lines and lines[0].lower().startswith("targets"):
        lines = lines[1:]
    rules = {}
    for ln in lines:
        tgt, expr = ln.split(",", 1)
        tgt, expr = tgt.strip(), expr.strip()
        inputs = sorted(set(re.findall(r"v_[A-Za-z0-9_]+", expr)))
        py = expr.replace("&", " and ").replace("|", " or ").replace("!", " not ").strip()
        # order: replace longest names first to avoid prefix collisions
        code = compile(py, "<bnet>", "eval")
        rules[tgt] = (code, inputs)
    return rules


def successor(state, rules, externals=None):
    env = dict(state)
    if externals:
        env.update(externals)
    nxt = {}
    for tgt, (code, _in) in rules.items():
        if tgt in (externals or {}):
            nxt[tgt] = externals[tgt]
            continue
        nxt[tgt] = int(bool(eval(code, {"__builtins__": {}}, env)))
    return nxt


def attractor_analysis(rules, externals=None):
    """Exact enumeration over all non-external variables."""
    fixed_vars = set((externals or {}).keys())
    free = [v for v in rules if v not in fixed_vars]
    seen, fixed_keys, cycle_keys = set(), set(), set()
    fixed_points, cycles = [], []
    for bits in itertools.product([0, 1], repeat=len(free)):
        s0 = dict(zip(free, bits))
        key0 = tuple(s0[v] for v in free)
        if key0 in seen:
            continue
        path, cur, pos = [], s0, {}
        while True:
            key = tuple(cur[v] for v in free)
            if key in seen:  # joins an earlier trajectory
                for s in path:
                    seen.add(tuple(s[v] for v in free))
                break
            if key in pos:
                cyc = path[pos[key]:]
                if len(cyc) == 1:
                    if key not in fixed_keys:
                        fixed_keys.add(key)
                        fixed_points.append(cyc[0])
                else:
                    ck = min(tuple(s[v] for v in free) for s in cyc)
                    if ck not in cycle_keys:
                        cycle_keys.add(ck)
                        cycles.append(cyc)
                for s in path:
                    seen.add(tuple(s[v] for v in free))
                break
            pos[key] = len(path)
            path.append(cur)
            cur = successor(cur, rules, externals)
    return {"fixed": fixed_points, "cycles": cycles,
            "n_states": 2 ** len(free), "free_vars": free}
