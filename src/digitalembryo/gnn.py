"""GNN over a regulatory graph: predict basin fate from initial state."""
from __future__ import annotations

import numpy as np


def graph_arrays(rules, externals):
    """Adjacency (targets x sources) + node list for the free variables."""
    free = [v for v in rules if v not in externals]
    idx = {v: i for i, v in enumerate(free)}
    A = np.zeros((len(free), len(free)))
    for tgt, (_c, inputs) in rules.items():
        if tgt not in idx:
            continue
        for src in inputs:
            if src in idx:
                A[idx[tgt], idx[src]] = 1.0
    return free, A


def basin_labels(rules, externals, n_samples=2000, seed=0, max_steps=100):
    """Random initial states -> attractor label (1 if wg-protein ON at fixpoint)."""
    from .bnet import successor

    free = [v for v in rules if v not in externals]
    rng = np.random.default_rng(seed)
    X, y = [], []
    for _ in range(n_samples):
        s = {v: int(rng.integers(0, 2)) for v in free}
        cur, seen = s, set()
        for _ in range(max_steps):
            key = tuple(cur[v] for v in free)
            if key in seen:
                break
            seen.add(key)
            nxt = successor(cur, rules, externals)
            if all(nxt[v] == cur[v] for v in free):
                break
            cur = nxt
        X.append([cur0 for cur0 in (s[v] for v in free)])
        y.append(int(cur.get("v_wg", 0)))
    return np.array(X, dtype=np.float32), np.array(y, dtype=np.float32), free


class GNNClassifier:
    def __init__(self, n_nodes, hidden=32, seed=0):
        import torch

        torch.manual_seed(seed)
        self.torch = torch
        self.net = torch.nn.Sequential(
            torch.nn.Linear(2 * hidden, hidden), torch.nn.ReLU(),
            torch.nn.Linear(hidden, 1))

        class _Emb(torch.nn.Module):
            def __init__(self):
                super().__init__()
                self.emb = torch.nn.Parameter(torch.randn(n_nodes, hidden) * 0.1)

        self.emb = _Emb()
        self.params = list(self.net.parameters()) + list(self.emb.parameters())

    def _forward(self, A_norm, x0):
        t = self.torch
        h = self.emb.emb + x0.unsqueeze(-1)  # (batch, nodes, hidden) broadcast-add init
        m = t.einsum("ij,bjh->bih", A_norm, h)   # aggregate from sources
        z = t.cat([h, m], dim=-1)
        out = self.net(z).mean(dim=1).squeeze(-1)  # graph readout
        return out

    def fit(self, A, X, y, epochs=200, lr=5e-3):
        t = self.torch
        A_n = t.tensor(A / (A.sum(1, keepdims=True) + 1e-6), dtype=t.float32)
        Xt, yt = t.tensor(X), t.tensor(y)
        opt = t.optim.Adam(self.params, lr=lr)
        lossf = t.nn.BCEWithLogitsLoss()
        for _ in range(epochs):
            opt.zero_grad()
            loss = lossf(self._forward(A_n, Xt), yt)
            loss.backward()
            opt.step()
        return float(loss.item())

    def predict(self, A, X):
        t = self.torch
        A_n = t.tensor(A / (A.sum(1, keepdims=True) + 1e-6), dtype=t.float32)
        with t.no_grad():
            return (t.sigmoid(self._forward(A_n, t.tensor(X))) > 0.5).float().numpy()


class MLPBaseline:
    def __init__(self, n_in, hidden=64, seed=0):
        import torch

        torch.manual_seed(seed)
        self.torch = torch
        self.net = torch.nn.Sequential(
            torch.nn.Linear(n_in, hidden), torch.nn.ReLU(),
            torch.nn.Linear(hidden, 1))

    def fit(self, X, y, epochs=200, lr=5e-3):
        t = self.torch
        opt = t.optim.Adam(self.net.parameters(), lr=lr)
        lossf = t.nn.BCEWithLogitsLoss()
        Xt, yt = t.tensor(X), t.tensor(y)
        for _ in range(epochs):
            opt.zero_grad()
            loss = lossf(self.net(Xt).squeeze(-1), yt)
            loss.backward()
            opt.step()

    def predict(self, X):
        t = self.torch
        with t.no_grad():
            return (t.sigmoid(self.net(t.tensor(X)).squeeze(-1)) > 0.5).float().numpy()
