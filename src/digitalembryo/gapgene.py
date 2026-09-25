"""Gap-gene immunofluorescence data (Liu et al. PNAS 2013, Zenodo 4942019).

1000-point AP-axis protein profiles (dorsal+ventral) for Hb, Kr, Gt, Kni and
Eve across fly lines with different Bcd dosages. Each embryo record is one
accession-level dataset.
"""
from __future__ import annotations

import pathlib

import numpy as np
import scipy.io

DATA = pathlib.Path(__file__).resolve().parents[2] / "data" / "raw" / "liu2013"
AXIS = np.linspace(0.0, 1.0, 1000)  # the 1000 points span the AP axis


def load_if_sessions(path=None):
    """Return list of session dicts: line, genes, embryos[{dorsal,ventral,..}]."""
    path = path or DATA / "Immunofluorescence.mat"
    m = scipy.io.loadmat(str(path), squeeze_me=True, struct_as_record=False)
    sessions = []
    for s in np.atleast_1d(m["Sessions"]):
        genes = [str(g) for g in np.atleast_1d(getattr(s, "GeneName"))]
        embryos = []
        for e in np.atleast_1d(getattr(s, "Embryos")):
            dorsal = np.asarray(getattr(e, "Dorsal"), dtype=float)
            ventral = np.asarray(getattr(e, "Ventral"), dtype=float)
            if dorsal.ndim == 1:
                dorsal = dorsal[:, None]
            if ventral.ndim == 1:
                ventral = ventral[:, None]
            embryos.append({
                "dorsal": dorsal, "ventral": ventral,
                "egg_length_um": float(getattr(e, "EggLength", np.nan)),
                "invagination_depth_um": float(getattr(e, "InvaginationDepth", np.nan)),
            })
        sessions.append({"line": str(getattr(s, "FlyLineName")), "genes": genes,
                         "embryos": embryos})
    return sessions


def boundary_position(profile, axis=AXIS, lo=0.2, hi=0.8):
    """Boundary = AP position where the (normalised) profile crosses its
    mid-range value, restricted to [lo, hi] EL. NaN if no clean crossing."""
    p = np.asarray(profile, dtype=float)
    ok = np.isfinite(p) & (axis >= lo) & (axis <= hi)
    if ok.sum() < 20:
        return np.nan
    x, y = axis[ok], p[ok]
    span = np.nanmax(y) - np.nanmin(y)
    if span <= 0:
        return np.nan
    mid = np.nanmin(y) + 0.5 * span
    s = np.sign(y - mid)
    crossings = np.nonzero(np.diff(s) != 0)[0]
    if len(crossings) != 1:
        return np.nan  # multi-peak profile: no single boundary
    i = crossings[0]
    x0, x1, y0, y1 = x[i], x[i + 1], y[i], y[i + 1]
    return float(x0 + (mid - y0) * (x1 - x0) / (y1 - y0))


def hb_boundaries(session):
    """Hb boundary per embryo (dorsal+ventral mean) for an Hb-containing session."""
    gi = session["genes"].index("Hb")
    out = []
    for e in session["embryos"]:
        if e["dorsal"].shape[1] <= gi:
            continue
        prof = np.nanmean(np.stack([e["dorsal"][:, gi], e["ventral"][:, gi]]), axis=0)
        out.append(boundary_position(prof))
    return np.array(out)
