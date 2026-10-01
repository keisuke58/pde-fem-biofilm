#!/usr/bin/env python3
"""Reference for the one-species coupling, independent of the Fortran.

Klempt et al. (2024) Eq. 36, ``alpha_dot = k_alpha * phi``, integrated the way
``growth_from_phi.f`` integrates it (explicit in phi). Used by
``tests/test_one_species_coupling.py`` and, once the partner-side deck runs on
IKMHIWI03, to predict ``SVAR(84)`` from the phi history that deck prints.

Kept deliberately small: if this file grows a cap, a clamp or a gate, it has
stopped being Klempt 2024.
"""
from __future__ import annotations

import numpy as np


def alpha_history(phi, k_alpha: float, dt: float) -> np.ndarray:
    """alpha after each step, from alpha(0) = 0, for a sequence of phi."""
    alpha = 0.0
    out = []
    for p in np.asarray(phi, dtype=float):
        alpha = alpha + k_alpha * p * dt
        out.append(alpha)
    return np.array(out)


def alpha_klempt(alpha_repo: float) -> float:
    """This repository's alpha in Klempt's convention, Fg = alpha_K * I."""
    return 1.0 + alpha_repo


# ---------------------------------------------------------------------------
# reading the diagnostic trace the call site writes
# (apdl/callsite/phi_mode_exec.inc). Comma-separated, one row per call:
#   elem, ip, ldstep, isubst, dtime, bio1, locbio1, bio2, locbio2,
#   phi_used, alpha_n, alpha_new
# ---------------------------------------------------------------------------

import csv
from collections import defaultdict

COLS = ("elem", "ip", "ldstep", "isubst", "dtime", "bio1", "locbio1",
        "bio2", "locbio2", "phi_used", "alpha_n", "alpha_new")


def read_trace(path):
    rows = []
    with open(path, newline="") as f:
        for r in csv.reader(f):
            if not r or r[0].strip().lower() == "elem":
                continue
            v = [x.strip() for x in r]
            rows.append({"elem": int(v[0]), "ip": int(v[1]),
                         "ldstep": int(v[2]), "isubst": int(v[3]),
                         **{k: float(x) for k, x in zip(COLS[4:], v[4:])}})
    return rows


def check_trace(rows, k_alpha, phi_col="phi_used", rtol=1e-9):
    """Return a dict of named pass/fail checks plus the evidence for each.

    once_per_increment -- alpha_n is identical on every call within one
        sub-step, i.e. ANSYS hands back the converged value each iteration and
        growth is not counted once per Newton iteration.
    eq36 -- on every call, alpha_new = alpha_n + k_alpha * phi * dtime.
    carried -- the next sub-step starts from the alpha the previous one
        converged to (last call of each sub-step).
    phi_* -- what each candidate phi column does, to decide which is the
        Gauss-point value: whether it differs between points, changes in
        time, and stays in [0, 1].
    """
    by_pt = defaultdict(list)
    for r in rows:
        by_pt[(r["elem"], r["ip"])].append(r)
    once, eq36, carried = True, True, True
    worst = 0.0
    for pts in by_pt.values():
        subs = defaultdict(list)
        for r in pts:
            subs[(r["ldstep"], r["isubst"])].append(r)
            want = r["alpha_n"] + k_alpha * r[phi_col] * r["dtime"]
            err = abs(r["alpha_new"] - want)
            worst = max(worst, err)
            if err > rtol * max(1.0, abs(want)):
                eq36 = False
        keys = sorted(subs)
        for k in keys:
            if len({r["alpha_n"] for r in subs[k]}) != 1:
                once = False
        for a, b in zip(keys[:-1], keys[1:]):
            last = subs[a][-1]["alpha_new"]
            if abs(subs[b][0]["alpha_n"] - last) > rtol * max(1.0, abs(last)):
                carried = False

    def phi_stats(col):
        vals = [r[col] for r in rows]
        per_t = defaultdict(set)
        per_pt = defaultdict(set)
        for r in rows:
            per_t[(r["ldstep"], r["isubst"])].add(round(r[col], 12))
            per_pt[(r["elem"], r["ip"])].add(round(r[col], 12))
        return {"varies_between_points": any(len(s) > 1 for s in per_t.values()),
                "changes_in_time": any(len(s) > 1 for s in per_pt.values()),
                "in_0_1": min(vals) >= -1e-12 and max(vals) <= 1 + 1e-12,
                "min": min(vals), "max": max(vals)}

    return {"once_per_increment": once, "eq36": eq36, "carried": carried,
            "eq36_worst_abs_error": worst,
            **{f"phi_{c}": phi_stats(c)
               for c in ("bio1", "locbio1", "bio2", "locbio2")}}


def k_alpha_for_target(alpha_target, phi_max, time_total):
    """k_alpha that keeps total growth near alpha_target over the deck.

    The stand-in while the deck's physical time unit is unknown: choose growth
    by its result (a small, safe alpha) rather than by a rate in physical
    units, and report it that way.
    """
    return alpha_target / (phi_max * time_total)


def k_alpha_from_trace(rows, alpha_target, phi_col="bio1"):
    """k_alpha that brings the fastest-growing traced point to alpha_target.

    Uses the phi the run actually produced (a k_alpha = 0 stage-1 trace),
    integrated over converged sub-steps -- the last call of each sub-step --
    rather than assuming phi_max = 1. In the 1 Oct stage 2 run phi peaked
    near 0.01, so the phi_max = 1 rule landed alpha at 1/200 of its target.
    """
    by_pt = defaultdict(lambda: defaultdict(list))
    for r in rows:
        by_pt[(r["elem"], r["ip"])][(r["ldstep"], r["isubst"])].append(r)
    best = 0.0
    for subs in by_pt.values():
        integ = sum(v[-1][phi_col] * v[-1]["dtime"] for v in subs.values())
        best = max(best, integ)
    if best <= 0.0:
        raise ValueError(f"{phi_col} never rises above 0 in this trace")
    return alpha_target / best
