"""Reference and checks for the composition mode (prop(28) = 7).

The 3D field (Klempt 2024, the partner's USSFin, unchanged) decides the
amount phi_3D; the point model (Klempt et al. 2026, a Table 1 case) decides
only the composition chi_i = phi_i / sum(phi) and psi. At each Gauss point
and substep: rescale the carried state to sum(phi_i) = phi_3D (phi0 =
1 - phi_3D), integrate dt with n_sub inner steps, rescale again.

reference(...) runs that scheme stand-alone for a prescribed phi_3D(t), the
exact target for a gradient-free region (check (b) of SPLIT_COUPLING.md);
check_comp_trace(...) judges a trace written by the fragment.
"""
from __future__ import annotations

import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent / "coupling"))

PHIMAX = 1.0 - 1.0e-6


def rescale(g, phi3):
    """Same operations, same order as the Fortran."""
    g = np.array(g, dtype=np.float64)
    s = g[0] + g[1]
    f = phi3 / s
    g[0] = g[0] * f
    g[1] = g[1] * f
    g[5] = 1.0 - phi3
    return g


def seed(chi1=0.5):
    g = np.zeros(12)
    g[0], g[1], g[6], g[7] = chi1, 1.0 - chi1, 0.999, 0.999
    return g


def step(g, phi3, theta, hp, dt, n_sub):
    """One substep: rescale in, integrate, rescale out."""
    import ecology_jax as eco
    g = rescale(g, phi3)
    gn, _ = eco.ecology_substeps(g, theta, dt, n_sub, 2, hp)
    return rescale(np.asarray(gn), phi3)


def reference(phi3_series, theta, hp, dt, dt_max=1.0e-4, chi1=0.5,
              phi_min=0.0, s=1.0, phi_cap=PHIMAX):
    """Stand-alone scheme for phi_3D(t) given per substep; returns states.
    s: the point model's clock, dt_pm = s * dt (prop(31)).
    phi_cap: the point model sees min(phi_3D, phi_cap) (prop(32))."""
    g = seed(chi1)
    dt_pm = s * dt
    n_sub = 1 if dt_pm <= dt_max else math.ceil(dt_pm / dt_max)
    out = []
    for p in phi3_series:
        p = min(max(p, 0.0), min(phi_cap, PHIMAX))
        if p >= phi_min and p > 0.0:
            g = step(g, p, theta, hp, dt_pm, n_sub)
        out.append(g.copy())
    return out


def read_comp_trace(path):
    rows = []
    with open(path, newline="") as f:
        for r in csv.reader(f):
            v = [x.strip() for x in r]
            if not v or not v[0].isdigit():
                continue
            g = [float(x) for x in v[12:36]]
            rows.append({"elem": int(v[0]), "ip": int(v[1]),
                         "ldstep": int(v[2]), "isubst": int(v[3]),
                         "nsub": int(v[4]), "hit": int(v[5]),
                         "dtime": float(v[6]), "dt_pm": float(v[7]),
                         "phi_used": float(v[8]),
                         "phi3": float(v[9]), "alpha_n": float(v[10]),
                         "alpha_new": float(v[11]),
                         "g_old": g[:12], "g_new": g[12:]})
    return rows


def check_comp_trace(rows, k_alpha, theta=None, hp=None, dt_max=1.0e-4):
    """once_per_increment, carried, amount (sum phi_i = phi_3D after the
    rescale), one_call (one server call per point and substep), eq36
    (alpha_new = alpha_n + k_alpha * phi_3D * dt) and, with theta/hp,
    replay (each server call recomputed: bit for bit)."""
    by_pt = defaultdict(lambda: defaultdict(list))
    for r in rows:
        by_pt[(r["elem"], r["ip"])][(r["ldstep"], r["isubst"])].append(r)
    once = carried = amount = one_call = eq36 = True
    for subs in by_pt.values():
        keys = sorted(subs)
        for k in keys:
            rs = subs[k]
            if len({tuple(r["g_old"]) for r in rs}) != 1 or \
               len({r["alpha_n"] for r in rs}) != 1:
                once = False
            if rs[0]["hit"] != 2 and [r["hit"] for r in rs].count(0) != 1:
                one_call = False
            for r in rs:
                if r["hit"] != 2 and abs(r["g_new"][0] + r["g_new"][1]
                                         - r["phi3"]) > 4e-16:
                    amount = False
                if r["alpha_new"] != r["alpha_n"] + k_alpha * r["phi_used"] \
                        * r["dtime"]:
                    eq36 = False
        for a, b in zip(keys[:-1], keys[1:]):
            la, fb = subs[a][-1], subs[b][0]
            if la["alpha_new"] != fb["alpha_n"]:
                carried = False
            want = la["g_new"] if fb["hit"] == 2 else \
                list(rescale(la["g_new"], fb["phi3"]))
            if want != fb["g_old"]:
                carried = False
    out = {"once_per_increment": once, "carried": carried,
           "amount": amount, "one_call": one_call, "eq36": eq36}
    if theta is not None:
        import ecology_jax as eco
        worst = 0.0
        for subs in by_pt.values():
            for rs in subs.values():
                r = rs[-1]
                if r["hit"] == 2:
                    continue
                gn, _ = eco.ecology_substeps(r["g_old"], theta, r["dt_pm"],
                                             r["nsub"], 2, hp)
                gn = rescale(np.asarray(gn), r["phi3"])
                worst = max(worst, float(np.max(np.abs(
                    gn - np.asarray(r["g_new"])))))
        out.update({"replay": worst == 0.0, "replay_worst": worst})
    return out


# --------------------------------------------------------------------------
# prop(28) = 8: the point model runs on its own (no scaling to the field)
# from the substep in which the field first reaches phi_min. Klempt et al.
# 2026 call their model a material point model and phi0 a numerical
# auxiliary, not a physical void; this mode keeps it that way.

def seed_age(phi1=0.2, phi2=0.2):
    g = np.zeros(12)
    g[0], g[1] = phi1, phi2
    g[5] = 1.0 - (phi1 + phi2)
    g[6] = g[7] = 0.999
    return g


def reference_age(phi3_series, theta, hp, dt, dt_max=1.0e-4, phi_min=0.0,
                  s=1.0, phi_init=(0.2, 0.2)):
    """Stand-alone scheme of prop(28) = 8 for phi_3D(t) given per substep.
    Returns the point-model state after each substep (zeros before the
    start)."""
    import ecology_jax as eco
    g = np.zeros(12)
    dt_pm = s * dt
    n_sub = 1 if dt_pm <= dt_max else math.ceil(dt_pm / dt_max)
    out = []
    for p in phi3_series:
        p = min(max(p, 0.0), PHIMAX)
        if g[6] == 0.0 and p < phi_min:
            out.append(g.copy())
            continue
        if g[6] == 0.0:
            g = seed_age(*phi_init)
        gn, _ = eco.ecology_substeps(g, theta, dt_pm, n_sub, 2, hp)
        g = np.asarray(gn)
        out.append(g.copy())
    return out


def check_age_trace(rows, k_alpha, theta=None, hp=None):
    """once_per_increment, carried (unscaled: g_new of one substep is g_old
    of the next, bit for bit), held (before the start the state is zero),
    one_call, eq36, and with theta/hp a bit-exact replay of every call."""
    by_pt = defaultdict(lambda: defaultdict(list))
    for r in rows:
        by_pt[(r["elem"], r["ip"])][(r["ldstep"], r["isubst"])].append(r)
    once = carried = held = one_call = eq36 = True
    for subs in by_pt.values():
        keys = sorted(subs)
        for k in keys:
            rs = subs[k]
            if len({tuple(r["g_old"]) for r in rs}) != 1 or \
               len({r["alpha_n"] for r in rs}) != 1:
                once = False
            if rs[0]["hit"] != 2 and [r["hit"] for r in rs].count(0) != 1:
                one_call = False
            for r in rs:
                if r["hit"] == 2 and any(x != 0.0 for x in r["g_new"]):
                    held = False
                if r["alpha_new"] != r["alpha_n"] + k_alpha * r["phi_used"] \
                        * r["dtime"]:
                    eq36 = False
        for a, b in zip(keys[:-1], keys[1:]):
            la, fb = subs[a][-1], subs[b][0]
            if la["alpha_new"] != fb["alpha_n"]:
                carried = False
            # the substep that starts the point model follows a held one
            # (state zero) and begins from the seed instead
            if la["hit"] != 2 and la["g_new"] != fb["g_old"]:
                carried = False
    out = {"once_per_increment": once, "carried": carried, "held": held,
           "one_call": one_call, "eq36": eq36}
    if theta is not None:
        import ecology_jax as eco
        worst = 0.0
        for subs in by_pt.values():
            for rs in subs.values():
                r = rs[-1]
                if r["hit"] == 2:
                    continue
                gn, _ = eco.ecology_substeps(r["g_old"], theta, r["dt_pm"],
                                             r["nsub"], 2, hp)
                worst = max(worst, float(np.max(np.abs(
                    np.asarray(gn) - np.asarray(r["g_new"])))))
        out.update({"replay": worst == 0.0, "replay_worst": worst})
    return out

