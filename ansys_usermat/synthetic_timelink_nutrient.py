#!/usr/bin/env python3
"""synthetic_timelink_nutrient.py -- the time link s and the local nutrient c*(x), tested with synthetic data only.

Point model: the four-species synthetic case of Fritsch et al. 2025 (Case II, submodel M3, true parameters of their
Fig. 16, fritsch2025_cases.CASES["fritsch_4sp_M3"]): phi_i(0) = 0.02 each, psi_i(0) = 0.999, c* = 25, alpha* = 0,
eta_i = 1, 750 steps of 1e-4. Their synthetic data cover normalised times tau = t / 0.075 in [0.05, 1]. The "truth" is
this free point model; no measured data enter.

Coupling: the one-way scheme of Chapter 4 (coupling (a)): every coupling step dT* = 0.025 the point model is rescaled
to the field amount p = min(phi_3D, phi_cap) (phi_cap = 0.9, phi_0 = 1 - p), advanced by s dT* in substeps of at most
1e-4, rescaled again. phi_3D(T*) is the field at a seed Gauss point of the Chapter 4 seed problem (16^3, the traced
point of results/2026-10-cal5/w16_CH_g1_s015_final_d1_pm.json, column phi_used; the field does not depend on the point
model, so it is the same for any point model). Start: the synthetic initial shares (1/4 each).

1. Time link (improvement 6). With synthetic data the time window is known, so s follows from it: T* = 1 spans the
   data window when s = 750 x 1e-4 = 0.075. The coupled composition at T* = tau is compared with the truth at tau over
   the data times (RMS of the shares), for s = 0.075 and a scan of s; scan minimum = the s the data would select.
2. Local nutrient (improvement 4). c* = 25 c(x), c from d lap c = g phi c on the partner's 8^3 seed problem with the
   nutrient held on one face (composition_local_nutrient.nutrient, Thiele number Lambda^2 = g L^2 / d), s = 0.075.

    python ansys_usermat/synthetic_timelink_nutrient.py   -> assets/fig_synthetic_timelink_nutrient.pdf/.png,
                                                            ansys_usermat/results_synthetic_timelink_nutrient.json
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import jax
import jax.numpy as jnp
import numpy as np

jax.config.update("jax_enable_x64", True)
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "JAXFEM"))
import composition_local_nutrient as L  # noqa: E402
import figstyle  # noqa: E402
from fritsch2025_cases import CASES  # noqa: E402

eco = L.T.eco
CASE = CASES["fritsch_4sp_M3"]
N_SP = 4
DT_PM, STEPS = CASE["dt"], CASE["steps"]
WINDOW = DT_PM * STEPS                         # 0.075, the synthetic data window in point-model time
DT_C, PHI_CAP = 0.025, 0.9
TAU_DATA = np.linspace(0.05, 1.0, 20)          # Fritsch Case I/II: 20 data points on [0.05, 1]
OUT = ROOT / "assets" / "fig_synthetic_timelink_nutrient"
JSON = HERE / "results_synthetic_timelink_nutrient.json"


def theta20():
    A, b = CASE["A"], CASE["b"]
    th = [0.0] * 20
    th[0], th[1], th[2], th[3], th[4] = A[0][0], A[0][1], A[1][1], b[0], b[1]
    th[5], th[6], th[7], th[8], th[9] = A[2][2], A[2][3], A[3][3], b[2], b[3]
    th[10], th[11], th[12], th[13] = A[0][2], A[0][3], A[1][2], A[1][3]
    th[14], th[15] = 1.0, 1.0                  # species 5 masked off
    return jnp.asarray(th, dtype=jnp.float64)


TH = theta20()
MASK = eco.active_mask(N_SP)


def g0(n=1):
    G = np.zeros((n, 12))
    G[:, :N_SP] = CASE["phi_init"]
    G[:, 5] = 1.0 - sum(CASE["phi_init"])
    G[:, 6:11] = CASE.get("psi0", 0.999)
    G[:, 6 + N_SP:11] = 0.999
    return G


def scan_fn(n_sub, dt_sub):
    steps = jnp.arange(n_sub)
    return jax.jit(jax.vmap(lambda g, c: eco._substep_scan(g, TH, dt_sub, steps, MASK, c, jnp.float64(0.0),
                                                            jnp.ones(5))))


def truth():
    """Free point model: states after every 1e-4 step (STEPS + 1 rows)."""
    G, traj = scan_fn(STEPS, DT_PM)(jnp.asarray(g0()), jnp.asarray([CASE["c_star"]]))
    tr = np.concatenate([g0(), np.asarray(traj)[0]], 0)
    return tr


def rescale(G, p):
    G = G.copy()
    G[:, :5] *= (p / G[:, :5].sum(1))[:, None]
    G[:, 5] = 1.0 - p
    return G


def phi3d_series():
    pm = json.load(open(HERE / "apdl/results/2026-10-cal5/w16_CH_g1_s015_final_d1_pm.json"))
    C = {c: i for i, c in enumerate(pm["columns"])}
    rows = [r for r in pm["seed_trace"] if r[C["elem"]] == pm["seed_trace"][0][C["elem"]] and r[C["ip"]] == "1"]
    by = {}
    for r in rows:
        by.setdefault(int(r[C["isubst"]]), float(r[C["phi_used"]]))
    phi = np.array([by[k] for k in sorted(by)])
    assert phi.size == round(1 / DT_C), phi.size
    return phi                                  # phi_3D at the start of each coupling step


PHI3D = phi3d_series()


def coupled(s, cvals, scheme="a"):
    """Shares after every coupling step, for each nutrient level (rows: T* = 0, 0.025, ..., 1).
    scheme "a": rescaled to the field (Chapter 4); "b": the point model keeps its own amount."""
    cv = np.asarray(cvals, float)
    dt = s * DT_C
    n_sub = max(1, math.ceil(dt / 1e-4 - 1e-9))
    run = scan_fn(n_sub, dt / n_sub)
    G = np.repeat(g0(), cv.size, 0)
    cst = jnp.asarray(CASE["c_star"] * cv)
    out = [G[:, :N_SP] / G[:, :N_SP].sum(1, keepdims=True)]
    for k in range(PHI3D.size):
        p = np.full(cv.size, min(PHI3D[k], PHI_CAP))
        if scheme == "a":
            G = rescale(G, p)
        G = np.asarray(run(jnp.asarray(G), cst)[0])
        if scheme == "a":
            G = rescale(G, p)
        out.append(G[:, :N_SP] / G[:, :N_SP].sum(1, keepdims=True))
    return np.array(out)                        # (41, n, 4)


def at_tau(series_T, tau):
    """Linear interpolation of a (41, 4) series on T* = 0..1 at the times tau."""
    T = np.linspace(0, 1, series_T.shape[0])
    return np.stack([np.interp(tau, T, series_T[:, i]) for i in range(series_T.shape[1])], 1)


def main():
    tr = truth()
    sh_tr = tr[:, :N_SP] / tr[:, :N_SP].sum(1, keepdims=True)
    tau_tr = np.arange(STEPS + 1) / STEPS
    tr_data = np.stack([np.interp(TAU_DATA, tau_tr, sh_tr[:, i]) for i in range(N_SP)], 1)
    res = {"case": "Fritsch et al. 2025, Case II M3, true parameters", "window_point_model_time": WINDOW,
           "truth_end_shares": sh_tr[-1].tolist(), "truth_end_amount": float(tr[-1, :N_SP].sum()),
           "phi3d_start_end": [float(PHI3D[0]), float(PHI3D[-1])]}

    # 1. time link
    def rms(s, scheme="a"):
        c = coupled(s, [1.0], scheme)[:, 0, :]
        return float(np.sqrt(np.mean((at_tau(c, TAU_DATA) - tr_data) ** 2))), c
    s_grid = [0.02, 0.03, 0.04, 0.05, 0.06, 0.075, 0.09, 0.11, 0.15, 0.2, 0.3]
    scan = {s: rms(s)[0] for s in s_grid}
    r_win, c_win = rms(WINDOW)
    r_b, c_b = rms(WINDOW, "b")
    r_015, c_015 = rms(0.15)
    s_best = min(scan, key=scan.get)
    res["time_link"] = dict(s_window=WINDOW, rms_window=r_win, rms_s015=r_015, rms_scheme_b_window=r_b,
                            scan={str(k): v for k, v in scan.items()}, s_best=s_best, rms_best=scan[s_best],
                            end_shares_window=c_win[-1].tolist(), end_shares_s015=c_015[-1].tolist(),
                            end_shares_scheme_b=c_b[-1].tolist())
    print("time link: RMS of the shares over the 20 data times")
    for k, v in scan.items():
        print(f"   s = {k:<6} RMS {v:.4f}")
    print(f"   s = window {WINDOW}: {r_win:.4f}; scheme (b) at the window: {r_b:.2e}; s = 0.15: {r_015:.4f}")
    print("   end shares truth", np.round(sh_tr[-1], 3), "coupled s=window", np.round(c_win[-1], 3),
          "s=0.15", np.round(c_015[-1], 3))

    # 2. local nutrient
    L.N, L.H = 8, 0.25
    S, D = L.grid_masks()
    res["nutrient"] = {}
    nut = {}
    for lam in (1.0, 3.0, 4.0):
        cs = L.nutrient(S, D, lam)[S]
        sh = coupled(WINDOW, np.concatenate([[1.0], cs]))[-1]
        nut[lam] = (cs, sh[1:], sh[0])
        lo, hi = sh[1:].min(0), sh[1:].max(0)
        res["nutrient"][str(lam)] = dict(c_min=float(cs.min()), c_max=float(cs.max()), uniform=sh[0].tolist(),
                                         share_min=lo.tolist(), share_max=hi.tolist())
        print(f"Lambda = {lam}: c {cs.min():.2f}-{cs.max():.2f}; shares " +
              " ".join(f"{a:.3f}-{b:.3f}" for a, b in zip(lo, hi)) + " | c = 1: " + " ".join(f"{x:.3f}" for x in sh[0]))
    JSON.write_text(json.dumps(res, indent=1) + "\n")

    # figure
    figstyle.apply(size=9)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(6.3, 2.3))
    T = np.linspace(0, 1, c_win.shape[0])
    for i in range(N_SP):
        l, = ax[0].plot(tau_tr, sh_tr[:, i], lw=1.2, label=f"species {i + 1}")
        ax[0].plot(T, c_win[:, i], ls="--", lw=1.0, color=l.get_color())
        ax[0].plot(TAU_DATA, tr_data[:, i], "o", ms=2, color=l.get_color())
    ax[0].set(xlabel=r"$T^*$ = $\tau$", ylabel=r"$\phi_i/\Sigma_j\phi_j$", xlim=(0, 1), ylim=(0.05, 0.42))
    ax[0].legend(fontsize=6.5, frameon=False, loc="lower left", ncol=2, handlelength=1.2, columnspacing=0.8)
    ax[1].plot(list(scan), list(scan.values()), "o-", ms=3, color="0.2", lw=1)
    ax[1].axvline(WINDOW, color="C3", lw=0.8, ls=":")
    ax[1].set(xscale="log", xlabel="time link $s$", ylabel="RMS of the shares")
    ax[1].set_xticks([0.02, 0.05, 0.1, 0.2]); ax[1].set_xticklabels(["0.02", "0.05", "0.1", "0.2"])
    ax[1].minorticks_off()
    cs, sh, ref = nut[4.0]
    o = np.argsort(cs)
    for i in range(N_SP):
        ax[2].plot(cs[o], sh[o, i], ".", ms=2.5, color=f"C{i}")
        ax[2].axhline(ref[i], color=f"C{i}", lw=0.6, ls=":")
    ax[2].set(xlabel="local nutrient $c$", ylabel=r"$\phi_i/\Sigma_j\phi_j$ at $T^*=1$", ylim=(0.05, 0.42))
    for a, k, x in zip(ax, "abc", (0.97, 0.05, 0.97)):
        a.text(x, 0.97, f"({k})", transform=a.transAxes, va="top", ha="right" if x > 0.5 else "left")
    fig.tight_layout(w_pad=0.8)
    fig.savefig(OUT.with_suffix(".pdf")); fig.savefig(OUT.with_suffix(".png"), dpi=300)
    print("wrote", OUT.with_suffix(".pdf"), JSON)


if __name__ == "__main__":
    main()
