#!/usr/bin/env python3
"""How much does it matter that the composition is not carried with the front?

In the coupled run every Gauss point has its own point model (Klempt et al.
2026). A point that the growing biofilm reaches for the first time starts
from the seed composition, not from that of the biofilm that grew into it.
This check runs the same composition scheme (A: rescaled to the field's
amount every coupling step) on a 3D growth field, two ways:

  V0 (as in the coupled run): a newly colonised point starts from the seed,
     phi_1 / (phi_1 + phi_2) = 0.5, psi_i = 0.999;
  V1 (inherit): it starts from the state of the nearest point that was
     already colonised, i.e. the biofilm that grew into it.

V1 is a first-order stand-in for transport: new biomass at the front is made
by the biofilm behind it. The difference V1 - V0 is the size of the
limitation. The field is Klempt 2024 test case 4.1 as fitted in
JAXFEM/klempt2024_timescale.py (Table 2, edge-line source, s = 1.5), on the
paper's 21^3 grid; the point model runs at every grid node with phi >= phi_min.
Cases 3 and 6 of Klempt et al. 2026 Table 1. Coupling step 0.1, s = 0.15,
phi_cap = 0.9, phi_min = 0.01 (as composition_figs.py; assumptions, not from
a paper).

    python ansys_usermat/composition_transport_check.py
        -> assets/fig_composition_transport.png (maps and per-age means)
           assets/fig_composition_transport_slide.png (per-age means only)

Result (2026-10-03), biomass-weighted mean of phi_1/(phi_1+phi_2) at T* = 1:
  case 3: seed restart 0.656, inherit 0.656, largest pointwise difference 0.007
  case 6: seed restart 0.026, inherit 0.027, largest pointwise difference 0.093
Small, because the diffusive edge of the field reaches phi_min long before
the dense front does: every point starts its point model early, at a small
amount, and has run for most of the time by T* = 1. The difference is largest
at the points colonised last (case 6, T* = 0.4-0.5). The composition depends
more on the local amount (the rescaling) than on where it came from.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import jax
import jax.numpy as jnp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "coupling"))
sys.path.insert(0, str(ROOT / "JAXFEM"))
jax.config.update("jax_enable_x64", True)
import ecology_jax as eco  # noqa: E402
import figstyle  # noqa: E402
import klempt2024_case1_bc as B  # noqa: E402
import klempt2024_quantitative as K  # noqa: E402
import material_server as ms  # noqa: E402

OUT = ROOT / "assets" / "fig_composition_transport.png"
OUT_SLIDE = ROOT / "assets" / "fig_composition_transport_slide.png"
DT, S, PHI_CAP, PHI_MIN, S_MACRO = 0.1, 0.15, 0.9, 1.0e-2, 1.5
CASES = {"2sp_case3": "case 3", "2sp_case6": "case 6"}


def macro_snapshots():
    """phi of the fitted 4.1 field at T* = 0, 0.1, ..., 1."""
    K.R, K.BETA, K.K_A = 100 * S_MACRO, 2 * S_MACRO, 1e-3 * S_MACRO
    phi0, _, _ = K.setup("fig4_edge")
    snaps = [phi0.copy()]
    for k in range(1, 11):
        rec, phi, _ = B.run_setup(snaps[-1], B.nutrient_mask("edge"), 1e8, "ic",
                                  "first_order", "abs", t_end=DT)
        snaps.append(phi)
    return snaps


def seed_state(n):
    g = np.zeros((n, 12))
    g[:, 0] = g[:, 1] = 0.5
    g[:, 6] = g[:, 7] = 0.999
    return g


def rescale(g, p):
    f = p / (g[:, 0] + g[:, 1])
    g = g.copy()
    g[:, 0] *= f
    g[:, 1] *= f
    g[:, 5] = 1.0 - p
    return g


def stepper(theta, hp):
    dt_pm = S * DT
    n_sub = max(1, math.ceil(dt_pm / 1e-4 - 1e-9))
    mask = eco.active_mask(2)
    extra = (jnp.float64(hp["c"]), jnp.float64(hp["alpha"]),
             jnp.asarray(hp["eta"], dtype=jnp.float64))
    th = jnp.asarray(theta, dtype=jnp.float64)
    steps = jnp.arange(n_sub)

    @jax.jit
    def run(G):
        return jax.vmap(lambda g: eco._substep_scan(g, th, dt_pm / n_sub, steps, mask,
                                                    *extra)[0])(G)
    return lambda G: np.asarray(run(jnp.asarray(G)))


def simulate(snaps, case, inherit):
    ms.set_case(case)
    adv = stepper(ms.ECOLOGY_CASE["theta"], ms.ECOLOGY_CASE["hp"])
    ms.set_case(None)
    xyz = np.stack([K.X.ravel(), K.Y.ravel(), K.Z.ravel()], 1)
    n = xyz.shape[0]
    G = np.zeros((n, 12))
    started = np.zeros(n, bool)
    born = np.full(n, np.nan)
    for k, phi in enumerate(snaps[1:], start=1):
        p = np.minimum(phi.ravel(), PHI_CAP)
        new = (~started) & (p >= PHI_MIN)
        if new.any():
            if inherit and started.any():
                _, j = cKDTree(xyz[started]).query(xyz[new])
                G[new] = G[np.flatnonzero(started)[j]]
            else:
                G[new] = seed_state(new.sum())
            started |= new
            born[new] = (k - 1) * DT
        act = started & (p >= PHI_MIN)
        G[act] = rescale(G[act], p[act])
        G[act] = rescale(adv(G[act]), p[act])
    chi = np.where(started, G[:, 0] / np.maximum(G[:, 0] + G[:, 1], 1e-300), np.nan)
    return chi.reshape(phi.shape), born.reshape(phi.shape)


def main():
    figstyle.apply()
    snaps = macro_snapshots()
    phi = snaps[-1]
    w = np.where(phi >= PHI_MIN, phi, 0.0)
    fig, ax = plt.subplots(2, 3, figsize=(13, 7.4))
    fs, axs = plt.subplots(1, 2, figsize=(9, 3.6))
    print("biomass-weighted mean of phi_1/(phi_1+phi_2) at T* = 1")
    for r, (case, name) in enumerate(CASES.items()):
        chi0, born = simulate(snaps, case, False)
        chi1, _ = simulate(snaps, case, True)
        m0 = np.nansum(w * chi0) / w.sum()
        m1 = np.nansum(w * chi1) / w.sum()
        d = np.abs(chi1 - chi0)
        print(f"  {name}: seed restart {m0:.3f}  inherit {m1:.3f}  "
              f"max |difference| {np.nanmax(d):.3f}")
        k = K.N // 2
        ext = (0, K.L, 0, K.L)
        for c, (field, title) in enumerate(((chi0, "seed restart (as run)"),
                                            (chi1, "inherit from neighbour"))):
            im = ax[r, c].imshow(field[:, :, k].T, origin="lower", extent=ext,
                                 vmin=0, vmax=1, cmap="PuOr_r")
            ax[r, c].set_title(f"{name}: {title}")
        fig.colorbar(im, ax=ax[r, :2], label=r"$\phi_1/(\phi_1+\phi_2)$", shrink=0.85)
        tb = born.ravel()
        for t0 in np.unique(tb[~np.isnan(tb)]):
            sel = (tb == t0) & (w.ravel() > 0)
            for a, ms_ in ((ax[r, 2], 5), (axs[r], 7)):
                a.plot(t0, np.nanmean(chi0.ravel()[sel]), "o", color="#1d4e89", ms=ms_,
                       label="starts from the seed (as run)" if t0 == 0 else None)
                a.plot(t0, np.nanmean(chi1.ravel()[sel]), "s", color="#b45309", ms=ms_ * 0.8,
                       label="starts from the neighbour" if t0 == 0 else None)
        axs[r].set_title(f"{name}: mean {m0:.3f} / {m1:.3f}")
        axs[r].set_xlabel(r"time the point is first reached, $T^*$")
        axs[r].set_ylim(0, 1)
        ax[r, 2].set_xlabel(r"time the point was colonised, $T^*$")
        ax[r, 2].set_ylabel(r"mean $\phi_1/(\phi_1+\phi_2)$ at $T^*=1$")
        ax[r, 2].set_ylim(0, 1)
        ax[r, 2].legend(loc="best", fontsize=9)
        for a in ax[r, :2]:
            a.set_xlabel(r"$x$ [$\mu$m]"); a.grid(False)
        ax[r, 0].set_ylabel(r"$y$ [$\mu$m]")
        ax[r, 1].set_yticklabels([])
    fig.suptitle(r"Composition at $T^*=1$, mid-plane $z = 10\ \mu$m; growth field of "
                 r"Klempt et al. 2024 test case 4.1, nutrient at the top right edge", fontsize=12)
    fig.savefig(OUT, dpi=200)
    axs[0].set_ylabel(r"$\phi_1/(\phi_1+\phi_2)$ at $T^*=1$")
    axs[1].legend(loc="upper left", fontsize=9)
    fs.savefig(OUT_SLIDE, dpi=200)
    print("wrote", OUT, OUT_SLIDE)


if __name__ == "__main__":
    main()
