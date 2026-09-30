#!/usr/bin/env python3
"""Does spatial structure bring the clinical conditions back?

In 0D the four conditions collapse. The condition reaches growth only through
alpha = k_alpha * integral(phi_tot), phi_tot = sum_i phi_i psi_i is a single
total, CLSM compositions are fractions summing to 1, and psi is seeded
species-independent -- so the differences cancel in that sum and the four come
out within 0.4% in stress. psi_spread_sensitivity.py measures one way out of
that (give the species different viabilities). This measures the other, and the
more interesting one, because it is the thesis's actual subject.

The supervisor's task is stated as: the ecology model is *pointwise* and does
not account for spatial variation; integrate it so that spatial effects can be
represented. The 0D result is therefore a statement about exactly the setting
the thesis exists to leave behind, and it should not be carried into the
spatial one without checking.

There is a concrete mechanism for it to break. core_hamilton_2d_nutrient.py
gives each species its own diffusivity (D_eff) and its own nutrient
consumption (g_consumption). Both are species-dependent, so in a nutrient
gradient the species need not stay mixed: they can separate in space. Once
they do, phi_tot is no longer a single number standing for the whole
composition, and two conditions with the same total can differ.

So: seed each condition's real Day-1 CLSM composition *uniformly*, let the
coupled system run, and see whether the spread in alpha across CS/CH/DS/DH
exceeds what 0D gives. Either answer is a result. If it grows, the conclusion
changes and the reason is spatial effects -- the assignment's own words. If it
does not, the conclusion survives the mechanism most likely to overturn it.

Not a reproduction of Klempt 2024 section 4.2: that paper holds the nutrient on
one face of the domain, while nutrient_step here holds it on all walls, so the
gradient runs edge-to-centre rather than bottom-to-top. The point being tested
is whether a gradient separates the species at all, which either geometry poses.

STATUS 2026-09-30: the question is not answered yet, and the reason is worth
recording because it is not a bug in this file.

At the ecology's own pace nothing spatial happens at all. Over the horizon the
0D work uses (t ~ 1e-2, where the condition spread peaks) the nutrient barely
moves -- c stays at 0.9996 of its boundary value and the species separation is
exactly 0 -- so control and full give identical answers to every digit. That is
not a null result about biology; it is a timescale mismatch. The Hamilton step
is dt_h = 1e-5..1e-4 while the nutrient field needs t ~ 3-10 to develop a
gradient against D_c = 0.01 and g ~ 0.3-1.0, a separation of roughly a
thousandfold. It is the same shape as the gap between the ecology ODE and a
mechanical deck's increment that n_sub exists to bridge, and it will need the
same kind of deliberate answer rather than a longer loop.

Reaching t ~ 3-10 means thousands of macro steps, and this container did not
get there: raising dt_h to its measured ceiling of 1e-4, or n_react_sub to
100, makes JAX's compilation fail with "Cannot allocate memory" while 15 GB
sit free, so it is an LLVM JIT limit rather than the machine. Run it somewhere
with room before reading anything into the ratio; the script prints a refusal
instead of a verdict when the species never separated.

What it did already establish is the control itself. A first version compared
the 2D spread against the 0D figure quoted from psi_spread_sensitivity.py and
appeared to show spatial structure widening the condition gap 2.4-fold. It
does not: run the same grid with the coupling switched off and the number is
identical to four decimals. The apparent effect was the two scripts using
different horizons, and nothing else. Hence spatial=False here rather than a
number carried across from another file.

Verified working: --nx 8 --macro 5 (defaults otherwise). Larger runs are
what the paragraph above is about.

    python JAXFEM/condition_spread_2d.py --nx 8 --macro 5
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
_REPO = _HERE.parent
for p in (_REPO, _HERE, _REPO / "ansys_usermat" / "coupling",
          _REPO / "ansys_usermat" / "apdl"):
    sys.path.insert(0, str(p))

import jax.numpy as jnp                                        # noqa: E402
import core_hamilton_2d_nutrient as H                          # noqa: E402
from jax_hamilton_0d_5species_demo import THETA_DEMO           # noqa: E402
from plot_heine_phi_psi import load                            # noqa: E402

XLSX = _REPO / "data" / "heine_species_distribution_biofilm.xlsx"
PSI0 = 0.999
K_ALPHA = 50.0
CONDITIONS = [("CS", "Commensal", "Static all cells"),
              ("CH", "Commensal", "HOBIC all cells"),
              ("DS", "Dysbiotic", "Static all cells"),
              ("DH", "Dysbiotic", "HOBIC all cells")]


def day1_compositions():
    import openpyxl
    wb = openpyxl.load_workbook(XLSX, data_only=True)
    out = {}
    for name, sheet, key in CONDITIONS:
        rows = load(wb[sheet])[key][1]
        phi = np.array([np.nanmean(np.array(r, dtype=float)) for r in rows]) / 100.0
        out[name] = phi * (0.999999 / phi.sum())
    return out


def uniform_state(cfg, phi):
    """Every node starts at the same measured composition.

    Uniform on purpose: the only difference between the four runs is then the
    composition itself, so any spatial structure that appears was produced by
    the dynamics rather than put in by hand.
    """
    G = jnp.zeros((cfg.Nx, cfg.Ny, 12), dtype=jnp.float64)
    G = G.at[:, :, :5].set(jnp.asarray(phi))
    G = G.at[:, :, 5].set(1.0 - float(phi.sum()))
    G = G.at[:, :, 6:11].set(PSI0)
    c = jnp.ones((cfg.Nx, cfg.Ny), dtype=jnp.float64) * cfg.c_boundary
    return G.reshape(cfg.Nx * cfg.Ny, 12), c


def run(theta, cfg, phi0, spatial=True):
    """Mirrors core_hamilton_2d_nutrient.run_simulation's macro loop, but from
    a CLSM-seeded uniform state and accumulating alpha as it goes.

    spatial=False is the control, and the only comparison worth making. It
    zeroes the species diffusivities and the nutrient consumption, so the
    nutrient stays uniform at its boundary value and nothing couples one node
    to another: every node is then an independent 0D system with the same
    initial condition, and the whole grid reproduces the 0D answer. Running it
    here rather than quoting the 0D figure from another script keeps the
    horizon, the Hamilton constants and the integration identical between the
    two, so a difference between them is spatial coupling and not bookkeeping.
    """
    A, b_diag = H.theta_to_matrices(jnp.asarray(theta, dtype=jnp.float64))
    params = {"dt_h": cfg.dt_h, "Kp1": cfg.Kp1, "Eta": jnp.ones(5),
              "EtaPhi": jnp.ones(5), "c": cfg.c_hamilton, "alpha": cfg.alpha,
              "K_hill": jnp.array(cfg.K_hill), "n_hill": jnp.array(cfg.n_hill),
              "A": A, "b_diag": b_diag,
              "active_mask": jnp.ones(5, dtype=jnp.int64)}
    react = H._make_reaction_step(cfg.n_react_sub, cfg.newton_iters)
    D_eff = jnp.array(cfg.D_eff) if spatial else jnp.zeros(5)
    G, c = uniform_state(cfg, phi0)
    alpha = jnp.zeros((cfg.Nx, cfg.Ny), dtype=jnp.float64)

    for _ in range(cfg.n_macro):
        G = react(G, params)
        g2 = G.reshape(cfg.Nx, cfg.Ny, 12)
        phi2 = H.diffusion_step_species_2d(g2[:, :, :5], D_eff, cfg.dt_macro,
                                           cfg.dx, cfg.dy)
        g2 = g2.at[:, :, :5].set(phi2)
        g2 = g2.at[:, :, 5].set(1.0 - jnp.sum(phi2, axis=-1))
        G = g2.reshape(cfg.Nx * cfg.Ny, 12)
        phi_tot = jnp.sum(phi2 * g2[:, :, 6:11], axis=-1)      # sum_i phi_i psi_i
        alpha = alpha + cfg.dt_macro * K_ALPHA * phi_tot
        if spatial:
            c = H.nutrient_step(c, phi2, cfg, cfg.dt_macro)

    phi2 = np.asarray(G.reshape(cfg.Nx, cfg.Ny, 12)[:, :, :5])
    return float(jnp.mean(alpha)), phi2, np.asarray(c)


def separation(phi2):
    """How far the species have come apart in space.

    Each species' spatial pattern is normalised by its own mean, so this sees
    shape and not abundance; 0 means every species has the same pattern and the
    field is just one profile scaled five ways, which is the case in which
    phi_tot can still stand for the whole composition.
    """
    shapes = []
    for i in range(5):
        f = phi2[:, :, i]
        m = f.mean()
        shapes.append(f / m if m > 0 else f)
    shapes = np.array(shapes)
    return float(np.max(np.std(shapes, axis=0)))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--nx", type=int, default=16)
    ap.add_argument("--macro", type=int, default=40)
    ap.add_argument("--react-sub", type=int, default=20)
    ap.add_argument("--dt-h", type=float, default=1e-5,
                    help="ecology step; 1e-4 is the measured ceiling")
    ap.add_argument("-o", "--out", type=Path,
                    default=Path("JAXFEM/klempt2024_results/condition_spread_2d.json"))
    a = ap.parse_args(argv)

    cfg = H.Config2D(Nx=a.nx, Ny=a.nx, n_macro=a.macro, dt_h=a.dt_h,
                     n_react_sub=a.react_sub, save_every=10**9)
    print(f"grid {cfg.Nx}x{cfg.Ny}  n_macro={cfg.n_macro}  "
          f"dt_macro={cfg.dt_macro:.2e}  total t={cfg.dt_macro*cfg.n_macro:.2e}")
    print(f"D_eff={list(cfg.D_eff)}  g_cons={list(cfg.g_consumption)}\n")

    comps = day1_compositions()
    out, spreads = {}, {}
    for spatial in (False, True):
        tag = "full" if spatial else "control"
        print(f"-- {tag} "
              f"({'species diffusion + nutrient consumption on' if spatial else 'no coupling between nodes'})")
        rec = {}
        for name, _, _ in CONDITIONS:
            al, phi2, c = run(THETA_DEMO, cfg, comps[name], spatial=spatial)
            rec[name] = {"alpha_mean": al, "separation": separation(phi2),
                         "c_mean": float(c.mean())}
            print(f"   {name}  alpha={al:.6e}  separation={rec[name]['separation']:.4f}"
                  f"  c_mean={c.mean():.4f}", flush=True)
        al = np.array([rec[n]["alpha_mean"] for n, _, _ in CONDITIONS])
        spreads[tag] = float((al.max() - al.min()) / al.mean())
        out[tag] = rec
        print(f"   condition spread: {spreads[tag]:.4%}\n")

    ratio = spreads["full"] / spreads["control"] if spreads["control"] else float("inf")
    sep = max(out["full"][n]["separation"] for n, _, _ in CONDITIONS)
    print(f"control {spreads['control']:.4%}  ->  full {spreads['full']:.4%}"
          f"   ({ratio:.2f}x)")
    print(f"largest species separation reached: {sep:.4f}")
    if sep < 1e-6:
        print("  -> the species never came apart in space, so this run says"
              " nothing yet:")
        print("     no gradient developed. Run longer or on a finer grid"
              " before reading the ratio.")
    elif ratio > 1.5:
        print("  -> spatial coupling widens the condition gap")
    else:
        print("  -> spatial coupling does not restore the condition effect")

    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(
        {"grid": cfg.Nx, "n_macro": cfg.n_macro, "dt_macro": cfg.dt_macro,
         "D_eff": list(map(float, cfg.D_eff)),
         "g_consumption": list(map(float, cfg.g_consumption)),
         "conditions": out, "spreads": spreads}, indent=1))
    print(f"wrote {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
