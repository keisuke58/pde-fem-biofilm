"""klempt2024_quantitative.py -- quantitative check of Klempt, Soleimani,
Wriggers & Junker, "A Hamilton principle-based model for diffusion-driven
biofilm growth", Biomech Model Mechanobiol 23:2091-2113 (2024),
doi:10.1007/s10237-024-01883-x.

Solves the paper's evolution equations exactly as printed (Eq. 34-36) with the
Table 2 parameters, in 3D on the paper's own mesh resolution (20 um cube, 1 um
spacing -> 21^3 nodes), over the paper's normalised time T* in [0, 1], and
compares the domain averages of phi and c with the curves the paper plots:

  Fig. 4  test case 1 (nutrient at one edge, sphere r = 5 um at the centre)
  Fig. 7  "biofilm on nutrients" (c = 1 on the bottom face, disk d = 5 um just
          above it), consumption g = 1e8 ("high") and g = 1e10 ("low")

    Eq. 34  phi_dot - beta lap(phi) - k_a alpha
                    + |grad phi| r c/(k+c) n_gradphi . n_gradc = 0
    Eq. 35  c_dot - d lap(c) + g phi = 0          (quasi-static: d >> growth)
    Eq. 36  alpha_dot - k_a phi = 0                (Fg = alpha I -> alpha(0) = 1)

Consumption is run two ways, because the printed form does not reproduce the
paper's own Fig. 7 (see README note in klempt2024_results/):
    "printed"   : g phi        (Eq. 24/35 as printed, zero order in c; c is
                                clipped at 0 where it would go negative)
    "first_order": g phi c     (first order in c)

This is an independent finite-difference implementation (upwind for the
growth/transport term, 7-point Laplacian, Neumann walls by reflection), not
the paper's AceGen FEM, so agreement is expected at the level of the curves'
shape and values, not digit for digit. The paper's curves are read off the
figures by eye (+-0.01..0.02).

RESULT (2026-09-29, klempt2024_results/summary.json): not reproduced.
  - As printed (growth = transport, consumption g phi): far off on every
    figure (Fig. 4 max |diff| phi 0.63, Fig. 7 "high" phi 0.98).
  - growth |grad phi . n_c| + consumption g phi c: Fig. 4 close in shape and
    within ~0.17 (phi) / 0.10 (c); Fig. 7 grows ~5-10x too slowly ("high"
    reaches phi 0.79 at T* = 1 where the paper fills the cube by 0.2), though
    its c(phi) relation tends to the paper's 0.49 plateau as phi -> 1, which
    supports first-order consumption.
  - klempt2024_sensitivity.py has since tested the guesses this file used to
    name. The one-node seed is cleared directly: six nodes only reaches
    phi(0.20) = 0.238 against the paper's 1.000, and it scales with the seed
    rather than changing the rate. Nutrient starvation is cleared too --
    consumption a hundredfold weaker moves 0.105 to 0.108. The upwind front
    speed was NOT tested directly (no grid-refinement study was run); what can
    be said is that the two terms below account for the gap on their own,
    which leaves little for it to explain.
  - Where it does live, one change at a time from this file's own baseline:
    dropping Eq. 34's n_gradphi . n_gradc projection is worth 3.2x (0.105 ->
    0.333 -- the projection is ~0 on the colony's sides, so it grows upward as
    a column instead of spreading, and filling a cube needs the spreading),
    and K_M = 0.01 in place of Table 2's 1.0 is worth another 2.7x (0.333 ->
    0.886, against the paper's 1.000). Table 2's K_M with the paper's own
    plotted c holds f = c/(K_M+c) below 0.5.
  - ROOT CAUSE, 2026-09-30, read off the paper itself (it is bundled at the
    repository root; see THIRD_PARTY.md). The two variants that fit best are
    the two the paper does not use, and the combination the paper does use is
    the worst fit here. That is the finding, and it is not a small discrepancy.
      * Growth. Eq. 34 is phi_dot - beta lap(phi) - k_a alpha
        + ||grad phi|| (r c)/(k+c) n_gradphi . n_gradc = 0, a plain dot
        product with no absolute value, and sec. 4.1 leans on exactly that:
        "the gradient of biofilm and the gradient of nutrients are almost
        perpendicular to each other resulting in a small value for the vector
        product and consequently in minimal to no growth", which is what gives
        Fig. 3 its egg shape. So growth="abs" is not a reading of the paper,
        it contradicts the mechanism the paper describes.
      * Consumption. Eq. 35 is c_dot - d lap(c) + g phi = 0 and Eq. 24 fixes
        it: g*_c = g_bar phi, "the simplest possible functional dependency, a
        linear relation". Zeroth order in c. So consumption="first_order" is
        not the paper either.
      * The constants are right, which removes the other suspicion. Table 2's
        values are post-division: the paper sets eta_phi = eta_c = 1e-10 and
        says "parameters which have been divided by their respective eta will
        lose their bar", so d = 1e10, beta = 2, k_a = 1e-3, k = 1, g = 1e8 and
        r = 100 are already the coefficients of Eqs. 34-36 and no eta enters
        separately. The eta_phi visible in Table 1's weak form is the
        pre-division form of the same thing.
      * One genuine inconsistency in the paper: Table 1's step 3a solves
        (alpha_n1 - alpha_n)/((1+alpha_n1) dt) - (k_a/alpha_n)(phi_n1 -
        phi_n)/dt = 0, which is driven by phi_dot, while Eq. 36 is
        alpha_dot = k_a phi, driven by phi. Which one was run is not
        recoverable from the text. It barely moves phi here (k_a alpha ~ 1e-3
        against a growth term of order 50), so it is not the factor of ten,
        but it does mean "the paper's alpha equation" is ambiguous.
  - What is left, and not yet separated: an explicit upwind finite-difference
    scheme against the paper's implicit Galerkin FEM with bisection to as many
    as 1e6 substeps; a fixed grid here against a domain that swells there,
    since Fg feeds back into the geometry the averages are taken over; and the
    clip of phi to [0,1] here, which the paper does not have.
  - So this is a localisation, NOT a reproduction, and the earlier "Fig. 4
    within 0.17/0.10" must not be quoted as agreement with Klempt 2024: it is
    the agreement of a variant the paper does not use. Either Fig. 7 came from
    something other than the
    literal Eq. 34 / Table 2, or c is normalised differently than its axis
    suggests -- a question for the authors. Treat the 2024 PDE as not
    independently reproduced.

Nothing here uses this repo's 5-species model: Klempt 2024 is a different
(single-species, interface-growth) PDE. Only Eq. 36, alpha_dot = k_a phi, is
shared with this repo's growth law.

    python JAXFEM/klempt2024_quantitative.py            # all cases, both variants
"""
import json
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

# ---- Table 2 ---------------------------------------------------------------
D = 1e10        # nutrient diffusivity   [um^2 / T*]
BETA = 2.0      # phase-field regularisation [um^2 / T*]
K_A = 1e-3      # growth factor          [1 / T*]
K_M = 1.0       # half-velocity constant [-]
G_TABLE = 1e8   # consumption            [1 / T*]
R = 100.0       # growth parameter       [1 / T*]  (um / T* as a front speed)

L = 20.0
N = 21
H = L / (N - 1)
DT = 1e-3       # the paper's nominal 1e3 substeps over T* in [0, 1]
T_END = 1.0

# ---- the paper's curves, read off the figures ------------------------------
PAPER = {
    "fig4_edge": {
        "t":   [0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00],
        "phi": [0.13, 0.21, 0.28, 0.325, 0.41, 0.47, 0.53, 0.58, 0.63, 0.67, 0.705, 0.74],
        "c":   [0.33, 0.24, 0.19, 0.165, 0.135, 0.12, 0.11, 0.105, 0.098, 0.093, 0.089, 0.085],
    },
    "fig7_high": {
        "t":   [0.01, 0.05, 0.10, 0.15, 0.20, 0.50, 1.00],
        "phi": [0.02, 0.12, 0.55, 0.90, 1.00, 1.00, 1.00],
        "c":   [1.00, 0.80, 0.55, 0.50, 0.49, 0.49, 0.49],
    },
    "fig7_low": {
        "t":   [0.01, 0.05, 0.10, 0.15, 0.20, 0.50, 1.00],
        "phi": [0.02, 0.07, 0.17, 0.22, 0.25, 0.29, 0.30],
        "c":   [0.80, 0.42, 0.15, 0.10, 0.09, 0.08, 0.075],
    },
}

_ax = np.arange(N) * H
X, Y, Z = np.meshgrid(_ax, _ax, _ax, indexing="ij")
_w1 = np.ones(N); _w1[0] = _w1[-1] = 0.5
W = _w1[:, None, None] * _w1[None, :, None] * _w1[None, None, :]
W = W / W.sum()        # trapezoid weights = domain average of the trilinear field


def avg(u):
    return float((W * u).sum())


def lap(u):
    p = np.pad(u, 1, mode="reflect")      # zero-flux walls
    return (p[2:, 1:-1, 1:-1] + p[:-2, 1:-1, 1:-1] + p[1:-1, 2:, 1:-1]
            + p[1:-1, :-2, 1:-1] + p[1:-1, 1:-1, 2:] + p[1:-1, 1:-1, :-2]
            - 6.0 * u) / H**2


def grad_c(c):
    p = np.pad(c, 1, mode="reflect")
    return ((p[2:, 1:-1, 1:-1] - p[:-2, 1:-1, 1:-1]) / (2 * H),
            (p[1:-1, 2:, 1:-1] - p[1:-1, :-2, 1:-1]) / (2 * H),
            (p[1:-1, 1:-1, 2:] - p[1:-1, 1:-1, :-2]) / (2 * H))


def upwind_dot(phi, v):
    """v . grad(phi), upwinded per axis."""
    p = np.pad(phi, 1, mode="reflect")
    out = np.zeros_like(phi)
    sl = [(slice(2, None), slice(1, -1), slice(1, -1)),
          (slice(1, -1), slice(2, None), slice(1, -1)),
          (slice(1, -1), slice(1, -1), slice(2, None))]
    sr = [(slice(None, -2), slice(1, -1), slice(1, -1)),
          (slice(1, -1), slice(None, -2), slice(1, -1)),
          (slice(1, -1), slice(1, -1), slice(None, -2))]
    for i in range(3):
        fwd = (p[sl[i]] - phi) / H
        bwd = (phi - p[sr[i]]) / H
        out += np.where(v[i] > 0, v[i] * bwd, v[i] * fwd)
    return out


def neumann_laplacian():
    idx = np.arange(N**3).reshape(N, N, N)
    rows, cols, vals = [], [], []
    for axis in range(3):
        for sgn in (-1, 1):
            nb = np.roll(idx, -sgn, axis=axis)
            # reflect at the walls: the missing neighbour is the interior one
            edge = [slice(None)] * 3
            edge[axis] = -1 if sgn == 1 else 0
            inner = [slice(None)] * 3
            inner[axis] = -2 if sgn == 1 else 1
            nb[tuple(edge)] = idx[tuple(inner)]
            rows.append(idx.ravel()); cols.append(nb.ravel())
            vals.append(np.full(N**3, 1.0 / H**2))
    rows.append(idx.ravel()); cols.append(idx.ravel())
    vals.append(np.full(N**3, -6.0 / H**2))
    return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows),
                          np.concatenate(cols))), shape=(N**3, N**3))


LAP = neumann_laplacian()


def solve_c(phi, g, dir_mask, variant):
    """Quasi-static Eq. 35: d lap(c) = g phi [c]; c = 1 on dir_mask."""
    free = ~dir_mask.ravel()
    q = (g / D) * phi.ravel()
    if variant == "first_order":
        A = -LAP + sp.diags(q)
        rhs = np.zeros(N**3)
    else:
        A = -LAP
        rhs = -q
    Aff = A[free][:, free]
    rhs_f = rhs[free] - A[free][:, ~free] @ np.ones((~free).sum())
    c = np.ones(N**3)
    c[free] = spla.spsolve(Aff.tocsc(), rhs_f)
    return np.clip(c, 0.0, 1.0).reshape(N, N, N)


def setup(case):
    if case == "fig4_edge":
        phi = ((X - 10)**2 + (Y - 10)**2 + (Z - 10)**2 <= 5.0**2).astype(float)
        mask = (np.isclose(X, L) & np.isclose(Y, L))       # one edge of the cube
        g = G_TABLE
    else:
        phi = (np.isclose(Z, H) & ((X - 10)**2 + (Y - 10)**2 <= 2.5**2)).astype(float)
        mask = np.isclose(Z, 0.0)                           # bottom face
        g = 1e8 if case == "fig7_high" else 1e10
    return phi, mask, g


def run(case, variant, growth="printed", dt=DT, t_end=T_END):
    """growth="printed": Eq. 34 as printed, |grad phi| n_phi.n_c = grad phi . n_c
    -- a transport of phi towards the nutrient (the back of the colony erodes).
    growth="abs": |grad phi . n_c| as a source -- growth on both faces aligned
    with the nutrient gradient, none on the faces perpendicular to it, which
    is what the paper's text and Fig. 3/Table 4 describe."""
    phi, mask, g = setup(case)
    alpha = np.ones_like(phi)
    c = solve_c(phi, g, mask, variant)
    n = int(round(t_end / dt))
    rec = {"t": [0.0], "phi": [avg(phi)], "c": [avg(c)]}
    for s in range(1, n + 1):
        gx, gy, gz = grad_c(c)
        mag = np.sqrt(gx**2 + gy**2 + gz**2)
        speed = np.where(mag > 1e-14, R * c / (K_M + c) / np.maximum(mag, 1e-14), 0.0)
        v = (speed * gx, speed * gy, speed * gz)
        if growth == "abs":
            src = np.abs(upwind_dot(phi, v))
        else:
            src = -upwind_dot(phi, v)
        phi = phi + dt * (BETA * lap(phi) + K_A * alpha + src)
        phi = np.clip(phi, 0.0, 1.0)
        alpha = alpha + dt * K_A * phi
        c = solve_c(phi, g, mask, variant)
        if s % max(1, n // 100) == 0:
            rec["t"].append(s * dt); rec["phi"].append(avg(phi)); rec["c"].append(avg(c))
    return rec


def compare(rec, paper):
    t = np.array(rec["t"])
    out = {}
    for k in ("phi", "c"):
        sim = np.interp(paper["t"], t, rec[k])
        out[k] = {"sim": [round(float(v), 4) for v in sim],
                  "paper": paper[k],
                  "max_abs_diff": round(float(np.max(np.abs(sim - np.array(paper[k])))), 4)}
    return out


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--growth", nargs="+", default=["printed", "abs"])
    ap.add_argument("--consumption", nargs="+", default=["printed", "first_order"])
    ap.add_argument("--cases", nargs="+", default=list(PAPER))
    a = ap.parse_args(argv)
    out = Path(__file__).resolve().parent / "klempt2024_results"
    out.mkdir(exist_ok=True)
    f = out / "summary.json"
    try:
        results = json.loads(f.read_text())
    except (OSError, ValueError):
        results = {}
    for case in a.cases:
        for growth in a.growth:
          for variant in a.consumption:
            t0 = time.time()
            rec = run(case, variant, growth)
            cmp = compare(rec, PAPER[case])
            results[f"{case}/growth={growth}/consumption={variant}"] = cmp
            print(f"\n== {case}  growth={growth}  consumption={variant}  "
                  f"({time.time() - t0:.0f}s)", flush=True)
            print("   t     " + "  ".join(f"{t:5.2f}" for t in PAPER[case]["t"]))
            for k in ("phi", "c"):
                print(f"   {k:3s} sim " + "  ".join(f"{v:5.3f}" for v in cmp[k]["sim"]))
                print(f"   {k:3s} pap " + "  ".join(f"{v:5.3f}" for v in cmp[k]["paper"]))
                print(f"   max |diff| {k}: {cmp[k]['max_abs_diff']}", flush=True)
            f.write_text(json.dumps(results, indent=1))
    print(f"\nwrote {f}")


if __name__ == "__main__":
    main()
