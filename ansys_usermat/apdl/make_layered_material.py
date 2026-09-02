#!/usr/bin/env python3
"""Turn a depth-resolved growth field alpha(x) into ANSYS layered materials.

Why this exists. The assignment behind this work is that the ecology model is
pointwise and should be integrated into FEM *so that spatial variation is
represented*. The ANSYS side reads the growth driver from a state slot that is
set per material (`TBDATA,10`), so it takes alpha piecewise-uniformly -- there
is no INISTATE path in use here. The Abaqus side takes a continuous field
through the temperature slot, so the two halves differ in what they can carry,
and the ANSYS route needs the field binned into layers.

That is a discretisation, and its error is controllable rather than unknown:
`layer_error` below reports what binning into N layers costs against the field
itself, so N is chosen from a number instead of a guess.

Input is the CSV written by `JAXFEM/extract_alpha_field_1d.py`
(`x_norm, x_mm, c_final, phi_total_final, alpha_x_monod, eps_growth_x_monod`),
or any CSV with a depth column and an alpha column. Both names are reported on
the first line of the output, so a wrong column is visible rather than silent.

    python ansys_usermat/apdl/make_layered_material.py alpha_field_1d.csv -n 8

Emits the `TB,USER` / `TB,STATE` blocks and the layer boundaries, ready to
paste into a deck built on the pattern of `t_growth_cylinder_shell_wrapper.dat`
(one swept volume per layer, each under its own `MAT,n`). It does not generate
geometry: the meshing in that deck was arrived at over several days of BC and
mesh work, and is not something to regenerate blindly.
"""
import argparse
import sys
from pathlib import Path

import numpy as np

# prop layout the wrapper harness unpacks (usermat_wrapper_v01_smoketest.f):
#   prop(1)=E  prop(2)=nu  prop(3)=eta  prop(4)=C01/C10  prop(5)=mtype
N_PROP = 5
N_STATE = 10          # ustatev(1:9)=Fv, ustatev(10)=alpha
TBDATA_MAX = 6        # APDL drops values past the sixth, silently


# Depth candidates are the same quantity in different units, so the order is a
# deliberate preference (mm before normalised) and picking the first is right.
# Alpha candidates are NOT interchangeable -- `alpha_x_monod` is the
# Monod-weighted field the extractor actually solves for, `alpha_x` the
# unweighted one -- so a file carrying both is ambiguous and is refused rather
# than resolved by list order. Choosing silently between two different physical
# fields is exactly the kind of failure this repository tries not to ship.
# Two producers feed this, and they do not agree on names:
#   extract_alpha_field_1d.py            x_mm / x_norm      alpha_x_monod
#   macro_eigenstrain_<condition>*.csv   depth_mm / _norm   alpha_monod
# The second is the per-condition field the results chapter needs, so both
# spellings have to be here; neither is a hypothetical.
DEPTH_COLS = ("x_mm", "depth_mm", "x_norm", "depth_norm", "x", "depth")
ALPHA_COLS = ("alpha_x_monod", "alpha_monod", "alpha_x", "alpha")


def _skip_comment_block(path):
    """Number of leading `#` lines before the header.

    The per-condition files carry a comment block, and genfromtxt cannot be
    left to handle it: with `names=True` it takes the first line as the header
    *after* stripping `#`, so a comment containing commas becomes the column
    names and the read fails with a column-count error that says nothing about
    the real cause. Counting the block and skipping it is unambiguous.
    """
    n = 0
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.lstrip().startswith("#"):
                break
            n += 1
    return n


def load_field(path, depth_col=None, alpha_col=None):
    rows = np.genfromtxt(path, delimiter=",", names=True,
                         skip_header=_skip_comment_block(path))
    names = list(rows.dtype.names)
    depth = depth_col or next((c for c in DEPTH_COLS if c in names), None)
    alpha = alpha_col
    if alpha is None:
        found = [c for c in ALPHA_COLS if c in names]
        if len(found) > 1:
            raise SystemExit(
                f"{path.name} carries more than one alpha column {found}; "
                f"they are different fields, so pass --alpha-col to say which")
        alpha = found[0] if found else None
    if depth is None or alpha is None:
        raise SystemExit(f"cannot find depth/alpha columns in {names}")
    x = np.asarray(rows[depth], dtype=float)
    a = np.asarray(rows[alpha], dtype=float)
    order = np.argsort(x)
    return x[order], a[order], depth, alpha


def bin_field(x, a, n):
    """Equal-width layers; each layer's alpha is the field's mean over it."""
    edges = np.linspace(x[0], x[-1], n + 1)
    idx = np.clip(np.digitize(x, edges[1:-1]), 0, n - 1)
    means = np.array([a[idx == k].mean() if np.any(idx == k) else np.nan
                      for k in range(n)])
    return edges, means


def layer_error(x, a, n):
    """Relative L2 error of the N-layer approximation against the field.

    Reported because the point of binning is to trade accuracy for a
    representation the solver can take, and an untracked trade is a guess.
    """
    edges, means = bin_field(x, a, n)
    idx = np.clip(np.digitize(x, edges[1:-1]), 0, n - 1)
    approx = means[idx]
    denom = np.linalg.norm(a)
    return float(np.linalg.norm(approx - a) / denom) if denom > 0 else 0.0


def material_block(mat, alpha, E, nu, eta, c01_ratio, mtype, comment=""):
    """One TB,USER + TB,STATE pair. Fv is initialised to the identity."""
    out = [f"! ---- material {mat}{(': ' + comment) if comment else ''}",
           f"TB,USER,{mat},1,{N_PROP}",
           f"TBDATA,1,{E:.6G},{nu:.6G},{eta:.6G},{c01_ratio:.6G},{mtype:.1f}",
           f"TB,STATE,{mat},,{N_STATE}",
           # Fv = I, split at six values because APDL takes no more per call
           "TBDATA,1,1.0,0.0,0.0, 0.0,1.0,0.0",
           "TBDATA,7,0.0,0.0,1.0",
           f"TBDATA,10,{alpha:.6G}                ! alpha (growth driver)"]
    return "\n".join(out)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("csv", type=Path)
    ap.add_argument("-n", "--layers", type=int, default=8)
    ap.add_argument("--E", type=float, default=0.9e-3, help="Young's modulus")
    ap.add_argument("--nu", type=float, default=0.125)
    ap.add_argument("--eta", type=float, default=0.0, help="0 = elastic")
    ap.add_argument("--c01-ratio", type=float, default=0.0)
    ap.add_argument("--mtype", type=float, default=0.0)
    ap.add_argument("--substrate", action="store_true",
                    help="emit material 1 as a non-growing substrate")
    ap.add_argument("--convergence", action="store_true",
                    help="report binning error against N and stop")
    ap.add_argument("--depth-col", default=None,
                    help=f"depth column name; default the first of {DEPTH_COLS} present")
    ap.add_argument("--alpha-col", default=None,
                    help=f"alpha column name; required if the CSV carries "
                         f"more than one of {ALPHA_COLS}")
    a = ap.parse_args(argv)

    x, alpha, dcol, acol = load_field(a.csv, a.depth_col, a.alpha_col)
    print(f"! generated from {a.csv.name} ({dcol} vs {acol}), "
          f"{len(x)} nodes, alpha in [{alpha.min():.4G}, {alpha.max():.4G}]")

    if a.convergence:
        print("! layers   relative L2 error of the binned field")
        for n in (1, 2, 4, 8, 16, 32):
            if n <= len(x):
                print(f"!   {n:3d}      {layer_error(x, alpha, n):.4G}")
        return 0

    edges, means = bin_field(x, alpha, a.layers)
    err = layer_error(x, alpha, a.layers)
    print(f"! {a.layers} layers, relative L2 error vs the field: {err:.3G}")
    print(f"! layer edges ({dcol}): "
          + ", ".join(f"{e:.5G}" for e in edges))
    print()

    mat = 1
    if a.substrate:
        print(material_block(mat, 0.0, a.E, a.nu, a.eta, a.c01_ratio, a.mtype,
                             "substrate, does not grow"))
        print()
        mat += 1
    for k, m in enumerate(means):
        lo, hi = edges[k], edges[k + 1]
        print(material_block(mat, float(m), a.E, a.nu, a.eta, a.c01_ratio,
                             a.mtype, f"layer {k+1}/{a.layers}, "
                             f"{dcol} {lo:.4G}..{hi:.4G}"))
        print()
        mat += 1

    print(f"! {mat-1} materials emitted. Assign one swept volume per layer with")
    print("! MAT,n before each VSWEEP, as in t_growth_cylinder_shell_wrapper.dat,")
    print("! then check the result with check_deck.py before reporting numbers.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
