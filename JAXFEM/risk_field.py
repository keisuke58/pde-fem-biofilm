#!/usr/bin/env python3
"""risk_field.py
===============
T2 / Fig4 — per-location uncertainty + risk fields.

Where `risk_metric.py` reports the *scalar* peak-Mises exceedance per condition,
this module lifts the same idea to a **spatial field**: given a posterior stress
stack over mesh nodes it produces the two Fig4 panels from
`research_goals_1_2.md`:

  (a) a stress **credible band** (mean + p05/p95) along a periodontal-pocket line;
  (b) a jaw-surface **risk map**  P[σ_node > τ]  per node.

Data contract (same layout as ``aggregate_di_credible.py``)
----------------------------------------------------------
A stack directory contains:
  sigma_stack.npy   (n_samples, n_nodes)   posterior peak/nodal von Mises [MPa]
  coords.npy        (n_nodes, 3)           node coordinates
  line_nodes.npy    (n_line,)  [optional]  ordered node indices of the pocket line

If ``line_nodes.npy`` is absent, a line is auto-extracted (a thin slab around the
principal axis of the point cloud, ordered by arc length) so the module runs on
any stack.

Because a clinically calibrated stress threshold does not yet exist (see
`VERIFICATION_SENSITIVITY_LIMITATIONS.md` L1), the risk map is a **relative**
(condition-to-condition) read-out; the threshold is a CLI parameter (kPa) and the
mean/credible band is the primary uncertainty representation.

Outputs
-------
  JAXFEM/_risk/
    risk_field_summary_{tag}.json     (tracked; global stats + band arrays)
    risk_field_{tag}.csv              (per-node x,y,z,mean,p05,p50,p95,risk)
    risk_map_{tag}.png                (gitignored; jaw-surface risk scatter)
    pocket_band_{tag}.png             (gitignored; credible band along the line)

Usage
-----
  python JAXFEM/risk_field.py --stack-dir _di_credible/commensal_hobic --tag CH
  python JAXFEM/risk_field.py --demo            # synthetic stack, end-to-end
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
_OUT = _HERE / "_risk"
MPA_PER_KPA = 1.0e-3
QTAGS = ("p05", "p50", "p95")
QVALS = (0.05, 0.50, 0.95)


# ── core field statistics ────────────────────────────────────────────────────
def risk_field(sigma_stack: np.ndarray, tau_mpa: float) -> dict:
    """Per-node statistics from a (n_samples, n_nodes) stress stack.

    Returns dict of (n_nodes,) arrays: mean, p05, p50, p95, risk = P[σ>τ].
    """
    s = np.asarray(sigma_stack, dtype=float)
    if s.ndim != 2:
        raise ValueError(f"sigma_stack must be 2-D (n_samples, n_nodes); got {s.shape}")
    q = np.quantile(s, QVALS, axis=0)  # (3, n_nodes)
    return {
        "mean": s.mean(axis=0),
        "p05": q[0],
        "p50": q[1],
        "p95": q[2],
        "risk": np.mean(s > tau_mpa, axis=0),
        "n_samples": s.shape[0],
        "n_nodes": s.shape[1],
    }


def _principal_axis(coords: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (centroid, unit principal direction) of a point cloud."""
    c = coords.mean(axis=0)
    X = coords - c
    # leading right-singular vector = direction of maximum variance
    _, _, vt = np.linalg.svd(X, full_matrices=False)
    return c, vt[0]


def auto_pocket_line(coords: np.ndarray, slab_frac: float = 0.15) -> np.ndarray:
    """Auto-extract an ordered node line: a thin slab hugging the principal axis,
    ordered by projection (arc length) along that axis.

    Nodes whose transverse distance from the axis is within ``slab_frac`` of the
    cloud's transverse spread are kept. Falls back to all nodes if the slab is
    empty. Deterministic (no RNG)."""
    coords = np.asarray(coords, dtype=float)
    c, u = _principal_axis(coords)
    X = coords - c
    t = X @ u                      # arc length along axis (n_nodes,)
    perp = X - np.outer(t, u)      # transverse component
    d = np.linalg.norm(perp, axis=1)
    scale = d.max() if d.max() > 0 else 1.0
    mask = d <= slab_frac * scale
    if mask.sum() < 2:
        mask = np.ones(len(coords), dtype=bool)
    idx = np.where(mask)[0]
    return idx[np.argsort(t[idx])]


def pocket_line_band(
    sigma_stack: np.ndarray,
    coords: np.ndarray,
    line_nodes: np.ndarray | None = None,
) -> dict:
    """Credible band (mean, p05, p50, p95) along an ordered node line vs arc length."""
    coords = np.asarray(coords, dtype=float)
    if line_nodes is None:
        line_nodes = auto_pocket_line(coords)
    line_nodes = np.asarray(line_nodes, dtype=int)

    pts = coords[line_nodes]
    seg = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    arclen = np.concatenate([[0.0], np.cumsum(seg)])  # (n_line,)

    s = np.asarray(sigma_stack, dtype=float)[:, line_nodes]  # (n_samp, n_line)
    q = np.quantile(s, QVALS, axis=0)
    return {
        "line_nodes": line_nodes,
        "arclen": arclen,
        "mean": s.mean(axis=0),
        "p05": q[0],
        "p50": q[1],
        "p95": q[2],
    }


# ── IO ───────────────────────────────────────────────────────────────────────
def load_stack(stack_dir: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
    """Load (sigma_stack, coords, line_nodes|None) from a stack directory."""
    stack_dir = Path(stack_dir)
    sig_p = stack_dir / "sigma_stack.npy"
    crd_p = stack_dir / "coords.npy"
    if not sig_p.exists() or not crd_p.exists():
        raise FileNotFoundError(
            f"expected {sig_p.name} + {crd_p.name} in {stack_dir} "
            f"(produced by the posterior FEM ensemble; see aggregate_di_credible.py)"
        )
    sigma = np.load(sig_p)
    coords = np.load(crd_p)
    line_p = stack_dir / "line_nodes.npy"
    line = np.load(line_p) if line_p.exists() else None
    return sigma, coords, line


def synthesize(n_samples: int = 40, nx: int = 12, seed: int = 0):
    """Deterministic synthetic stack for tests/--demo.

    A 3-D grid whose stress peaks in a pocket band (mid-plane), with node-wise
    posterior scatter that grows toward the peak — a stand-in for the real
    ensemble so the plumbing is exercised end-to-end.
    """
    rng = np.random.default_rng(seed)
    lin = np.linspace(0.0, 1.0, nx)
    gx, gy, gz = np.meshgrid(lin, lin, lin, indexing="ij")
    coords = np.column_stack([gx.ravel(), gy.ravel(), gz.ravel()])
    # mean field (MPa): Gaussian ridge along x at the y,z mid-plane (the "pocket")
    r2 = (coords[:, 1] - 0.5) ** 2 + (coords[:, 2] - 0.5) ** 2
    mean = 2.0e-3 + 12.0e-3 * np.exp(-r2 / 0.02)      # 2 → 14 kPa
    scatter = 0.15 * mean + 1.0e-3
    stack = mean[None, :] + rng.normal(0.0, 1.0, (n_samples, coords.shape[0])) * scatter[None, :]
    stack = np.clip(stack, 0.0, None)
    return stack, coords


# ── plotting ─────────────────────────────────────────────────────────────────
def plot_pocket_band(band: dict, out_path: Path, tag: str) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    a = band["arclen"]
    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    ax.fill_between(a, band["p05"] / MPA_PER_KPA, band["p95"] / MPA_PER_KPA,
                    alpha=0.30, color="#4c78a8", label="90% credible band")
    ax.plot(a, band["mean"] / MPA_PER_KPA, color="#1f4e79", lw=2.0, label="mean")
    ax.plot(a, band["p50"] / MPA_PER_KPA, color="#1f4e79", lw=1.0, ls="--", label="median")
    ax.set_xlabel("arc length along pocket line")
    ax.set_ylabel("von Mises σ (kPa)")
    ax.set_title(f"Pocket-line stress credible band — {tag}")
    ax.legend(frameon=False)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_risk_map(coords: np.ndarray, risk: np.ndarray, out_path: Path, tag: str,
                  tau_kpa: float) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # project onto the two highest-variance axes for a 2-D surface view
    c, u = _principal_axis(coords)
    X = coords - c
    _, _, vt = np.linalg.svd(X, full_matrices=False)
    p = X @ vt[0]
    q = X @ vt[1]
    fig, ax = plt.subplots(figsize=(5.6, 4.4))
    sc = ax.scatter(p, q, c=risk, cmap="inferno", vmin=0.0, vmax=1.0, s=14)
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label(f"P[σ > {tau_kpa:g} kPa]")
    ax.set_xlabel("principal axis 1")
    ax.set_ylabel("principal axis 2")
    ax.set_title(f"Growth-stress risk map — {tag}")
    ax.set_aspect("equal", "box")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def write_field_csv(coords: np.ndarray, fld: dict, out_path: Path) -> None:
    import csv

    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["x", "y", "z", "mean_kpa", "p05_kpa", "p50_kpa", "p95_kpa", "risk"])
        for i in range(coords.shape[0]):
            w.writerow([
                f"{coords[i,0]:.6g}", f"{coords[i,1]:.6g}", f"{coords[i,2]:.6g}",
                f"{fld['mean'][i]/MPA_PER_KPA:.6g}", f"{fld['p05'][i]/MPA_PER_KPA:.6g}",
                f"{fld['p50'][i]/MPA_PER_KPA:.6g}", f"{fld['p95'][i]/MPA_PER_KPA:.6g}",
                f"{fld['risk'][i]:.4f}",
            ])


# ── driver ───────────────────────────────────────────────────────────────────
def build_fig4(
    sigma_stack: np.ndarray,
    coords: np.ndarray,
    tag: str,
    threshold_kpa: float,
    line_nodes: np.ndarray | None,
    out_dir: Path,
    make_plots: bool = True,
) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    tau = threshold_kpa * MPA_PER_KPA

    fld = risk_field(sigma_stack, tau)
    band = pocket_line_band(sigma_stack, coords, line_nodes)

    summary = {
        "tag": tag,
        "threshold_kpa": threshold_kpa,
        "n_samples": int(fld["n_samples"]),
        "n_nodes": int(fld["n_nodes"]),
        "risk_field": {
            "min": round(float(fld["risk"].min()), 4),
            "mean": round(float(fld["risk"].mean()), 4),
            "max": round(float(fld["risk"].max()), 4),
            "frac_nodes_risk_gt_0.5": round(float(np.mean(fld["risk"] > 0.5)), 4),
        },
        "peak_stress_kpa": {
            "p50_of_nodal_max": round(float(np.median(np.max(sigma_stack, axis=1)) / MPA_PER_KPA), 4),
            "max_mean_node": round(float(fld["mean"].max() / MPA_PER_KPA), 4),
        },
        "pocket_line": {
            "n_line_nodes": int(len(band["line_nodes"])),
            "arclen": [round(float(x), 5) for x in band["arclen"]],
            "mean_kpa": [round(float(x) / MPA_PER_KPA, 4) for x in band["mean"]],
            "p05_kpa": [round(float(x) / MPA_PER_KPA, 4) for x in band["p05"]],
            "p95_kpa": [round(float(x) / MPA_PER_KPA, 4) for x in band["p95"]],
        },
        "note": (
            "Relative (condition-to-condition) risk map; threshold is a modelling "
            "choice, not a calibrated constant. Band is the primary uncertainty "
            "read-out. Data contract: sigma_stack (n_samples,n_nodes) + coords (n_nodes,3)."
        ),
    }

    summary_path = out_dir / f"risk_field_summary_{tag}.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    write_field_csv(coords, fld, out_dir / f"risk_field_{tag}.csv")

    if make_plots:
        try:
            plot_pocket_band(band, out_dir / f"pocket_band_{tag}.png", tag)
            plot_risk_map(coords, fld["risk"], out_dir / f"risk_map_{tag}.png", tag, threshold_kpa)
        except Exception as exc:  # figures are non-essential
            print(f"  [plot] skipped: {exc}")

    return summary


def _print_summary(summary: dict) -> None:
    s = summary
    print(f"\nFig4 risk field — {s['tag']}   τ = {s['threshold_kpa']:g} kPa")
    print(f"  nodes={s['n_nodes']}  samples={s['n_samples']}")
    rf = s["risk_field"]
    print(f"  risk field  P[σ>τ]:  min={rf['min']:.2f}  mean={rf['mean']:.2f}  "
          f"max={rf['max']:.2f}   nodes>0.5: {rf['frac_nodes_risk_gt_0.5']*100:.0f}%")
    ps = s["peak_stress_kpa"]
    print(f"  peak stress: median nodal-max={ps['p50_of_nodal_max']:.2f} kPa  "
          f"max mean-node={ps['max_mean_node']:.2f} kPa")
    print(f"  pocket line: {s['pocket_line']['n_line_nodes']} nodes")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Fig4 per-location risk fields (T2)")
    ap.add_argument("--stack-dir", help="dir with sigma_stack.npy + coords.npy (+ line_nodes.npy)")
    ap.add_argument("--tag", default="field", help="output tag / condition label")
    ap.add_argument("--threshold-kpa", type=float, default=5.0)
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--no-plot", action="store_true")
    ap.add_argument("--demo", action="store_true", help="run on a synthetic stack end-to-end")
    args = ap.parse_args(argv)

    out_dir = Path(args.out_dir) if args.out_dir else _OUT

    if args.demo or not args.stack_dir:
        if not args.demo:
            print("no --stack-dir given; running --demo on a synthetic stack\n")
        sigma, coords = synthesize()
        line = None
        tag = args.tag if args.tag != "field" else "demo"
    else:
        sigma, coords, line = load_stack(args.stack_dir)
        tag = args.tag

    summary = build_fig4(
        sigma, coords, tag, args.threshold_kpa, line, out_dir,
        make_plots=not args.no_plot,
    )
    _print_summary(summary)
    print(f"\n  summary → {out_dir / f'risk_field_summary_{tag}.json'}")
    print(f"  field   → {out_dir / f'risk_field_{tag}.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
