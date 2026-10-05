#!/usr/bin/env python3
"""Which isotropic share w of the front term (KLEMPT2024_REPRODUCTION.md sec. 15)
gives Klempt 2024's shapes? The mean curves do not decide w (0.35-0.7 within
0.003); here the paper's contour panels do.

The panels (Table 3: 4.1, diagonal cut, phi column; Table 4: 4.2 "high" and
"low", vertical mid-plane) are read back into phi by their banded colours
(figstyle.KLEMPT_PANEL_COLORS, 7 bands of 1/7; grid lines and anti-aliased
pixels, farther than 80 in RGB from every band colour, are left out) on a
cell grid, and compared with the finite-difference reproduction
(klempt2024_case1_bc.run_setup, growth "blend<w>", consumption g phi c, the
time scale per run that fits the mean curves best) on the same cut:
  IoU   overlap of the regions phi >= 4/7 (paper band 5-7), per panel;
  MAE   mean |phi_sim - band centre| over the cells.
Averaged over the panels of a case. A DIAGNOSTIC of a hypothesis about what
the paper computed, not the paper's equations.

    python JAXFEM/klempt2024_shape_scan.py   -> JAXFEM/klempt2024_results/shape_scan.json
"""
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from PIL import Image
from scipy.interpolate import RegularGridInterpolator

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path[:0] = [str(HERE), str(ROOT / "ansys_usermat")]
import figstyle  # noqa: E402
import klempt2024_case1_bc as B  # noqa: E402
import klempt2024_case2 as C  # noqa: E402
import klempt2024_quantitative as K  # noqa: E402

PDF = ROOT / "references" / "Klempt2024_Hamilton_biofilm_growth_BMMB.pdf"
OUT = HERE / "klempt2024_results" / "shape_scan.json"
SCAN = HERE / "klempt2024_results" / "blend_scan.json"
CASES = {
    "4.1": {"page": 11, "cols": (957, 1287), "rows": [(291, 524), (536, 770), (782, 1015), (1028, 1262)],
            "times": [0.05, 0.25, 0.45, 1.00], "grid": (30, 42)},
    "4.2 high": {"page": 13, "cols": (1306, 1653),
                 "rows": [(344, 691), (706, 1053), (1068, 1415), (1431, 1778), (1793, 2139), (2156, 2501)],
                 "times": [0.01, 0.02, 0.04, 0.07, 0.13, 1.00], "grid": (40, 40)},
    "4.2 low": {"page": 13, "cols": (1801, 2148),
                "rows": [(344, 691), (706, 1053), (1068, 1415), (1431, 1778), (1793, 2139), (2156, 2501)],
                "times": [0.01, 0.02, 0.04, 0.07, 0.13, 1.00], "grid": (40, 40)},
}
# best time scale per w from the curves (blend_scan.json; w = 0 and 0.25 from sec. 15's table)
EXTRA_S = {0.0: {"4.1": 2.0, "4.2 high": 10.0, "4.2 low": 4.0},
           0.25: {"4.1": 1.5, "4.2 high": 6.0, "4.2 low": 4.0}}
WS = [0.0, 0.25, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7]
COLORS = np.array([[int(c[i:i + 2], 16) for i in (1, 3, 5)] for c in figstyle.KLEMPT_PANEL_COLORS], float)


def paper_fields():
    out = {}
    with tempfile.TemporaryDirectory() as tmp:
        for page in (11, 13):
            subprocess.run(["pdftoppm", "-f", str(page), "-l", str(page), "-r", "300", "-png", str(PDF),
                            str(Path(tmp) / f"p{page}")], check=True)
        for case, d in CASES.items():
            im = np.asarray(Image.open(next(Path(tmp).glob(f"p{d['page']}*.png"))).convert("RGB"), float)
            ny, nx = d["grid"]
            for r, (y0, y1) in enumerate(d["rows"]):
                pan = im[y0:y1, d["cols"][0]:d["cols"][1]][::-1]          # origin at the bottom
                dist = np.linalg.norm(pan[:, :, None, :] - COLORS[None, None], axis=-1)
                band = dist.argmin(-1)
                ok = dist.min(-1) < 80
                val = np.full((ny, nx), np.nan)
                H, W = band.shape
                for j in range(ny):
                    for i in range(nx):
                        sl = (slice(j * H // ny, (j + 1) * H // ny), slice(i * W // nx, (i + 1) * W // nx))
                        b = band[sl][ok[sl]]
                        if b.size:
                            val[j, i] = (np.median(b) + 0.5) / 7
                out[(case, r)] = val
    return out


def cut(case, phi):
    n = phi.shape[0]
    if case == "4.1":                                  # diagonal plane x = z, vertical y
        u = np.array([[phi[i, j, i] for i in range(n)] for j in range(n)])
        lx = K.L * np.sqrt(2)
    else:                                              # vertical mid-plane y = 10, vertical z
        u = phi[:, n // 2, :].T
        lx = K.L
    ny, nx = CASES[case]["grid"]
    f = RegularGridInterpolator((np.linspace(0, K.L, n), np.linspace(0, lx, n)), u)
    yy, xx = np.meshgrid((np.arange(ny) + 0.5) / ny * K.L, (np.arange(nx) + 0.5) / nx * lx, indexing="ij")
    return f(np.stack([yy.ravel(), xx.ravel()], -1)).reshape(ny, nx)


def run(job):
    w, case, s = job
    K.R, K.BETA, K.K_A = 100 * s, 2 * s, 1e-3 * s
    if case == "4.1":
        phi, _, _ = K.setup("fig4_edge")
        mask = (K.X >= K.L - 2 - 1e-9) & (K.Y >= K.L - 2 - 1e-9) & (K.Z >= K.L - 2 - 1e-9)
        g = 1e8
    else:
        phi, mask = C.setup(1)
        g = C.G[case.split()[1]]
    snaps, t = [], 0.0
    for tn in CASES[case]["times"]:
        _, phi, _ = B.run_setup(phi, mask, g, "ic", "first_order", f"blend{w}", t_end=round(tn - t, 6))
        snaps.append(cut(case, phi))
        t = tn
    return w, case, s, snaps


def best_s(w, case, scan):
    if w in EXTRA_S:
        return EXTRA_S[w][case]
    v = scan[f"{w:g}"][case]
    return float(min(v, key=lambda k: v[k]["score"]))


def main():
    scan = json.loads(SCAN.read_text())
    P = paper_fields()
    jobs = [(w, c, best_s(w, c, scan)) for w in WS for c in CASES]
    res = {}
    with Pool(4) as pool:
        for w, case, s, snaps in pool.imap_unordered(run, jobs):
            iou, mae = [], []
            for r, u in enumerate(snaps):
                p = P[(case, r)]
                ok = ~np.isnan(p)
                a, b = (u >= 4 / 7) & ok, (p >= 4 / 7) & ok
                iou.append(float((a & b).sum() / max((a | b).sum(), 1)))
                mae.append(float(np.abs(u - p)[ok].mean()))
            res.setdefault(f"{w:g}", {})[case] = {"s": s, "iou": iou, "mae": mae,
                                                  "iou_mean": float(np.mean(iou)), "mae_mean": float(np.mean(mae))}
            print(f"w={w:<4g} {case:9s} s={s:<4g} IoU {np.mean(iou):.3f} ({' '.join(f'{x:.2f}' for x in iou)})"
                  f"  MAE {np.mean(mae):.3f}", flush=True)
            OUT.write_text(json.dumps(res, indent=1))
    print("\n  w     4.1 IoU/MAE    4.2 high       4.2 low        mean IoU")
    for w in sorted(res, key=float):
        r = res[w]
        print(f"{w:>5s}  " + "  ".join(f"{r[c]['iou_mean']:.3f}/{r[c]['mae_mean']:.3f}" for c in CASES)
              + f"   {np.mean([r[c]['iou_mean'] for c in CASES]):.3f}")


if __name__ == "__main__":
    main()
