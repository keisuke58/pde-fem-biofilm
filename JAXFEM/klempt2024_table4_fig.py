#!/usr/bin/env python3
"""Klempt et al. 2024 Table 4 (test case 4.2, biofilm on nutrients, "high"
and "low"), redrawn from this repository's reproduction next to the paper's
panels.

Field: Table 2 unchanged, growth on both faces, consumption g phi c, seed a
5 um disk one node above the nutrient face, time scale 10 ("high") and 4
("low") (klempt2024_timescale.py). Cut: the vertical mid-plane y = 10 um,
nutrient face at the bottom. Colours: the paper's banded legend
(figstyle.klempt_cmap). The paper's panels are cut from the PDF (CC BY 4.0).

    python JAXFEM/klempt2024_table4_fig.py   -> assets/fig_klempt2024_table4.png
"""
from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from PIL import Image  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "ansys_usermat"))
import figstyle  # noqa: E402
import klempt2024_case1_bc as B  # noqa: E402
import klempt2024_case2 as C  # noqa: E402
import klempt2024_quantitative as K  # noqa: E402

OUT = ROOT / "assets" / "fig_klempt2024_table4.png"
PDF = ROOT / "references" / "Klempt2024_Hamilton_biofilm_growth_BMMB.pdf"
TIMES = [0.01, 0.02, 0.04, 0.07, 0.13, 1.00]
SCALE = {"high": 10.0, "low": 4.0}
# Table 4 on page 13 at 300 dpi
PAPER_COLS = {"high": (1306, 1653), "low": (1801, 2148)}
PAPER_ROWS = [(344, 691), (706, 1053), (1068, 1415), (1431, 1778), (1793, 2139), (2156, 2501)]


def fields(sub, growth="abs"):
    s = SCALE[sub]
    K.R, K.BETA, K.K_A = 100 * s, 2 * s, 1e-3 * s
    phi, mask = C.setup(1)
    out, t = {}, 0.0
    for t_next in TIMES:
        _, phi, _ = B.run_setup(phi, mask, C.G[sub], "ic", "first_order", growth,
                                t_end=round(t_next - t, 6))
        out[t_next] = phi.copy()
        t = t_next
    return out


def paper_panels():
    with tempfile.TemporaryDirectory() as tmp:
        subprocess.run(["pdftoppm", "-f", "13", "-l", "13", "-r", "300", "-png", str(PDF),
                        str(Path(tmp) / "p")], check=True)
        im = Image.open(next(Path(tmp).glob("p*.png"))).convert("RGB")
        return {(k, r): np.asarray(im.crop((x0, y0, x1, y1)))
                for k, (x0, x1) in PAPER_COLS.items() for r, (y0, y1) in enumerate(PAPER_ROWS)}


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--growth", default="abs", help='"abs" or e.g. "blend0.5" (diagnostic)')
    ap.add_argument("--s-high", type=float, default=SCALE["high"])
    ap.add_argument("--s-low", type=float, default=SCALE["low"])
    ap.add_argument("--out", default=str(OUT))
    args = ap.parse_args(argv)
    SCALE["high"], SCALE["low"] = args.s_high, args.s_low
    figstyle.apply(size=11)
    cmap, norm = figstyle.klempt_cmap()
    F = {sub: fields(sub, args.growth) for sub in SCALE}
    P = paper_panels()
    L, ky = K.L, K.N // 2
    fig, ax = plt.subplots(len(TIMES), 4, figsize=(9.6, 13.2))
    heads = ['"high", this work', '"high", paper', '"low", this work', '"low", paper']
    for r, t in enumerate(TIMES):
        for col, (sub, kind) in enumerate((("high", "sim"), ("high", "pap"),
                                           ("low", "sim"), ("low", "pap"))):
            a = ax[r, col]
            if kind == "sim":
                a.imshow(F[sub][t][:, ky, :].T, origin="lower", extent=(0, L, 0, L),
                         cmap=cmap, norm=norm, interpolation="bilinear")
                figstyle.element_grid(a, np.arange(0, 21), np.arange(0, 21), alpha=0.3)
            else:
                a.imshow(P[(sub, r)], extent=(0, L, 0, L))
            a.set_xticks([]); a.set_yticks([]); a.grid(False)
            for sp in a.spines.values():
                sp.set_visible(True)
            if r == 0:
                a.set_title(heads[col], fontsize=10.5)
        ax[r, 0].set_ylabel(f"{round(t * 100)} %")
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    fig.colorbar(sm, ax=ax, shrink=0.4, pad=0.02, label=r"$\phi$ [-]")
    fig.suptitle("Klempt et al. 2024, test case 4.2 (Table 4): vertical mid-plane, nutrient face at "
                 "the bottom.\nThis work: Python reproduction (Table 2, time scale "
                 f"{SCALE['high']:g} / {SCALE['low']:g}, growth {args.growth}, consumption $g\\phi c$); "
                 "paper: Table 4 (CC BY 4.0)", fontsize=10.5)
    fig.savefig(args.out, dpi=200)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
