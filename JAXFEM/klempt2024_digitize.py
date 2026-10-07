#!/usr/bin/env python3
"""Read the curves of Klempt et al. 2024, Fig. 4 and Fig. 7, off the PDF.

The reproduction scripts compared against values read off the figures by
eye. Read by colour instead: render the page at 300 dpi
(references/Klempt2024_Hamilton_biofilm_growth_BMMB.pdf, pdftoppm), find the
grey tick marks inside the plot frame to calibrate both axes, and take the
red (phi) and blue (c) pixels column by column, the legend box masked. In
Fig. 7 two curves share each colour ("high" dashed, "low" solid); the "high"
curve lies above the "low" one for both phi and c at every time, so at each
column the upper cluster is "high" and the lower "low"; where the dashed curve
has a gap, the single cluster goes to the curve it continues. Accuracy is
about one pixel, 0.001 in value.

    python JAXFEM/klempt2024_digitize.py
        -> JAXFEM/klempt2024_results/paper_curves_digitized.json

Checked against the eye readings (2026-10-02): Fig. 7 "low" and both c curves
agree within 0.02; Fig. 7 "high" phi was read far too low early on (0.12 /
0.55 at T* = 0.05 / 0.10, the figure gives 0.35 / 0.73).
"""
from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
PDF = ROOT / "references" / "Klempt2024_Hamilton_biofilm_growth_BMMB.pdf"
OUT = Path(__file__).resolve().parent / "klempt2024_results" / "paper_curves_digitized.json"
T = np.round(np.arange(0.01, 1.0001, 0.01), 2)


def render(page, tmp):
    stem = Path(tmp) / f"p{page}"
    subprocess.run(["pdftoppm", "-f", str(page), "-l", str(page), "-r", "300", "-png",
                    str(PDF), str(stem)], check=True)
    f = next(Path(tmp).glob(f"p{page}*.png"))
    return np.array(Image.open(f).convert("RGB")).astype(int)


def groups(a, gap=3):
    if len(a) == 0:
        return []
    out, cur = [], [a[0]]
    for v in a[1:]:
        if v - cur[-1] <= gap:
            cur.append(v)
        else:
            out.append(float(np.mean(cur))); cur = [v]
    out.append(float(np.mean(cur)))
    return out


def frame(im, y_from):
    dark = im.sum(2) < 250
    H, W = dark.shape
    rows = [y for y in range(y_from, H) if dark[y, int(W * 0.45):int(W * 0.9)].mean() > 0.6]
    top, bottom = rows[0], rows[-1]
    mid = (top + bottom) // 2
    cols = [x for x in range(W) if dark[mid - 100:mid + 100, x].mean() > 0.9]
    return top, bottom, cols[0], cols[-1]


def ticks(im, top, bottom, left, right):
    grey = (im.sum(2) < 600) & (np.abs(im[:, :, 0] - im[:, :, 2]) < 30)
    bx = grey[bottom - 16:bottom - 4, left + 3:right - 3].mean(0)
    xs = groups([left + 3 + i for i, v in enumerate(bx) if v > 0.7])
    by = grey[top + 3:bottom - 3, left + 4:left + 16].mean(1)
    ys = groups([top + 3 + i for i, v in enumerate(by) if v > 0.7])
    return xs, ys


def read(im, box, x_ticks, x_vals, y_ticks, y_vals, legend, two):
    top, bottom, left, right = box
    R, G, B = im[:, :, 0], im[:, :, 1], im[:, :, 2]
    masks = {"phi": (R > 180) & (G < 90) & (B < 90), "c": (B > 180) & (R < 90) & (G < 90)}
    lx0, ly0, lx1, ly1 = legend
    for m in masks.values():
        m[ly0:ly1, lx0:lx1] = False
    ax = np.polyfit(x_vals, x_ticks, 1)                  # value -> pixel
    ay = np.polyfit(y_ticks, y_vals, 1)                  # pixel -> value
    out = {}
    for name, m in masks.items():
        hi = lo = None
        H, L = [], []
        for t in T:
            x = int(round(np.polyval(ax, t)))
            cs = []
            for dx in (0, -1, 1, -2, 2, -3, 3):
                ys = np.nonzero(m[top + 2:bottom - 2, x + dx])[0] + top + 2
                cs = [float(np.polyval(ay, y)) for y in groups(list(ys))]
                if len(cs) >= (2 if two else 1):
                    break
            if not two:
                H.append(round(float(np.mean(cs)), 4) if cs else None)
                continue
            if len(cs) >= 2:
                h, l = max(cs), min(cs)
            elif len(cs) == 1:
                v = cs[0]
                if hi is not None and lo is not None and abs(v - lo) < abs(v - hi):
                    h, l = hi, v
                else:
                    h, l = v, lo
            else:
                h, l = hi, lo
            hi, lo = h, l
            H.append(None if h is None else round(h, 4))
            L.append(None if l is None else round(l, 4))
        if two:
            out[f"{name}_high"], out[f"{name}_low"] = H, L
        else:
            out[name] = H
    return out


def main():
    res = {"t": [float(t) for t in T], "source": str(PDF.relative_to(ROOT))}
    with tempfile.TemporaryDirectory() as tmp:
        im = render(11, tmp)                              # Fig. 4, page 2101
        box = frame(im, int(im.shape[0] * 0.42))
        xs, ys = ticks(im, *box)
        xv = [round(0.1 * i, 1) for i in range(len(xs))]
        yv = [round(0.8 - 0.1 * i, 1) for i in range(len(ys))]
        top, bottom, left, right = box
        legend = (left, top, left + int(0.42 * (right - left)), top + int(0.13 * (bottom - top)))
        res["fig4"] = read(im, box, xs, xv, ys, yv, legend, two=False)
        res["fig4"]["calibration"] = {"x_ticks_px": xs, "y_ticks_px": ys, "frame": box}

        im = render(14, tmp)                              # Fig. 7, page 2104
        box = frame(im, int(im.shape[0] * 0.5))
        xs, ys = ticks(im, *box)
        xv = [round(0.1 * i, 1) for i in range(len(xs))]
        yv = [round(1.0 - 0.2 * i, 1) for i in range(len(ys))]
        top, bottom, left, right = box
        legend = (left + int(0.30 * (right - left)), top + int(0.17 * (bottom - top)),
                  right, top + int(0.42 * (bottom - top)))
        res["fig7"] = read(im, box, xs, xv, ys, yv, legend, two=True)
        res["fig7"]["calibration"] = {"x_ticks_px": xs, "y_ticks_px": ys, "frame": box}
    OUT.write_text(json.dumps(res, indent=1))
    for fig, keys in (("fig4", ("phi", "c")), ("fig7", ("phi_high", "c_high", "phi_low", "c_low"))):
        print(fig, "ticks", len(res[fig]["calibration"]["x_ticks_px"]), "x,",
              len(res[fig]["calibration"]["y_ticks_px"]), "y")
        for i in (0, 4, 9, 14, 19, 29, 49, 99):
            print(f"  t={T[i]:.2f} " + "  ".join(f"{k}={res[fig][k][i]}" for k in keys))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
