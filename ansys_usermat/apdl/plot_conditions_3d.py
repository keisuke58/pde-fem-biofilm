#!/usr/bin/env python3
"""The four-condition shell in space, coloured by the stress ANSYS returned.

Companion to plot_condition_comparison.py: it puts the same four means on the
geometry they were measured on, so the 2x2 arrangement is visible as the thing
it physically is -- four quadrants of a single curved tile.

Geometry is not read from a listing; the run printed stresses and state, not
nodal coordinates. It is reconstructed from the deck's own parameters, which
define it exactly:

    CYLIND,4.0,4.3,0,0.15,0,5      substrate, does not grow
    CYLIND,4.3,4.4,0,0.15,0,5      growth layer, split 2x2 at ARC/2 and LEN/2

so the drawing is the deck's geometry at true proportions, with no axis scaled.
What is measured -- the mean von Mises of each quadrant -- is what the colour
carries.

The two high-arc quadrants come out warm here. An earlier version of this file
read that as the arc; plot_condition_paired.py shows it follows the conditions
that happen to sit there, the arc itself being worth about a tenth of a
percent.

    python ansys_usermat/apdl/plot_conditions_3d.py -o assets/v222_conditions_3d.png
"""
import argparse
import statistics as st
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_condition_comparison import LAYOUT, parse  # noqa: E402

R_IN, R_MID, R_OUT = 4.0, 4.3, 4.4
LEN, ARC = 0.15, 5.0          # Z extent, and arc in degrees


def _solid(ax, r0, r1, th0, th1, z0, z1, **kw):
    """A curved hexahedral block, all six faces.

    Drawing only the outer surface makes the two layers read as two detached
    plates: at true proportions the 0.1 radial gap is a quarter of the tile's
    width, so the separation is real and has to be closed by actually drawing
    the sides.
    """
    import numpy as np
    t = np.radians(np.linspace(th0, th1, 40))
    r = np.linspace(r0, r1, 2)
    z = np.linspace(z0, z1, 2)

    def surf(X, Y, Z):
        ax.plot_surface(X, Y, Z, shade=False, edgecolor="none", **kw)

    for rr in (r0, r1):                      # inner / outer
        T, Z = np.meshgrid(t, z)
        surf(rr * np.cos(T), rr * np.sin(T), Z)
    for tt in (np.radians(th0), np.radians(th1)):   # the two arc ends
        R, Z = np.meshgrid(r, z)
        surf(R * np.cos(tt), R * np.sin(tt), Z)
    for zz in (z0, z1):                      # bottom / top
        T, R = np.meshgrid(t, r)
        surf(R * np.cos(T), R * np.sin(T), np.full_like(T, zz))


def _trim(path, pad=12):
    """Crop uniform white margins, leaving a small pad."""
    from PIL import Image, ImageChops
    im = Image.open(path).convert("RGB")
    bg = Image.new("RGB", im.size, (255, 255, 255))
    box = ImageChops.difference(im, bg).getbbox()
    if not box:
        return
    l, t, r, b = box
    im.crop((max(l - pad, 0), max(t - pad, 0),
             min(r + pad, im.width), min(b + pad, im.height))).save(path)


def plot(mean, out, caption=True):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    # Without the caption the canvas can be narrower, which removes the dead
    # margin the caption's own width was forcing. Cropping the wide version
    # instead would cut the caption text, which spans the full width.
    fig = plt.figure(figsize=(10.6, 4.8) if caption else (7.4, 3.1))
    ax = fig.add_subplot(111, projection="3d")

    lo, hi = min(mean.values()), max(mean.values())
    cmap = plt.get_cmap("YlOrRd")
    labels = []

    # substrate as a solid block, so the tile reads as two bonded layers
    _solid(ax, R_IN, R_MID, 0, ARC, 0, LEN, color="0.78")

    for m, (code, full, _short, th_hi, z_hi) in LAYOUT.items():
        if m not in mean:
            continue
        t0, t1 = (ARC / 2, ARC) if th_hi else (0, ARC / 2)
        z0, z1 = (LEN / 2, LEN) if z_hi else (0, LEN / 2)
        f = (mean[m] - lo) / (hi - lo) if hi > lo else 0.5
        _solid(ax, R_MID, R_OUT, t0, t1, z0, z1,
               color=cmap(0.25 + 0.55 * f))
        tm = np.radians((t0 + t1) / 2)
        labels.append((R_OUT * 1.002 * np.cos(tm), R_OUT * 1.002 * np.sin(tm),
                       (z0 + z1) / 2, code, mean[m]))

    # True proportions come from the data's own extents. The tile is a 5 deg
    # arc at r ~ 4.3, so it spans only ~0.017 in x against 0.383 in y: setting
    # the box aspect from anything but the real ranges squashes it into a
    # sliver, which is what a hand-picked aspect did on the first attempt.
    th = np.radians(np.linspace(0, ARC, 64))
    xs = R_OUT * np.cos(th)
    ys = R_OUT * np.sin(th)
    dx = max(xs.max() - R_IN * np.cos(th).min(), 1e-9)
    dy, dz = ys.max() - ys.min(), LEN
    ax.set_xlim(R_IN * np.cos(th).min(), xs.max())
    ax.set_ylim(ys.min(), ys.max())
    ax.set_zlim(0, LEN)
    ax.set_box_aspect((dx, dy, dz))
    # Look along -x, i.e. radially inward at the outer face: that is the face
    # the four quadrants are on. Viewing down the arc end instead shows mostly
    # substrate and hides the thing the figure is about.
    ax.view_init(elev=11, azim=6)
    ax.set_axis_off()

    sm = plt.cm.ScalarMappable(cmap=cmap,
                               norm=plt.Normalize(lo * 1e6, hi * 1e6))
    cb = fig.colorbar(sm, ax=ax, fraction=0.026, pad=0.10)
    cb.set_label(r"mean von Mises per region  [$\times 10^{-6}$]", fontsize=9)

    # The 3D box is very flat in x, so matplotlib centres a bounding box whose
    # visible content sits well to one side, leaving a wide dead margin. Force
    # the axes to fill the canvas; do it after the colorbar, which shrinks it.
    ax.set_position([-0.06, 0.02 if caption else 0.0,
                     0.92, 0.90 if caption else 1.0])

    # mplot3d does not depth-sort text against surfaces, so a 3D label sits
    # behind whichever face it lands on and zorder does not help. Project each
    # quadrant centre to display coordinates and draw the label on the figure
    # instead, where it is unconditionally on top.
    from mpl_toolkits.mplot3d import proj3d
    fig.canvas.draw()
    for x, y, z, code, v in labels:
        px, py, _ = proj3d.proj_transform(x, y, z, ax.get_proj())
        fx, fy = fig.transFigure.inverted().transform(
            ax.transData.transform((px, py)))
        fig.text(fx, fy, f"{code}\n{v:.3e}", ha="center", va="center",
                 fontsize=9 if caption else 8, weight="bold", color="0.12",
                 bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none",
                           alpha=0.75))


    if caption:
        # On a slide the frame title already says this.
        fig.suptitle("The four clinical conditions are four quadrants of "
                     "one tile", fontsize=11.5, y=0.95)
    if caption:
        fig.text(0.5, 0.045,
                 f"Growth layer, r = {R_MID}--{R_OUT}, over the non-growing "
                 f"substrate (grey, r = {R_IN}--{R_MID}); arc {ARC}°, "
                 f"Z = 0--{LEN}. True proportions, no axis scaled.\n"
                 "Colour is measured; geometry is reconstructed from the deck's own "
                 "CYLIND parameters, since the run printed stress and state but no "
                 "nodal coordinates.\n"
                 "The two high-arc quadrants are the warm ones in this run; the "
                 "paired runs show that follows the conditions placed\nthere, "
                 "not the arc.",
                 ha="center", va="bottom", fontsize=8, color="0.35")
    fig.savefig(out, dpi=200)
    if not caption:
        # Safe only without the caption: the caption spans the full width, so
        # cropping the wide version would cut its text. Here the dead margin
        # is genuinely dead -- the flat 3D box leaves it and no axes position
        # reclaims it.
        _trim(out)
    print(f"wrote {out}")


def main(argv=None):
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser()
    ap.add_argument("--listing", type=Path, default=here /
                    "growth_result_cylinder_ecology_4region_all_real_clsm.txt")
    ap.add_argument("-o", "--out", type=Path,
                    default=Path("assets/v222_conditions_3d.png"))
    ap.add_argument("--no-caption", action="store_true",
                    help="omit the explanatory block; for slides, where the "
                         "speaker says it and the margin is worth more")
    a = ap.parse_args(argv)

    seqv, mats = parse(a.listing)
    mean = {m: st.mean(seqv[e] for e in mats[m]) for m in mats}
    for m in sorted(mean):
        t = "high" if LAYOUT[m][3] else "low "
        z = "high" if LAYOUT[m][4] else "low "
        print(f"  {LAYOUT[m][0]}  arc {t}  Z {z}  mean SEQV {mean[m]:.4e}")
    a.out.parent.mkdir(parents=True, exist_ok=True)
    plot(mean, a.out, caption=not a.no_caption)
    return 0


if __name__ == "__main__":
    sys.exit(main())
