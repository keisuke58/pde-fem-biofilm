"""compare_blend.py -- Klempt 2024's test cases 4.1 / 4.2 in Abaqus with the
closest-reproduction model of KLEMPT2024_REPRODUCTION.md sec. 15 (growth on
every face with share w, consumption g phi c, a time scale s per run), next
to the finite-difference reproduction with the same settings
(JAXFEM/klempt2024_case1_bc.run_setup, growth "blend<w>") and the paper's
digitised curves.

    python abaqus_composition/compare_blend.py JOB.dat --case fig4_corner|fig7_high|fig7_low --w 0.5 --scale s

Prints mean phi and c at the paper's times and the RMS differences
(phi and c averaged, as klempt2024_variant_search.score).
"""
import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "JAXFEM"), str(Path(__file__).resolve().parent)]
import klempt2024_case1_bc as B  # noqa: E402
import klempt2024_case2 as C  # noqa: E402
import klempt2024_quantitative as K  # noqa: E402
import klempt2024_variant_search as VS  # noqa: E402
from compare_klempt import series  # noqa: E402


def reference(case, w, s):
    K.R, K.BETA, K.K_A = 100 * s, 2 * s, 1e-3 * s
    if case == "fig4_corner":
        p, _, _ = K.setup("fig4_edge")
        mask = (K.X >= K.L - 2 - 1e-9) & (K.Y >= K.L - 2 - 1e-9) & (K.Z >= K.L - 2 - 1e-9)
        rec, _, _ = B.run_setup(p, mask, 1e8, "ic", "first_order", f"blend{w}")
    else:
        p, m = C.setup(1)
        rec, _, _ = B.run_setup(p, m, C.G[case.split("_")[1]], "ic", "first_order", f"blend{w}")
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dat")
    ap.add_argument("--case", required=True, choices=("fig4_corner", "fig7_high", "fig7_low"))
    ap.add_argument("--w", type=float, default=0.5)
    ap.add_argument("--scale", type=float, default=1.0)
    a = ap.parse_args()
    t, phi, c = series(a.dat)
    ab = {"t": list(t), "phi": list(phi), "c": list(c)}
    fd = reference(a.case, a.w, a.scale)
    key, ts = ("4.1", VS.T41) if a.case == "fig4_corner" else ("4.2 " + a.case.split("_")[1], C.T_CMP)
    paper = VS.digitized()[key]
    ts = [x for x in ts if x <= t[-1] + 1e-9]
    at = lambda rec, k: np.interp(ts, rec["t"], rec[k])  # noqa: E731
    print(f"{key}, w = {a.w}, time scale {a.scale}; mean phi / mean c")
    print(f"{'T*':>5s} {'Abaqus':>15s} {'FD':>15s} {'paper':>15s}")
    for i, x in enumerate(ts):
        pp, pc = paper["phi"][i], paper["c"][i]
        print(f"{x:5.2f} {at(ab, 'phi')[i]:7.3f} {at(ab, 'c')[i]:7.3f} {at(fd, 'phi')[i]:7.3f} {at(fd, 'c')[i]:7.3f} "
              f"{pp:7.3f} {pc:7.3f}")
    n = len(ts)
    rms = lambda x, y: float(np.sqrt(np.nanmean((np.asarray(x) - np.asarray(y)) ** 2)))  # noqa: E731
    for name, (x1, x2) in {"Abaqus - paper": (ab, None), "FD - paper": (fd, None), "Abaqus - FD": (ab, fd)}.items():
        r = []
        for k in ("phi", "c"):
            y = at(x2, k) if x2 is not None else np.asarray(paper[k][:n])
            r.append(rms(at(x1, k), y))
        print(f"RMS {name:15s} phi {r[0]:.3f}  c {r[1]:.3f}  mean {0.5 * (r[0] + r[1]):.3f}")


if __name__ == "__main__":
    main()
