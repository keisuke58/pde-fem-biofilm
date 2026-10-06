"""compare_paper.py -- an Abaqus run of make_klempt_inp.py against Klempt 2024's
digitised curves (Fig. 4: test case 4.1; Fig. 7: 4.2 "high" / "low") only, for
runs without a finite-difference counterpart (e.g. the set-up Felix Klempt
described on 5 Oct 2026: Eq. 34 as printed, consumption g phi, T* = 0..1).

    python abaqus_composition/compare_paper.py JOB.dat --case fig4_corner|fig7_high|fig7_low

Prints domain mean phi and c (centroid means, as compare_klempt.py), the largest
and smallest element phi, at the paper's times, and the RMS differences to the
digitised curves (phi and c, and their mean, as klempt2024_variant_search.score).
"""
import argparse
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "JAXFEM"))
import klempt2024_case2 as C  # noqa: E402
import klempt2024_variant_search as VS  # noqa: E402


def tables(dat):
    """(t, element phi array, element c array) for every printed increment."""
    out, cur, rows, tab_t = [], None, {}, None
    for l in Path(dat).read_text(errors="replace").splitlines():
        m = re.search(r"STEP TIME COMPLETED\s+([0-9.Ee+-]+)", l)
        if m:
            cur = float(m.group(1).rstrip(","))
            continue
        if "ELEMENT  FOOT-" in l:
            if rows and tab_t is not None:
                v = np.array(list(rows.values()))
                out.append((tab_t, v[:, 0], v[:, 1]))
            rows, tab_t = {}, cur
            continue
        m = re.match(r"^\s+(\d+)\s+((?:[-+]?\d+\.\d*(?:E[-+]\d+)?\s*){2})$", l)
        if m:
            rows[int(m.group(1))] = [float(x) for x in m.group(2).split()]
    if rows and tab_t is not None:
        v = np.array(list(rows.values()))
        out.append((tab_t, v[:, 0], v[:, 1]))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dat")
    ap.add_argument("--case", required=True, choices=("fig4_corner", "fig7_high", "fig7_low"))
    a = ap.parse_args()
    tb = tables(a.dat)
    t = np.array([x[0] for x in tb])
    phi = np.array([x[1].mean() for x in tb])
    c = np.array([x[2].mean() for x in tb])
    pmax = np.array([x[1].max() for x in tb])
    pmin = np.array([x[1].min() for x in tb])
    cmin = np.array([x[2].min() for x in tb])
    key, ts = ("4.1", VS.T41) if a.case == "fig4_corner" else ("4.2 " + a.case.split("_")[1], C.T_CMP)
    paper = VS.digitized()[key]
    n_all = len(ts)
    ts = [x for x in ts if x <= t[-1] + 1e-9]
    at = lambda y: np.interp(ts, t, y)  # noqa: E731
    print(f"{key}: Abaqus (mean phi, mean c, max / min element phi, min element c) against the paper (mean phi, mean c)")
    print(f"{'T*':>5s} {'phi':>7s} {'c':>7s} {'max phi':>8s} {'min phi':>8s} {'min c':>8s} {'paper phi':>10s} {'paper c':>8s}")
    for i, x in enumerate(ts):
        print(f"{x:5.2f} {at(phi)[i]:7.3f} {at(c)[i]:7.3f} {at(pmax)[i]:8.3f} {at(pmin)[i]:8.3f} {at(cmin)[i]:8.3f} "
              f"{paper['phi'][i]:10.3f} {paper['c'][i]:8.3f}")
    n = len(ts)
    rms = lambda x, y: float(np.sqrt(np.nanmean((np.asarray(x) - np.asarray(y)) ** 2)))  # noqa: E731
    rp, rc = rms(at(phi), paper["phi"][:n]), rms(at(c), paper["c"][:n])
    print(f"RMS Abaqus - paper: phi {rp:.3f}  c {rc:.3f}  mean {0.5 * (rp + rc):.3f}"
          + ("" if n == n_all else f"  (only T* <= {t[-1]:.3f}: {n} of {n_all} points)"))


if __name__ == "__main__":
    main()
