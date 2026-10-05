"""summarize_implant.py -- depth profiles of a make_implant_inp.py run.

    python abaqus_composition/summarize_implant.py JOB.dat JOB.json

Last increment, element centroids, averaged over theta, for the ring next to
the titanium (inner) and the outer ring: phi, the share phi_1/(phi_1+phi_2)
(where the point model ran), c, alpha - 1, and the stress in cylindrical
components (Pa): sigma_rr and sigma_rz on the titanium are the normal and
shear load on the bond, sigma_tt and sigma_zz the in-layer stresses.
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_ansys import last_tables  # noqa: E402

MPA_TO_PA = 1e6


def main(dat, js):
    m = json.loads(Path(js).read_text())
    t = last_tables(dat)
    S = next(x for x in t if "S11" in x["head"])
    V = next(x for x in t if "SDV84" in x["head"])
    el = np.array(sorted(V["rows"]))
    cen = np.array([m["cen"][str(e)] for e in el])
    th, ir, kz = cen[:, 1], cen[:, 3].astype(int), cen[:, 5].astype(int)
    col = lambda tab, name: np.array([tab["rows"][e][tab["head"].index(name)] for e in el])  # noqa: E731
    sxx, syy, szz, sxy, sxz, syz = (col(S, h) * MPA_TO_PA for h in ("S11", "S22", "S33", "S12", "S13", "S23"))
    c_, s_ = np.cos(th), np.sin(th)
    srr = c_**2 * sxx + s_**2 * syy + 2 * c_ * s_ * sxy
    stt = s_**2 * sxx + c_**2 * syy - 2 * c_ * s_ * sxy
    srz = c_ * sxz + s_ * syz
    vm = np.sqrt(0.5 * ((sxx - syy)**2 + (syy - szz)**2 + (szz - sxx)**2) + 3 * (sxy**2 + sxz**2 + syz**2))
    phi, a1, cn = col(V, "TEMP"), col(V, "SDV84"), col(V, "SDV51")
    p1, p2 = col(V, "SDV72"), col(V, "SDV73")
    s = p1 + p2
    ran = (s > 1e-6) & (s < 0.95)
    share = np.where(ran, p1 / np.where(ran, s, 1), np.nan)
    nz, nr, h = m["nz"], m["nr"], m["height"]
    spread = {"theta spread of sigma_tt (inner, max over z) %":
              max(np.ptp(stt[(ir == 0) & (kz == k)]) / max(abs(stt[(ir == 0) & (kz == k)].mean()), 1e-30) * 100
                  for k in range(nz))}
    for ring, name in ((0, "inner (on the titanium)"), (nr - 1, "outer")):
        print(f"\n{name} ring, theta means; z from the base (0) to the gingival margin ({h})")
        print(f"{'z mm':>6s} {'phi':>6s} {'share':>6s} {'c':>6s} {'alpha-1':>9s} "
              f"{'s_rr':>10s} {'s_rz':>10s} {'s_tt':>10s} {'s_zz':>10s} {'vM':>10s}")
        for k in range(nz):
            q = (ir == ring) & (kz == k)
            print(f"{(k + 0.5) * h / nz:6.2f} {phi[q].mean():6.3f} {np.nanmean(share[q]) if ran[q].any() else float('nan'):6.3f} "
                  f"{cn[q].mean():6.3f} {a1[q].mean():9.2e} {srr[q].mean():10.3e} {srz[q].mean():10.3e} "
                  f"{stt[q].mean():10.3e} {szz[q].mean():10.3e} {vm[q].mean():10.3e}")
    for k_, v in spread.items():
        print(f"\n{k_}: {v:.2e}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
