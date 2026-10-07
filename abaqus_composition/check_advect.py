"""check_advect.py -- a well-posed check of the front term of Eq. 34 in Abaqus.

    python abaqus_composition/check_advect.py JOB.dat [n]

make_klempt_inp.py --case advect: 20 um cube, c held at 1 on the bottom and 0 on
the top, no consumption, so c = 1 - z/20 and n_c = grad c/|grad c| = (0, 0, -1)
everywhere; phi = 1 on a ball of radius 3 um at the centre. The front term then
carries phi down with v = r c/(k+c) (r = 100 um/T*, k = 1), diffusion
beta = 2 spreads it. Compared over time: the height of the phi centroid,
  - Abaqus (element centroid phi, from the .dat),
  - the finite-difference scheme of JAXFEM/klempt2024_quantitative.py (upwind,
    explicit, phi clipped to [0, 1]) on the same grid with the same c,
  - the ODE dz/dt = -v(c(z)) for a point starting at z = 10 um.
"""
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "JAXFEM"))
import klempt2024_quantitative as k  # noqa: E402

R, KM = 100.0, 1.0


def abaqus(dat, n=20):
    t, zc, cur, rows = [], [], None, {}
    tab_t = [None]                                   # time of the table being read

    def flush():
        cur_t = tab_t[0]
        if rows and cur_t is not None:
            e = np.array(sorted(rows))
            ph = np.array([rows[i] for i in e])
            z = ((e - 1) // n // n + 0.5) * (20.0 / n)
            t.append(cur_t); zc.append(float((ph * z).sum() / ph.sum()))

    for l in Path(dat).read_text(errors="replace").splitlines():
        m = re.search(r"STEP TIME COMPLETED\s+([0-9.Ee+-]+)", l)
        if m:
            cur = float(m.group(1).rstrip(","))
            continue
        if "ELEMENT  FOOT-" in l:
            flush(); rows.clear(); tab_t[0] = cur      # the increment summary comes before its table
            continue
        m = re.match(r"^\s+(\d+)\s+([-+]?\d+\.\d*(?:E[-+]\d+)?)\s+([-+]?\d+\.\d*(?:E[-+]\d+)?)\s*$", l)
        if m:
            rows[int(m.group(1))] = float(m.group(2))
    flush()
    return np.array(t), np.array(zc)


def fd(t_end, dt=1e-3):
    phi = ((k.X - 10) ** 2 + (k.Y - 10) ** 2 + (k.Z - 10) ** 2 <= 9.0 + 1e-9).astype(float)
    c = 1.0 - k.Z / 20.0
    v = (np.zeros_like(c), np.zeros_like(c), -R * c / (KM + c))
    alpha = np.ones_like(phi)
    t, zc = [0.0], [float((k.W * phi * k.Z).sum() / (k.W * phi).sum())]
    for s in range(1, int(round(t_end / dt)) + 1):
        phi = phi + dt * (k.BETA * k.lap(phi) + k.K_A * alpha - k.upwind_dot(phi, v))
        phi = np.clip(phi, 0.0, 1.0)
        alpha = alpha + dt * k.K_A * phi
        t.append(s * dt); zc.append(float((k.W * phi * k.Z).sum() / (k.W * phi).sum()))
    return np.array(t), np.array(zc)


def ode(ts, z0=10.0, dt=1e-5):
    z, out, s = z0, [], 0.0
    for tt in ts:
        while s < tt - 1e-12:
            c = 1 - z / 20
            z -= dt * R * c / (KM + c)
            s += dt
        out.append(z)
    return np.array(out)


def main(dat, n=20):
    ta, za = abaqus(dat, n)
    tf, zf = fd(float(ta[-1]))
    zo = ode(ta)
    print(f"{'T*':>6s} {'z Abaqus':>9s} {'z FD':>8s} {'z ODE':>8s}   (centroid height, um; start 10)")
    for t, z, o in zip(ta, za, zo):
        i = int(np.argmin(np.abs(tf - t)))
        print(f"{t:6.3f} {z:9.3f} {zf[i]:8.3f} {o:8.3f}")


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 20)
