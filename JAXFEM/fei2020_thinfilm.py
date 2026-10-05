#!/usr/bin/env python3
"""Expansion of a V. cholerae biofilm on agar (Fei et al. 2020, PNAS 117:7622),
re-implemented in JAX, and the same mechanics with the growth law of this work.

Fei et al. model a flat, growing biofilm as an axisymmetric plane-stress film
on a substrate with viscous friction (their Eq. 6, before wrinkling), solved in
the Lagrangian frame on R in [0, 1] (initial radius R_b0 = 2 mm, time in units
of tau0 = 1/k_g at the rim). Unknowns at the nodes: radial displacement u_r,
in-plane growth stretch lambda, nutrient c. Their FEniCS script
(github.com/f-chenyi/biofilm-mechanics-theory, biofilm_morphogenesis_circle_
CLEANED.py) is followed term by term, with wrinkling off:

  friction  xi (u - u0) r F_rr w + dt (gamma s_tt F_rr w + gamma s_rr w' r) = 0
  growth    lambda - lambda0 = dt k(c) lambda          (their law, exponential)
  nutrient  (c - c0) F_rr w - c' (u - u0) w + dt D c' w'/F_rr
            + dt q(c) F_rr w / J_a = 0                  (all times r)
  s_rr = A_rr^2 - gamma^2, s_tt = A_tt^2 - gamma^2, A = F/lambda,
  gamma = 1/(A_rr A_tt) (incompressible plane stress), r = R + u.
  k(c) = 1.5 (1 - k_r) c/(c + 1/2) + k_r, q(c) = Q0 c/(c + 1/2),
  xi = 36, D = 1/2, Q0 = 2 (a_c = 0.5), k_r = 0.15 as in their script.

Growth laws compared (only the growth line changes):
  fei       d lambda/dt = k(c) lambda            (Fei et al.)
  eq36      d alpha/dt  = k_alpha phi, phi = 1   (Klempt 2024 Eq. 36; no nutrient)
  eq36c     d alpha/dt  = k_alpha phi c/(c+1/2)  (Eq. 36 limited by the nutrient)
For eq36 and eq36c the growth acts in plane only, as in Fei et al., so that
only the growth law differs. Each law has one rate (k_g or k_alpha) set by the
time unit tau0, fitted to the measured velocity profiles at 1, 6 and 16 h on
0.7 % agar (their Fig. 2D; data from the repository above, not bundled here).

    git clone --depth 1 https://github.com/f-chenyi/biofilm-mechanics-theory /tmp/fei
    python JAXFEM/fei2020_thinfilm.py --data /tmp/fei
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import jax
import jax.numpy as jnp
import numpy as np

jax.config.update("jax_enable_x64", True)

XI, DN, Q0, KR, CHALF = 36.0, 0.5, 2.0, 0.15, 0.5
RB0_MM = 2.0
HOURS = (1, 6, 16)


def mesh(n=160):
    s = np.linspace(0.0, 1.0, n + 1)
    return jnp.asarray(1.0 - (1.0 - s) ** 1.6)          # finer near R = 1


def make_step(R, law, dt, xi=XI):
    ne = R.size - 1
    g = jnp.array([-1.0, 1.0]) / jnp.sqrt(3.0)
    wq = jnp.array([1.0, 1.0])
    Rl, Rr = R[:-1], R[1:]
    h = Rr - Rl
    # shape functions at the two Gauss points of every element
    N0 = (1 - g) / 2
    N1 = (1 + g) / 2
    Rq = Rl[:, None] * N0 + Rr[:, None] * N1                  # (ne, 2)
    dN = jnp.stack([-1.0 / h, 1.0 / h], 1)                     # (ne, 2)
    jac = h[:, None] / 2 * wq                                  # dR weights
    n = R.size

    def interp(v):
        return v[:-1, None] * N0 + v[1:, None] * N1, (v[1:] - v[:-1])[:, None] / h[:, None]

    phi = lambda c: c / (c + CHALF)

    def rate(c):
        if law == "fei":
            return 1.5 * (1 - KR) * phi(c) + KR
        if law == "eq36":
            return jnp.ones_like(c)
        return 1.5 * phi(c)                                    # eq36c: 1 at c = 1

    def residual(x, x0):
        u, lam, c = x[:n], x[n:2 * n], x[2 * n:]
        u0, lam0, c0 = x0[:n], x0[n:2 * n], x0[2 * n:]
        uq, du = interp(u)
        lq, _ = interp(lam)
        cq, dc = interp(c)
        u0q, _ = interp(u0)
        c0q, _ = interp(c0)
        Rc = jnp.maximum(Rq, 1e-12)
        Frr = 1 + du
        Ftt = 1 + uq / Rc
        Arr, Att = Frr / lq, Ftt / lq
        Ja = Arr * Att
        gam = 1.0 / Ja
        srr = Arr ** 2 - gam ** 2
        stt = Att ** 2 - gam ** 2
        r = Rq + uq
        # element vectors for the two local nodes, (ne, 2)
        Nloc = jnp.stack([jnp.broadcast_to(N0, uq.shape), jnp.broadcast_to(N1, uq.shape)], -1)
        dNloc = dN[:, None, :]
        f1 = (xi * (uq - u0q) * r * Frr)[..., None] * Nloc + dt * (
            (gam * stt * Frr)[..., None] * Nloc + (gam * srr * r)[..., None] * dNloc)
        f3 = ((cq - c0q) * Frr - dc * (uq - u0q) + dt * Q0 * phi(cq) * Frr / Ja)[..., None] * Nloc * r[..., None] \
            + (dt * DN * dc / Frr * r)[..., None] * dNloc
        r1 = (f1 * jac[..., None]).sum(1)
        r3 = (f3 * jac[..., None]).sum(1)
        R1 = jnp.zeros(n).at[:-1].add(r1[:, 0]).at[1:].add(r1[:, 1])
        R3 = jnp.zeros(n).at[:-1].add(r3[:, 0]).at[1:].add(r3[:, 1])
        R1 = R1.at[0].set(u[0])                                # u(0) = 0
        R3 = R3.at[-1].set(c[-1] - 1.0)                        # c(1) = 1
        # growth, pointwise at the nodes (backward Euler)
        if law == "fei":
            R2 = lam - lam0 - dt * rate(c) * lam
        else:
            R2 = lam - lam0 - dt * rate(c)
        return jnp.concatenate([R1, R2, R3])

    jacf = jax.jit(jax.jacfwd(residual))
    resf = jax.jit(residual)

    def step(x0):
        x = x0
        for _ in range(30):
            r = resf(x, x0)
            dx = jnp.linalg.solve(jacf(x, x0), -r)
            x = x + dx
            if float(jnp.max(jnp.abs(dx))) < 1e-11:
                return x
        raise RuntimeError("Newton did not converge")
    return step


def simulate(law, T, dt=0.01, n=160, xi=XI):
    """Return times (tau), current radius of every node and velocity, per output."""
    R = mesh(n)
    nn = R.size
    x = jnp.concatenate([jnp.zeros(nn), jnp.ones(nn), jnp.ones(nn)])
    step = make_step(R, law, dt, xi)
    res = []
    steps = int(round(T / dt))
    for k in range(1, steps + 1):
        x0 = x
        x = step(x0)
        u, u0 = np.asarray(x[:nn]), np.asarray(x0[:nn])
        res.append((k * dt, np.asarray(R) + u, (u - u0) / dt))
    return res


def load_data(root: Path):
    import scipy.io as sio
    d = {}
    for hr in HOURS:
        m = sio.loadmat(root / "paper data" / "Fig 2" / "velocity profile" / f"velocity-profile-0.7-{hr}h.mat")
        d[hr] = (m["r"].ravel(), m["v_avg"].ravel(), m["v_std"].ravel())   # mm, um/min
    k = sio.loadmat(root / "paper data" / "Fig 1" / "kymograph" / "kymo-0.7.mat")
    t, rb = k["tAll"].ravel(), k["Rout"].ravel() / 1000.0                     # h, mm
    d["radius"] = {hr: float(rb[np.argmin(np.abs(t - hr))]) for hr in HOURS}
    return d


def model_profiles(sim, tau0_h):
    """Velocity profiles in mm and um/min at the measured hours for time unit tau0 (h)."""
    t = np.array([s[0] for s in sim])
    out = {}
    for hr in HOURS:
        k = int(np.argmin(np.abs(t - hr / tau0_h)))
        r = sim[k][1] * RB0_MM
        v = sim[k][2] * RB0_MM * 1000 / (tau0_h * 60)          # mm/tau -> um/min
        out[hr] = (r, v, t[k] * tau0_h)
    return out


def misfit(sim, data, tau0_h):
    if max(HOURS) / tau0_h > sim[-1][0] + 1e-9:
        return np.inf                                         # beyond the simulated time
    prof = model_profiles(sim, tau0_h)
    e = []
    for hr in HOURS:
        r, v, _ = prof[hr]
        rd, vd, sd = data[hr]
        if rd.max() > r.max() * 1.02:
            return np.inf                                     # model biofilm too small
        vm = np.interp(rd, r, v)
        e.append(((vm - vd) / 3.0) ** 2)                      # V0 = 3 um/min as in their MSD
        e.append(np.array([((r.max() - data["radius"][hr]) / 5.0) ** 2]))   # L0 = 5 mm, radius
    return float(np.mean(np.concatenate(e)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True, help="clone of f-chenyi/biofilm-mechanics-theory")
    ap.add_argument("--laws", nargs="+", default=["fei", "eq36", "eq36c"])
    ap.add_argument("--T", type=float, default=6.0)
    ap.add_argument("--xi", type=float, nargs="+", default=[XI], help="friction values to scan")
    ap.add_argument("--out", default=None, help="figure path")
    args = ap.parse_args()
    data = load_data(Path(args.data))
    fits = {}
    taus = np.linspace(2.0, 30.0, 281)
    for law in args.laws:
        best = None
        for xi in args.xi:
            sim = simulate(law, args.T, xi=xi)
            m = [misfit(sim, data, t) for t in taus]
            i = int(np.argmin(m))
            print(f"{law:6s} xi = {xi:5.1f}  tau0 = {taus[i]:5.2f} h  normalised MSD = {m[i]:.4f}  "
                  f"radius at 16 h = {model_profiles(sim, taus[i])[16][0].max():.2f} mm "
                  f"(measured {data['radius'][16]:.2f})", flush=True)
            if best is None or m[i] < best[2]:
                best = (sim, taus[i], m[i], xi)
        fits[law] = best
    if args.out:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "ansys_usermat"))
        import figstyle
        figstyle.apply(11)
        fig, ax = plt.subplots(1, len(fits), figsize=(4.2 * len(fits), 3.6), sharey=True)
        ax = np.atleast_1d(ax)
        names = {"fei": "Fei et al. law, $\\dot\\lambda=k(c)\\lambda$",
                 "eq36": "Eq. 36, $\\dot\\alpha=k_\\alpha\\phi$",
                 "eq36c": "Eq. 36 with nutrient, $\\dot\\alpha=k_\\alpha\\phi\\,c/(c+K)$"}
        for a, (law, (sim, tau0, ms, xi)) in zip(ax, fits.items()):
            prof = model_profiles(sim, tau0)
            for hr, col in zip(HOURS, ("C0", "C1", "C3")):
                rd, vd, sd = data[hr]
                a.errorbar(rd, vd, sd, fmt="o", ms=3, color=col, alpha=0.6, label=f"{hr} h, measured")
                r, v, _ = prof[hr]
                a.plot(r, v, "-", color=col, lw=1.6)
            a.set_title(f"{names[law]}\n$\\xi$ = {xi:g}, $\\tau_0$ = {tau0:.1f} h, MSD = {ms:.3f}", fontsize=10)
            a.set_xlabel("$r$ [mm]")
        ax[0].set_ylabel("radial velocity [$\\mu$m/min]")
        ax[0].legend(frameon=False, fontsize=8)
        fig.tight_layout()
        fig.savefig(args.out, dpi=200)
        print("wrote", args.out)


if __name__ == "__main__":
    main()
