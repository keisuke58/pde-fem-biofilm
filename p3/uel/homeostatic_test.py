"""homeostatic_test.py -- the homeostatic-pressure growth law in the unified user element (unified_uel.f90,
PROPS(19) = p_ref, PROPS(20) = form), checked on one element against closed forms written without the
element's source, in the limit of Klempt et al. 2024 (klempt2024_limit_test.py: one species, living fraction
frozen near 1, c = 1 and phi = 1 fixed, lambda = 0, so the growth reaction is Eq. 36, alpha_dot = k_alpha phi).

The law (keio_wp2/README.ja.md, stages 2B and 4-7; an assumption of this work, not from a paper):
    alpha_dot = k_alpha phi g,   g = min(1, max(0, 1 - p / p_h)),   p = -tr(sigma)/3 (positive in compression),
with p_h = p_ref (form 0) or p_ref (phi^2 + f) (form 1); form 2 is g = 1/(1 + max(p, 0)/p_ref) (the ANSYS cal5gp
build). The element takes p of the start of the increment (explicit gate; see the module header), or, with
PROPS(20) = form + 10 (`implicit=True` here), p of the end of the increment (the implicit gate: a local solve at every
Gauss point and the consistent tangent); the checks (b)-(d) are repeated for it: the constrained cube is then backward
Euler on the ODE (first order as well), p reaches p_h from below without overshoot, and the tangent includes the
implicit function's derivative d xi / dU.

  (a) free growth (supports 3-2-1): p = 0, so g = 1 and alpha = 1 + k_alpha t as Eq. 36, with any p_ref;
  (b) fully constrained (all u = 0): the hyperelastic law gives sigma = (phi^2 + f) K (alpha^-3 - 1), i.e.
      p(alpha) = (phi^2 + f) K (1 - alpha^-3), and the law is the ODE alpha_dot = k_alpha g(p(alpha)), solved
      here with scipy (rtol 1e-12). The element's explicit gate is forward Euler on this ODE: the element
      reproduces that Euler sequence to round-off, and its error against the ODE halves with dt (first order);
  (c) homeostatic limit (form 0, long run): growth stops where p reaches p_h, alpha_h = (1 - p_h/((phi^2+f)K))^(-1/3);
      the explicit gate lets p overshoot p_h by at most one increment of growth, dt k_alpha dp/dalpha;
  (d) the element tangent with the gate active (0 < g < 1) against central differences of the residual.
E = 10 Pa, nu = 0.49, f = 1e-3 (Klempt et al. 2024, Table 2).

    python homeostatic_test.py            (writes the eco file itself when BIOFILM_ECO_CASE is not set)
"""
import os
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
if not os.environ.get("BIOFILM_ECO_CASE"):
    from make_column_inp import eco_file
    _p = Path(tempfile.mkdtemp()) / "uel_eco.txt"; eco_file(_p); os.environ["BIOFILM_ECO_CASE"] = str(_p)
from driver import NDN, NGP, NSV, props_base  # noqa: E402
from emulate_abaqus import NSVT, F, increment, uel  # noqa: E402
from klempt2024_limit_test import cube  # noqa: E402

E, NU, FV = 10.0, 0.49, 1e-3
MU, K = E / (2 * (1 + NU)), E / (3 * (1 - 2 * NU))
LIV = 0.999 * 0.999                       # living fraction of the frozen one-species state; L0 divides it out


def props(k_alpha, p_ref=0.0, form=0, implicit=False):
    """element properties of the Klempt 2024 limit (L0 = k_alpha / (liv ln 10), so the element's rate is k_alpha phi
    as in seed_frozen_reference.py); 20 entries when the law is on (17, 18: antibiotic off); implicit (or form >= 10):
    PROPS(20) = form + 10, the implicit gate"""
    P = props_base(E=E, nu=NU, f=FV, lam=0.0, beta=0.0, L0=k_alpha / (LIV * np.log(10.0)), Y=1.0, mu_b=np.log(0.1),
                   s=1e-14, phi_min=2.0)
    if p_ref > 0:
        P = np.r_[P, 1e30, 0.0, p_ref, float(form + (10 if implicit and form < 10 else 0))]
    return F(P)


def frozen_svars():
    sv = np.zeros((NGP, NSV)); sv[:, 0] = 0.999; sv[:, 5] = 0.001; sv[:, 6] = 0.999; sv[:, 13] = 1.0
    return F(np.r_[sv.ravel(), np.zeros(NSVT - NSV * NGP)])


def setup(support):
    nodes, elems = cube()
    U = np.zeros(NDN * len(nodes)); U[3::NDN] = 1.0; U[4::NDN] = 1.0
    fixed = np.zeros(len(U), bool); fixed[3::NDN] = fixed[4::NDN] = True
    if support == "all":
        for d in range(3):
            fixed[d::NDN] = True
    else:                                                      # 3-2-1
        n0 = np.where((nodes == [0, 0, 0]).all(1))[0][0]; n1 = np.where((nodes == [1, 0, 0]).all(1))[0][0]
        n2 = np.where((nodes == [0, 1, 0]).all(1))[0][0]
        fixed[NDN * n0:NDN * n0 + 3] = True; fixed[NDN * n1 + 1:NDN * n1 + 3] = True; fixed[NDN * n2 + 2] = True
    return nodes, elems, U, fixed


def run(k_alpha, T, steps, support, p_ref=0.0, form=0, md=False, implicit=False, iters=None):
    """one element, `steps` increments of T/steps; returns alpha - 1 and p (Gauss-point means) after every increment;
    iters: a list that receives the Newton iterations of every increment"""
    nodes, elems, U, fixed = setup(support)
    P = props(k_alpha, p_ref, form + (10 if implicit else 0)); SV = [frozen_svars()]
    xi, p = [], []
    dt = T / steps
    for k in range(steps):
        U, SV, it = increment(U.copy(), U, fixed, nodes, elems, SV, P, dt, k + 1, ttot=k * dt, finite=True, md=md)
        xi.append(SV[0][:NSV * NGP].reshape(NGP, NSV)[:, 12].mean()); p.append(SV[0][NSV * NGP:].reshape(NGP, 2)[:, 1].mean())
        if iters is not None:
            iters.append(it)
    return np.array(xi), np.array(p), U, nodes


# ---- references written without the element -------------------------------------------------------------------
def p_closed(alpha, phi=1.0):
    """pressure of the fully constrained cube, F = I, F_e = I / alpha: sigma = (phi^2 + f) K (alpha^-3 - 1)"""
    return (phi ** 2 + FV) * K * (1.0 - alpha ** -3)


def gate(p, p_ref, form, phi=1.0):
    if p_ref <= 0:
        return 1.0
    if form == 2:
        return 1.0 / (1.0 + max(p, 0.0) / p_ref)
    ph = p_ref * (phi ** 2 + FV) if form == 1 else p_ref
    return min(1.0, max(0.0, 1.0 - p / ph))


def ode_constrained(k_alpha, T, p_ref, form, n_out=1):
    """alpha(T) of alpha_dot = k_alpha g(p(alpha)), alpha(0) = 1 (scipy, rtol 1e-12)"""
    from scipy.integrate import solve_ivp
    kk = k_alpha                                               # k_alpha phi with phi = 1
    s = solve_ivp(lambda t, a: [kk * gate(p_closed(a[0]), p_ref, form)], (0.0, T), [1.0], method="DOP853",
                  rtol=1e-12, atol=1e-14, t_eval=np.linspace(0, T, n_out + 1)[1:])
    return s.y[0]


def euler_constrained(k_alpha, T, steps, p_ref, form):
    """the explicit gate as a scalar forward-Euler sequence on the same ODE"""
    kk = k_alpha; a = 1.0; dt = T / steps; out = []
    for _ in range(steps):
        a = a + dt * kk * gate(p_closed(a), p_ref, form); out.append(a)
    return np.array(out)


def backward_euler_constrained(k_alpha, T, steps, p_ref, form, tol=1e-15):
    """the implicit gate as a scalar backward-Euler sequence on the same ODE, a = a_n + dt k g(p(a)), each step solved
    by bisection on the monotone equation (the root lies in [a_n, a_n + dt k])"""
    kk = k_alpha; a = 1.0; dt = T / steps; out = []
    for _ in range(steps):
        lo, hi = a, a + dt * kk
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if mid - a - dt * kk * gate(p_closed(mid), p_ref, form) > 0:
                hi = mid
            else:
                lo = mid
            if hi - lo < tol:
                break
        a = 0.5 * (lo + hi); out.append(a)
    return np.array(out)


def alpha_h(p_ref):
    return (1.0 - p_ref / ((1.0 + FV) * K)) ** (-1.0 / 3.0)


# ---- tangent with the gate active ---------------------------------------------------------------------------------
def tangent_check(p_ref=1.0, form=0, h=1e-6, implicit=False, md=False, small=False):
    """AMATRX of UEL against central differences of -RHS on a distorted element whose start-of-increment pressure puts
    the gate at g = 0.5 (form 0: p_n = p_ref/2); returns (relative difference, g, whether the gate changed the tangent).
    implicit: the implicit gate (p_ref chosen so that 0 < g < 1 at the end of the increment; g is then the mean of the
    Gauss-point gates, read from the element's state); md: mean dilatation; small: small strain"""
    nodes, elems, U, fixed = setup("free")
    P = props(0.5, p_ref, form + (10 if implicit else 0))
    sv = frozen_svars(); sv[NSV * NGP + 1::2] = 0.5 * p_ref       # pressures of the start of the increment
    sv[12::NSV][:] = 0.0
    sv[NSV * NGP:] = np.tile([0.0, 0.5 * p_ref], NGP)
    rng = np.random.default_rng(1)
    Ud = U.copy(); Ud[0::NDN] += 1e-2 * rng.standard_normal(8); Ud[1::NDN] += 1e-2 * rng.standard_normal(8)
    Ud[2::NDN] += 1e-2 * rng.standard_normal(8); Ud[3::NDN] = 0.9 + 0.05 * rng.standard_normal(8); Ud[4::NDN] = 0.8
    DU = (Ud - U)[:, None]
    co = F(nodes[elems[0]].T); dt = 0.05
    kw = dict(finite=not small, md=md)
    R0, A, svo, _ = uel(F(Ud), F(DU), co, sv, P, dt, 2, 0.1, **kw)
    g = 1.0 - 0.5 * p_ref / p_ref if form == 0 else None
    if implicit:                                    # g at the end of the increment: xi - xi_n = dt xidot g
        xi = svo[12::NSV][:8]; sv_n = sv[12::NSV][:8]
        P0 = props(0.5); _, _, sv0, _ = uel(F(Ud), F(DU), co, sv, P0, dt, 2, 0.1, **kw)
        g = float(np.mean((xi - sv_n) / (sv0[12::NSV][:8] - sv_n)))
    Kfd = np.zeros_like(A)
    for m in range(len(Ud)):
        Up = Ud.copy(); Up[m] += h; Um = Ud.copy(); Um[m] -= h
        Rp = uel(F(Up), F((Up - U)[:, None]), co, sv, P, dt, 2, 0.1, **kw)[0]
        Rm = uel(F(Um), F((Um - U)[:, None]), co, sv, P, dt, 2, 0.1, **kw)[0]
        Kfd[:, m] = (Rp - Rm) / (2 * h)
    rel = np.abs(A - Kfd).max() / np.abs(Kfd).max()
    # the same element with the law off: the growth columns (phi, c -> mechanics) must differ
    A0 = uel(F(Ud), F(DU), co, sv, props(0.5), dt, 2, 0.1, **kw)[1]
    changed = np.abs(A - A0).max() / np.abs(A).max()
    return rel, g, changed


if __name__ == "__main__":
    print(f"E = {E}, nu = {NU}, f = {FV}; K = {K:.4f}")
    # (a) free growth: the gate is 1 at p = 0 for every form (Klempt et al. 2024 Table 2 rate; the 1e-6 above
    # k_alpha T is the drift of the frozen point model, as in klempt2024_limit_test.py)
    for form in (0, 1, 2):
        xi, p, U, X = run(1e-3, 1.1, 11, "free", p_ref=1e-3, form=form)
        print(f"(a) free, form {form}: alpha-1 {xi[-1]:.6e} (Eq. 36: {1e-3 * 1.1:.6e}), |p| {np.abs(p).max():.1e}")
    # (b) constrained: Euler sequence to round-off, first order against the ODE. p_h = 0.25 Pa is reached at
    # alpha - 1 = 5.0e-4, half the free growth over T = 1.1; dt k_alpha (dp/dalpha) / p_h = 2.2 / steps
    k, T, p_ref = 1e-3, 1.1, 0.25
    a_ode = ode_constrained(k, T, p_ref, 0)[-1]
    print(f"(b) constrained, form 0, p_ref {p_ref}: ODE alpha-1 {a_ode - 1:.8e}, alpha_h-1 {alpha_h(p_ref) - 1:.8e}")
    errs = []
    for steps in (11, 22, 44, 88):
        xi, p, _, _ = run(k, T, steps, "all", p_ref=p_ref, form=0)
        eul = euler_constrained(k, T, steps, p_ref, 0)
        errs.append(abs(xi[-1] + 1 - a_ode))
        print(f"    {steps:3d} steps: alpha-1 {xi[-1]:.8e}, |element - Euler| {np.abs(xi + 1 - eul).max():.1e}, "
              f"|element - ODE| {errs[-1]:.2e}, p {p[-1]:.6f} (closed form {p_closed(1 + xi[-1]):.6f})")
    print("    error ratios between successive halvings of dt:", " ".join(f"{e1 / e2:.3f}" for e1, e2 in zip(errs, errs[1:])))
    for form in (1, 2):
        xi, p, _, _ = run(k, T, 44, "all", p_ref=p_ref, form=form)
        print(f"    form {form}, 44 steps: alpha-1 {xi[-1]:.8e}, |element - Euler| {np.abs(xi + 1 - euler_constrained(k, T, 44, p_ref, form)).max():.1e}, "
              f"ODE {ode_constrained(k, T, p_ref, form)[-1] - 1:.8e}")
    # (c) homeostatic limit: long run (T = 4, free growth would give 4e-3), growth stops at p_h
    for steps in (40, 160):
        xi, p, _, _ = run(k, 4.0, steps, "all", p_ref=p_ref, form=0)
        dpda = 3 * (1 + FV) * K * (1 + xi[-1]) ** -4
        print(f"(c) T = 4, {steps} steps: alpha-1 {xi[-1]:.8e} (alpha_h-1 {alpha_h(p_ref) - 1:.8e}), p {p[-1]:.6f} (p_h {p_ref}), "
              f"overshoot bound dt k dp/dalpha {4.0 / steps * k * dpda:.4f}, last increment of alpha-1 {xi[-1] - xi[-2]:.1e} "
              f"(first {xi[0]:.1e})")
    xi, p, _, _ = run(k, 4.0, 40, "all", p_ref=p_ref, form=2)
    print(f"    form 2 (no homeostatic value): alpha-1 {xi[-1]:.8e}, p {p[-1]:.4f}, last increment {xi[-1] - xi[-2]:.1e}")
    # (d) tangent
    rel, g, changed = tangent_check()
    print(f"(d) tangent with the gate at g = {g}: closed form vs central differences {rel:.1e}; gate changed the tangent by {changed:.1e} relative")
    # ---- the implicit gate (PROPS(20) = form + 10) ----
    print("implicit gate:")
    errs = []
    for steps in (11, 22, 44, 88):
        xi, p, _, _ = run(k, T, steps, "all", p_ref=p_ref, form=0, implicit=True)
        be = backward_euler_constrained(k, T, steps, p_ref, 0)
        errs.append(abs(xi[-1] + 1 - a_ode))
        print(f"(b) {steps:3d} steps: alpha-1 {xi[-1]:.8e}, |element - backward Euler| {np.abs(xi + 1 - be).max():.1e}, "
              f"|element - ODE| {errs[-1]:.2e}, largest p/p_h {p.max() / p_ref:.6f}")
    print("    error ratios between successive halvings of dt:", " ".join(f"{e1 / e2:.3f}" for e1, e2 in zip(errs, errs[1:])))
    for steps, TT in ((40, 4.0), (160, 4.0), (40, 16.0)):
        xi, p, _, _ = run(k, TT, steps, "all", p_ref=p_ref, form=0, implicit=True)
        print(f"(c) T = {TT}, {steps} steps: alpha-1 {xi[-1]:.8e} (alpha_h-1 {alpha_h(p_ref) - 1:.8e}), p {p[-1]:.8f}, largest p/p_h "
              f"{p.max() / p_ref:.8f}, 1 - p/p_h {1 - p[-1] / p_ref:.1e}, last increment {xi[-1] - xi[-2]:.1e}")
    for md in (False, True):
        for h in (1e-5, 1e-6):
            rel, g, changed = tangent_check(p_ref=0.3, implicit=True, md=md, h=h)
            print(f"(d) {'mean dilatation' if md else 'no mean dilatation'}, h = {h:.0e}: tangent vs central differences {rel:.1e} "
                  f"(g = {g:.3f}; the implicit terms change the tangent by {changed:.1f} relative)")
