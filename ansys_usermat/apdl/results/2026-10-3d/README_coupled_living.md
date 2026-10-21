# Python runs of the agreed coupling (`coupled_living_consumption*.json`): which results to use

**Current main results (21 Oct, late): `_front1_phimin0.001_..._a10_gab1_gc1_tab0.333333`** (antibiotic from T* = 1/3,
phi_min = 0.001 converged); the same without `_tab0.333333` and with `_phimin0.1` were the step before (antibiotic 0.5).
Earlier: names with `_front1` ... `_a10_gab1_gc1`:
front term r of Klempt 2024 Table 2 (10 mm/T*), d and g of Table 2 (quasi-static, g l^2/d = 1), alpha*_0 = 10
(Fritsch case I, Klempt 2026 Table 1), antibiotic with the same diffusion and uptake as the nutrient (Oliver: the
same equation as c), antibiotic from T* = 0.5 (Klempt 2026 Table 2, case 4). phi_min = 0.1 has no paper value.
Everything else (r/10, alpha*_0 = 50, g_a = 0.25, g = 4) is superseded; `_gc1` with `_front0.1` was an intermediate step.

Earlier note: quasi-static (no `_transc`) or `_transc1e+08`
which agree to < 0.005: d of Klempt 2024 Table 2 converted to the 2 mm cube is 1e8 mm^2/T*, so the time derivative
of Eq. 35 does not matter. `_transc1`, `_transc0.3`, `_transc3` (d = 1, 0.3, 3 mm^2/T*, Oliver's example deck value
and around it) are a **slow-diffusion sensitivity** only. Never mix them in one table without saying which is which.

| tag | meaning |
|---|---|
| `_transc1` (`_transc0.3`, `_transc3`) | with dc/dt, da/dt, d = 1 (0.3, 3) mm^2/T* -> **main** (d sensitivity) |
| no `_transc` | quasi-static -> comparison with the quasi-static ANSYS runs `u1_` |
| `_ansysbc` | boundary conditions of the ANSYS decks (layers held) -> use for comparisons with ANSYS |
| `_startpos` | new points start by position (the ANSYS rule); without it: from the neighbours' state |
| `_mixeq` | mixed start with equal amounts (the agreed example); without it in a mixed run: Fritsch's start (0.25, 0.30) |
| `_nf16`, `_nf48` | mesh (default 32^3) |
| `_v2` | the v2 set-up (B in corner blocks, phi_min 0.1, antibiotic on all outer faces, alpha*_0 = 50) |
| files without `_v2` or with `_fortran` missing | earlier development runs (before 20 Oct evening): do not use |
