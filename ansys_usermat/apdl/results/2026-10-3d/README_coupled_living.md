# Python runs of the agreed coupling (`coupled_living_consumption*.json`): which results to use

**Current main results (21 Oct, evening): names with `_gc1`** (nutrient uptake g l^2/d = 1, Klempt 2024 Table 2 converted),
quasi-static. Everything without `_gc1` used g l^2/d = 4 (Oliver's example deck) and is superseded.

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
