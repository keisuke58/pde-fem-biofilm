# Python runs of the agreed coupling (`coupled_living_consumption*.json`): which results to use

The file name carries the settings. **Main results: names with `_transc1`** (Eq. 35 and the antibiotic field with
their time derivative, d = 1 mm^2/T*, as in Oliver's sketch). Names without `_transc` are **quasi-static**
(0 = d lap c - ...): comparison only. Never mix the two in one table without saying which is which.

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
