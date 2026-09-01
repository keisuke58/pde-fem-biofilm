# Hand-off to Felix — biofilm material model → ANSYS FEM

Goal: replace the phenomenological constitutive model in Felix/IKM's existing
ANSYS FE model with the paper's calibrated biofilm growth / viscoelastic law,
called at each Gauss point (eventually via a Python material model). This sheet
lists **what to hand over, in what order, and what is needed back** to run it.

Status of the deliverable: the ANSYS `USERMAT` **runs and converges in ANSYS
MAPDL 2022 R2 (v222)** — interface arguments, `keycut`/`cutFactor` substepping
and the `dsdePl` Jacobian validated on a `SOLID185` / `NLGEOM,ON` uniaxial
benchmark — and its constitutive core is **bit-identical to the verified Abaqus
UMAT** (0 ULP over 8017 states). Still to check in-solver: the **growth** path
(`α > 0`), for which `apdl/` carries a closed-form single-element case.

---

## ① Core — hand over first (self-contained, compiles today)

- [ ] **`ansys_usermat/usermat_biofilm.f`** — the ANSYS `USERMAT`. Felix swaps his
      phenomenological law for this. Growth `Fg=(1+α)I` + Neo-Hookean/`D1` +
      backward-Euler viscoelasticity + F-perturbation consistent tangent.
- [ ] **`ansys_usermat/README.md`** — read first: Abaqus↔ANSYS mapping table,
      `prop`/`ustatev` layout, `TB,USER`/`TB,STATE` setup, build & syntax-check.

Minimum viable hand-off = ① only. With it Felix can compile and drop the law
into a single-element ANSYS model.

## ② Evidence — so the law is trusted

- [ ] **`umat_biofilm_visco.f`**, **`umat_biofilm_visco_phase2.f`** — the verified
      Abaqus UMATs the port mirrors (the reference implementation).
- [ ] **`phase2_patch_test.py`** — Python replica; patch tests 13/13, tangent
      vs FD ~2.4e-8. (The ANSYS port reproduces the same algebra → 2.97e-8.)
- [ ] **`VERIFICATION_SENSITIVITY_LIMITATIONS.md`** — what is verified vs assumed.
- [ ] *(optional)* `umat_flow/` — UMAT algorithm-flow figures.
- [ ] **Notebooks (bilingual EN/JA, reproduce all of the above live)**:
      [`crosscheck/crosscheck_analysis.ipynb`](crosscheck/crosscheck_analysis.ipynb),
      [`coupling/python_core_vs_fortran_verification.ipynb`](coupling/python_core_vs_fortran_verification.ipynb),
      [`apdl/growth_closed_form_verification.ipynb`](apdl/growth_closed_form_verification.ipynb),
      [`../umat_tangent_test/abaqus_1elem/abaqus_visco_1elem_analysis.ipynb`](../umat_tangent_test/abaqus_1elem/abaqus_visco_1elem_analysis.ipynb)
      (a real single-element Abaqus 2024 run of the production UMAT, not just
      the crosscheck driver — the closest existing analog to the "single-
      element smoke test" asked for below, on Abaqus rather than ANSYS).

## ③ Coupling plan — the Python-at-Gauss-point step (the thesis work)

**Code-side wiring is done** — `kUsePy=1` calls the Python material hook end
to end through the real `usermat()` subroutine (`MAP6` Voigt reindex,
fallback to the inline core on failure) and is verified against `kUsePy=0`
in `tests/test_usermat_kusepy_e2e.py`. Still open: an actual ANSYS build/link
of this path (`apdl/RUNBOOK.md` has the new build-dependency note) and
swapping in the calibrated JAX model in place of the current NumPy mirror.

- [ ] **`ch5_flow/flow_impl_architecture.png`** — Python(JAX) ↔ ISO_C_BINDING /
      socket ↔ USERMAT ↔ ANSYS Gauss-point loop.
- [ ] **`ch5_flow/flow_python_material_hook`** — the wire protocol and
      per-Gauss-point payload, grounded in the actual `protocol.py`/
      `material_server.py` implementation.
- [ ] **`ch5_flow/flow_growth_kinematics.png`** — `F=Fe·Fg` kinematics (the model's basis).
- [ ] **`material_models.py`**, **`JAXFEM/`** — the Python material model that the
      `PYTHON MATERIAL HOOK` in the USERMAT will call (currently a NumPy mirror
      of the verified Fortran law, not yet the calibrated JAX model).

---

## Needed back from Felix (blocks an actual ANSYS run)

> **Partly answered, 2026-09-01.** Oliver sent the Workbench archive
> `BiofilmImplementation.wbpz`, which answers item 1 (it targets **ANSYS 2024
> R2**, not the v222 this repo is verified on) and shows a working, converging
> `SOLID185`/`NLGEOM,ON` model. It does **not** answer item 2 or 3, and it
> arrived **without its USERMAT Fortran source**, so it cannot be run as
> received. It is also a *nonlocal* material routine — a different scope from
> this local law. Full findings and the follow-up questions:
> [`OLIVER_MODEL_NOTES.md`](OLIVER_MODEL_NOTES.md).

1. [ ] **Target ANSYS version + exact `usermat` argument list** (`var1..var8`,
       `tsstif`, `epsZZ`, `cutFactor` vary by release), and whether **`dsdePl`**
       is expected as the Jaumann material Jacobian (may need rotation `rotateM`
       / symmetrisation).
2. [ ] **Material parameters** `C10, C01, D1, eta` (units: stress / 1-stress /
       stress·time) — and the **growth `α` field**: how the JAXFEM α-field is
       delivered to each integration point (`TB,STATE` init vs user field).
3. [ ] **A single-element smoke test** in ANSYS (uniaxial / simple shear) to
       confirm interface wiring before touching the full model. (The
       constitutive law is already proven **bit-identical to the Abaqus UMAT**
       — `crosscheck/`, 20/20 cases — so any mismatch would be interface, not
       physics.)

## Suggested sequence

1. Hand over ①; Felix confirms the `usermat` signature for his ANSYS version.
2. Single-element smoke test → compare stress/tangent against the Abaqus UMAT
   (same F, same params) — should match to machine precision.
3. Plug into the full model with the α-field (item 2 above), phenomenological
   law removed.
4. ~~Wire the `PYTHON MATERIAL HOOK` (③)~~ — the thesis' core contribution;
   **code-side done**, verified against the inline Fortran core (kept as the
   verification reference / fallback) via `tests/test_usermat_kusepy_e2e.py`.
   Still needed: an ANSYS build/link smoke test and the JAX model swap-in.

---

*One-line framing for Felix:* "A USERMAT that implements our verified biofilm
growth/viscoelastic law (matches the Abaqus UMAT to machine precision) with an
explicit hook to later call the Python material model per Gauss point — needs a
version/interface check and a single-element smoke test before going into the
full model."
