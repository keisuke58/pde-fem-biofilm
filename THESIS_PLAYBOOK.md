# Thesis Playbook — how to drive the thesis to submission, check it, and continue

A single systematised view of **what exists, how to advance it, how to verify it,
and where it goes next.** It ties together the detailed docs rather than repeating
them:

- **scope as set by the supervisor** → [`THESIS_ASSIGNMENT.md`](THESIS_ASSIGNMENT.md) (roles, deliverables, open questions)
- rigor audit → [`VERIFICATION_SENSITIVITY_LIMITATIONS.md`](VERIFICATION_SENSITIVITY_LIMITATIONS.md) (read first)
- roadmap → [`PLAN_NEXT.md`](PLAN_NEXT.md) · research ladder → [`research_goals_1_2.md`](research_goals_1_2.md)
- methods → [`methods_supplement_fem.md`](methods_supplement_fem.md) · code map → [`REPO_MAP.md`](REPO_MAP.md) · doc index → [`DOCS.md`](DOCS.md)

---

## 0. TL;DR (one screen)

- **State:** the constitutive core is verified and the repo is a citable,
  CI-tested, reproducible artifact. The thesis contribution is **complete and
  defensible**.
- **Now → submission:** **FREEZE** (no new physics). Do the 3 must-dos in §3,
  then write, checking every claim against §2.
- **After (Keio M2, 2027-04):** the science frontier — first paper, the
  Python↔Fortran coupling (skeleton in place), and closing the measurement gaps.

---

## 1. 今までの体系化 — what is built & verified (inventory)

| Layer | Component | Status | Evidence / file |
|---|---|---|---|
| Physics | Growth kinematics `F=Fe·Fg`, `Fg=(1+α)I` | 🟢 verified | production UMAT; `VERIFICATION…` V1 |
| Physics | Consistent tangent (DDSDDE) | 🟢 verified | analytic + F-perturbation vs FD 2.4–2.9e-8; `phase2_patch_test.py` 13/13 |
| Physics | Dual-solver port (Abaqus↔ANSYS) | 🟢 bit-identical | `ansys_usermat/crosscheck/` — 0 ULP over 8017 cases |
| Physics | USERMAT **in ANSYS 2022 R2**, mechanical branch | 🟢 runs & converges | `SOLID185`/`NLGEOM,ON` uniaxial benchmark; interface args, `keycut`/`cutFactor`, `dsdePl` validated in-solver |
| Physics | USERMAT **in ANSYS 2022 R2**, growth branch (`Fg=(1+α)I`, `α≠0`) | 🟢 verified 2026-08-19/20 | closed-form checks, `ansys_usermat/apdl/`: constrained cube (all 4 `reference_values.json` cases + `KEYOPT` sweep) **and** the complementary free/traction-free cube (`t_growth_free.dat`, stress≡0 confirmed) both match exactly; see `RUNBOOK.md` |
| Physics | Mesh convergence | 🟢 verified | `VERIFICATION…` V-series |
| Physics | Single-element, partial-constraint (non-closed-form) smoke test on both solvers | 🟢 verified 2026-08-20 | real Abaqus 2024 run (`umat_tangent_test/abaqus_1elem/`, completed successfully) and real ANSYS run (`ansys_usermat/apdl/t_growth_baseclamped.dat`, 0 errors) — base fixed, top free, growth+viscosity, qualitatively sensible on both; prep ahead of Felix's UserElement code per `HANDOFF.md` |
| Physics | Geometric realism beyond the unit cube (curved, bonded two-layer shell) | 🟡 partial, 2026-08-20 | `ansys_usermat/apdl/t_growth_cylinder_shell.dat` — two real input bugs found & fixed (VGLUE attribute loss, undersized mesh in the thin layer); **converges cleanly at α=0.01** (target α=0.05 does not; threshold between 0.01 and 0.015, four targeted BC/mesh attempts to push past it all reverted). At converged α, outer surface shows a two-lobe displacement pattern, consistent with early buckling but not confirmed against a real eigenvalue analysis; see `cylinder_shell_bulge_analysis.ipynb` |
| Result | Headline `σ_CH/σ_DH ≈ 6.44×` (early) | 🟢 frozen | `tests/test_golden_stress.py`; `JAXFEM/_posterior_ci/` |
| Result | Model ↔ experiment (dysbiotic/static) | 🟢 validated | `validate_composition.py` — MAE 4.2 pp, TVD 0.11 |
| Input | Composition φ (CLSM-measured) | 🟢 anchored | TMCMC calibrates interaction matrix A, not φ |
| Input | Per-species stiffness `E_SPEC` | 🟡 assumed | reported as a band (3.7–12×); needs species AFM |
| Input | Per-condition growth `α` | 🟡 magnitude-anchored | not calibrated per condition |
| Engineering | CI, reproducibility, citability, site | 🟢 done | `ci.yml`, `reproduce.sh`, `LICENSE`/`CITATION.cff`, `docs/` |
| Continuation | Python↔Fortran coupling | 🟢 skeleton | `ansys_usermat/coupling/` (runnable, tested) |

**Two analysis lineages** (keep them distinct in the write-up):
1. **Klempt growth-stress** — the verified core above → the thesis headline.
2. **DI-bridge FEM** — an alternative bridging variable (Dysbiosis Index → E(DI));
   ratios are scale-invariant. See [`FEM_README.md`](FEM_README.md).

---

## 2. 修論の進め方 — chapter-driven, evidence-backed workflow

**Principle: every quoted number traces to a committed script.** Write the thesis
as claims, and attach the evidence file to each. Suggested claim → evidence map:

| Thesis claim | Backed by |
|---|---|
| Constitutive law is Klempt-faithful & verified | UMAT sources; `VERIFICATION…` V1–V2; `crosscheck/` |
| Growth field from calibrated ecology | `JAXFEM/` (PDE); `methods_supplement_fem.md` (TMCMC→Monod) |
| σ_CH/σ_DH ≈ 6.4× early, ~2× mature | `test_golden_stress.py`; σ(t) figure |
| Result robust to depth model, sensitive to E_SPEC | `VERIFICATION…` sensitivity; README bands |
| Model matches experiment (composition) | `validate_composition.py` |
| Every figure regenerates from code | `reproduce.sh`, `audit_all.py` |

**Figure/number pipeline:** regenerate with `./reproduce.sh`; confirm coverage with
`python JAXFEM/audit_all.py --strict-env` on the full workspace.

**Write-up order (suggested):** Methods (model → calibration → FEM) → Verification
→ Results (σ(t), the ratio, validation) → Sensitivity & Limitations (verbatim from
the rigor audit) → Discussion → Future work (reference, don't include, the held
extensions).

---

## 3. チェック方法 — the verification / QA loop

**Automated (run any time, green = safe):**
```bash
pytest tests/                                   # unit + validation + golden guard
pytest ansys_usermat/crosscheck/crosscheck.py ansys_usermat/crosscheck/adversarial.py
python JAXFEM/audit_all.py --quick              # runnable-subset audit
```
CI runs all of the above on every push.

**Manual rigor (the source of truth):** `VERIFICATION_SENSITIVITY_LIMITATIONS.md`
— what is Verified vs Sensitive vs a disclosed Limitation. Keep the thesis claims
inside what it supports.

### Recorded state — full `./reproduce.sh` run (2026-08, clean checkout)

| Check | Result |
|---|---|
| `./reproduce.sh` | **exit 0** — all steps completed |
| Figures regenerated | `heine_species_composition`, `heine_phi_psi_joint`, `validation_composition_dysbiotic` — **pixel-identical** to the committed PNGs (max pixel diff 0) |
| `audit_all.py --quick` | **ALL CLEAR** (runnable subset); 4 sections SKIPPED — external inputs absent (Abaqus extracts / sibling repos / author workspace) |
| Test suite | **142 collected → 139 passed, 3 xfailed** (the 3 xfails are intentional trackers: `E_di` bounds + the two stale-DS artifact entries) |
| Model↔experiment validation | MAE **4.216 pp**, TVD **0.1054**, worst species 10.54 pp (dysbiotic/static) |
| Coupling equivalence | Python core ≡ Fortran core, worst **6.8e-14** relative over 28 states; live-reproducible in `ansys_usermat/coupling/python_core_vs_fortran_verification.ipynb` |

**Reproducibility note.** The figures regenerate *pixel-identically*, but the PNG
bytes also carry matplotlib's version string, so byte-level identity additionally
requires the pinned toolchain (`pip install -r requirements.txt` — matplotlib
3.11.0). Running with a different patch release leaves `git status` showing the
figures as modified even though the rendered content is unchanged; verify with a
pixel comparison, not `git diff`, before assuming a real change.

**Still gated on the full workspace** (cannot be closed here without further
work): the 4 SKIPPED audit sections and the stale DS artifact both need an
Abaqus run — see the checklist below.

> **Correction, 2026-08-20:** earlier text here (and elsewhere this session)
> assumed IKMHIWI03 has no Abaqus at all, on the strength of "this machine is
> the primary ANSYS environment" in CLAUDE.md. That's wrong — **Abaqus 2024
> is actually installed and licensed on IKMHIWI03**
> (`C:\SIMULIA\Commands\abaqus.bat`, confirmed working:
> `abaqus information=release` completes, valid Site ID). What's genuinely
> missing here is the *prior run output* — no `.odb`/`.sta`/`.msg`/`.dat`
> scratch files exist anywhere on this machine (checked the whole `C:` drive)
> — not Abaqus itself. So a **fresh** regeneration of the DS artifact and a
> fresh cost-timing run are technically possible on IKMHIWI03 now, *if* the
> job-generation pipeline (`tier2b_real/`, `configs/`) and its required input
> data are also available here — that has not yet been checked.

**Pre-submission checklist (the 3 must-dos + closure):**
- [ ] **Regenerate the stale DS artifact** (`tooth_klempt_comparison.json` `_flat_golden`)
      via the Abaqus run → clears the `test_golden_stress.py` xfail tracker.
- [ ] **`audit_all.py --strict-env` = ALL CLEAR** on the complete workspace (every
      thesis figure regenerates; external inputs present).
- [ ] **Citations correct** — Klempt, Geisler, **Soleimani** et al. (2026),
      *Arch. Appl. Mech.* 96, 164, doi:10.1007/s00419-026-03160-y; confirm any
      referenced equation numbers match this published version. *(Soleimani is a
      co-author — get this exactly right.)*
- [ ] **Capture the computational cost** — run `extract_abaqus_cost.py` over the
      Abaqus scratch directories (the `.sta`/`.msg`/`.dat` files are still on the
      workstation) and commit `runs/abaqus_cost.json`. Currently the repo records
      stress results but no timings, so the cost section has no measured basis.
- [ ] Every quoted number traces to a committed script (grep the thesis, check each).
- [ ] Limitations section carries `E_SPEC`-assumed / `α`-magnitude-anchored /
      composition-CLSM-anchored verbatim in intent.

---

## 3.5 Blocking questions for the supervisors

Two items from [`THESIS_ASSIGNMENT.md`](THESIS_ASSIGNMENT.md) sit on the critical
path and cannot be resolved by working harder here — **status as of
2026-08-19, narrowed but not settled:**

- **UserElement or UserMat?** The assignment names *UserElement/UserMat*; this
  repo has a `USERMAT` only. Whether the spatial field is solved as extra DOFs in
  a UserElement or precomputed and passed in per integration point determines the
  entire coupling architecture. **Still open** — but no longer a blind guess:
  [`ansys_usermat/USERELEM_NOTES.md`](ansys_usermat/USERELEM_NOTES.md) works
  through ANSYS's own `UserElem.F` example and confirms a UserElement can call
  our *existing, verified* USERMAT unchanged for the mechanical response
  (`KEYANSMAT=1` → `ElemGetMat`) — so the constitutive-core verification work
  above transfers either way. The genuinely new, un-templated work is the
  field-DOF residual/stiffness for spatial ecology transport itself. A
  deadline-based default assumption (extra DOFs in a UserElement) has been
  proposed to Oliver/Meisam; proceeding on it if no objection lands.
- **Felix's final implementation** is the stated starting point and is not in
  this repo. **Update, 2026-08-20:** Felix replied — he has left IKM and can
  only give informal input going forward, not a detailed USERMAT review (that
  now routes to Oliver/Meisam/Hendrik instead). He also confirmed his own
  approach has no viscoelastic split (`F_v`); this repo's model provably
  reduces to exactly his `F=F_eF_g` at `η=0` (proven live over a 50-step
  chained load history, see the appendix in
  `python_core_vs_fortran_verification.ipynb`), and that `η=0` slice is
  already inside the closed three-way verification chain, not a hypothetical
  reduction. [`ansys_usermat/crosscheck/crosscheck.py`](ansys_usermat/crosscheck/crosscheck.py)
  still supports diffing against a third source (`--right-src` etc.) once his
  actual code arrives (template: `xcheck_driver_template.f`). **Code itself
  still not received; timeline now uncertain given his departure.**

Separately, **the build/toolchain blocker that RUNBOOK.md flagged Timo for is
resolved** (not classic-ifort, as first suspected — ANSYS's shipped
`ansys.lrf` was just missing one `-defaultlib` line) — see `RUNBOOK.md` for
the full story, including a real APDL deck bug (`TBDATA`'s 6-value-per-call
limit) found and fixed along the way.

## 4. 今後の方針 — go-forward (systematised)

```
  now ──────────────► thesis submission ──────────► Keio M2 (2027-04) ──────────►
  FREEZE the core        3 must-dos (§3)             T2 first paper (jaw-level twin)
  no new physics         + write-up                  Python↔Fortran coupling (skeleton ready)
                                                      close E_SPEC / α measurement gaps
```

- **T1 — thesis:** freeze; ship. Highest priority, shortest horizon.
- **T2 — first paper:** jaw-level mechano-ecological digital twin (Level 2) — the
  publishable next step; builds on the verified core.
- **T3 — continuation (Keio):** wire the coupling (`ansys_usermat/coupling/` →
  C shim → single-element smoke test → JAX model), viscoelastic/phase-field/
  multiscale extensions. The IKM↔Keio bridge.

Do **not** interleave T1/T2/T3 on one branch. See `PLAN_NEXT.md` §1 for sequencing.

---

## 5. リスクと守り

| Risk | Mitigation (already in place / to do) |
|---|---|
| Result sensitive to `E_SPEC` | reported as a **band**, not a point; disclosed as a limitation |
| DS composition bug re-creeping in | fixed end-to-end; `test_golden_stress.py` **guards** it; stale artifact tracked by xfail |
| A number silently changes | golden-value + validation regression tests fail CI |
| Citation wrong before a co-author | §3 checklist; correct Klempt 2026 recorded |
| Scope creep before submission | FREEZE principle; continuation work is scaffolded, not started |
| IKMHIWI03's `C:` drive hitting 0 bytes free mid-build | hit once, 2026-08-20 (killed an ANSYS run outright). **Resolved same day**: all ANSYS/Abaqus work (build + run + scratch) moved to `F:\` (~3.7 TB free, essentially unused), via `run_apdl.ps1`/`run_abaqus.ps1` which default there and refuse to run below a free-space threshold. `C:`'s underlying VSS/System-Restore retention issue itself is still unresolved (needs admin) but is now moot for this workflow — flagged to Timo separately |

---

*This playbook is a living index — update the checkboxes in §3 as the submission
must-dos are cleared, and move §4 items as they start.*
