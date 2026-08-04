# Thesis Playbook — how to drive the thesis to submission, check it, and continue

A single systematised view of **what exists, how to advance it, how to verify it,
and where it goes next.** It ties together the detailed docs rather than repeating
them:

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
| Physics | Mesh convergence | 🟢 verified | `VERIFICATION…` V-series |
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

**Pre-submission checklist (the 3 must-dos + closure):**
- [ ] **Regenerate the stale DS artifact** (`tooth_klempt_comparison.json` `_flat_golden`)
      via the Abaqus run → clears the `test_golden_stress.py` xfail tracker.
- [ ] **`audit_all.py --strict-env` = ALL CLEAR** on the complete workspace (every
      thesis figure regenerates; external inputs present).
- [ ] **Citations correct** — Klempt, Geisler, **Soleimani** et al. (2026),
      *Arch. Appl. Mech.* 96, 164, doi:10.1007/s00419-026-03160-y; confirm any
      referenced equation numbers match this published version. *(Soleimani is a
      co-author — get this exactly right.)*
- [ ] Every quoted number traces to a committed script (grep the thesis, check each).
- [ ] Limitations section carries `E_SPEC`-assumed / `α`-magnitude-anchored /
      composition-CLSM-anchored verbatim in intent.

---

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

---

*This playbook is a living index — update the checkboxes in §3 as the submission
must-dos are cleared, and move §4 items as they start.*
