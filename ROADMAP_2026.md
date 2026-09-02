# Roadmap — submission November 2026, defence December 2026

Set 2026-09-01. Supersedes the "early December" date in
[`ansys_usermat/INTEGRATION_PLAN.md`](ansys_usermat/INTEGRATION_PLAN.md).
Process companion: [`THESIS_PLAYBOOK.md`](THESIS_PLAYBOOK.md).

## Fixed points

| | |
|---|---|
| **Submission** | November 2026 — plan to **week 12 (Nov 23–27)**, not Nov 30, so the last week is buffer |
| **Defence** | December 2026 |
| **Weeks from now to submission** | **12** |

Confirm the exact submission date with M. Soleimani at the first check-in
(§4); everything below is anchored to late November and shifts as a block.

## 1. The Chapter 5 decision

**Ch5 = the ANSYS material-law contribution** (this repo's `ansys_usermat/`,
the Oliver-framework integration). The JAXFEM-based Ch5 draft that already
exists is **held back for Keio**, not submitted here.

This is not a change of direction — it is the split
[`THESIS_PLAYBOOK.md`](THESIS_PLAYBOOK.md) §4 already draws, made concrete:

```
T1  thesis (LUH, now)      ANSYS USERMAT + growth law + the handover     → Ch5
T3  continuation (Keio)    JAXFEM, coupling, phase-field, multiscale     → held back
```

Two reasons it is the right way round, worth stating in the chapter itself:

- **It answers the assignment.** [`THESIS_ASSIGNMENT.md`](THESIS_ASSIGNMENT.md)
  §1 sets the task as *integrating the pointwise TMCMC model into an existing,
  working FEM framework so spatial effects are represented.* Oliver's framework
  **is** that existing working framework. JAXFEM is a framework of our own —
  a good one, but answering a different question.
- **It is the part that is verified.** The core is 0 ULP against the Abaqus
  UMAT over 8017 states; the wrapper inherits that (PR #46). A chapter can be
  written on it without any new physics — which is what the FREEZE principle
  asks for.

## 2. What Ch5 can already claim — today

This is the load-bearing paragraph of the whole plan. Write the chapter so
that **everything below Tier A is already true**, and nothing in it waits on
anyone else's calendar.

### Tier A — done, verifiable, in the repo

| Claim | Evidence |
|---|---|
| A growth + viscoelastic constitutive law implemented as an ANSYS USERMAT | `ansys_usermat/usermat_biofilm.f` |
| It is numerically identical to the independent Abaqus UMAT | `crosscheck/` + `adversarial.py`, 0 ULP over 8017 states |
| It is packaged to drop into the partner framework at their own call site | `biofilm_material_v01.f`, `BIOFILM_GROWTH_VISCO_V01` |
| The repackaging changes no physics | `test_wrapper_is_only_an_adapter`, rtol=0 atol=0 |
| The deliverable is ANSYS-release-independent | no ANSYS includes / common block / UPF calls, same as their `AceGenNeoHookV04` |
| The interface constraints are characterised, not assumed | `Fv` must be 9 slots (measured asymmetry ~6e-5); `dt` must resolve `η/(2·C10)` (σ₁₁ flips sign at `dt/τ≈0.5`) |

That is a complete methods-and-verification chapter. It is defensible in the
viva on its own.

### Tier B — the results section, if the runs land

Coupled ANSYS runs producing von Mises fields under the growth law. **This is
where the risk is**, and §3 is about removing our dependence on other people
for it.

## 3. Critical path — and the one move that shortens it

The integration plan has four steps. Their dependency structure matters more
than their order:

| Step | Status | Depends on |
|---|---|---|
| 2. Wrap the core to their convention | ✅ done (PR #46) | us |
| 3. Verify the wrapper | ✅ done (PR #46) | us |
| 1. Port to ANSYS v222 locally | ✅ **done 2026-09-02**, ahead of the Sep 7 target — see [`ansys_usermat/apdl/V222_PORT_INSTRUCTIONS.md`](ansys_usermat/apdl/V222_PORT_INSTRUCTIONS.md) | us + this machine |
| 4. They wire it into their framework | open | **Oliver's calendar** |

**Step 1 is now the highest-value item, and its priority has changed.** It was
scoped as "a realistic test bed, not a prerequisite." With a November deadline
it is more than that: a working local v222 build lets us run our own ANSYS jobs
with the wrapper and produce Tier B results **without waiting for step 4 at
all**. It converts the one dependency we do not control into a bonus.

So: **Tier B must be reachable through step 1 alone.** Step 4 landing before
submission would be excellent and belongs in the chapter if it happens — but
the chapter must not be written assuming it will.

## 4. Communication with M. Soleimani and Oliver

Regular, short, online. The point is that neither of them is ever surprised,
and that step 4 has been asked for early enough to be possible.

| Cadence | With | Form |
|---|---|---|
| **Every 2 weeks** | Meisam | Short written update — what landed, what is next, anything blocking. 5 lines. |
| **Monthly** | Meisam | Online call — scope, chapter structure, results as they appear |
| **On a concrete milestone** | Oliver | Written, one topic at a time — the handover package, then the interface notes, then the wiring |
| **Once, early** | Oliver | Online call when handing over the routine, so the interface notes get discussed rather than just read |

Three things to raise early rather than late:

1. **Confirm the exact submission date** with Meisam (week 1).
2. **Send Oliver the handover package** (week 2–3) — the routine plus the two
   interface constraints. He needs lead time if step 4 is to happen at all.
3. **Agree the Ch5 scope with Meisam explicitly** (week 1–2), including that
   the JAXFEM material is held back. Better as a decision he signs off on now
   than as a surprise in November.

## 5. Week by week

| Weeks | Focus | Done means |
|---|---|---|
| **1** (Sep 7) | v222 build on IKMHIWI03 (step 1) — ✅ done 2026-09-02, ahead of schedule. Still open: confirm date + Ch5 scope with Meisam. | custom `ANSYS.exe` runs a deck with the wrapper linked — ✅ exact closed-form match |
| **2–3** (Sep 14, 21) | Handover package to Oliver + call. First smoke runs with the wrapper on v222. | Oliver has the routine; a single-element run gives sane stress |
| **4–5** (Sep 28, Oct 5) | Real geometry runs — coupon, then tooth/implant. **Coupon-scale piece done 2026-09-02** (re-ran the existing two-layer curved-shell case through the wrapper — identical SEQV to the original build across all 12240 elements, see `V222_PORT_INSTRUCTIONS.md`). Tooth/implant scale still open. | von Mises fields out, `dt` inside the stable range and shown to be |
| **6–8** (Oct 12, 19, 26) | Condition comparison, the actual Ch5 results. Figures. | the comparison table and figures Ch5 needs exist |
| **9–10** (Nov 2, 9) | **Write Ch5.** FREEZE — no new physics from here. | full draft to Meisam |
| **11** (Nov 16) | Revise on his comments. Citations, numbers, consistency checks. | `audit_all.py` clean, no loose claims |
| **12** (Nov 23) | **Submit.** | submitted |
| **13** (Nov 30) | Buffer. | — |
| **14–15** (Dec 7, 14) | Defence prep — slides from the Ch5 figures. | rehearsed |

## 6. Risks

| Risk | What we do about it |
|---|---|
| Step 4 does not happen before submission | Tier A/B split (§2) + step 1 as the independent route (§3). The chapter never depends on it. |
| v222 build fails | `apdl/V222_PORT_INSTRUCTIONS.md` has the pre-flight findings; Timo supports toolchain issues; fallback is Abaqus, which is installed and licensed here |
| Production runs use `dt` past `η/(2·C10)` | Now a known, pinned constraint — **check it on the first real run**, not after the results are in |
| Scope creep back into JAXFEM | It is held for Keio by decision, not by accident. Any JAXFEM work before submission needs a reason to exist in Ch5. |
| Silent regression in the verified core | CI: golden values, crosscheck, adversarial hunt |

## Later, not now

Two continuation ideas are written up rather than started, so the weeks to
November stay on the ANSYS contribution:

- [`PINN_DESIGN.md`](PINN_DESIGN.md) — a physics-informed surrogate for the
  biofilm field. Lands in the JAXFEM material held back for Keio.
- Option (D) in [`OLIVER_MODEL_NOTES.md`](ansys_usermat/OLIVER_MODEL_NOTES.md)
  — borrowing the NEM derivative operator to strengthen `JAXFEM/`.

They compete for the same slot; §6 of the PINN note says why, and the choice
belongs at Keio, not now.
