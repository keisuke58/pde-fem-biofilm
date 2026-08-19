# Thesis assignment — TMCMC↔FEM coupling in ANSYS

The scope as set by the supervisor, and what it means for this repository.
Companion documents: [`THESIS_PLAYBOOK.md`](THESIS_PLAYBOOK.md) (process),
[`RESEARCH_MODEL.md`](RESEARCH_MODEL.md) (science),
[`ansys_usermat/`](ansys_usermat/README.md) (the code this concerns).

---

## 1. The task in one sentence

> The Python TMCMC model is **pointwise** — it captures the temporal evolution at
> a material point but does not account for spatial variation. The task is to
> integrate it into an existing, working FEM framework **so that spatial effects
> can also be represented.**

That sentence is the thesis abstract in miniature; the Methods chapter should
open from it. Concretely: 0-D ecology ODE → spatially resolved field → FEM.

## 2. People and roles

| Who | Role |
|---|---|
| **Meisam Soleimani** | Supervisor. Sets the scope; co-author on the Klempt constitutive paper and on the TMCMC paper. |
| **Oliver** | **Co-supervisor** (newly asked). Holds the UserElement/UserMat implementation inherited from Felix. Technical support on UPFs. Contributes to any resulting publication. |
| **Felix Klempt** | Original author of the UserElement/UserMat implementation and of the continuum growth model. |
| **Timo** | Software installation support (ANSYS, compilers, VS Code) — as he did for Oliver and Felix. |

The thesis is described as *"a miniature version of the work underlying Oliver's
PhD thesis"* — so his PhD is the closest methodological precedent, and worth
reading early.

## 3. Immediate plan, as stated

1. Install the required software — ANSYS, the necessary **compilers**, VS Code.
   **Timo** supports this.
2. Visit **Oliver's** office for the basic procedure for user-programmable
   routines in ANSYS — **how to compile, link and invoke** the UserElement/UserMat
   routines.
3. **Oliver provides the final version** of the UserElement/UserMat
   implementation developed by Felix and subsequently transferred to him.

---

## 4. What this changes for this repository

Three things follow that are not currently reflected in the code, in decreasing
order of importance.

### 4.1 🔴 "UserElement/UserMat" — the repo has only USERMAT

The assignment says **UserElement/UserMat** throughout. This repo implements a
`USERMAT` only (`ansys_usermat/usermat_biofilm.f`). These are different
user-programmable features with different reach:

| | Sees | Can express |
|---|---|---|
| `USERMAT` | one integration point: `F`, state, `dt` | a **constitutive law** — stress and tangent from local deformation |
| `USERELEM` | the whole element: nodal DOFs, residual, element stiffness | extra **field DOFs** and their coupling — e.g. φ or α transported as nodal unknowns |

This matters because of §1. A constitutive law at a Gauss point cannot, by
itself, transport anything spatially — it only reacts to the deformation it is
handed. Representing spatial variation of the *ecology* (rather than just
receiving a precomputed α field) is naturally an element-level problem, which is
presumably why Felix's implementation is a UserElement **and** a UserMat.

**Open question to resolve early with Oliver and Meisam:** is the spatial field
to be (a) solved as extra DOFs inside a UserElement, or (b) precomputed
externally and passed in as a state variable per integration point? This repo
currently does (b) — the offline hand-off path in
`ch5_flow/flow_impl_architecture`. Option (a) is a materially larger and more
interesting piece of work, and the wording of the assignment leans that way.
**Do not start implementing until this is settled**; it determines the whole
architecture.

### 4.2 🟡 Felix's implementation is an incoming dependency this repo does not have

`ansys_usermat/usermat_biofilm.f` is **our own port** of the verified Abaqus UMAT
— not Felix's original. The assignment says the starting point is *Felix's final
version, as transferred to Oliver.*

So there are two codebases that need reconciling, and the answer is not obvious:

- **ours** — verified end to end (0 ULP vs the Abaqus UMAT over 8017 states,
  6.8e-14 vs the Python core), runs and converges in ANSYS 2022 R2, but is a
  `USERMAT` only;
- **Felix's** — the framework the assignment names, with the UserElement, but
  unseen here.

The sensible move once the code arrives: diff the constitutive cores. If they
agree, our verification evidence transfers to Felix's framework, which is
worth a lot — it means the growth law is trusted *inside the framework the thesis
is actually built on*. If they disagree, that discrepancy is itself a finding and
must be resolved before any result is quoted.

**Ask Oliver for it early** — it is on the critical path and nothing else about
the coupling design can be finalised without seeing it.

### 4.3 🟢 Two open items in the runbook now have owners

[`ansys_usermat/apdl/RUNBOOK.md`](ansys_usermat/apdl/RUNBOOK.md) flags two things
as unverified. Both are now answerable by a person rather than by guesswork:

| Runbook item | Ask |
|---|---|
| Step 0b — Intel Fortran / Visual Studio present and version-matched? | **Timo** |
| Step 1 — the actual compile/link procedure for UPFs in v222 on Windows | **Oliver** |

When Oliver walks through compile/link/invoke, **write down exactly what works**
and correct the runbook in the same sitting. That procedure is the single piece
of knowledge this repo most lacks, and it is currently held only in people's
heads — which is how it got lost between Felix and Oliver in the first place.

---

## 5. Consequences to keep in view

- **Authorship.** Oliver contributes to any resulting publication. Record this
  now so the author list is not reconstructed from memory later; `CITATION.cff`
  will need updating when the scope of his contribution is clear.
- **The verification chain is an asset, not overhead.** Walking into a shared
  framework with "this constitutive core is bit-identical across Abaqus, ANSYS
  and Python" is the strongest possible position from which to modify someone
  else's code — every later disagreement can be localised to the parts that
  changed.
- **This does not invalidate existing results.** The Abaqus tooth-stress results
  and the σ_CH/σ_DH ratio stand on their own. The ANSYS work is the coupling
  vehicle, not a replacement for them — keep the two lineages distinct in the
  write-up, as [`RESEARCH_MODEL.md`](RESEARCH_MODEL.md) §6 already insists.

---

## 6. Source

Email from **Meisam Soleimani**, subject *"Master thesis of Keisuke — TMCMC-FEM
coupling in ANSYS: support and supervision"*, addressed to Oliver, cc Keisuke,
2026-08-19. Quoted here because it is the authoritative statement of scope.

> As you know, Keisuke (cc'd here) has been on our team since last year and will
> begin his master's thesis with us today. His project can be regarded as a
> miniature version of the work underlying your PhD thesis.
>
> Keisuke has developed a multispecies material model in Python using
> TMCMC-based Bayesian updating. At present, it is a pointwise model: it captures
> the temporal evolution at a material point but does not account for spatial
> variation. His task will be to integrate this Python model into an existing,
> functioning FEM framework so that spatial effects can also be represented.
>
> You guessed correctly: he will work with the UserElement/UserMat implementation
> that you inherited from Felix. In practical terms, he will couple his Python
> code with the existing UserElement/UserMat routines.
>
> I would therefore like to ask you to co-supervise his thesis with me. I believe
> this would also be valuable for you, as his work follows a logic similar to
> yours and provides a useful hands-on application. Furthermore, you would
> contribute to any future publication resulting from this work. From my side,
> your involvement would be a great help in providing Keisuke with the necessary
> technical support.
>
> The immediate plan would be as follows:
>
> I have asked Keisuke to install all the required software, including ANSYS, the
> necessary compilers and Visual Studio Code. Timo will support him in this
> regard as he did greatly for you and Felix. He may visit your office for
> guidance on the basic procedure for using user-programmable routines in
> ANSYS—for example, how to compile, link and invoke the UserElement/UserMat
> routines.
>
> Please provide him with the final version of the UserElement/UserMat
> implementation developed by Felix and subsequently transferred to you.

> **Note on visibility.** This repository is **public**. The quotation above is
> private correspondence naming colleagues. Nothing in it is confidential and the
> names already appear elsewhere in the repo, but it was reproduced without
> asking the sender. If that is not wanted, replace §6 with a paraphrase — §§1–5
> carry all the information the thesis actually needs.
