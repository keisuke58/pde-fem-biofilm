# Research content, systematised — the model chain end to end

The **scientific** content of the project in one place: what is modelled, with
which equations, and — for every quantity — whether it is **measured**,
**calibrated**, or **assumed**. This is the skeleton of a Methods chapter.

Companion documents: process & submission workflow → [`THESIS_PLAYBOOK.md`](THESIS_PLAYBOOK.md) ·
rigor audit → [`VERIFICATION_SENSITIVITY_LIMITATIONS.md`](VERIFICATION_SENSITIVITY_LIMITATIONS.md) ·
methods detail → [`methods_supplement_fem.md`](methods_supplement_fem.md) ·
references → [`CITATION_AUDIT.md`](CITATION_AUDIT.md).

---

## 0. The chain in one line

```
CLSM composition φ ──▶ ecology (TMCMC-calibrated) ──▶ growth field α(x)
                                                          │
                          stiffness E(φ) ──────────────────┤
                                                          ▼
                    FEM with growth UMAT  F = Fe·Fg,  Fg = (1+α)I
                                                          ▼
                       max von Mises σ  ──▶  ratio σ_CH / σ_DH
```

Four clinical conditions: **CH** commensal-HOBIC · **DH** dysbiotic-HOBIC ·
**CS** commensal-static · **DS** dysbiotic-static.

---

## 1. Ecology layer — who is there, and how it changes

**Five species** (early coloniser → pathogen):
*S. oralis, A. naeslundii, Veillonella, F. nucleatum, P. gingivalis*, with Σφᵢ = 1.

**Model.** A 0-D ODE on the volume fractions,

$$\dot{\boldsymbol\varphi} = f(\boldsymbol\varphi, \boldsymbol\psi; \boldsymbol\theta)$$

with **20 parameters** θ: the species **interaction matrix A**, viability ψᵢ,
decay bᵢ, and Hill-function gating for the *F. nucleatum → P. gingivalis*
pathway (the dysbiotic cascade).

**Calibration.** TMCMC (β: 0→1, 8 stages, 150 particles) against the measured
time series (6 timepoints × 5 species) → MAP + 95 % CI.

> **Provenance, stated precisely:** TMCMC calibrates **A and the rates** — *not*
> the composition. φ is **CLSM-measured** and enters as an input. The full 15-D
> inverse problem is under-identified; this separation is a deliberate, audited
> choice (see the rigor audit).

**For the FEM**, the TMCMC growth parameters map to a Monod rate,
`r ≈ maxᵢ(aᵢᵢ) / t_scale` (derivation in `methods_supplement_fem.md` §3).

---

## 2. Growth layer — composition to a mechanical driver

The growth variable α is obtained from a reaction–diffusion PDE (Klempt et al.,
published version — see `CITATION_AUDIT.md` before quoting equation numbers):

$$\frac{\partial \alpha}{\partial t} = \nabla\!\cdot\!\left(D\,\nabla\alpha\right) + k_\alpha\,\varphi_{\text{tot}}(\mathbf{x},t)$$

Implemented as a JAX testbed in [`JAXFEM/`](JAXFEM/README.md); the resulting
α-field is mapped to the integration points of the FE model.

---

## 3. Mechanics layer — growth to stress

**Kinematics — multiplicative growth split:**

$$\mathbf{F} = \mathbf{F}_e\,\mathbf{F}_g, \qquad \mathbf{F}_g = (1+\alpha)\,\mathbf{I}$$

Growth is **isotropic**. The viscoelastic extension inserts a viscous factor,
$\mathbf{F} = \mathbf{F}_e\mathbf{F}_v\mathbf{F}_g$, integrated by backward Euler.

**Elastic response — Mooney–Rivlin with a volumetric term:**

$$\Psi = C_{10}(\bar I_1 - 3) + C_{01}(\bar I_2 - 3) + \tfrac{1}{D_1}(J_e-1)^2$$

(`mtype = 0` reduces it to neo-Hookean, `C01 = 0`.)

**Composition → stiffness.** Two distinct scales — keep them apart in the text:

| Mode | Range | Meaning |
|---|---|---|
| **Substrate** | E(φ) ∈ [0.5, 10] GPa | effective stiffness of the biofilm-covered **periodontal attachment** (enamel/cementum/interface), mm scale |
| **Biofilm / EPS** | E ≈ 10–10³ Pa | the **EPS matrix itself** (Billings et al. 2015; the scale used by Klempt et al.) |

**Solve.** FEM (Abaqus/Standard) with the growth UMAT → max von Mises σ per
condition → the reported ratio.

---

## 4. What is measured / calibrated / assumed

The single most important table for the Limitations section.

| Quantity | Status | Basis |
|---|---|---|
| Composition **φ** | 🟢 **measured** | CLSM, 5 species × 6 timepoints |
| Interaction matrix **A**, rates | 🟢 **calibrated** | TMCMC, MAP + 95 % CI |
| Growth kinematics, tangent, mesh | 🟢 **verified** | consistent tangent vs FD ≈ 2.4–3.0e-8; patch tests 13/13; mesh convergence |
| Constitutive implementation | 🟢 **verified across 3 codes** | Abaqus UMAT ≡ ANSYS USERMAT (0 ULP) ≡ Python core (6.8e-14 rel.) |
| Growth magnitude **α** | 🟡 **magnitude-anchored** | thickness 1.5–3.5× → α ≈ 0.5–2.5; **not** calibrated per condition |
| Per-species stiffness **E_SPEC** | 🟡 **assumed** | order-of-magnitude only; needs species-level AFM |
| Depth / nutrient model | 🟢 **result robust to it** | 5.3–6.6× via 3-D reaction–diffusion |

---

## 5. Result and its uncertainty

**Headline:** σ(commensal-HOBIC) / σ(dysbiotic-HOBIC) ≈ **6.44×** for early
biofilm, decaying to ≈ 2× at maturity. Reported as a **full σ(t) trajectory with
bands**, never as a single number.

| Sensitivity | Band | Reading |
|---|---|---|
| depth / nutrient model | **5.3 – 6.6×** | 🟢 robust |
| per-species stiffness `E_SPEC` | **3.7 – 12×** | 🟡 the dominant uncertainty |

**Independent validation.** The pipeline's dysbiotic-static composition versus
the Heine CLSM measurement: **MAE 4.22 pp, TVD 0.105**, Veillonella-dominant
structure reproduced; largest single miss *F. nucleatum* (10.5 pp).
(`validate_composition.py`, guarded by a regression test.)

---

## 6. Two analysis lineages — do not mix them in the write-up

1. **Klempt growth-stress** (the thesis headline) — everything above:
   `Fg=(1+α)I`, the verified UMAT, `JAXFEM/` α-field.
2. **DI-bridge FEM** (alternative bridging variable) — Dysbiosis Index → E(DI)
   with transverse isotropy: `E = E_max(1-r)^n + E_min·r`, `r = clip(DI/s, 0, 1)`.
   Absolute moduli are on a *nominal* GPa scale; the **ratios are
   scale-invariant**. See [`FEM_README.md`](FEM_README.md).

> Known property of lineage 2, already covered by a test: the blend is convex and
> dips ~1 % below `E_min` near `r ≈ 1 − E_min/(2E_max)`, i.e. E(DI) is *not*
> bounded below by `E_min` (`tests/test_di_e_mapping.py`, strict xfail).

---

## 7. Open scientific questions (the continuation)

1. **`E_SPEC` by measurement** — species-level AFM would convert the dominant
   sensitivity band into a data-backed number.
2. **Per-condition α calibration** — currently one magnitude for all conditions.
3. **Python-in-the-loop material model** — the equivalence is already proven
   (§4), so the calibrated JAX model can be swapped in at the Gauss point without
   changing the physics: `ansys_usermat/coupling/`.
4. **Jaw-level, uncertainty-aware pipeline** (Level 2) — the first-paper target.
