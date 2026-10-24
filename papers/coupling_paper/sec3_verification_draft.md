# Section 3 draft: implementations and verification (IKMHIWI03, 23-24 Oct 2026)

Draft text for Section 3 of `OUTLINE.md`. Numbers from `ansys_usermat/apdl/IKMHIWI03_TO_CLOUD.ja.md` (A-46 to A-75) and
`ansys_usermat/apdl/CONFIRMED_2026-10-22.ja.md` (items 23-40). Example throughout: the separate example with antibiotic
(two species of Fritsch et al. 2025 case I, front term on, alpha*_0 = 100 from T* = 1/3), seed means at T* = 1 unless
stated. Values in [..] wait for the last runs (ANSYS 24^3 with the upwind gradients). Wording is neutral ("this work"); the
author list is open.

## 3.1 Three implementations

The coupled model of Section 2 is solved in three codes that share the equations and the parameters but not the
discretisation:

| | Python | ANSYS (user element) | Abaqus (user element) |
|---|---|---|---|
| phi_3D (Eq. 34) | finite volumes on cells, RK4 | neighboured element method (MLS operators on the Gauss points, 30 neighbours), explicit Euler | Galerkin hexahedra (8 nodes, 2x2x2 Gauss points), implicit Euler, Newton |
| c, a (Eq. 35) | 7-point Laplacian, sparse direct | NEM Laplacian, PARDISO | Galerkin, same Newton system |
| mechanics | hexahedra, selective reduced integration | SOLID185 with the growth USERMAT | B-bar hexahedra in the same element |
| front term | upwind (TVD, minmod) | central MLS gradient + stabilisation (below) | Galerkin + stabilisation (below) |
| point model | the published one (Klempt et al. 2026) through the same Fortran routine in all three |

The point model is the same compiled routine in all three, so differences come from the field and the mechanics only.

## 3.2 Code verification

- **Point model in the field.** With a uniform field the coupled run reproduces the point model alone to 1e-15 (Python,
  two and four species) and 1.8e-12 (ANSYS); the hand-over of sum_i phibar_i to the consumption is exact.
- **Manufactured solutions.** Sources for phi and c built with JAX automatic differentiation and checked against
  sympy (agreement 1e-14). The Python and Abaqus solutions converge at second order in L2 (32^3: phi error 2.0e-4 and
  2.9e-4) for the stabilised equation. For the streamline-upwind form (Section 3.4) the order is about 1.5; the
  largest error sits in a thin layer next to the wall opposite the supply, where grad c -> 0 and the front direction
  grad c/|grad c| changes over a few cells; halving the time step does not change it.
- **Tangent.** The analytic tangent of the Abaqus element agrees with central differences to 5e-11, with and
  without the streamline term.
- **The element next to the walls.** With a manufactured solution the element does not converge from 8^3 to 16^3
  (phi error 2.2e-3 to 1.8e-3, c 8.2e-3 to 1.2e-2): the MLS stencils of the Gauss points next to a face include the
  surface points and are ill-conditioned (error of the Laplacian of the manufactured c: 0.88 on 8^3, 65 on 16^3, against
  |lap c| = 3.5). Rebuilding these stencils with the neighbours mirrored across the faces (zero flux built in) makes
  the operators consistent up to the walls, and with the weight width proportional to h the element converges (phi
  9.5e-4 to 3.4e-4, order 1.5; c 9.0e-4 to 2.2e-4, order 2.0). In the coupled example the mirrored operators run on
  8^3 but become unstable at a corner of the 16^3 mesh after the antibiotic onset; they are not used for the results.
- **Mechanics alone.** The last growth field of an Abaqus run, solved by the Python mechanics, gives the seed von
  Mises stress to 1.1 % (no obstacles) and 0.5 % (obstacles with frictionless contact). Differences in the stress
  between implementations therefore come from the growth field, not from the mechanics.

## 3.3 Solution verification

- **Time step.** Abaqus, 16^3, Delta t = 0.005, 0.00125, 0.000625: observed order 1.1; at Delta t = 0.00125 the error
  against the extrapolated value is 0.8 % in alpha - 1 and 2.2 % in von Mises. All mesh studies use 0.00125
  (ANSYS 24^3 with the constant stabilisation: 0.0005, its explicit stability limit).
- **Mesh.** Grid convergence index after Celik et al. (2008), Fs = 1.25 for three grids in the asymptotic range,
  Fs = 3 for two grids. h = 2 mm / n.

  Seed alpha - 1 [1e-4]:

  | scheme | 16^3 | 24^3 | 32^3 | 48^3 | 64^3 | extrapolated, GCI |
  |---|---|---|---|---|---|---|
  | Python, constant r h/4 | 6.27 | 7.26 | 7.72 | 8.20 | 8.43 | 9.0 (8.5 %) |
  | Python, TVD | 8.25 | | 8.82 | | 9.17 | 9.8 (8.1 %), 128^3: 9.33 |
  | Python, LW + sensor | 7.95 | 8.52 | 8.80 | 9.08 | 9.21 | 9.5 (3.8 %) |
  | Abaqus, constant r h/4 | 5.79 | 7.03 | 7.60 | 8.10 | | 8.8 (10 %) |
  | Abaqus, SUPG | 8.55 | 8.77 | 8.94 | 9.11 | | 9.45 (4.6 %) |
  | ANSYS, constant r h/4 | 5.73 | 6.99 | | | | 9.5 (two grids, order 1 assumed) |
  | ANSYS, LW + sensor | 8.19 | unstable | | | | |
  | ANSYS, upwind MLS gradients | 8.01 | [..] | | | | 8^3: 6.38 |

  Seed von Mises [1e-4 Pa]:

  | scheme | 16^3 | 24^3 | 32^3 | 48^3 | 64^3 | remark |
  |---|---|---|---|---|---|---|
  | Python, constant r h/4 | 1.73 | 2.78 | 3.40 | 4.29 | 4.80 | 96^3: 5.45, not in the asymptotic range |
  | Python, TVD | 6.51 | | 5.90 | | 5.98 | settled from 32^3 |
  | Python, LW + sensor | 4.20 | 5.31 | 5.85 | 6.28 | 6.26 | settled from 48^3 |
  | Abaqus, constant r h/4 | 1.52 | 2.89 | 3.75 | 4.62 | | not in the asymptotic range |
  | Abaqus, SUPG | 6.17 | 6.25 | 6.27 | 6.22 | | GCI 1.2 % at 48^3 |
  | ANSYS, constant r h/4 | 1.33 | 2.52 | | | | not in the asymptotic range |
  | ANSYS, LW + sensor | 4.20 | unstable | | | | at the corners of the supply face |
  | ANSYS, upwind MLS gradients | 4.56 | [..] | | | | 8^3: 2.14 |

  The growth variable converges for every scheme (the implementations agree within their GCI bands). The stress does
  not converge with a constant added diffusion: r h/4 smears the seed edge over a few cells, and the stress, which is
  set by the jump of alpha across that edge, follows the smearing. At 16^3 it is a quarter of the converged value
  (5.9 to 6.3e-4 Pa from the two second-order schemes of two codes).
- **Code-to-code at equal resolution.** With the same added diffusion (beta = 0.176 mm^2/T*), ANSYS on 16^3 elements
  (its Gauss points are h/2 apart) and Python on 32^3 cells differ by 1 % in alpha - 1 and 5 % in von Mises; with the
  sensor-switched form ANSYS 16^3 matches Python 16^3 of the same form (8.19 vs 7.95e-4, 4.20 vs 4.20e-4 Pa): the MLS
  operators over 30 neighbours give the element a resolution of about one element, not one Gauss-point spacing.

## 3.4 Stabilisation of the front term

The front term of Eq. 34 is an advection with speed v = r c/(K + c) along grad c / |grad c| (r = 10 mm/T* on the 2 mm
cube) against a small diffusion (beta = 0.02 mm^2/T*): the cell Peclet number r h/(2 beta) is about 31 on 16^3, so a
central discretisation oscillates. Four treatments were compared:

1. constant added diffusion r h/4 (first order; the element as delivered needs it to stay stable);
2. TVD upwinding (MUSCL, minmod; Python);
3. SUPG (Brooks and Hughes 1982) with tau = h/(2|v|)(coth Pe - 1/Pe) on the strong residual (Abaqus; consistent, no
   added diffusion away from the residual);
4. for the explicit element: the Lax-Wendroff term (Delta t/2) v . grad(v . grad phi) of the front plus an added
   diffusion switched by an extremum sensor s = (|sum_j L_ij (phi_j - phi_i)| / sum_j |L_ij (phi_j - phi_i)|)^2 on the
   conservative Laplacian (face value max(s_i, s_j); s = 1 at points whose stencil reaches the boundary). Without the
   sensor the Lax-Wendroff form oscillates (phi from -0.74 to 1.88 on 16^3); with it phi stays in [0, 1] up to
   -0.03/1.01 on 64^3.

None of these changes Eqs. 34-36. Treatments 2 to 4 give the same stress within a few per cent from 32^3 on; treatment
1 does not converge in the stress up to 96^3. For the stress the treatment of the front term matters more than the
mesh.

## 3.5 Where this leaves the examples

- Stresses are quoted from the second-order schemes (Python TVD, Abaqus SUPG) with their GCI; the element's values
  are given with their mesh error (16^3 and 24^3); the LW + sensor form is unstable on 24^3, the upwind MLS gradients
  are stable on 16^3 [24^3 pending].
- The growth variable and the composition (sum_i phibar_i, the shares) are insensitive to the stabilisation at the
  level of the GCI and can be quoted from any of the three codes.
- The case 4.2 of Klempt et al. (2024) ("biofilm on nutrients", high consumption) has a cell Peclet number of about 31
  on 16^3 and about 25 on the paper's own mesh; the paper does not state a stabilisation. In the element it either
  oscillates (no stabilisation) or the thin seed at the wall diffuses away (r h/4, and the LW + sensor form, which uses
  the full added diffusion next to walls). Not reproduced; not used.

- The maze of Klempt et al. (2024, Sect. 4.4) in Python and in the Abaqus element (4 um cells, T* = 60): without the
  antibiotic the colonised area agrees within 3 to 5 % and the living fraction within 0.006. With the antibiotic the
  result depends on the rule for newly reached points: when they take the state of their neighbours the biofilm
  recovers (living fraction 0.12 to 0.35 in Python), when they start from the start state it does not (0.04 to 0.07
  in both codes). This rule is a modelling assumption and is stated with each example.

## To finish
- ANSYS 24^3 with the upwind gradients (running); the maze in Abaqus with the inherit rule (running).
- Decide which code carries which figure (proposal: maps from Python TVD at 64^3 or 128^3, the element for the
  verification table only).
- Figure 2 of the outline (verification): the two tables above as a two-panel plot (alpha - 1 and von Mises against
  1/n, one line per scheme and code).
