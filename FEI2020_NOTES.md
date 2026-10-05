# Fei et al. (2020): growth, friction and wrinkling of V. cholerae biofilms

> Fei C, Mao S, Yan J, Alert R, Stone HA, Bassler BL, Wingreen NS, Košmrlj A.
> **Non-uniform growth and surface friction determine bacterial biofilm
> morphology on soft substrates.** *PNAS* **117**(14):7622-7632 (2020).
> doi:[10.1073/pnas.1919607117](https://doi.org/10.1073/pnas.1919607117)
> (free at PMC7148565; arXiv:1910.14645). Data and code:
> <https://github.com/f-chenyi/biofilm-mechanics-theory>

Read 2026-10-05. **The PDF is not bundled**: it is "published under the PNAS
license", not CC BY, and this repository is public. Link and notes only.

## Why it is relevant

The nearest published case in which **growth-induced stress of a biofilm is
modelled with the same kinematics as here and compared quantitatively with
experiment**. A candidate for validating the mechanics, which no oral-biofilm
data allow (`ROADMAP_TWO_WAY.md`, literature section).

## Model, set against this work

| | Fei et al. 2020 | this work (Klempt 2024 + growth UMAT) |
|---|---|---|
| kinematics | F = Fe · Fg | F = Fe Fg |
| growth | Fg = diag(g I_∥, 1): isotropic in plane, none in thickness | Fg = α I, isotropic |
| growth law | ∂g/∂t = k_g(c) g, Monod in the nutrient | α̇ = k_α φ (Eq. 36) |
| nutrient | ∂c/∂t = D∇²c − Q0 Je⁻¹ c/(K + c) (Eq. 1), Monod uptake | ċ = d∇²c − gφ (Eq. 35) |
| material | nearly incompressible neo-Hooke; results insensitive to the rheology | compressible neo-Hooke, ν = 0.49 |
| mechanics | plane-stress thin film, force balance ∇·(hσ) − ηv = 0 with substrate friction (Eq. 2) | 3D, div σ = 0, rigid-body constraints |
| solver | FEniCS, Lagrangian | ANSYS |
| species | one (V. cholerae) | one and two |

## What they compare with experiment

- Radial expansion **velocity profiles** v(r, t) on 0.7 % agar: three
  kinematic stages, matched quantitatively (Fig. 2D).
- Biofilm **leading angle** against friction, several agar concentrations
  (Fig. 2E).
- Where **wrinkles** first appear: at the rim on soft agar, at the centre on
  stiff agar (Fig. 1, 4).
- All parameters from experiment except the friction coefficient η and the
  uptake rate Q0, which are fitted to the velocity profiles.

## Points to use

1. **Precedent for Monod-limited growth driving the stress.** Growth confined
   to a nutrient-rich rim about 1 mm wide produces circumferential compression
   at the rim and isotropic compression in the centre. The same mechanism as
   the seed under compression in this work, at a larger scale.
2. **A quantitative validation target for the continuation at Keio**: the
   velocity profiles and wrinkle onset, with their data and code public. It
   would need a thin-film geometry and substrate friction, which the present
   models do not have.
3. **Friction with the substrate** is the constraint that creates the stress
   there; in the example model here the constraint is the soft surroundings of
   the seed. Worth one sentence in the thesis outlook, not more.

## Their data and code (looked at 5 Oct)

`git clone --depth 1 https://github.com/f-chenyi/biofilm-mechanics-theory` (91 MB).
**No licence file**, so nothing of it is copied into this repository; clone it
where needed.

| file | content |
|---|---|
| `paper data/Fig 2/velocity profile/velocity-profile-0.7-{1,6,16}h.mat` | measured radial velocity on 0.7 % agar: `r` (23 points, up to about 4.7 mm), `v_avg`, `v_std` (µm/min). **The quantitative target of the paper** |
| `paper data/Fig 1/kymograph/kymo-0.{4,7}.mat` | biofilm radius `Rout`, radii of the radial and zigzag patterns over time `tAll` (241 points, h) |
| `paper data/Fig 1/surface profiling/line-profile-0.6.mat` | height profile `z(dr1)`: leading angle, thickness, wrinkle wavelength |
| `paper data/Fig 2/leading angle/` | leading angles (3 biological replicates) and the fitted scaled friction |
| `biofilm_morphogenesis_circle_CLEANED.py` | the FEniCS model (377 lines), axisymmetric disc, dimensionless (Rc = 1, T = 6) |

What a comparison here would need: the radius `Rout(t)` and the velocity
profiles are the observables; matching them requires a thin film on a
substrate with friction (their Eq. 2), which neither the ANSYS example model
nor the Python field has. A task for the continuation, not for the thesis.
