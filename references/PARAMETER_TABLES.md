# Parameter tables of the papers, for quick reference (21 Oct 2026)

Transcribed from the text extractions in this folder (`*.txt`). Superscripts are lost in the extractions
("10 10" = 10^10, "1013" = 10^13); every value marked (check) should be read off the PDF before it goes into the
thesis. Only parameters are copied here, no derivations (Felix's confidential material is not in this folder).

## Klempt et al. 2024 (BMMB), the growth model this work follows

### Table 2: parameters of all simulations unless stated otherwise (test case 1: 20 x 20 x 20 um cube, 1 um mesh)

| parameter | symbol | value | units |
|---|---|---|---|
| shear modulus (E = 10 Pa, nu = 0.49) | mu | 3.3557 | Pa |
| nutrient diffusivity | d | 10^10 | um^2 / T* |
| phase-field regularisation | beta | 2 | um^2 / T* |
| growth factor | k_alpha | 10^-3 | 1 / T* |
| half-velocity constant | k | 1 | - |
| consumption parameter | g | 10^8 | 1 / T* |
| growth parameter (front term) | r | 100 | 1 / T* as printed (used as um / T* in this work) |

Equations (as printed): Eq. 34 `phi_dot - beta lap phi - k_alpha alpha + |grad phi| r c/(k+c) n_gradphi . n_gradc = 0`,
Eq. 35 `c_dot - d lap c + g phi = 0` (g phi is a misprint for g phi c), Eq. 36 `alpha_dot - k_alpha phi = 0`.
Computations: eta_c = 1e-10 ("steady state for the nutrients"). Time step: 10^3 substeps nominal, bisection up to 10^6.

### Table 6: maze simulation (all other parameters as Table 2)

| parameter | value | units |
|---|---|---|
| nutrient diffusivity d | 10^13 (check) | um^2 / T* |
| beta | 10 | um^2 / T* |
| k_alpha | 10^-6 | 1 / T* |
| r | 200 | 1 / T* |

### Converted to the 2 mm cube of this work (lengths x 100, T* unchanged, d/L^2 and beta/L^2 kept)

| quantity | Table 2 (20 um cube) | 2 mm cube | used in this work |
|---|---|---|---|
| beta | 2 um^2/T* | 0.02 mm^2/T* | 0.02 mm^2/T* (decided 5 Oct) |
| d | 10^10 um^2/T* | **10^8 mm^2/T*** | quasi-static solve, exact limit (21 Oct); d = 1 only as slow-diffusion sensitivity |
| relaxation time of c over the cube, L^2/d | 4e-8 T* | 4e-8 T* | (quasi-steady) |
| g | 10^8 1/T* | 10^8 1/T* | via the Thiele number below |
| Thiele number g L^2/d, L = cube edge | 4 | 4 | |
| Thiele number with the half edge (l = 1 mm here, 10 um there) | 1 | **1** | **1** since 21 Oct evening (the user: as the paper); earlier runs used 4 (Oliver's example deck), superseded |
| r | 100 | 10 mm/T* (as a speed) | **10 mm/T* since 21 Oct evening (paper values only)**; r/10 earlier, superseded |
| k_alpha | 10^-3 1/T* | 10^-3 1/T* | 10^-3 |
| E, nu | 10 Pa, 0.49 | same | same |

## Klempt et al. 2026 (AAM), the point model (used as published)

### Table 1: two species, c* = 100, alpha* = 10 constant

| | Case 1 | Case 2 | Case 3 | Case 4 | Case 5 | Case 6 |
|---|---|---|---|---|---|---|
| a11 | 2 | 1 | 1 | 1 | 1 | 1 |
| a12 | 0 | 0 | 1 | 0 | 0 | -1 |
| a22 | 1 | 1 | 1 | 1 | 1 | 1 |
| b1 | 0 | 0 | 0 | 1 | 1 | 0 |
| b2 | 0 | 0 | 0 | 2 | 2 | 0 |
| eta1 [kg/(m Time)] | 1 | 1 | 1 | 1 | 1 | 1 |
| eta2 | 1 | 2 | 2 | 2 | 2 | 2 |
| initial phi1 | 0.2 | 0.2 | 0.2 | 0.2 | 0.25 | 0.2 |
| initial phi2 | 0.2 | 0.2 | 0.2 | 0.3 | 0.3 | 0.2 |

### Table 2: four species; A and eta fixed (Eq. 19, 20)

A = 1/2 [[1, 5, 5, 5], [5, 1, 3, 3], [5, 3, 1, 2], [5, 3, 2, 1]] (check, extracted layout), eta = diag(0.8, 1.0, 1.5, 2.0) kg/(m Time).

| | Case 1 | Case 2 | Case 3 | Case 4 |
|---|---|---|---|---|
| b1..b4 | 0.4, 0.3, 0.2, 0.1 | 0.4, 0.3, 0.2, 0.1 | 0.4, 0.3, 0.2, 0.1 | 10, 2, 1, 0.01 |
| initial phi1..phi4 | 0.02 each | 0.02, 0.02, 0.02, 0.2 | 0.02 each | 0.02 each |
| c* | 100 | 100 | 50 + 50 sin(500 t) | 100 |
| alpha* | 10 | 10 | 10 | 0, then 100 for t > 500 |

### Table 3 (appendix): two species, varied viscosities, c* = 100, alpha* = 10

| | Case 4 | Case 5 | Case 4 | Case 5 |
|---|---|---|---|---|
| a11, a12, a22 | 1, 0, 1 | 1, 0, 1 | 1, 0, 1 | 1, 0, 1 |
| b1, b2 | 1, 2 | 1, 2 | 1, 2 | 1, 2 |
| eta1, eta2 | 2, 1 | 2, 1 | 1, 1 | 1, 1 |
| initial phi1, phi2 | 0.2, 0.3 | 0.25, 0.3 | 0.2, 0.3 | 0.25, 0.3 |

## Fritsch et al. 2025, synthetic case I (the data of the coupled examples)

From `JAXFEM/fritsch2025_cases.py` (transcribed earlier from the paper's Tables 1-2; the paper is not in this folder):
two species, A = [[1, 0.1], [0.1, 1]], b = (1, 2), eta = (1, 2), initial phi = (0.25, 0.30), c* = 100,
alpha* = 10 (this work's examples use alpha*_0 = 50, an assumption), dt = 1e-4, 1000 steps (T* = 1 spans them, s = 0.1),
theta* = [1, 0.1, 1, 1, 2]. Four-species case II: `fritsch_4sp_M1..M3` in the same file.

## Soleimani et al. 2023 (Sci. Rep.), co-aggregation, Table 1 (2D, 100 x 100 um, mesh 0.5 um)

| description | symbol | value | units |
|---|---|---|---|
| nutrient diffusivity | D_n | 10^-6 (check) | um^2 / Time |
| bacterial density threshold | phi_cri | 0.1 | ug / um^3 |
| consumption rate | R_c | 50 | um^3 / (ug Time) |
| co-aggregation coefficient | R_alpha | 15 | um / Time |
| co-aggregation coefficient | tau | 1 | - |
| growth coefficients | R_s1, R_s2 | 500, 500 | um^3 / (ug Time) |
| phase-field coefficient for alpha | M | 10^-4 | 1 / Time |
| phase-field coefficient for phi_1, phi_2, c | M | 0 | 1 / Time |
| phase-field diffusivity | eps^2 (D_b) | 10^-8 (check) | um^2 / Time |

## Feng et al. 2021 (BMB), symbiotic two-species biofilm, Table 1 (2D, 300 x 300 um)

| quantity | symbol | value | units | their source |
|---|---|---|---|---|
| saliva concentration at the top | s1 | 10 | kg/m^3 | estimated |
| diffusion coefficient of saliva | D1 | 5e-10 | m^2/s | Corbin et al. 2011 |
| diffusion coefficient of lactic acid | D2 | 1e-9 | m^2/s | Ribeiro et al. 2005 |
| densities (S. gordonii, Veillonella) | rho_s, rho_v | 1100 | kg/m^3 | from Duddu et al. 2009 |
| yields | Y_s1, Y_v1, Y_v2 | 0.1 | - | estimated |
| lactic acid production | lambda | 0.5 | - | estimated |
| maximum growth rates | mu1, mu2 | 3e-5, 8e-5 | 1/s | Rath et al. 2017, estimated |
| Monod constants | k11, k21, k22, kl | 0.878, 0.878, 0.0878, 1 | kg/m^3 | Rath et al. 2017, estimated |

(Further rows of their Table 1 continue in the text extraction after line 600.)

## What the coupled examples use (21 Oct 2026 evening, the user: "paper values only")

| value | used | source |
|---|---|---|
| beta, d, g, k_alpha, k, r, E, nu | Table 2 of Klempt 2024, converted (table above) | paper |
| A, b, eta, start phi, c*_0 = 100, alpha*_0 = 10, s = 0.1 | Fritsch case I (= Klempt 2026 Table 1 style) | paper |
| antibiotic from **T* = 1/3** (decided 21 Oct, the user) | Klempt 2026 case 4 inserts it after one third ("t > 500" of 1500 steps); 0.5 (earlier runs) as sensitivity | paper |
| antibiotic diffusion and uptake | same as the nutrient | Oliver (20 Oct): the second nutrient field, same equation |
| antibiotic source (all outer faces), nutrient source (layer y = -1 mm) | set-up of the examples | Oliver's sketch ("defined sources") |
| **phi_min = 0.001** (start of a point model; decided 21 Oct: the converged value) | 0.1 / 0.03 / 0.01 / 0.003 / 0.001 change the seed values by at most 0.018 (c), 0.003 to 0.001 by 0.001 | no paper value; converged |
| mixed start with equal amounts | Oliver (20 Oct) | meeting |
