# Coupling (c) with synthetic data (20 Oct 2026)

`scheme_c_synthetic.py`, numbers in `scheme_c_synthetic.json`, figure `assets/fig_scheme_c_synthetic.pdf`.
Not part of the thesis. Seed problem of Chapter 4 (32^3 field, 16^3 stress, E = 10 Pa, nu = 0.49). Point model:
Fritsch et al. 2025 Case II, true parameters, own void, start 0.02 per species, c* = 25 c, s = 0.15 (T* = 1 = their
1500-step validation run), antibiotic alpha* = 50 from T* = 0.5 as in their Fig. 18. Growth (c):
d alpha/dt = kappa phi sum phi_i psi_i, kappa = 1.72e-3 /T* (equal seed growth at T* = 1 for uniform c).

| run | seed alpha-1, T* = 0.5 | T* = 1 | seed vM T* = 1 (Pa) | living fraction in the seed, T* = 1 |
|---|---|---|---|---|
| (a) Eq. 36 | 3.73e-4 | 6.56e-4 | 6.07e-4 | |
| (c) uniform c | 2.01e-4 | 6.56e-4 | 6.82e-4 (+12 %) | 0.96 |
| (c) antibiotic from T* = 0.5 | 2.01e-4 | 5.74e-4 (-12 %) | 5.93e-4 (-13 %) | 0.91 (dip to 0.49) |
| (c) local nutrient, Lambda = 4 | 1.09e-4 | 2.78e-4 (-58 %) | 1.93e-4 (-68 %) | 0.35-0.83 |

- Same final growth by the choice of kappa, but the synthetic biofilm starts sparse (0.08) and fills, so the
  growth comes late; the seed stress is 12 % higher than (a). In the calibrated runs (ANSYS, CH) the living
  fraction falls over time instead, and the stress was 7 % lower: the sign follows the time course of the
  living biomass.
- The antibiotic kills about half the living biomass at once; the point model recovers within 0.2 T*, so the
  growth is 12 % lower at the end (the calibrated ANSYS runs: 32 %, a longer exposure without regrowth).
- With the local nutrient the living biomass and the growth vary across the seed (more growth on the nutrient
  side); growth and stress drop to about a third, with the same kappa.

## Mesh check and kappa = k_alpha (20 Oct, later)

Field grid 48^3 instead of 32^3 (`--nf 48`, `scheme_c_synthetic_nf48.json`): seed alpha - 1 changes by 0.9 % in every
run, kappa by 0.3 %; the ratios do not change: living biomass vs Eq. 36 von Mises +11.5 % (32^3: +12.4 %); against the
uniform-c run of (c), antibiotic -12.4 % / -13.0 % (alpha - 1 / von Mises), local nutrient -57.6 % / -71.6 %
(32^3: -12.4 / -13.1, -57.6 / -71.7).

kappa = k_alpha of Klempt 2024 (`--kappa-ka`, nothing fitted, `scheme_c_synthetic_kappa_ka.json`): growth 42 % and
von Mises 35 % below Eq. 36 for uniform c; the antibiotic and nutrient ratios are the same as with the matched kappa
(to 1e-8: the growth is linear in kappa).
