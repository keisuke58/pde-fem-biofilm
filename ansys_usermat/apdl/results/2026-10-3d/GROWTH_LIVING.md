# Growth from the living biomass (Keio continuation; summary of 15 Nov 2026)

Not part of the LUH thesis (the model there is frozen; the thesis names this in its outlook and limitations).
Runs of 29 Oct to 2 Nov on fifa (`CLOUD_TO_FIFA.ja.md`, `KEIO_CHAIN_LOG.md`), collected here; time histories of 15 Nov.

## The two forms

Klempt 2024 Eq. 36 lets every biomass grow, living or dead: $\dot\alpha=k_\alpha\,\phi_{3D}$. Two variants in which only
living biomass grows (`coupled_living_consumption.py`):

- **A, from the energy** (`--growth-living-energy`): the coupling term $-k_\alpha\,\alpha\,\phi$ of the energy (Eq. 20) becomes
  $-k_\alpha\,\alpha\,\phi\sum_i\bar\phi_i$, with $\sum_i\bar\phi_i$ a given field of the point model (not varied). The
  stationarity conditions then give $\dot\alpha=k_\alpha\,\phi_{3D}\sum_i\bar\phi_i$ and the source $k_\alpha\,\alpha\sum_i\bar\phi_i$
  in Eq. 34.
- **B, Eq. 36 only** (`--growth-living`): $\dot\alpha=k_\alpha\,\phi_{3D}\sum_i\bar\phi_i$, Eq. 34 unchanged (a viscosity
  $\eta_\alpha/\sum_i\bar\phi_i$). Growth of the active biomass only is the standard assumption of biofilm models since
  Wanner and Gujer (1986).

A and B agree within 0.1 % in every run (the source $k_\alpha\alpha$ of Eq. 34 is small against the front term
$-r\,\frac{c}{k+c}\,\mathbf n_c\cdot\nabla\phi_{3D}$). A is the one to carry forward, because it is derived.
$k_\alpha$ is kept at the value of Klempt 2024; nothing is fitted.

## Seed means at $T^*=1$

Two species (Fritsch case I), antibiotic from $T^*=1/3$ at $\alpha^*_0=100$, minmod, outer step 0.0025, $\dot c$ kept.
$\alpha-1$ in $10^{-4}$, von Mises in $10^{-4}$ Pa. Eq. 36 on $128^3$ with antibiotic: values of the thesis table.

| example | grid | Eq. 36: $\alpha-1$, vM | living (A): $\alpha-1$, vM |
|---|---|---|---|
| separate | 32^3 | 8.729, 7.440 | 7.963, 6.944 |
| separate | 64^3 | 9.014, 6.730 | 8.226, 6.254 |
| separate | 128^3 | 9.110, 6.585 | 8.313, 6.117 |
| separate, antibiotic | 32^3 | 8.455, 5.529 | 4.355, 4.417 |
| separate, antibiotic | 64^3 | 8.933, 6.003 | 4.600, 4.925 |
| separate, antibiotic | 128^3 | 9.112, 6.284 | 4.688, 5.251 |
| mixed | 32^3 | 7.461, 2.129 | 6.414, 1.871 |
| mixed | 64^3 | 8.080, 3.079 | 6.992, 2.734 |
| mixed | 128^3 | 8.350, 3.649 | 7.243, 3.266 |
| mixed, antibiotic | 32^3 | 7.422, 1.998 | 2.761, 1.315 |
| mixed, antibiotic | 64^3 | 8.041, 2.816 | 2.841, 1.930 |
| mixed, antibiotic | 128^3 | 8.316, 3.298 | 2.866, 2.309 |

## What the antibiotic does (with against without, same grid and law)

| example | grid | Eq. 36: $\alpha-1$, vM | living (A): $\alpha-1$, vM |
|---|---|---|---|
| separate | 32^3 | -3.1 %, -25.7 % | -45.3 %, -36.4 % |
| separate | 64^3 | -0.9 %, -10.8 % | -44.1 %, -21.3 % |
| separate | 128^3 | 0.0 %, -4.6 % | -43.6 %, -14.2 % |
| mixed | 32^3 | -0.5 %, -6.2 % | -57.0 %, -29.7 % |
| mixed | 64^3 | -0.5 %, -8.5 % | -59.4 %, -29.4 % |
| mixed | 128^3 | -0.4 %, -9.6 % | -60.4 %, -29.3 % |

- With Eq. 36 the antibiotic does not change the growth; with growth from the living biomass it lowers the growth of the
  seed by 44 % (separate) and 60 % (mixed).
- The stress falls by less than the growth: 14 % (separate) and 29 % (mixed) on $128^3$. The seed stress follows the
  difference in growth between the seed and its surroundings, and the antibiotic slows both.
- In the separate example the size of the stress effect still depends on the grid (36, 21, 14 %), as it does with Eq. 36
  (26, 11, 5 %); the difference between the two laws is steady at about 10 points. In the mixed example the effect is steady
  at 29 to 30 %, about 20 points more than with Eq. 36.
- Without antibiotic the living law gives 9 % (separate) and 13 % (mixed) less growth and 7 and 11 % less stress than Eq. 36,
  because the living volume fraction starts near 0.6 and is below 1 throughout. The same $k_\alpha$ therefore means a
  slightly slower growth; a comparison at equal final growth would need $k_\alpha$ rescaled (not done).

## Over time (`assets/fig_growth_living_over_time.png`, `fig_growth_living_over_time.py`)

$64^3$, outer step 0.025, snapshots every 0.025. After the onset the living volume fraction of the seed falls to 0.15
(separate) and 0.04 (mixed). With Eq. 36 $\alpha-1$ keeps rising at the same rate; with the living law it bends over
(separate) or stops (mixed, from about $T^*=0.4$). The stress of the two laws separates before the onset already and
stays apart by a nearly constant amount afterwards. Seed values at $T^*=1$ of these runs: separate $\alpha-1$ 8.95 against
4.69, von Mises 5.69 against 4.63; mixed 8.30 against 2.96, 3.23 against 2.15.

## Open

- Which data could fix the size of the effect: species-resolved survival under chlorhexidine
  (`references/CHLORHEXIDINE_VS_POINT_MODEL.md`).
- Combination with a stiffness that depends on composition (`--stiff-species`, runs of 2 Nov) and with the cervical window.
- The antibiotic effect on the stress in the separate example is not grid converged; no further grid runs are planned.
