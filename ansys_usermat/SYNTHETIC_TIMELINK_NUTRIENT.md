# Time link and local nutrient with synthetic data (20 Oct 2026)

Script: `synthetic_timelink_nutrient.py`; numbers: `results_synthetic_timelink_nutrient.json`;
figure: `assets/fig_synthetic_timelink_nutrient.pdf`. Not part of the thesis.

Point model: Fritsch et al. 2025, Case II M3, true parameters (four species, start 0.02 each, c* = 25,
750 steps of 1e-4). The "truth" is this free point model; it agrees with the verified reproduction
(`klempt2026_reproduction.simulate`) to four digits (end shares 0.262, 0.277, 0.322, 0.139; amount 0.753).
Coupling: the one-way scheme of Chapter 4 (rescaled to the field amount min(phi_3D, 0.9) every
dT* = 0.025), phi_3D(T*) from the seed point of the Chapter 4 seed problem (16^3; the field does not
depend on the point model).

## Time link s

The synthetic data window is 0.075 in point-model time, so T* = 1 spans it when **s = 0.075** (instead of
the assumed 0.15 of the thesis, which was set from the 1500-step horizon of the Klempt et al. 2026 cases).

| | RMS of the shares over the 20 data times | species 4 at the end (truth 0.139) |
|---|---|---|
| point model keeps its own amount, s = 0.075 | 7e-7 | 0.139 |
| rescaled to the field, s = 0.075 | 0.019 | 0.108 |
| rescaled to the field, s = 0.15 | 0.049 | 0.072 |
| rescaled to the field, best s = 0.05 | 0.005 | |

- With its own amount the coupled scheme reproduces the truth exactly, so the clock is right.
- The remaining error at s = 0.075 comes from the amount: the field holds the point model at 0.9 to 0.5
  from the start, while the free point model fills from 0.08 to 0.75. A fuller point model competes
  faster, so species 4 loses more.
- Fitting s to the data would pick 0.05, a third lower than the window: a fitted s absorbs the amount
  error. s should be set by the time window, and the amount handled separately.

## Local nutrient c* = 25 c(x)

8^3 seed problem, nutrient held on one face, s = 0.075. Range of the shares in the seed at T* = 1:

| Lambda | c in the seed | species 3 | species 4 | (c = 1: 0.337 / 0.108) |
|---|---|---|---|---|
| 1 | 0.86-0.91 | 0.330-0.332 | 0.117-0.122 | |
| 3 | 0.38-0.57 | 0.294-0.310 | 0.159-0.187 | |
| 4 | 0.23-0.44 | 0.280-0.300 | 0.177-0.210 | |

A lower local nutrient lowers c* and slows the competition: species 4 keeps up to 0.10 more share where
the nutrient is low, no new parameter.

## What follows

Clock (s from the window) and nutrient (c* = 25 c) can be linked without new parameters. The void cannot:
under the one-way scheme the point model's amount comes from the field (diffusion and k_alpha alpha), not
from its own growth. Making the two amounts agree is the step to coupling (c) (growth from the point
model), the Keio continuation.
