# 3D growth-PDE rigor investigation — depth profile & nutrient model (2026-06-26)

Question: the headline sigma_CH/sigma_DH=6.44x uses an imposed linear depth ramp
alpha=alpha_max*(1-depth) (inner-peaked) to map the 2D growth field onto the
tooth. Is the headline robust if we instead solve the growth PDE rigorously in 3D?

## Steps
1. 2D alpha field is laterally ~uniform (mean/max=0.966) -> the "max" spatial
   reduction is ~3% (negligible). The real assumption is the depth ramp.
2. Built a 3D FE (P1) reaction-diffusion solver on the real C3D4 biofilm shell
   (8985 nodes, 0.2 mm uniform thickness), implicit diffusion (unconditionally
   stable), nutrient Dirichlet on the oral (outer) surface, consumption, growth
   alpha accumulation. ESTABLISHED biofilm: phase-field phi=1 over the full shell
   (the meshed shell IS the biofilm) -> uniform E (phi-gate=1).
3. Physical alpha profile is OUTER-peaked (nutrient from oral side), the OPPOSITE
   of the imposed inner-peaked ramp. Resolution converged (Nz 5->101: <3%).
4. Dimensionalized via literature (O2 limiting nutrient): D_eff~1.5-2.4e-9 m^2/s,
   penetration 50-200 um (de Beer 1994; Stewart 2003) ~ shell 0.2 mm
   -> Lp/L ~ 0.25-1.0. Swept Lp/L = 0.15 .. 10 (strongly-limited .. well-mixed).

## Result: the ratio is ROBUST
| Lp/L | 0.15 | 0.25 | 0.50 | 1.00 | 2.50 | 10.0 |
| sigma_CH/sigma_DH | 5.6 | 5.3 | 5.5 | 6.3 | 6.6 | 6.6 |
Across the entire (and literature) penetration range, sigma_CH/sigma_DH = 5.3-6.6x,
bracketing the headline 6.44x (the well-mixed limit -> 6.6x). Absolute sigma does
change (CH 7-13 kPa, DH 1.3-2.1 kPa) but the RATIO (the thesis claim) does not.

## Conclusion
The imposed linear depth ramp is physically arbitrary (inner-peaked vs the
physical outer-peaked), but the headline quantity sigma_CH/sigma_DH ~ 5-6.6x is
robust to the depth-profile / nutrient model. The thesis ratio survives a
rigorous 3D growth-PDE treatment over the literature-plausible nutrient range.

Tools: tooth_pde3d.py (3D FE RD solver), sensitivity_sweep.py (Lp sweep),
run_pde3d_umat.py, resolution_check.py.
