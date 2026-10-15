"""calibrate_growth.py -- the growth rate L0 of the unified element chosen so that the seed's mean
alpha - 1 at Day 21 equals the measured volume increase of CH: V(21)/V(3) = 0.400/0.222 = 1.80
(Heine et al., aggregated volumes), read as J = alpha^3 = 1.80 from Day 1, alpha - 1 = 0.216. Assumptions:
the whole volume increase is growth (no detachment, no medium change), and the growth from Day 1 to
Day 3 is small. Then the stresses for measured stiffness values.

    BIOFILM_ECO_CASE=.../uel_eco.txt python calibrate_growth.py [n] [Y]

Result (15 Oct 2026, 8^3, finite strain, nu = 0.49, mean dilatation): with the yield Y = 1 the seed's growth has a
ceiling, alpha - 1 = 0.107 near L0 = 13 (faster growth eats the nutrient and starves the seed; L0 = 40 gives 0.052),
so the measured 1.80 cannot be reached. With a smaller yield (less nutrient per growth; c_stop kept at 0.1) it can:
Y = 0.25, L0 = 53.1: alpha - 1 = 0.2154, von Mises / E = 0.0794; Y = 0.1, L0 = 113.9: 0.2165, 0.0773. The stress
does not depend on the choice of Y within 3 %. Y and L0 are not from a paper.
"""
import sys

import emulate_abaqus as em

TARGET = 1.80 ** (1 / 3) - 1          # 0.2164


def run(L0, n, Y=1.0):
    import contextlib
    import io
    with contextlib.redirect_stdout(io.StringIO()):
        r = em.seed(n, finite=True, nu=0.49, md=True, L0=L0, Y=Y)
    return r


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    Y = float(sys.argv[2]) if len(sys.argv) > 2 else 0.25
    hist = []
    L0 = 20.0
    for it in range(8):
        r = run(L0, n, Y)
        hist.append((L0, r["alpha_seed"]))
        print(f"L0 {L0:.4g}: seed alpha-1 {r['alpha_seed']:.4f}, vM/E mean {r['vm_mean']:.4e}, max {r['vm_max']:.4e}, "
              f"p/E {r['p_mean']:.4e}, edge {r['edge']:.4f}", flush=True)
        if abs(r["alpha_seed"] / TARGET - 1) < 0.005:
            break
        if len(hist) == 1:
            L0 = L0 * TARGET / r["alpha_seed"]
        else:                                          # secant on log(alpha-1) against log(L0)
            import math
            (a0, y0), (a1, y1) = hist[-2], hist[-1]
            k = (math.log(y1) - math.log(y0)) / (math.log(a1) - math.log(a0))
            L0 = a1 * (TARGET / y1) ** (1 / k)
    ratio = r["vm_mean"] / r["alpha_seed"]
    print(f"\nL0 = {L0 if abs(r['alpha_seed'] / TARGET - 1) < 0.005 else hist[-1][0]:.4g}: seed alpha-1 {r['alpha_seed']:.4f} "
          f"(target {TARGET:.4f}), vM/E {r['vm_mean']:.4e}, vM/(E (alpha-1)) {ratio:.3f}")
    for E, src in ((10.0, "Klempt et al. 2024"), (600.0, "0.6 kPa"), (41000.0, "Pattem et al. 2018, upper end")):
        print(f"  E = {E:g} Pa ({src}): seed von Mises mean {r['vm_mean'] * E:.3g} Pa, max {r['vm_max'] * E:.3g} Pa, "
              f"pressure {r['p_mean'] * E:.3g} Pa")
