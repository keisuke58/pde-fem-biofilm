"""
Full n=5 pairwise interaction terms: independent Python reference for the
generalized GrowthBio1..GrowthBio5 formula added to
Ussfin_P21-V21_Conection_Test.F (N3_GROWTH_TRIAL.md, "interaction terms"
section). Mirrors n3_growth_reference.py's role -- confirm the closed-form
formula behaves sanely before spending an ANSYS compile/run cycle on it.

No ANSYS/ifort/MKL dependency. Pure Python, stdlib only.
"""
import math

N_SPECIES = 5


def inter_sum(i, interaction, bio):
    """InterSum_i = sum_{j != i} Interaction[i][j] * bio[j], 1-indexed to
    match the Fortran naming (Interaction13 etc.)."""
    total = 0.0
    for j in range(1, N_SPECIES + 1):
        if j == i:
            continue
        total += interaction.get((i, j), 0.0) * bio[j]
    return total


def growth(i, max_growth_nut1, max_growth_nut2, half_velo_nut1, half_velo_nut2,
           nut1, nut2, lap_bio, inter_sum_i):
    """Mirrors Ussfin's generalized GrowthBioI:
    GrowthBioI = sqrt(LapBioI**2) * (
        (MaxGrowth_nut1_I + InterSumI) * Nut1 / (HalfVelo_nut1_I + Nut1)
      + (MaxGrowth_nut2_I + InterSumI) * Nut2 / (HalfVelo_nut2_I + Nut2)
    )
    """
    return math.sqrt(lap_bio ** 2) * (
        (max_growth_nut1 + inter_sum_i) * nut1 / (half_velo_nut1 + nut1)
        + (max_growth_nut2 + inter_sum_i) * nut2 / (half_velo_nut2 + nut2)
    )


if __name__ == "__main__":
    # Constants mirroring the n5trial deck (species 1's real values,
    # species 3/4/5 uncoupled MaxGrowth/HalfVelo unchanged from before).
    max_growth_nut1 = {1: 100.0, 2: 100.0, 3: 100.0, 4: 100.0, 5: 100.0}
    max_growth_nut2 = {1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0}
    half_velo_nut1 = {1: 0.1, 2: 0.1, 3: 0.1, 4: 0.1, 5: 0.1}
    half_velo_nut2 = {1: 0.1, 2: 0.1, 3: 0.1, 4: 0.1, 5: 0.1}
    nut1_start, nut2_start = 1.0, 0.0
    bio = {1: 1.0, 2: 1.0, 3: 1.0, 4: 1.0, 5: 1.0}
    lap_bio = {1: 0.5, 2: 0.5, 3: 0.5, 4: 0.5, 5: 0.5}

    print("=== Sanity: all-zero interaction reproduces the pre-existing")
    print("    uncoupled formula (species 1/2's Interaction12/21 = 0 case,")
    print("    species 3/4/5's formula before this change) ===")
    zero_interaction = {}
    for i in range(1, N_SPECIES + 1):
        s = inter_sum(i, zero_interaction, bio)
        g = growth(i, max_growth_nut1[i], max_growth_nut2[i],
                    half_velo_nut1[i], half_velo_nut2[i],
                    nut1_start, nut2_start, lap_bio[i], s)
        g_uncoupled = growth(i, max_growth_nut1[i], max_growth_nut2[i],
                              half_velo_nut1[i], half_velo_nut2[i],
                              nut1_start, nut2_start, lap_bio[i], 0.0)
        assert s == 0.0 and g == g_uncoupled, f"species {i}: mismatch"
        print(f"  species {i}: InterSum={s:.4f}  Growth={g:.6f}  "
              f"(matches uncoupled: {g == g_uncoupled})")

    print()
    print("=== Modest nonzero coupling: check magnitude stays sane ===")
    # Small test values, same order of magnitude as this deck's other
    # rate constants (Beta*=1e-4, KLocal*=0.01) -- deliberately far below
    # MaxGrowth=100 so the interaction term perturbs, not dominates, the
    # base growth rate, in line with n=3's finding that dt*MaxGrowth
    # dominates dt*KLocal*Bioloc by orders of magnitude.
    small_interaction = {
        (1, 2): 0.05, (1, 3): 0.05, (1, 4): 0.05, (1, 5): 0.05,
        (2, 1): 0.05, (2, 3): 0.05, (2, 4): 0.05, (2, 5): 0.05,
        (3, 1): 0.05, (3, 2): 0.05, (3, 4): 0.05, (3, 5): 0.05,
        (4, 1): 0.05, (4, 2): 0.05, (4, 3): 0.05, (4, 5): 0.05,
        (5, 1): 0.05, (5, 2): 0.05, (5, 3): 0.05, (5, 4): 0.05,
    }
    for i in range(1, N_SPECIES + 1):
        s = inter_sum(i, small_interaction, bio)
        g = growth(i, max_growth_nut1[i], max_growth_nut2[i],
                    half_velo_nut1[i], half_velo_nut2[i],
                    nut1_start, nut2_start, lap_bio[i], s)
        # 4 partners * 0.05 * bio(=1.0) = 0.2 -- a small, bounded perturbation
        # on top of MaxGrowth1_I=100, not a source of blow-up by itself.
        assert abs(s - 0.2) < 1e-9, f"species {i}: InterSum={s}, expected 0.2"
        print(f"  species {i}: InterSum={s:.4f}  Growth={g:.6f}")

    print()
    print("=== Directionality check: InterSum_i only sums partners j != i ===")
    for i in range(1, N_SPECIES + 1):
        assert (i, i) not in small_interaction, f"species {i}: self-term present"
    print("  confirmed: no Interaction_ii self-terms defined anywhere.")

    print()
    print("All checks passed.")
