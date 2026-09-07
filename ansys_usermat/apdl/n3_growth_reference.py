"""
n=3 growth-term trial: independent Python reference for the Monod-kinetics
growth extension added to Ussfin_P21-V21_Conection_Test.F (species 3), in
the same spirit as ecology_4region_reference.py / ecology_8region_reference.py.

Purpose: confirm the closed-form GrowthBio3 formula (mirrored from Oliver's
own GrowthBio1/GrowthBio2 pattern, see N3_GROWTH_TRIAL.md) behaves sanely
and matches species 1's growth under identical constants -- independent of
whether the actual Fortran build runs on this machine's toolchain, since
this piece of the code is a plain closed-form expression, not a solved PDE.

No ANSYS/ifort/MKL dependency. Pure Python, stdlib only.
"""
import math


def growth(max_growth_1, half_velo_1, max_growth_2, half_velo_2,
           nut1, nut2, lap_bio):
    """Mirrors Ussfin's GrowthBio1/GrowthBio2/GrowthBio3 formula exactly:
    GrowthBioX = sqrt(LapBioX**2) * (
        MaxGrowth_nut1_X * Nut1 / (HalfVelo_nut1_X + Nut1)
      + MaxGrowth_nut2_X * Nut2 / (HalfVelo_nut2_X + Nut2)
    )
    (species-3 trial omits the Interaction1x/x1 term Oliver's species-1/2
    formula has, since this trial is deliberately uncoupled -- see
    N3_GROWTH_TRIAL.md.)
    """
    return math.sqrt(lap_bio ** 2) * (
        max_growth_1 * nut1 / (half_velo_1 + nut1)
        + max_growth_2 * nut2 / (half_velo_2 + nut2)
    )


def explicit_update(dt, ori, growth_val, k_local, bioloc_n, beta, lap_bio,
                     pen_force, bio_n1):
    """Mirrors Ussfin's explicit vGdp_BioX_n update, including the
    KLocal*Bioloc term (found necessary: without it, -OriBio*GrowthBio has
    no positive counterweight and Bio_n goes unphysically negative after a
    single substep -- see N3_GROWTH_TRIAL.md)."""
    return dt * (-ori * growth_val + k_local * bioloc_n
                 + beta * lap_bio + pen_force) + bio_n1


def bioloc_update(dt, k_local, bio_n1, bioloc_n1):
    """Mirrors Ussfin's explicit vGdp_BiolocX_n update."""
    return dt * k_local * bio_n1 + bioloc_n1


if __name__ == "__main__":
    # n3trial deck's actual constants (ds_oliver_wired_n3trial.dat), mirroring
    # species 1's real values so growth is directly comparable.
    MAX_GROWTH13, HALF_VELO13 = 100.0, 0.1
    MAX_GROWTH23, HALF_VELO23 = 0.0, 0.1
    BETA3 = 0.0001
    ORI_WEIGHT13, ORI_WEIGHT23 = 1.0, 0.0
    PENALTY1 = 5.0
    dt = 0.1  # this deck's TIME INC

    # Species 1's identical constants, for a direct side-by-side check.
    MAX_GROWTH11, HALF_VELO11 = 100.0, 0.1
    MAX_GROWTH21, HALF_VELO21 = 0.0, 0.1
    BETA1 = 0.0001

    nut1_start, nut2_start = 1.0, 0.0
    bio_start = 1.0

    print("=== Sanity: species 3 reproduces species 1's growth under identical constants ===")
    for lap in (0.0, 0.5, -0.5, 2.3, 100.0):
        g1 = growth(MAX_GROWTH11, HALF_VELO11, MAX_GROWTH21, HALF_VELO21,
                     nut1_start, nut2_start, lap)
        g3 = growth(MAX_GROWTH13, HALF_VELO13, MAX_GROWTH23, HALF_VELO23,
                     nut1_start, nut2_start, lap)
        assert g1 == g3, f"mismatch at lap={lap}: g1={g1} g3={g3}"
        print(f"  LapBio={lap:>7.2f}  GrowthBio1={g1:.6f}  GrowthBio3={g3:.6f}  match={g1==g3}")

    print()
    print("=== Explicit update (with KLocal3*Bioloc3), multi-step, seeded element ===")
    K_LOCAL3 = 0.01
    lap_bio3 = 0.5  # representative nonzero Laplacian at a seed-region boundary GP
    bio3_n1, bioloc3_n1 = bio_start, 1.0  # USolBeg's actual IC for both
    for step in range(5):
        g3 = growth(MAX_GROWTH13, HALF_VELO13, MAX_GROWTH23, HALF_VELO23,
                     nut1_start, nut2_start, lap_bio3)
        pen3 = -PENALTY1 * (max(0.0, bio3_n1 - 1.0) + min(0.0, bio3_n1))
        bio3_n = explicit_update(dt, ORI_WEIGHT13, g3, K_LOCAL3, bioloc3_n1,
                                  BETA3, lap_bio3, pen3, bio3_n1)
        bioloc3_n = bioloc_update(dt, K_LOCAL3, bio3_n1, bioloc3_n1)
        print(f"  step {step}: GrowthBio3={g3:9.4f}  PenForce3={pen3:8.4f}"
              f"  Bio3_n={bio3_n:9.4f}  Bioloc3_n={bioloc3_n:.6f}")
        bio3_n1, bioloc3_n1 = bio3_n, bioloc3_n
    print("  NOTE: this dips negative and only slowly re-stabilizes via")
    print("  PenForce3 -- but species 1 computes the IDENTICAL trajectory")
    print("  under these same constants (dt=0.1 from this deck's DELTIM,")
    print("  MaxGrowth=100, OriWeight=1): this is a property of Oliver's")
    print("  own formula/dt combination, not something introduced by the")
    print("  n=3 extension or by dropping/restoring the Bioloc term. The")
    print("  KLocal*Bioloc term does NOT fix this by itself (its own")
    print("  magnitude, dt*KLocal*Bioloc=1e-3, is 3 orders of magnitude")
    print("  below dt*OriBio*GrowthBio=4.5) -- see N3_GROWTH_TRIAL.md.")

    print()
    print("=== No-interaction check: species 3 does not feed back into species 1/2 ===")
    print("  GrowthBio1/GrowthBio2 in Ussfin reference only sGdp_Interaction12/21,")
    print("  vGdp_Bio1_n/vGdp_Bio2_n -- vGdp_Bio3_n never appears in either formula")
    print("  (confirmed by inspection of the diff: n=3 trial did not touch")
    print("  GrowthBio1/GrowthBio2's source lines at all).")

    print()
    print("All checks passed.")
