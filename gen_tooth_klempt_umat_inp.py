"""gen_tooth_klempt_umat_inp.py
============================================================
TMCMC + Klempt multi-species PDE → conformal tooth UMAT inp

Modes (--mode):
  A  : umat_klempt_voigt.f   (Option A+B: Voigt E + condition-specific α)
       PREDEF 1=α_total, 2..6=φ_i
  C  : umat_klempt2025.f     (Option C: per-species α_s + Voigt E)
       PREDEF 1..5=α_So..Pg, 6..10=φ_i

Pre-requisites:
  python JAXFEM/klempt_pde_multispecies.py   (generates α per condition)

Output
------
  p23_klempt_{mode}_{condition}.inp  (4 conditions × 2 modes = 8 files)

Run Abaqus
----------
  bash run_tooth_klempt.sh --mode A
  bash run_tooth_klempt.sh --mode C
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

import numpy as np

# ── paths ─────────────────────────────────────────────────────────────────────
HERE   = Path(__file__).resolve().parent
JAXFEM = HERE / "JAXFEM"
MSCL   = HERE / "_multiscale_2d_results"
TMPL   = HERE / "p23_biofilm_commensal_static.inp"
UMAT_DIR = (HERE.parent / "nife/masterarbeit_ansys_fem/coupling_prototype/abaqus")

UMAT = {
    "A": UMAT_DIR / "umat_klempt_voigt.f",
    "C": UMAT_DIR / "umat_klempt2025.f",
}

# ── mesh constants ─────────────────────────────────────────────────────────────
N_NODES           = 8985
N_VERTS_PER_LAYER = 1797   # 8985 / 5
N_LAYERS          = 4

# ── species ────────────────────────────────────────────────────────────────────
SPECIES    = ["So", "An", "Vd", "Fn", "Pg"]
E_SPEC_MPa = np.array([1e-3, 8e-4, 6e-4, 2e-4, 1e-5])   # [MPa] = Pa × 1e-6

# ── conditions ─────────────────────────────────────────────────────────────────
CONDITIONS = [
    "commensal_static",
    "commensal_hobic",
    "dysbiotic_static",
    "dysbiotic_hobic",
]

NU_KLEMPT = 0.49


def load_tmcmc() -> dict[str, np.ndarray]:
    """Load TMCMC φ_i (5-vector) per condition."""
    out = {}
    for cond in CONDITIONS:
        d = json.load(open(MSCL / f"ref_0d_{cond}.json"))
        out[cond] = np.array(d["phi_final"])
    return out


def load_alpha_condition(cond: str) -> np.ndarray:
    """Load condition-specific Klempt α_final from multi-species PDE (Option B)."""
    p = JAXFEM / f"klempt_alpha_final_{cond}.npy"
    if not p.exists():
        raise FileNotFoundError(
            f"{p} not found. Run: python JAXFEM/klempt_pde_multispecies.py")
    return np.load(p)


def load_phi_condition(cond: str) -> np.ndarray:
    """Load condition-specific Klempt ϕ_final (biofilm density field) from PDE."""
    p = JAXFEM / f"klempt_phi_final_{cond}.npy"
    if not p.exists():
        raise FileNotFoundError(
            f"{p} not found. Run: python JAXFEM/klempt_pde_multispecies.py")
    return np.load(p)


def depth_norm(node_id_0: int) -> float:
    """depth_norm ∈ [0,1]: 0=inner(tooth), 1=outer(oral). From layer index."""
    layer = node_id_0 // N_VERTS_PER_LAYER   # 0..4
    return layer / N_LAYERS                   # 0.0..1.0


def alpha_at_depth(depth: float, alpha2d: np.ndarray) -> float:
    """depth=0 (inner, attachment) → α_max; depth=1 (outer, planktonic) → 0."""
    return float(alpha2d.max() * (1.0 - depth))


def phi_local_at_depth(depth: float, phi2d: np.ndarray) -> float:
    """
    Local biofilm density ϕ at depth_norm.
    Used for Klempt Eq.20 ϕ²-gated stiffness: E_gated = (ϕ_local/ϕ_total)² × E_Voigt.
    peak of phi2d ≈ phi_total_condition (center of biofilm blob).
    depth=0 → full biofilm density (ϕ_max); depth=1 → 0 (planktonic).
    """
    return float(phi2d.max() * (1.0 - depth))


def parse_template(inp_path: Path) -> dict:
    """Parse template inp → node/element/nset blocks."""
    blocks = {"nodes": [], "elements": [], "nsets": [], "elsets": []}
    cur = None
    with open(inp_path) as f:
        for line in f:
            ls = line.strip().upper()
            if ls.startswith("*NODE") and "OUTPUT" not in ls:
                cur = "nodes"; blocks["nodes"].append(line.rstrip()); continue
            if ls.startswith("*ELEMENT,") or ls == "*ELEMENT":
                cur = "elements"; blocks["elements"].append(line.rstrip()); continue
            if ls.startswith("*NSET"):
                cur = "nsets"; blocks["nsets"].append(line.rstrip()); continue
            if ls.startswith("*ELSET"):
                cur = "elsets"; blocks["elsets"].append(line.rstrip()); continue
            if ls.startswith("*") and not ls.startswith("**"):
                cur = None
            if cur:
                blocks[cur].append(line.rstrip())
    return blocks


def write_mode_A(out_path: Path, blocks: dict,
                 alpha_nodes: np.ndarray, phi_nodes: np.ndarray,
                 phi_vec: np.ndarray, cond: str) -> None:
    """Mode A: umat_klempt_voigt.f — PREDEF 1=α (ramped), 2..6=φ_i, 7=ϕ_local (fixed)."""
    L = []
    ap = L.append
    E_voigt_mean = float(np.dot(phi_vec, E_SPEC_MPa))

    ap("*HEADING")
    ap(f" Klempt+Voigt UMAT (Option A+B) — condition: {cond}")
    ap(f" E_Voigt={E_voigt_mean:.3e} MPa  nu={NU_KLEMPT}  alpha_max={alpha_nodes.max():.4f}")
    ap(f" UMAT: umat_klempt_voigt.f  |  alpha ramped via *Field+AMPLITUDE")
    ap("**")

    for line in blocks["nodes"]:    ap(line)
    ap("**")
    for line in blocks["elements"]: ap(line)
    ap("**")
    for line in blocks["nsets"]:    ap(line)
    ap("**")

    _write_all_elems(ap, blocks)
    ap("**")

    # Material
    ap("*Solid Section, elset=ALL_ELEMS, material=KLEMPT_VOIGT")
    ap(",")
    ap("*Material, name=KLEMPT_VOIGT")
    ap("*User Material, constants=1")
    ap(f" {NU_KLEMPT:.4f}")
    ap("*Depvar")
    ap(" 5")
    ap("**  SDV1=s=1+alpha  SDV2=Je  SDV3=alpha  SDV4=E_gated[MPa]  SDV5=phi_gate")
    ap("**")

    # Linear amplitude ramp 0→1 for α (avoids single-increment full eigenstrain)
    ap("*Amplitude, name=RAMP_ALPHA, definition=TABULAR")
    ap(" 0.0, 0.0")
    ap(" 1.0, 1.0")
    ap("**")

    # alpha_0 = 0 (initial, before ramp)
    ap("*Initial Conditions, type=field, variable=1")
    ap(" ALL_NODES, 0.0")
    ap("**")

    # phi_i: fixed composition fractions (uniform per condition)
    for var_idx, (sp, phi_val) in enumerate(zip(SPECIES, phi_vec), start=2):
        ap(f"*Initial Conditions, type=field, variable={var_idx}")
        ap(f" ALL_NODES, {phi_val:.8e}")
        ap("**")

    # phi_local (PREDEF 7): depth-varying total biofilm density — Klempt Eq.20 gate
    ap("*Initial Conditions, type=field, variable=7")
    for nid0, pl in enumerate(phi_nodes):
        ap(f" {nid0+1}, {pl:.8e}")
    ap("**")

    _write_step_ramped(ap, cond, alpha_nodes, var_start=1, var_end=1)
    out_path.write_text("\n".join(L) + "\n")
    print(f"  [A] {out_path.name}  ({len(L)} lines)  E_Voigt={E_voigt_mean:.3e} MPa")


def write_mode_C(out_path: Path, blocks: dict,
                 alpha_nodes_per_species: dict[str, np.ndarray],
                 phi_nodes: np.ndarray, phi_vec: np.ndarray, cond: str) -> None:
    """Mode C: umat_klempt2025.f — PREDEF 1..5=α_s, 6..10=φ_i, 11=ϕ_local."""
    L = []
    ap = L.append
    E_voigt_mean = float(np.dot(phi_vec, E_SPEC_MPa))
    alpha_total_max = sum(a.max() for a in alpha_nodes_per_species.values())

    ap("*HEADING")
    ap(f" Klempt 2025 multi-species UMAT (Option C) — condition: {cond}")
    ap(f" E_Voigt={E_voigt_mean:.3e} MPa  alpha_total_max={alpha_total_max:.4f}")
    ap(f" UMAT: umat_klempt2025.f  |  alpha_s ramped via *Field+AMPLITUDE")
    ap("**")

    for line in blocks["nodes"]:    ap(line)
    ap("**")
    for line in blocks["elements"]: ap(line)
    ap("**")
    for line in blocks["nsets"]:    ap(line)
    ap("**")
    _write_all_elems(ap, blocks)
    ap("**")

    # Material
    ap("*Solid Section, elset=ALL_ELEMS, material=KLEMPT_2025")
    ap(",")
    ap("*Material, name=KLEMPT_2025")
    ap("*User Material, constants=1")
    ap(f" {NU_KLEMPT:.4f}")
    ap("*Depvar")
    ap(" 7")
    ap("**  SDV1=alpha_total  SDV2=Je  SDV3=E_gated  SDV4=s  SDV5=alpha_So  SDV6=alpha_Vd  SDV7=phi_gate")
    ap("**")

    # Linear amplitude ramp for all per-species α fields
    ap("*Amplitude, name=RAMP_ALPHA, definition=TABULAR")
    ap(" 0.0, 0.0")
    ap(" 1.0, 1.0")
    ap("**")

    # alpha_s: initial = 0 for all 5 species fields
    for var_idx in range(1, 6):
        ap(f"*Initial Conditions, type=field, variable={var_idx}")
        ap(f" ALL_NODES, 0.0")
        ap("**")

    # phi_i: fixed composition fractions (uniform per condition)
    for var_idx, (sp, phi_val) in enumerate(zip(SPECIES, phi_vec), start=6):
        ap(f"*Initial Conditions, type=field, variable={var_idx}")
        ap(f" ALL_NODES, {phi_val:.8e}")
        ap("**")

    # phi_local (PREDEF 11): depth-varying total biofilm density — Klempt Eq.20 gate
    ap("*Initial Conditions, type=field, variable=11")
    for nid0, pl in enumerate(phi_nodes):
        ap(f" {nid0+1}, {pl:.8e}")
    ap("**")

    _write_step_ramped_multispecies(ap, cond, alpha_nodes_per_species)
    out_path.write_text("\n".join(L) + "\n")
    print(f"  [C] {out_path.name}  ({len(L)} lines)  α_total_max={alpha_total_max:.4f}")


def _write_all_elems(ap, blocks):
    in_all = False
    for line in blocks["elsets"]:
        ls = line.strip().upper()
        if "ALL_ELEMS" in ls: in_all = True
        elif ls.startswith("*ELSET") and "ALL_ELEMS" not in ls: in_all = False
        if in_all: ap(line)


def _write_step_ramped(ap, cond, alpha_nodes, var_start=1, var_end=1):
    """Write step that ramps α field via *Field + RAMP_ALPHA amplitude.

    Ramped application avoids convergence failure from full-eigenstrain
    single-increment load (α_max can exceed 1.0 for commensal conditions).
    *Static, stabilize adds viscous dissipation to handle contact singularities.
    """
    ap("*Boundary")
    ap(" INNER_FACE, ENCASTRE")
    ap("**")
    ap(f"*Step, name=GROWTH, nlgeom=YES")
    ap(f" Klempt multi-species growth stress — {cond}")
    ap("*Static, stabilize=1e-4, allsdtol=0.05, continue=NO")
    ap(" 5e-3, 1.0, 1e-7, 0.05")
    ap("**")
    # Prescribe α field (ramped from initial 0 to final α_nodes values)
    for var_idx in range(var_start, var_end + 1):
        ap(f"*Field, variable={var_idx}, amplitude=RAMP_ALPHA")
        for nid0, a in enumerate(alpha_nodes):
            ap(f" {nid0+1}, {a:.8e}")
        ap("**")
    ap("*Output, field")
    ap("*Node Output")
    ap(" U, RF")
    ap("*Element Output")
    ap(" S, E, MISES, SDV")
    ap("**")
    ap("*End Step")


def _write_step_ramped_multispecies(ap, cond, alpha_nodes_per_species):
    """Write step that ramps 5 per-species α fields for mode C."""
    ap("*Boundary")
    ap(" INNER_FACE, ENCASTRE")
    ap("**")
    ap(f"*Step, name=GROWTH, nlgeom=YES")
    ap(f" Klempt 2025 multi-species growth stress — {cond}")
    ap("*Static, stabilize=1e-4, allsdtol=0.05, continue=NO")
    ap(" 5e-3, 1.0, 1e-7, 0.05")
    ap("**")
    for var_idx, sp in enumerate(SPECIES, start=1):
        alpha_sp = alpha_nodes_per_species[sp]
        ap(f"*Field, variable={var_idx}, amplitude=RAMP_ALPHA")
        for nid0, a in enumerate(alpha_sp):
            ap(f" {nid0+1}, {a:.8e}")
        ap("**")
    ap("*Output, field")
    ap("*Node Output")
    ap(" U, RF")
    ap("*Element Output")
    ap(" S, E, MISES, SDV")
    ap("**")
    ap("*End Step")


def main(mode: str):
    print("=" * 65)
    print(f"gen_tooth_klempt_umat_inp.py  mode={mode}")
    print("=" * 65)

    tmcmc = load_tmcmc()
    print(f"\n[1] TMCMC loaded. E_Voigt per condition [MPa]:")
    for cond, phi in tmcmc.items():
        e = np.dot(phi, E_SPEC_MPa)
        print(f"  {cond:25s}  E_Voigt={e:.3e} MPa  "
              f"({'stiff' if e > 5e-4 else 'soft'})")

    print(f"\n[2] Parsing conformal mesh template...")
    blocks = parse_template(TMPL)

    print(f"\n[3] Loading condition-specific α and ϕ fields (klempt_pde_multispecies)...")
    alpha_fields = {}
    phi_fields   = {}
    for cond in CONDITIONS:
        alpha_fields[cond] = load_alpha_condition(cond)
        phi_fields[cond]   = load_phi_condition(cond)
        print(f"  {cond:25s}  α_max={alpha_fields[cond].max():.4f}  ϕ_max={phi_fields[cond].max():.4f}")

    print(f"\n[4] Generating {mode} inp files...")

    for cond in CONDITIONS:
        phi_vec  = tmcmc[cond]
        alpha2d  = alpha_fields[cond]
        phi2d    = phi_fields[cond]

        # phi_local per node: depth-varying ϕ for Klempt Eq.20 gating
        phi_nodes = np.array([
            phi_local_at_depth(depth_norm(nid0), phi2d)
            for nid0 in range(N_NODES)
        ])

        if mode == "A":
            alpha_nodes = np.array([
                alpha_at_depth(depth_norm(nid0), alpha2d)
                for nid0 in range(N_NODES)
            ])
            out = HERE / f"p23_klempt_A_{cond}.inp"
            write_mode_A(out, blocks, alpha_nodes, phi_nodes, phi_vec, cond)

        elif mode == "C":
            phi_total = phi_vec.sum()
            alpha_nodes_per_sp = {}
            alpha_total_nodes = np.array([
                alpha_at_depth(depth_norm(nid0), alpha2d)
                for nid0 in range(N_NODES)
            ])
            for sp_idx, sp in enumerate(SPECIES):
                sp_frac = phi_vec[sp_idx] / (phi_total + 1e-12)
                alpha_nodes_per_sp[sp] = alpha_total_nodes * sp_frac
            out = HERE / f"p23_klempt_C_{cond}.inp"
            write_mode_C(out, blocks, alpha_nodes_per_sp, phi_nodes, phi_vec, cond)

    print(f"\n{'='*65}")
    print(f"Done. Run with:")
    print(f"  abaqus job=p23_klempt_{mode}_{{cond}} user={UMAT[mode].name} cpus=1 ask_delete=OFF")
    print(f"{'='*65}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["A", "C"], default="A",
                        help="A=Voigt+total_alpha, C=Klempt2025_per_species")
    args = parser.parse_args()
    main(args.mode)
