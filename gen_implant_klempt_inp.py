"""gen_implant_klempt_inp.py
===========================
Generate Abaqus .inp for the 3-D helical implant screw with Klempt UMAT v3.

Geometry: same as gen_implant_umat_3d_inp.py (C3D8 hex, helical thread).
UMAT:     umat_klempt_voigt.f  (mode A, phi^2-gated E, PREDEF 1-7)
          umat_klempt2025.f    (mode C, 5-species, PREDEF 1-11) [--mode C]
Data:     phi_i from ultimate_10000p MAP (ref_0d_{cond}.json)
          alpha_max from PDE (klempt_alpha_final_{cond}.npy, max value)
Depth:    linear radial — inner face j=0 (bonded to Ti) → full growth/phi
                          outer face j=Nr → zero growth/phi

Usage:
  python gen_implant_klempt_inp.py [--cond DH] [--mode A] [--out implant_klempt_A_DH.inp]
"""

import argparse
import json
from pathlib import Path

import numpy as np

# ── paths ──────────────────────────────────────────────────────────────────────
HERE  = Path(__file__).resolve().parent
JAXF  = HERE / "JAXFEM"
MSCL  = HERE / "_multiscale_2d_results"
NIFE_ABAQUS = Path("/home/nishioka/IKM_Hiwi/nife/masterarbeit_ansys_fem/coupling_prototype/abaqus")

# ── geometry parameters (match nife gen_implant_umat_3d_inp.py) ────────────────
R0, HF, P = 2.0, 1.0, 1.0      # core radius, biofilm thickness, thread pitch [mm]
AMP    = 0.18                    # thread amplitude [mm]
NTURNS = 3.0
Nt, Nr, Nz = 48, 4, 60          # circumferential, radial, axial divisions

# ── species names for mode C alpha split ──────────────────────────────────────
SPECIES  = ["So", "An", "Vd", "Fn", "Pg"]
K_ALPHA  = np.array([1.0, 0.8, 0.4, 0.6, 0.3])   # species growth accumulation rates


def r_surf(theta, z):
    return R0 + AMP * 0.5 * (1.0 - np.cos(2.0 * np.pi * z / P - theta))


def nid(i, j, k):
    """Node ID (1-based): i=circ (periodic Nt), j=radial (0=inner), k=axial."""
    return k * (Nt * (Nr + 1)) + j * Nt + (i % Nt) + 1


def depth_norm_of(j):
    """Normalised depth in [0,1]: 0 = inner face (Ti surface), 1 = outer (planktonic)."""
    return j / Nr


def write_inp(cond: str, mode: str, out_path: Path):
    # ── load TMCMC ultimate_10000p MAP ──────────────────────────────────────────
    ref = json.load(open(MSCL / f"ref_0d_{cond}.json"))
    phi_i    = np.array(ref["phi_final"])         # shape (5,), ultimate_10000p MAP
    phi_sum  = float(phi_i.sum())
    phi_So, phi_An, phi_Vd, phi_Fn, phi_Pg = phi_i

    # ── load PDE alpha_max ──────────────────────────────────────────────────────
    alpha_2d   = np.load(JAXF / f"klempt_alpha_final_{cond}.npy")
    alpha_max  = float(alpha_2d.max())

    # ── per-species alpha split for mode C (Σ alpha_s = alpha_total) ──────────
    k_alpha_eff = float(phi_i @ K_ALPHA)
    # weight_s proportional to phi_i * k_alpha_s, normalised to sum=1
    weights = (phi_i * K_ALPHA) / (k_alpha_eff + 1e-15)
    alpha_s_max = alpha_max * weights          # shape (5,)

    # ── UMAT file ───────────────────────────────────────────────────────────────
    if mode == "A":
        umat_file = NIFE_ABAQUS / "umat_klempt_voigt.f"
        n_depvar  = 5
        n_predef  = 7
    else:
        umat_file = NIFE_ABAQUS / "umat_klempt2025.f"
        n_depvar  = 7
        n_predef  = 11

    # ── build all nodes ─────────────────────────────────────────────────────────
    HZ     = NTURNS * P
    thetas = np.linspace(0, 2 * np.pi, Nt, endpoint=False)
    zs     = np.linspace(0, HZ, Nz + 1)

    # Collect node coords (for phi_local / alpha profile)
    # Also collect node_list for sorted output
    node_js = {}  # nid → j (radial index)
    L = []
    ap = L.append

    # ── header ──────────────────────────────────────────────────────────────────
    ap("*HEADING")
    ap(f" Klempt implant screw mode {mode} cond={cond} | phi^2-gated E | ultimate_10000p MAP")
    ap(f" UMAT: {umat_file.name}  PREDEF {n_predef}  DEPVAR {n_depvar}")
    ap(f" alpha_max={alpha_max:.4f}  phi_sum={phi_sum:.4f}  k_alpha_eff={k_alpha_eff:.4f}")

    # ── nodes ──────────────────────────────────────────────────────────────────
    ap("*NODE")
    for k in range(Nz + 1):
        z = zs[k]
        for j in range(Nr + 1):
            for i in range(Nt):
                th = thetas[i]
                r  = r_surf(th, z) + (j / Nr) * HF
                n  = nid(i, j, k)
                ap(f" {n}, {r*np.cos(th):.6f}, {r*np.sin(th):.6f}, {z:.6f}")
                node_js[n] = j

    # ── elements (C3D8) ────────────────────────────────────────────────────────
    ap("*ELEMENT, TYPE=C3D8, ELSET=FILM")
    eid = 0
    for k in range(Nz):
        for j in range(Nr):
            for i in range(Nt):
                eid += 1
                ns = [nid(i,   j,   k), nid(i,   j+1, k),
                      nid(i+1, j+1, k), nid(i+1, j,   k),
                      nid(i,   j,   k+1), nid(i,   j+1, k+1),
                      nid(i+1, j+1, k+1), nid(i+1, j,   k+1)]
                ap(f" {eid}, " + ", ".join(str(x) for x in ns))

    # ── node sets ──────────────────────────────────────────────────────────────
    all_ids   = sorted(node_js.keys())
    inner_ids = [n for n, j in node_js.items() if j == 0]
    zbot_ids  = [nid(i, j, 0)  for j in range(Nr+1) for i in range(Nt)]
    ztop_ids  = [nid(i, j, Nz) for j in range(Nr+1) for i in range(Nt)]

    for name, ids in (("ALL_NODES", all_ids), ("INNER", inner_ids),
                      ("ZBOT", zbot_ids), ("ZTOP", ztop_ids)):
        ap(f"*NSET, NSET={name}")
        ids_s = sorted(set(ids))
        for kk in range(0, len(ids_s), 16):
            ap(" " + ",".join(str(x) for x in ids_s[kk:kk+16]))

    # ── material ────────────────────────────────────────────────────────────────
    ap("*SOLID SECTION, ELSET=FILM, MATERIAL=KLEMPT_BIOFILM")
    ap("*MATERIAL, NAME=KLEMPT_BIOFILM")
    ap("*USER MATERIAL, CONSTANTS=1")
    ap(" 0.49")
    ap("*DEPVAR")
    ap(f" {n_depvar}")
    ap(f"**  Mode {mode} SDVs: " + (
        "SDV1=s  SDV2=Je  SDV3=alpha  SDV4=E_gated[MPa]  SDV5=phi_gate"
        if mode == "A" else
        "SDV1=alpha_total  SDV2=Je  SDV3=E_gated  SDV4=s  SDV5=alpha_So  SDV6=alpha_Vd  SDV7=phi_gate"
    ))

    # ── amplitude ──────────────────────────────────────────────────────────────
    ap("*Amplitude, name=RAMP_ALPHA, definition=TABULAR")
    ap(" 0.0, 0.0")
    ap(" 1.0, 1.0")

    # ── initial conditions ─────────────────────────────────────────────────────
    # Variable 1 (alpha_total or alpha_So) — initial = 0, ramped in step via *Field
    ap("**")
    if mode == "A":
        ap("*Initial Conditions, type=field, variable=1")
        ap(" ALL_NODES, 0.0")
        # phi_i (uniform per condition)
        for vidx, val in enumerate([phi_So, phi_An, phi_Vd, phi_Fn, phi_Pg], start=2):
            ap(f"*Initial Conditions, type=field, variable={vidx}")
            ap(f" ALL_NODES, {val:.10e}")
        # phi_local (depth-varying, variable 7)
        ap("*Initial Conditions, type=field, variable=7")
        for n in all_ids:
            d = depth_norm_of(node_js[n])
            ap(f" {n}, {(1.0 - d):.10e}")
    else:  # mode C
        for sp in range(5):
            ap(f"*Initial Conditions, type=field, variable={sp+1}")
            ap(f" ALL_NODES, 0.0")
        # phi_i (uniform, variables 6..10)
        for vidx, val in enumerate([phi_So, phi_An, phi_Vd, phi_Fn, phi_Pg], start=6):
            ap(f"*Initial Conditions, type=field, variable={vidx}")
            ap(f" ALL_NODES, {val:.10e}")
        # phi_local (variable 11)
        ap("*Initial Conditions, type=field, variable=11")
        for n in all_ids:
            d = depth_norm_of(node_js[n])
            ap(f" {n}, {(1.0 - d):.10e}")

    # ── BCs ────────────────────────────────────────────────────────────────────
    ap("*BOUNDARY")
    ap(" INNER, 1, 3")   # bonded to Ti surface: u_x=u_y=u_z=0
    ap(" ZBOT,  3, 3")
    ap(" ZTOP,  3, 3")

    # ── step: ramp alpha via *Field ────────────────────────────────────────────
    ap("*STEP, NLGEOM=YES, INC=200")
    ap(f" Klempt implant {mode} {cond}: phi^2-gated growth stress")
    ap("*Static, stabilize=1e-4, allsdtol=0.05, continue=NO")
    ap(" 5e-3, 1.0, 1e-7, 0.05")
    ap("**")
    if mode == "A":
        ap("*Field, variable=1, amplitude=RAMP_ALPHA")
        for n in all_ids:
            d = depth_norm_of(node_js[n])
            ap(f" {n}, {alpha_max * (1.0 - d):.10e}")
    else:
        for sp in range(5):
            ap(f"*Field, variable={sp+1}, amplitude=RAMP_ALPHA")
            for n in all_ids:
                d = depth_norm_of(node_js[n])
                ap(f" {n}, {alpha_s_max[sp] * (1.0 - d):.10e}")

    # ── output ─────────────────────────────────────────────────────────────────
    ap("*OUTPUT, FIELD")
    ap("*NODE OUTPUT")
    ap(" U, COORD")
    ap("*ELEMENT OUTPUT, POSITION=CENTROID")
    ap(" S, SDV, COORD")
    ap("*END STEP")

    out_path.write_text("\n".join(L) + "\n")
    n_nodes = len(all_ids)
    n_elems = Nt * Nr * Nz
    print(f"Wrote: {out_path}")
    print(f"  Mesh: {Nt}x{Nr}x{Nz} C3D8  ({n_nodes} nodes, {n_elems} elems, {NTURNS:.0f} turns)")
    print(f"  Cond: {cond}  phi=[{phi_So:.3f},{phi_An:.3f},{phi_Vd:.3f},{phi_Fn:.3f},{phi_Pg:.3f}]")
    print(f"  alpha_max={alpha_max:.4f}  k_alpha_eff={k_alpha_eff:.4f}")
    print(f"  UMAT: {umat_file.name}  PREDEF={n_predef}  DEPVAR={n_depvar}")


def main():
    pa = argparse.ArgumentParser()
    pa.add_argument("--cond", default="dysbiotic_hobic",
                    choices=["commensal_static","commensal_hobic",
                             "dysbiotic_static","dysbiotic_hobic"])
    pa.add_argument("--mode", default="A", choices=["A","C"])
    pa.add_argument("--out", default=None)
    args = pa.parse_args()

    out = Path(args.out) if args.out else HERE / f"p23imp_klempt_{args.mode}_{args.cond}.inp"
    write_inp(args.cond, args.mode, out)


if __name__ == "__main__":
    main()
