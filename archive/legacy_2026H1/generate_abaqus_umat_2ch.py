"""
generate_abaqus_umat_2ch.py
===========================
Abaqus .inp generator for the 2-channel Prony UMAT (umat_biofilm_visco_2ch.f).

Fixes vs v1:
  - C3D8  (full integration) instead of C3D8R — avoids hourglass stiffness issue
  - *STATIC instead of *VISCO — correct for UMAT-managed viscoelasticity
  - *NODE OUTPUT  instead of *NODE HISTORY — Abaqus 2024 syntax
  - NSET data lines capped at 16 node IDs
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np

# ─────────────────────────────────────────────────────────────────────────────
# Physical defaults
# ─────────────────────────────────────────────────────────────────────────────
L_M = 200.0e-6   # film thickness [m]  (0.2 mm — matches existing 1D bar INP)
L_MM = L_M * 1e3  # [mm]
W_MM = 1.0e-3     # slab width [mm]  (1 element in X and Y)

NZ_DEFAULT = 29
T_PERIOD   = 1000.0    # [s]  covers ~12*tau2; shortens run, full relaxation captured
DT_INIT    = 0.5       # << tau1=5s so backward-Euler Prony update stays stable
DT_MIN     = 1.0e-4
DT_MAX     = 5.0       # = tau1; keeps dt/tau1 <= 1
FREQOUT    = 20


# ─────────────────────────────────────────────────────────────────────────────
# Eigenstrain profile
# ─────────────────────────────────────────────────────────────────────────────
def read_profile_from_inp(inp_path: Path):
    txt = inp_path.read_text()
    node_block = re.search(r'\*NODE\n(.*?)\n\*', txt, re.DOTALL)
    nodes = {}
    if node_block:
        for line in node_block.group(1).splitlines():
            parts = line.split(',')
            if len(parts) >= 4:
                try:
                    nodes[int(parts[0])] = float(parts[3])
                except ValueError:
                    pass
    temp_block = re.search(r'\*TEMPERATURE\n(.*?)\n\*', txt, re.DOTALL)
    temps = {}
    if temp_block:
        for line in temp_block.group(1).splitlines():
            parts = line.split(',')
            if len(parts) >= 2:
                try:
                    temps[int(parts[0])] = float(parts[1])
                except ValueError:
                    pass
    if nodes and temps:
        nids = sorted(nodes.keys())
        return (np.array([nodes[n] for n in nids]),
                np.array([temps[n] for n in nids]))
    return None, None


# ─────────────────────────────────────────────────────────────────────────────
# UMAT Prony parameters
# ─────────────────────────────────────────────────────────────────────────────
def compute_props(G0_Pa, nu, e_inf, e1, tau1_s, tau2_s):
    e2    = 1.0 - e_inf - e1
    C10   = G0_Pa * 0.5e-6        # G = 2*C10  [MPa]
    G_MPa = G0_Pa * 1e-6
    K_MPa = 2 * G_MPa * (1 + nu) / (3 * (1 - 2 * nu))
    D1    = 2.0 / K_MPa
    eta1  = 2.0 * C10 * e1 * tau1_s
    eta2  = 2.0 * C10 * e2 * tau2_s
    return dict(C10=C10, C01=0.0, D1=D1, eta1=eta1, eta2=eta2,
                a1=e1, a2=e2, mtype=0.0,
                G0_MPa=G_MPa, G_inf_MPa=G_MPa * e_inf,
                K_MPa=K_MPa, nu=nu, tau1_s=tau1_s, tau2_s=tau2_s)


# ─────────────────────────────────────────────────────────────────────────────
# Mesh
# ─────────────────────────────────────────────────────────────────────────────
def make_mesh(nz, L_mm, W_mm):
    x_vals = [0.0, W_mm]
    y_vals = [0.0, W_mm]
    z_vals = [L_mm * iz / nz for iz in range(nz + 1)]
    nodes = {}
    node_idx = {}
    nid = 1
    for iz in range(nz + 1):
        for iy in range(2):
            for ix in range(2):
                nodes[nid] = (x_vals[ix], y_vals[iy], z_vals[iz])
                node_idx[(ix, iy, iz)] = nid
                nid += 1
    elements = []
    for iz in range(nz):
        n = [node_idx[(ix, iy, jz)]
             for jz in (iz, iz + 1)
             for iy in range(2) for ix in range(2)]
        # C3D8 connectivity: bottom face n1-4, top face n5-8
        n1 = node_idx[(0, 0, iz)];  n2 = node_idx[(1, 0, iz)]
        n3 = node_idx[(1, 1, iz)];  n4 = node_idx[(0, 1, iz)]
        n5 = node_idx[(0, 0, iz+1)]; n6 = node_idx[(1, 0, iz+1)]
        n7 = node_idx[(1, 1, iz+1)]; n8 = node_idx[(0, 1, iz+1)]
        elements.append((iz + 1, n1, n2, n3, n4, n5, n6, n7, n8))
    return nodes, node_idx, elements


def write_nset(f, name, nids, generate=False):
    if generate:
        lo, hi = min(nids), max(nids)
        f.write(f"*NSET, NSET={name}, GENERATE\n")
        f.write(f"  {lo}, {hi}, 1\n")
    else:
        f.write(f"*NSET, NSET={name}\n")
        nids = sorted(nids)
        for i in range(0, len(nids), 16):
            chunk = nids[i:i + 16]
            f.write("  " + ", ".join(str(n) for n in chunk) + "\n")


# ─────────────────────────────────────────────────────────────────────────────
# INP writer
# ─────────────────────────────────────────────────────────────────────────────
def write_inp(filepath, case, props, nz, alpha_profile):
    nodes, node_idx, elements = make_mesh(nz, L_MM, W_MM)
    C = props

    bottom_nodes = [node_idx[(ix, iy, 0)]  for ix in range(2) for iy in range(2)]
    x0_nodes     = list({node_idx[(0, iy, iz)] for iy in range(2) for iz in range(nz+1)})
    y0_nodes     = list({node_idx[(ix, 0, iz)] for ix in range(2) for iz in range(nz+1)})

    case_str = case.upper()
    mat_name = f"BIOFILM_{case_str}_2CH"

    with open(filepath, 'w') as f:
        f.write("**\n")
        f.write(f"** ABAQUS INPUT — BIOFILM VISCO 2-CH UMAT  ({case_str})\n")
        f.write(f"** UMAT: umat_biofilm_visco_2ch.f  NSTATV=18  NPROPS=8\n")
        f.write(f"** Element: C3D8  Step: *STATIC  Growth: TEMP=alpha_g(z) × RAMP(t)\n")
        f.write(f"** G0={C['G0_MPa']*1e6:.0f}Pa  nu={C['nu']:.2f}  tau1={C['tau1_s']:.2f}s  tau2={C['tau2_s']:.1f}s\n")
        f.write("**\n")
        f.write("*HEADING\n")
        f.write(f"Biofilm viscoelastic 2-ch UMAT  ({case_str})\n")
        f.write("SI-like units: MPa, mm, N, s\n")
        f.write("**\n")

        # Nodes
        f.write("*NODE\n")
        for nid, (x, y, z) in nodes.items():
            f.write(f" {nid:6d}, {x:.8f}, {y:.8f}, {z:.10f}\n")
        f.write("**\n")

        # Elements
        f.write("*ELEMENT, TYPE=C3D8, ELSET=BIOFILM_SLAB\n")
        for row in elements:
            eid, *conn = row
            f.write(f" {eid:5d},  " + ",  ".join(str(n) for n in conn) + "\n")
        f.write("**\n")

        # Node sets
        write_nset(f, "BOTTOM",    bottom_nodes)
        write_nset(f, "XFIX",      x0_nodes)
        write_nset(f, "YFIX",      y0_nodes)
        write_nset(f, "ALL_NODES", list(nodes.keys()), generate=True)
        f.write("**\n")

        # Material
        f.write(f"*MATERIAL, NAME={mat_name}\n")
        f.write("*USER MATERIAL, CONSTANTS=8\n")
        f.write("** C10[MPa], C01[MPa], D1[1/MPa], eta1[MPa*s], eta2[MPa*s], a1, a2, mtype\n")
        # All 8 constants on ONE line — Abaqus 2024 rejects multi-line for CONSTANTS<=8
        f.write(f" {C['C10']:.6e}, {C['C01']:.6e}, {C['D1']:.6e}, "
                f"{C['eta1']:.6e}, {C['eta2']:.6e}, "
                f"{C['a1']:.6e}, {C['a2']:.6e}, {C['mtype']:.1f}\n")
        f.write("*DEPVAR\n")
        f.write("  18\n")
        f.write("** SDV1-9=Fv_1  SDV10-18=Fv_2  (3x3 row-major each)\n")
        f.write("**\n")

        # Section
        f.write(f"*SOLID SECTION, ELSET=BIOFILM_SLAB, MATERIAL={mat_name}\n")
        f.write("**\n")

        # BCs — full lateral constraint: all nodes Ux=Uy=0 (1D uniaxial bar in Z)
        f.write("*BOUNDARY\n")
        f.write("BOTTOM,    3, 3, 0.0\n")
        f.write("ALL_NODES, 1, 2, 0.0\n")
        f.write("**\n")

        # Initial temperature = 0
        f.write("*INITIAL CONDITIONS, TYPE=TEMPERATURE\n")
        for iz in range(nz + 1):
            for iy in range(2):
                for ix in range(2):
                    f.write(f"  {node_idx[(ix,iy,iz)]:6d},  0.0\n")
        f.write("**\n")

        # Amplitude
        f.write("*AMPLITUDE, NAME=GROWTH_RAMP, TIME=TOTAL TIME, VALUE=RELATIVE\n")
        f.write(f"  0.0, 0.0,  {T_PERIOD:.1f}, 1.0\n")
        f.write("**\n")

        # Step: STATIC (correct for UMAT with internal viscous state variables)
        f.write("*STEP, NLGEOM=YES, NAME=GROWTH_STATIC, INC=50000\n")
        f.write("*STATIC\n")
        f.write(f"  {DT_INIT:.4f},  {T_PERIOD:.1f},  {DT_MIN:.2e},  {DT_MAX:.4f}\n")
        f.write("*CONTROLS, PARAMETERS=TIME INCREMENTATION\n")
        # I4=max eq-iters per attempt(16), I5=max attempts(20), I6=max cutbacks(10)
        f.write(" , , , 16, 20, 10\n")
        f.write("**\n")

        # Temperature field
        f.write("** TEMP(node,t) = alpha_g(z) * amplitude(t)\n")
        f.write("*TEMPERATURE, AMPLITUDE=GROWTH_RAMP\n")
        for iz in range(nz + 1):
            av = float(alpha_profile[iz])
            for iy in range(2):
                for ix in range(2):
                    f.write(f"  {node_idx[(ix,iy,iz)]:6d},  {av:.10f}\n")
        f.write("**\n")

        # Output
        f.write(f"*OUTPUT, FIELD, FREQUENCY={FREQOUT}\n")
        f.write("*NODE OUTPUT\n")
        f.write("U, NT\n")
        f.write("*ELEMENT OUTPUT, ELSET=BIOFILM_SLAB\n")
        f.write("S, E, EE, SDV\n")
        f.write("**\n")
        f.write("*OUTPUT, HISTORY, FREQUENCY=1\n")
        f.write("*NODE OUTPUT, NSET=BOTTOM\n")
        f.write("U3\n")
        f.write("**\n")
        f.write("*END STEP\n")
        f.write(f"** Submit: abaqus job={filepath.stem} user=umat_biofilm_visco_2ch.f\n")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case",  choices=["ch", "dh", "both"], default="both")
    parser.add_argument("--nz",    type=int,   default=NZ_DEFAULT)
    parser.add_argument("--g0",    type=float, default=2000.0)
    parser.add_argument("--nu",    type=float, default=0.45)
    parser.add_argument("--einf",  type=float, default=0.1977)
    parser.add_argument("--e1",    type=float, default=0.4553)
    parser.add_argument("--tau1",  type=float, default=5.01)
    parser.add_argument("--tau2",  type=float, default=80.72)
    parser.add_argument("--inp-dir",  default="~/IKM_Hiwi/FEM/_abaqus_input")
    parser.add_argument("--outdir",   default="~/IKM_Hiwi/FEM/_abaqus_input")
    parser.add_argument("--k-ratio",  type=float, default=6.44**(1.0/2.68),
                        help="k_eff^b(CH)/k_eff^b(DH); DH profile = CH_profile/k_ratio")
    args = parser.parse_args()

    K_RATIO = args.k_ratio

    props = compute_props(args.g0, args.nu, args.einf, args.e1, args.tau1, args.tau2)
    nz = args.nz
    z_nodes_mm = np.array([L_MM * iz / nz for iz in range(nz + 1)])

    inp_dir = Path(args.inp_dir).expanduser()
    out_dir = Path(args.outdir).expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)

    cases = ["ch", "dh"] if args.case == "both" else [args.case]

    print(f"G0={args.g0:.0f}Pa  C10={props['C10']:.4e}MPa  D1={props['D1']:.4e}/MPa")
    print(f"a1={props['a1']:.4f} tau1={args.tau1:.2f}s  a2={props['a2']:.4f} tau2={args.tau2:.2f}s")
    print(f"k_ratio={K_RATIO:.3f}  → alpha_DH = alpha_CH / {K_RATIO:.3f}")

    # Always derive CH profile from existing INP or fallback
    existing_ch = inp_dir / "biofilm_1d_bar_commensal_hobic.inp"
    if existing_ch.exists():
        z_ref, a_ref = read_profile_from_inp(existing_ch)
        alpha_ch = np.interp(z_nodes_mm, z_ref, a_ref) if z_ref is not None else None
    else:
        alpha_ch = None

    if alpha_ch is None:
        z_norm = np.clip(z_nodes_mm / L_MM, 0, 1)
        alpha_ch = 0.1399 * z_norm ** 0.8

    alpha_dh = alpha_ch / K_RATIO   # DH grows slower by k_ratio

    profiles = {"ch": alpha_ch, "dh": alpha_dh}

    for case in cases:
        alpha_profile = profiles[case]

        out = out_dir / f"biofilm_visco_2ch_{case}.inp"
        write_inp(out, case, props, nz, alpha_profile)
        print(f"  [{case.upper()}] alpha=[{alpha_profile.min():.4f}, {alpha_profile.max():.4f}] → {out.name}")


if __name__ == "__main__":
    main()
