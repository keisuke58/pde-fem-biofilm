"""felix_abaqus_pressure.py
Staggered JAXFEM → Abaqus pipeline for Felix Klempt 2024 Fig.8 reproduction.

Stage 1: Run JAXFEM PDE (felix_complete_reproduction.py) → save α(x,y) at 4 times
Stage 2: Generate Abaqus .inp per snapshot (59×59 C3D8, FIELD=α)
Stage 3: Run Abaqus with umat_klempt_alpha.f (Felix F=FeFg + neo-Hookean)
Stage 4: Extract hydrostatic pressure PRESS(x,y) per snapshot
Stage 5: Plot 4-panel figure → felix_abaqus_pressure.png

Goal: reproduce Felix Fig.8 tension ring (missing from 2D FD JAXFEM).

Run:
    PATH=.../texlive/2025/.../bin:$PATH python felix_abaqus_pressure.py
    # Abaqus must be on PATH: /home/nishioka/DassaultSystemes/SIMULIA/Commands/
"""
from __future__ import annotations
import os
import sys
import subprocess
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

# ── paths ────────────────────────────────────────────────────────────────────
HERE       = Path(__file__).resolve().parent
UMAT_F     = Path("/home/nishioka/IKM_Hiwi/nife/masterarbeit_ansys_fem/"
                  "coupling_prototype/abaqus/umat_klempt_alpha.f")
EXTRACT_PY = HERE / "extract_felix_press.py"
WORKDIR    = HERE / "_felix_abaqus"
ABAQUS_CMD = "/home/nishioka/DassaultSystemes/SIMULIA/Commands/abaqus"

sys.path.insert(0, str(HERE))

# ── import JAXFEM ─────────────────────────────────────────────────────────────
import felix_complete_reproduction as F

# Snapshot steps to send to Abaqus (subset of SAVE_AT)
SNAP_STEPS = [2000, 5000, 8000, 12000]

# ── mesh constants (must match JAXFEM grid) ──────────────────────────────────
NX, NY = F.NX, F.NY          # 60, 60  (nodes per axis)
LX, LY = F.LX, F.LY          # 20 μm,  20 μm
DX      = LX / (NX - 1)      # node spacing in x (μm)
DY      = LY / (NY - 1)      # node spacing in y
DZ      = DX                  # quasi-2D: 1 element thick in z

NNX = NX        # nodes in x: 60
NNY = NY        # nodes in y: 60
NNZ = 2         # 2 z-layers → 1 element in z
NEX = NX - 1    # elements in x: 59
NEY = NY - 1    # elements in y: 59

E_MOD = F.E_MOD   # 10.0 Pa  (Felix Table 2)
NU    = F.NU       # 0.49

# ── helper: linear node id (1-based Abaqus) ──────────────────────────────────
def nid(ix, iy, iz):
    """1-based Abaqus node id for grid position (ix, iy, iz)."""
    return iz * (NNX * NNY) + iy * NNX + ix + 1


# ── generate Abaqus .inp ──────────────────────────────────────────────────────
def gen_inp(jobname: str, alpha: np.ndarray, workdir: Path) -> Path:
    """
    Write a C3D8 Abaqus input for quasi-2D growth problem.
    alpha: shape (NY, NX) = (60, 60) — JAXFEM accumulated growth field.
    Field variable 1 = alpha is prescribed at each node via
    *INITIAL CONDITIONS, TYPE=FIELD, VARIABLE=1.
    """
    lines = []

    # --- heading ---
    lines.append(f"*HEADING")
    lines.append(f" Felix Klempt 2024 staggered: {jobname}")
    lines.append(f" JAXFEM alpha -> Abaqus UMAT F=FeFg, neo-Hookean")
    lines.append(f" Mesh: {NEX}x{NEY}x1 C3D8 on {NNX}x{NNY}x{NNZ} nodes")
    lines.append(f" Domain: {LX}x{LY}x{DZ:.4f} um")

    # --- nodes ---
    lines.append("**")
    lines.append("*NODE")
    for iz in range(NNZ):
        for iy in range(NNY):
            for ix in range(NNX):
                x = ix * DX
                y = iy * DY
                z = iz * DZ
                nid_ = nid(ix, iy, iz)
                lines.append(f" {nid_}, {x:.6f}, {y:.6f}, {z:.6f}")

    # --- elements (C3D8) ---
    lines.append("**")
    lines.append("*ELEMENT, TYPE=C3D8, ELSET=ALL")
    eid = 1
    for iy in range(NEY):
        for ix in range(NEX):
            # Bottom face (z=0): ix,iy ; ix+1,iy ; ix+1,iy+1 ; ix,iy+1
            n1 = nid(ix,   iy,   0)
            n2 = nid(ix+1, iy,   0)
            n3 = nid(ix+1, iy+1, 0)
            n4 = nid(ix,   iy+1, 0)
            # Top face (z=1)
            n5 = nid(ix,   iy,   1)
            n6 = nid(ix+1, iy,   1)
            n7 = nid(ix+1, iy+1, 1)
            n8 = nid(ix,   iy+1, 1)
            lines.append(f" {eid}, {n1},{n2},{n3},{n4},{n5},{n6},{n7},{n8}")
            eid += 1

    lines.append("*SOLID SECTION, ELSET=ALL, MATERIAL=BIOFILM")
    lines.append("**")

    # --- material ---
    lines.append("*MATERIAL, NAME=BIOFILM")
    lines.append("*USER MATERIAL, CONSTANTS=2")
    lines.append(f"** E,       nu")
    lines.append(f" {E_MOD},  {NU}")
    lines.append("*DEPVAR")
    lines.append(" 3")
    lines.append("**")

    # --- boundary conditions: prevent rigid body modes ---
    # Fix z=0 face: UZ=0 everywhere on front face
    lines.append("*NSET, NSET=ZFRONT, GENERATE")
    lines.append(f" 1, {NNX*NNY}, 1")
    # x=0 edge at z=0: fix UX (ix=0)
    lines.append("*NSET, NSET=XLEFT")
    xleft = [nid(0, iy, iz) for iz in range(NNZ) for iy in range(NNY)]
    for chunk_start in range(0, len(xleft), 16):
        lines.append(" " + ",".join(str(n) for n in xleft[chunk_start:chunk_start+16]))
    # y=0 edge: fix UY (iy=0)
    lines.append("*NSET, NSET=YBOTTOM")
    ybot = [nid(ix, 0, iz) for iz in range(NNZ) for ix in range(NNX)]
    for chunk_start in range(0, len(ybot), 16):
        lines.append(" " + ",".join(str(n) for n in ybot[chunk_start:chunk_start+16]))

    lines.append("*BOUNDARY")
    lines.append(" ZFRONT, 3, 3")      # UZ=0 at z=0
    lines.append(" XLEFT, 1, 1")       # UX=0 at x=0 (symmetry)
    lines.append(" YBOTTOM, 2, 2")     # UY=0 at y=0 (substrate anchor)
    lines.append("**")

    # --- initial alpha field as field variable 1 ---
    lines.append("*INITIAL CONDITIONS, TYPE=FIELD, VARIABLE=1")
    # alpha[iy, ix] corresponds to node at (ix*DX, iy*DY) in JAXFEM
    for iz in range(NNZ):
        for iy in range(NNY):
            for ix in range(NNX):
                nid_ = nid(ix, iy, iz)
                a = float(np.clip(alpha[iy, ix], 0.0, None))
                lines.append(f" {nid_}, {a:.8f}")
    lines.append("**")

    # --- step: single static increment (alpha prescribed, solve mechanics) ---
    lines.append("*STEP, NLGEOM=YES, INC=1")
    lines.append(f" {jobname}: quasi-static growth mechanics")
    lines.append("*STATIC")
    lines.append(" 1.0, 1.0, 1.0e-8, 1.0")
    lines.append("**")
    lines.append("*OUTPUT, FIELD, FREQ=1")
    lines.append("*NODE OUTPUT")
    lines.append(" U, COORD")
    lines.append("*ELEMENT OUTPUT, POSITION=CENTROID")
    lines.append(" S, SDV, COORD")
    lines.append("*END STEP")

    inp_path = workdir / f"{jobname}.inp"
    inp_path.write_text("\n".join(lines))
    return inp_path


# ── Abaqus extraction script (written next to the job) ───────────────────────
EXTRACT_SCRIPT = '''"""Abaqus Python ODB extraction for felix_klempt snapshots.
Run: abaqus python extract_felix_press.py <jobname>
Writes: <jobname>_press.npy  -- shape (N_elem, 4) = [x, y, PRESS, alpha]
"""
from __future__ import print_function
import sys, os
import numpy as np
from odbAccess import openOdb

job = sys.argv[1] if len(sys.argv) > 1 else 'felix_klempt_t0'
odb = openOdb(job + '.odb')
step = list(odb.steps.values())[-1]
fr   = step.frames[-1]

S     = fr.fieldOutputs['S']
COORD = fr.fieldOutputs['COORD']

try:
    SDV3 = fr.fieldOutputs['SDV3']
    sdv3 = {v.elementLabel: v.data for v in SDV3.values}
except KeyError:
    sdv3 = {}

coord_map = {v.elementLabel: v.data for v in COORD.values}

rows = []
for v in S.values:
    el   = v.elementLabel
    sd   = v.data          # (S11, S22, S33, S12, S13, S23) or (S11,S22,S33,S12)
    if len(sd) >= 3:
        press = -(sd[0] + sd[1] + sd[2]) / 3.0
    else:
        press = 0.0
    xy = coord_map.get(el, [0.0, 0.0, 0.0])
    alpha = sdv3.get(el, 0.0)
    if hasattr(alpha, '__len__'):
        alpha = alpha[0] if len(alpha) > 0 else 0.0
    rows.append([xy[0], xy[1], press, alpha])

odb.close()
arr = np.array(rows, dtype=float)
out = job + '_press.npy'
np.save(out, arr)
print('Saved', out, arr.shape)
'''


# ── run one Abaqus job ────────────────────────────────────────────────────────
def run_abaqus_job(jobname: str, workdir: Path) -> bool:
    """Run Abaqus interactively. Returns True if .odb was created."""
    odb_path = workdir / f"{jobname}.odb"
    if odb_path.exists():
        print(f"  [skip] {jobname}.odb already exists")
        return True
    cmd = [ABAQUS_CMD, f"job={jobname}", f"user={UMAT_F}", "interactive"]
    print(f"  Running: {' '.join(str(c) for c in cmd)}")
    result = subprocess.run(cmd, cwd=str(workdir),
                            capture_output=True, text=True, timeout=3600)
    if result.returncode != 0 or not odb_path.exists():
        print(f"  ERROR: Abaqus failed for {jobname}")
        print(result.stdout[-2000:])
        print(result.stderr[-1000:])
        return False
    print(f"  OK: {jobname}.odb created")
    return True


# ── extract PRESS from .odb ───────────────────────────────────────────────────
def extract_press(jobname: str, workdir: Path) -> np.ndarray | None:
    """
    Call 'abaqus python extract_felix_press.py <jobname>' in workdir.
    Returns numpy array of shape (N_elem, 4) = [x, y, PRESS, alpha].
    """
    npy_path = workdir / f"{jobname}_press.npy"
    if npy_path.exists():
        print(f"  [skip] {jobname}_press.npy already exists")
        return np.load(str(npy_path))
    extract_script = workdir / "extract_felix_press.py"
    extract_script.write_text(EXTRACT_SCRIPT)
    cmd = [ABAQUS_CMD, "python", str(extract_script), jobname]
    result = subprocess.run(cmd, cwd=str(workdir),
                            capture_output=True, text=True, timeout=120)
    if not npy_path.exists():
        print(f"  ERROR extracting {jobname}")
        print(result.stdout)
        print(result.stderr)
        return None
    return np.load(str(npy_path))


# ── scatter → regular grid (inverse-distance or nearest-neighbour) ────────────
def to_grid(data: np.ndarray, N: int = 59) -> np.ndarray:
    """Map centroid (x,y,val) array to NxN regular grid (0..LX, 0..LY)."""
    from scipy.interpolate import griddata
    xs = data[:, 0]
    ys = data[:, 1]
    vs = data[:, 2]
    xi = np.linspace(xs.min(), xs.max(), N)
    yi = np.linspace(ys.min(), ys.max(), N)
    XX, YY = np.meshgrid(xi, yi)
    grid = griddata((xs, ys), vs, (XX, YY), method='nearest')
    return grid, xi, yi


# ── main ─────────────────────────────────────────────────────────────────────
def main():
    WORKDIR.mkdir(exist_ok=True)

    # ── Stage 1: JAXFEM PDE ────────────────────────────────────────────────
    print("=== Stage 1: JAXFEM B_hi PDE ===")
    phi0_b = F.sphere_phi(F.LX/2, 2.5, r=5.0)
    c_bc_b, c_bv_b = F.bc_bottom()
    snaps, avgs = F.run("B_hi", phi0_b, c_bc_b, c_bv_b,
                        g_ratio=F.G_HI/F.D_STAR)
    print("  Available snaps:", sorted(snaps.keys()))

    snap_list = [(s, snaps[s]) for s in SNAP_STEPS if s in snaps]
    if not snap_list:
        print("ERROR: No requested snapshots available")
        sys.exit(1)

    # ── Stage 2: Generate .inp ─────────────────────────────────────────────
    print("\n=== Stage 2: Generate Abaqus .inp ===")
    jobnames = []
    alphas_jax = []
    phis_jax   = []
    times_t    = []
    for i, (step, (phi, c, alpha)) in enumerate(snap_list):
        t_star = step * F.DT
        jobname = f"felix_klempt_t{i}"
        jobnames.append(jobname)
        alphas_jax.append(alpha)
        phis_jax.append(phi)
        times_t.append(t_star)
        inp = gen_inp(jobname, alpha, WORKDIR)
        print(f"  Written: {inp.name}  (step={step}, T*={t_star:.2f},"
              f" alpha_max={alpha.max():.5f})")

    # ── Stage 3: Run Abaqus ────────────────────────────────────────────────
    print("\n=== Stage 3: Abaqus runs ===")
    success = []
    for jn in jobnames:
        ok = run_abaqus_job(jn, WORKDIR)
        success.append(ok)

    if not any(success):
        print("All Abaqus jobs failed. Exiting.")
        sys.exit(1)

    # ── Stage 4: Extract PRESS ─────────────────────────────────────────────
    print("\n=== Stage 4: Extract PRESS ===")
    press_grids  = []
    alpha_grids  = []
    for jn, ok in zip(jobnames, success):
        if not ok:
            press_grids.append(None)
            alpha_grids.append(None)
            continue
        data = extract_press(jn, WORKDIR)
        if data is None:
            press_grids.append(None)
            alpha_grids.append(None)
        else:
            pg, xi, yi = to_grid(data, N=59)
            ag, _,  _  = to_grid(data[:, [0, 1, 3, 2]], N=59)  # col 3=alpha
            press_grids.append(pg)
            alpha_grids.append(ag)

    # ── Stage 5: Plot ──────────────────────────────────────────────────────
    print("\n=== Stage 5: Plot ===")

    sys.path.insert(0, str(HERE))
    try:
        from thesis_style import use as thesis_use
        thesis_use()
    except Exception:
        pass

    n_panels = sum(1 for pg in press_grids if pg is not None)
    if n_panels == 0:
        print("No data to plot.")
        sys.exit(0)

    fig, axes = plt.subplots(2, n_panels, figsize=(3.5*n_panels, 7))
    if n_panels == 1:
        axes = [[axes[0]], [axes[1]]]

    col = 0
    for i, (jn, pg, ag, phi_jax, alpha_jax, t_star, ok) in enumerate(
            zip(jobnames, press_grids, alpha_grids,
                phis_jax, alphas_jax, times_t, success)):
        if not ok or pg is None:
            continue

        xi = np.linspace(0, LX, 59)
        yi = np.linspace(0, LY, 59)

        # Row 0: hydrostatic pressure (Abaqus PRESS = -(S11+S22+S33)/3)
        ax0 = axes[0][col]
        vlim = max(abs(pg[~np.isnan(pg)]).max(), 0.1)
        im0 = ax0.pcolormesh(xi, yi, pg,
                              cmap='RdBu_r', vmin=-vlim, vmax=vlim,
                              shading='auto')
        # Overlay phi=0.5 contour from JAXFEM
        xx, yy = np.meshgrid(
            np.linspace(0, LX, F.NX),
            np.linspace(0, LY, F.NY))
        ax0.contour(xx, yy, phi_jax, levels=[0.5],
                    colors='k', linewidths=0.8)
        ax0.set_aspect('equal')
        ax0.set_title(f"$T^*={t_star:.1f}$", fontsize=9)
        if col == 0:
            ax0.set_ylabel("Pressure $p$ [Pa]", fontsize=8)
        plt.colorbar(im0, ax=ax0, shrink=0.8, label="PRESS [Pa]")

        # Row 1: alpha field from Abaqus SDV3 vs JAXFEM alpha
        ax1 = axes[1][col]
        im1 = ax1.pcolormesh(xi, yi, ag,
                              cmap='viridis', vmin=0,
                              vmax=max(alpha_jax.max(), 1e-8),
                              shading='auto')
        ax1.contour(xx, yy, phi_jax, levels=[0.5],
                    colors='w', linewidths=0.8)
        ax1.set_aspect('equal')
        ax1.set_title(f"$\\alpha$ (Abaqus SDV3)", fontsize=9)
        if col == 0:
            ax1.set_ylabel(r"$\alpha$ field", fontsize=8)
        plt.colorbar(im1, ax=ax1, shrink=0.8, label=r"$\alpha$")

        col += 1

    fig.suptitle("Felix Klempt 2024 Fig.8 analog — Abaqus staggered\n"
                 r"Pressure $p=-(S_{11}+S_{22}+S_{33})/3$ [Pa],"
                 r" $\phi{=}0.5$ contour", fontsize=10)
    fig.tight_layout()
    out = HERE / "felix_abaqus_pressure.png"
    fig.savefig(str(out), dpi=180, bbox_inches='tight')
    plt.close(fig)
    print(f"\nSaved: {out}")
    print("\nExpected physics:")
    print("  - Blue interior (phi>0.5): compression from growth (PRESS > 0)")
    print("  - Red ring at phi~0.5 boundary: tension ring (PRESS < 0)")
    print("  - This is the missing piece from JAXFEM v4 (KAPPA proxy only)")


if __name__ == "__main__":
    main()
