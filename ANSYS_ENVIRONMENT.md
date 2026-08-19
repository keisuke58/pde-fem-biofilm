# ANSYS Environment Report

- **Host:** IKMHIWI03 (Dell Pro Tower Plus QBT1250)
- **OS:** Windows 11 Pro, 64-bit (build 10.0.26200)
- **CPU:** Intel(R) Core(TM) Ultra 5 235 — 14 cores / 14 logical processors
- **RAM:** 31.5 GB
- **GPU:** NVIDIA RTX A400 (+ Intel integrated graphics)
- **Report generated:** 2026-08-19

## Installed Version

**ANSYS 2022 R2** (v222), installed 2025-07-09
Install directory: `C:\Program Files\ANSYS Inc\v222`
Total install footprint: ~50.3 GB

## Environment Variables

| Variable | Value |
|---|---|
| `AWP_ROOT222` | `C:\Program Files\ANSYS Inc\v222` |
| `ANSYS222_DIR` | `C:\Program Files\ANSYS Inc\v222\ANSYS` |
| `CADOE_LIBDIR222` | `C:\Program Files\ANSYS Inc\v222\CommonFiles\Language\de` |
| `AWP_LOCALE222` | `en-us` |

No `ANSYSLMD_LICENSE_FILE` / `ANSYSLI_SERVERS` variables are set in the user/system environment — licensing is resolved via the client config file instead (see below). No PATH entries reference ANSYS directly (products are launched via Start Menu / AnsysEDT shortcuts rather than a bare CLI on PATH).

## License Configuration

Config file: `C:\Program Files\ANSYS Inc\Shared Files\Licensing\ansyslmd.ini`

```
SERVER=1055@ansys-lic.rrzn.uni-hannover.de
ANSYSLI_SERVERS=2325@ansys-lic.rrzn.uni-hannover.de
```

License server is the **RRZN (Leibniz Universität Hannover)** floating license server — this machine is a client, not a license host. Network/VPN access to `ansys-lic.rrzn.uni-hannover.de` is required for licenses to check out.

## Installed Products / Toolset (2022 R2)

Identified from Start Menu shortcuts and program folders:

**Structural / Mechanical**
- Mechanical, Mechanical APDL (+ Product Launcher)
- ACP (Composite PrepPost)
- LS-DYNA / LS-PrePost 4.8.29 / LS-Run
- Statistics on Structures (+ Viewer)

**Fluids**
- Fluent
- CFX / CFD-Post
- Polyflow
- ICEM CFD
- TurboGrid
- Icepak

**Electronics**
- Ansys Electronics (folder present: `Electronics`)

**Geometry / CAD**
- Discovery
- SpaceClaim
- CAD Configuration Manager / File Association / Product & CAD Configuration

**Marine / Offshore**
- Aqwa, AqwaGS, AqwaWave

**Systems / Workflow**
- Workbench
- System Coupling
- Design Point Service (DPS)
- DC Evaluator (DCE)
- optiSLang
- Remote Solve Manager (RSM) — Cluster Monitoring, Job Monitoring, Configuration
- DPF (Data Processing Framework)

**Photonics**
- Lumerical

**Visualization / Post**
- EnSight, EnSight Launcher, EnVe, EnVideo, EnVision
- CFD-Post
- Nexus Launcher

**Admin/Utility**
- ANS_ADMIN
- Ansys Client Licensing Settings
- ARC / DCG Configuration
- Ansys Help / Help Configuration

This is effectively the **full ANSYS 2022 R2 simulation suite** (structural, fluids, electronics, photonics, systems/workflow tools), not a single-product install.

## Notes / Things to Watch

- Only **v222 (2022 R2)** is installed — no newer or older release coexists on this machine.
- Licensing depends on network access to the Uni Hannover RRZN license server; if working off-campus, VPN is likely required.
- GPU is an NVIDIA RTX A400 (entry-level workstation GPU, 4 GB) — fine for Discovery/SpaceClaim viewport and light GPU-accelerated solves, but a constraint for large GPU-solver workloads (e.g. Fluent GPU solver, LS-DYNA GPU) compared to higher-tier RTX Ada cards.
- 31.5 GB RAM / 14 logical cores is a reasonable mid-range workstation spec — worth keeping in mind when sizing mesh counts or parallel domain decomposition for Fluent/Mechanical solves.
- **Update 2026-08-19 (custom UMAT build toolchain, verified):** Intel Fortran
  and Visual Studio are both present, contrary to the original report above.
  - `ifort.exe`: `C:\Program Files (x86)\Intel\oneAPI\compiler\2025.3\bin\ifort.exe`
    (Intel Fortran Compiler 2025.3.3, "ifort (IFX)" build 20260319).
  - Visual Studio: version **18** (2026 Developer), at
    `C:\Program Files\Microsoft Visual Studio\18\Community\`. This is newer
    than what ANSYS 2022 R2's platform-support table expects (historically
    VS2019-era) — a version mismatch is possible at link time; not yet
    confirmed either way since the build itself could not be completed
    non-interactively (see `ansys_usermat/apdl/RUNBOOK.md` Step 1).
  - `C:\Program Files (x86)\Intel\oneAPI\setvars.bat` (run with no args, or
    `intel64 vs2022`) does **not** reliably put `ifort` on `PATH`: it shells
    out to a bare `vswhere.exe` which isn't on `PATH` either
    (`C:\Program Files (x86)\Microsoft Visual Studio\Installer\vswhere.exe`),
    causing the VS init to warn/fail, which in turn makes the per-component
    `env\vars.bat` calls for `compiler`/`mpi`/`umf` fail with "command not
    found". Working init sequence instead calls the two env scripts directly:
    ```bat
    call "C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat"
    call "C:\Program Files (x86)\Intel\oneAPI\compiler\2025.3\env\vars.bat"
    ```
  - `ansysli_util.exe` is not at the path RUNBOOK.md originally guessed
    (`Shared Files\Licensing\winx64\`); it's at
    `C:\Program Files\ANSYS Inc\v222\licensingclient\winx64\ansysli_util.exe`.
    License checkout from this machine succeeded and resolved locally
    (`server=55206@ikmhiwi03...`, academic `Ansys Mechanical Enterprise`),
    not via the RRZN server address in `ansyslmd.ini` — worth re-checking
    that config file if a checkout ever fails unexpectedly.
  - `ANSCUST.BAT` (v222's UPF link script) is interactive and its Y/N prompts
    are read via a bundled `ASK.EXE` that reads the real console directly,
    not stdin — it cannot be driven by piping answers into a redirected
    `cmd.exe`. It must be run from an actual interactive terminal by a human
    sitting at the machine.
