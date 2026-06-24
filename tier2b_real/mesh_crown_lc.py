"""h-convergence variant of mesh_crown.py: accepts the crown tet edge LC from sys.argv[1] so the
headline mesh_crown.py (LC = 0.40) is preserved untouched.

Imports mesh_crown.main() after monkey-patching the module-level LC constant — this matches the
headline geometry / OCC steps / cache filename byte-for-byte; only the meshSize knob changes.

Run in gmsh_env with LD_LIBRARY_PATH=$CONDA_PREFIX/lib:

    python mesh_crown_lc.py 0.30          # → cache_crown.npz at LC=0.30
    python mesh_crown_lc.py 0.55          # → cache_crown.npz at LC=0.55

Use bash run_crown_hconvergence.sh to sweep the standard h-convergence ladder.
"""
import sys

import mesh_crown

if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("usage: python mesh_crown_lc.py <LC_mm>")
    lc = float(sys.argv[1])
    mesh_crown.LC = lc
    print(f"[mesh_crown_lc] overriding LC = {lc} (headline default {0.40})")
    mesh_crown.main()
