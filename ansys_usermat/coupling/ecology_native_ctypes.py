"""ecology_native_ctypes.py -- the Fortran point model of the ANSYS runs (ecology_native.f), called from Python
through ecology_native_capi.f90 and ctypes. Built on first use with gfortran into this folder's
libecology_native.so (not committed). The constants come from the file BIOFILM_ECO_CASE, read by the Fortran
code itself, as in ANSYS.

    import ecology_native_ctypes as F
    F.init("abaqus_composition/eco_fritsch_2sp_a50.txt")     # sets BIOFILM_ECO_CASE, builds and loads
    Gnew, ok = F.step_lin(G, dt, n_sub, c0, c1, a0, a1)     # G: (n, 12); dt = s dT* per increment

ECO_THREADS=k (environment, default 1): step_lin splits the points into k chunks and runs them in k threads (ctypes
releases the GIL; the Fortran only reads its module constants after init, locals on the stack via -frecursive).
Each point is computed alone, so the result is the same to the last bit for any k.
"""
import ctypes
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
LIB = HERE / "libecology_native.so"
_lib = None
CONST = {}


def build(force=False):
    src = [HERE / "ecology_native.f", HERE / "ecology_native_capi.f90"]
    if LIB.exists() and not force and LIB.stat().st_mtime > max(p.stat().st_mtime for p in src + [Path(__file__)]):
        return LIB
    tmp = HERE / "_build_ecology_native"
    tmp.mkdir(exist_ok=True)
    subprocess.run(["gfortran", "-O2", "-frecursive", "-fPIC", "-shared", "-J", str(tmp), *map(str, src), "-o", str(LIB)],
                   check=True, cwd=tmp)
    return LIB


def init(case_file):
    """Point the Fortran code at the point-model file, build and load the library, return its constants."""
    global _lib
    os.environ["BIOFILM_ECO_CASE"] = str(Path(case_file).resolve())
    build()
    _lib = ctypes.CDLL(str(LIB))
    D, I = ctypes.c_double, ctypes.c_int
    nact, nw = I(), I()
    cs, als, tol = D(), D(), D()
    th = np.zeros(20)
    _lib.eco_capi_init(ctypes.byref(nact), ctypes.byref(cs), ctypes.byref(als), ctypes.byref(nw),
                       ctypes.byref(tol), th.ctypes.data_as(ctypes.POINTER(D)))
    CONST.update(nact=nact.value, c_star=cs.value, alpha_star=als.value, newton_max=nw.value,
                 newton_tol=tol.value, theta=th)
    return dict(CONST)


THREADS = max(1, int(os.environ.get("ECO_THREADS", "1")))
_pool = None


def _ptr(a, off, T=ctypes.c_double):
    return ctypes.cast(ctypes.c_void_p(a.ctypes.data + int(off)), ctypes.POINTER(T))


def _lin(G, th, dt, n_sub, arr, out, ok, lo, hi):
    m = int(hi - lo)
    if m <= 0:
        return
    _lib.eco_capi_lin(ctypes.c_int(m), _ptr(G, 96 * lo), _ptr(th, 0), ctypes.c_double(dt), ctypes.c_int(int(n_sub)),
                      *[_ptr(x, 8 * lo) for x in arr], _ptr(out, 96 * lo), _ptr(ok, 4 * lo, ctypes.c_int))


def step_lin(G, dt, n_sub, c0, c1, a0, a1):
    """One increment for n points: ANSYS's biofilm_ecology_hook_lin per point. Returns (Gnew (n, 12), ok (n,))."""
    global _pool
    G = np.asfortranarray(np.asarray(G, float).T)            # (12, n)
    n = G.shape[1]
    arr = [np.ascontiguousarray(np.broadcast_to(np.asarray(x, float), (n,))) for x in (c0, c1, a0, a1)]
    out = np.zeros((12, n), order="F")
    ok = np.zeros(n, dtype=np.int32)
    th = np.ascontiguousarray(CONST["theta"])
    k = THREADS if n >= 256 * THREADS else 1
    if k == 1:
        _lin(G, th, dt, n_sub, arr, out, ok, 0, n)
    else:
        if _pool is None:
            _pool = ThreadPoolExecutor(THREADS)
        cut = np.linspace(0, n, k + 1).astype(int)
        list(_pool.map(lambda j: _lin(G, th, dt, n_sub, arr, out, ok, cut[j], cut[j + 1]), range(k)))
    return out.T.copy(), ok
