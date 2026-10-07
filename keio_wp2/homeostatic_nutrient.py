#!/usr/bin/env python3
"""Keio WP2, step 7: the nutrient and the front term on the implant collar
and the tooth crown, with the homeostatic-pressure growth law.

Geometry, mesh and mechanics as homeostatic_tooth.py (bottom u_z = 0). The
nutrient enters through the top face (z = 2 mm, the gingival margin) for the
implant, as abaqus_composition/make_implant_inp.py, and through the outer
face of the layer (saliva) for the tooth, as its --nut outer. Values as the
Abaqus front-term runs (abaqus_composition/README.md): Klempt 2024 Table 2
converted to the 2 mm cube (KLEMPT2024_REPRODUCTION.md sec. 8):

    0 = d lap(c) - g phi c,  c = 1 on the nutrient face        (Eq. 35 quasi-static, consumption g phi c)
    phi_dot = beta lap(phi) + k_alpha alpha
              + r c/(k+c) |grad phi| [w + (1 - w) |n_phi . n_c|]  (front term of Eq. 34)
    alpha_dot = k_alpha phi max(0, 1 - p / p_h)

d = 1 mm^2/T*, g = 1/T*, r = 10 mm/T*, k = 1, beta = 0.02 mm^2/T*. w = 0.5
(growth on every face of the colony) is my modification (KLEMPT2024_REPRODUCTION.md
sec. 15); w = 0 is the form with growth only along the nutrient gradient. phi
is clipped to [0, 1]. The front term is evaluated per element from the
element gradients (not upwinded) and spread to the nodes like the alpha
source; explicit sub-steps keep r dt / h below 0.4.

    python keio_wp2/homeostatic_nutrient.py -> keio_wp2/results_nutrient.json
    python keio_wp2/homeostatic_nutrient.py --bc -> keio_wp2/results_nutrient_bc.json
      bottom face free, u_z = 0 or clamped, without feedback and with P_h = 0.1
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import homeostatic_tooth as T  # noqa: E402

D_C, R_FRONT, K_MONOD, BETA = 1.0, 10.0, 1.0, 0.02


class NMesh(T.Mesh):
    def __init__(self, bulge, nut, bottom="uz"):
        super().__init__(bulge, bottom=bottom)
        nnz = self.nz + 1
        n = np.arange(self.nn)
        I, K = n // nnz, n % nnz
        self.cface = np.where(K == self.nz)[0] if nut == "top" else np.where(I == self.nr)[0]
        self.cfree = np.setdiff1d(n, self.cface)
        # element gradient operator at the centre (ne, 4 nodes, 2)
        ne = len(self.conn)
        self.G0 = np.zeros((ne, 4, 2))
        for e in range(ne):
            X = self.xn[self.conn[e]]
            _, dN = T.shape(0.0, 0.0)
            self.G0[e] = dN @ np.linalg.inv(dN.T @ X).T
        self.h_min = min(T.T_LAYER / self.nr, T.HZ / self.nz)

    def grad(self, fn):
        return np.einsum("eij,ei->ej", self.G0, fn[self.conn])

    def nutrient(self, phin, g):
        A = (D_C * self.Kphi + sp.diags(g * phin * self.Mphi)).tocsr()
        c = np.zeros(self.nn); c[self.cface] = 1.0
        rhs = -A[:, self.cface] @ c[self.cface]
        c[self.cfree] = spla.spsolve(A[self.cfree][:, self.cfree].tocsc(), rhs[self.cfree])
        return c


def run(bulge, nut, P_h, local=True, w=0.5, g=1.0, dt=0.01, t_end=1.0, bottom="uz"):
    M = NMesh(bulge, nut, bottom)
    phin = M.seed_n.copy()
    alpha = np.ones(len(M.conn))
    nsub = int(np.ceil(dt / min(1.8 / (BETA * M.lam), 0.4 * M.h_min / R_FRONT)))
    V = M.vol / M.vol.sum()
    hist, t_fill = [], None
    for k in range(int(round(t_end / dt))):
        phi = M.phi_e(phin)
        p, vm, stt, tn, ts = T.mechanics(M, phi, alpha - 1)
        if P_h is None:
            gate = np.ones_like(p)
        else:
            p_h = P_h * T.E0 * T.K_ALPHA * ((phi ** 2 + T.FLOOR) if local else 1.0)
            gate = np.clip(1 - p / p_h, 0.0, 1.0)
        c = M.nutrient(phin, g)
        ce = M.phi_e(c)
        gc = M.grad(c)
        nc = gc / np.maximum(np.linalg.norm(gc, axis=1), 1e-12)[:, None]
        src_a = T.K_ALPHA * M.source(alpha)
        for _ in range(nsub):
            gp = M.grad(phin)
            ng = np.linalg.norm(gp, axis=1)
            dot = np.abs(np.einsum("ej,ej->e", gp / np.maximum(ng, 1e-12)[:, None], nc))
            front = R_FRONT * ce / (K_MONOD + ce) * ng * (w + (1 - w) * dot)
            rate = -BETA * (M.Kphi @ phin) / M.Mphi + src_a + M.source(front)
            phin = np.clip(phin + dt / nsub * rate, 0, 1)
        phi1 = M.phi_e(phin)
        alpha = alpha + dt * T.K_ALPHA * phi1 * gate
        pm = float(V @ phi1)
        if t_fill is None and pm > 0.9:
            t_fill = (k + 1) * dt
        hist.append(((k + 1) * dt, float(V @ (alpha - 1)), pm, float(p.max()), float(vm.max()),
                     float(stt.min()), float(tn.max()), float(np.abs(ts).max()), float(V @ gate),
                     float(c.min()), float(V @ ce)))
    zc = (np.arange(M.nz) + 0.5) * T.HZ / M.nz
    prof = dict(z=zc.tolist(), alpha=(alpha - 1).reshape(M.nr, M.nz)[0].tolist(), stt=stt.tolist(),
                tn=tn.tolist(), ts=ts.tolist(), phi=M.phi_e(phin).reshape(M.nr, M.nz).tolist(),
                c=ce.reshape(M.nr, M.nz).tolist())
    return dict(bulge=bulge, nut=nut, P_h=P_h, local=local, w=w, g=g, bottom=bottom, t_fill=t_fill, hist=hist, field=prof)


def main():
    t0 = time.time()
    cases = []
    for bulge, nut in ((0.0, "top"), (1.0, "outer")):
        cases += [(bulge, nut, None, True, 0.5, 1.0), (bulge, nut, 0.1, False, 0.5, 1.0),
                  (bulge, nut, 0.1, True, 0.5, 1.0), (bulge, nut, None, True, 0.0, 1.0),
                  (bulge, nut, None, True, 0.5, 6.0)]
    out = {"runs": []}
    for bulge, nut, P_h, local, w, g in cases:
        r = run(bulge, nut, P_h, local, w, g)
        out["runs"].append(r)
        t, d, ph, px, vm, stt, tn, ts, gt, cmin, cm = r["hist"][-1]
        mid = -np.array(r["field"]["stt"])[len(r["field"]["stt"]) // 2]
        form = "none " if P_h is None else ("local" if local else "E    ")
        print(f"{'tooth  ' if bulge else 'implant'} nut={nut:5s} w={w} g={g:g} P_h={P_h} {form}: alpha-1 {d:.3e}, "
              f"phi {ph:.3f} (>0.9 at T*={r['t_fill']}), c min {cmin:.2f} mean {cm:.2f}, p max {px:.2e}, "
              f"hoop max {-stt:.2e} mid {mid:.2e}, normal {tn:.2e}, shear {ts:.2e} Pa, gate {gt:.2f} "
              f"[{time.time() - t0:.0f} s]", flush=True)
    (HERE / "results_nutrient.json").write_text(json.dumps(out))


def bc_sweep():
    t0 = time.time()
    out = {"runs": []}
    for bottom in ("free", "uz", "clamped"):
        for bulge, nut in ((0.0, "top"), (1.0, "outer")):
            for P_h in (None, 0.1):
                r = run(bulge, nut, P_h, bottom=bottom)
                out["runs"].append(r)
                t, d, ph, px, vm, stt, tn, ts, gt, cmin, cm = r["hist"][-1]
                f = r["field"]
                a = -np.array(f["stt"])
                print(f"bottom={bottom:7s} {'tooth  ' if bulge else 'implant'} P_h={P_h}: alpha-1 {d:.3e}, "
                      f"p max {px:.2e}, hoop max {a.max():.2e} at z={f['z'][int(a.argmax())]:.2f} "
                      f"mid {a[len(a) // 2]:.2e}, normal {tn:.2e}, shear {ts:.2e} Pa [{time.time() - t0:.0f} s]",
                      flush=True)
    (HERE / "results_nutrient_bc.json").write_text(json.dumps(out))


if __name__ == "__main__":
    bc_sweep() if "--bc" in sys.argv else main()
