#!/usr/bin/env python3
"""Seed stress with phi and alpha at the 2 x 2 x 2 Gauss points, as the
partner element carries them (wp3_stress_reference.py uses one value per
element).

The reference field (wp3_reference.solve on an NF^3 finite-volume grid) is
interpolated to the Gauss points of each hexahedral mesh; stiffness
E(phi) = (phi^2 + f) E_bio and eigenstrain (alpha - 1) I act per Gauss point;
mean-dilatation B-bar as SOLID185 with KEYOPT(2) = 0. The element value is the
mean of the eight Gauss-point stresses. On 8^3 this explains the gap between
the element's runs and the element-mean solve: inside a 0.25 mm element the
Gauss-point alpha smooths the step at the seed surface. The element-mean
series approaches the limit from above and this one from below, so the two
bracket the converged seed stress.

    python ansys_usermat/apdl/wp3_stress_gauss.py [NF] [meshes, e.g. 8,16,24,48]
"""
import sys, json, time
from pathlib import Path
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla
from scipy.interpolate import RegularGridInterpolator
R = Path(__file__).resolve().parent; sys.path.insert(0, str(R))
import mesh_study_seed as S
import wp3_reference as W
import wp3_stress_reference as SR

E_BIO, NU, FLOOR = 10.0, 0.49, 1e-3
cells8 = W.seed_cells_8(R / 'results/2026-10-wp3_fix/wp3c_n8_dt00125.json')
NF = int(sys.argv[1]) if len(sys.argv) > 1 else 96
axes, am1_f, phi_f = W.solve(NF, cells8, 0.15, verbose=False)
fa = RegularGridInterpolator(axes, am1_f, bounds_error=False, fill_value=None)
fp = RegularGridInterpolator(axes, phi_f, bounds_error=False, fill_value=None)

def gp_mats(h, nu=NU):
    g = np.array([-1, 1]) / np.sqrt(3)
    corners = np.array([[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],[-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]], float)
    def bmat(xi):
        dN = np.array([[cx*(1+xi[1]*cy)*(1+xi[2]*cz), cy*(1+xi[0]*cx)*(1+xi[2]*cz), cz*(1+xi[0]*cx)*(1+xi[1]*cy)]
                       for cx, cy, cz in corners]) / 8 * 2 / h
        b = np.zeros((6, 24))
        for a in range(8):
            dx, dy, dz = dN[a]; i = 3*a
            b[0,i]=dx; b[1,i+1]=dy; b[2,i+2]=dz; b[3,i]=dy; b[3,i+1]=dx; b[4,i+1]=dz; b[4,i+2]=dy; b[5,i]=dz; b[5,i+2]=dx
        return b
    lam = nu/((1+nu)*(1-2*nu)); mu = 1/(2*(1+nu))
    D = np.diag([2*mu]*3+[mu]*3); D[:3,:3] += lam
    m = np.array([1,1,1,0,0,0.]); P = np.outer(m, m)/3
    xis = [(a, b, c) for a in g for b in g for c in g]
    Bs = [bmat(x) for x in xis]
    Bvol = sum(P @ B for B in Bs) / 8                     # mean dilatation
    Bbar = [B - P @ B + Bvol for B in Bs]
    w = (h/2)**3
    K = np.array([w * Bb.T @ D @ Bb for Bb in Bbar]); f = np.array([w * Bb.T @ D @ m for Bb in Bbar])
    return np.array(xis), np.array(Bbar), K, f, D, m

def solve(n, gp=True):
    h = 2.0/n
    xis, Bbar, Kg, fg, D, m = gp_mats(h)
    idx = np.arange((n+1)**3).reshape(n+1, n+1, n+1)
    off = [(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)]
    conn = np.stack([idx[i:i+n, j:j+n, k:k+n].ravel() for i, j, k in off], 1)
    dofs = np.concatenate([3*conn[:, [a]] + np.arange(3) for a in range(8)], 1)
    ax = -1 + h*(np.arange(n)+0.5)
    cx, cy, cz = np.meshgrid(ax, ax, ax, indexing='ij'); cen = np.stack([cx.ravel(), cy.ravel(), cz.ravel()], 1)
    if gp:
        pts = cen[:, None, :] + xis[None] * h/2               # (ne, 8, 3)
        phi = fp(pts.reshape(-1, 3)).reshape(-1, 8); a = fa(pts.reshape(-1, 3)).reshape(-1, 8)
    else:                                                      # element means, as before
        f = NF // n
        phi = np.repeat(SR_block(phi_f, n).ravel()[:, None], 8, 1); a = np.repeat(SR_block(am1_f, n).ravel()[:, None], 8, 1)
    E = E_BIO * (np.clip(phi, 0, 1)**2 + FLOOR)
    Ke = np.einsum('eg,gij->eij', E, Kg); Fe = np.einsum('eg,gi->ei', E*a, fg)
    ndof = 3*(n+1)**3
    K = sp.csr_matrix((Ke.ravel(), (np.repeat(dofs, 24, 1).ravel(), np.tile(dofs, (1, 24)).ravel())), shape=(ndof, ndof))
    F = np.zeros(ndof); np.add.at(F, dofs.ravel(), Fe.ravel())
    node = lambda p: idx[tuple(int(round((c+1)/h)) for c in p)]
    fix = [3*node(S.CORNERS[0])+d for d in (0,1,2)] + [3*node(S.CORNERS[1])+d for d in (1,2)] + [3*node(S.CORNERS[2])+2]
    free = np.setdiff1d(np.arange(ndof), fix); u = np.zeros(ndof)
    Kff = K[free][:, free].tocsc()
    if n < 32:
        u[free] = spla.spsolve(Kff, F[free])
    else:
        dinv = 1.0/Kff.diagonal(); M = spla.LinearOperator(Kff.shape, matvec=lambda x: dinv*x)
        x, info = spla.cg(Kff, F[free], M=M, rtol=1e-12, maxiter=100000)
        print('   cg info', info, 'residual', np.linalg.norm(Kff@x-F[free])/np.linalg.norm(F[free]), flush=True)
        u[free] = x
    eps = np.einsum('gij,ej->egi', Bbar, u[dofs]) - a[:, :, None]*m           # (ne, 8, 6)
    sig = E[:, :, None] * np.einsum('ij,egj->egi', D, eps)
    def vm(s):
        p = s[..., :3].mean(-1); d = s[..., :3] - p[..., None]
        return np.sqrt(1.5*((d**2).sum(-1) + 2*(s[..., 3:]**2).sum(-1))), p
    vm_avg_comp, p_e = vm(sig.mean(1))          # components averaged over the 8 points, then vM
    vm_gp, _ = vm(sig); vm_mean_gp = vm_gp.mean(1)
    seed = W.seed_mask(n, cells8).ravel()
    return vm_avg_comp[seed].mean(), vm_mean_gp[seed].mean(), p_e[seed].mean()

def SR_block(field, n):
    return W.block_mean(field, n)

nem = SR.nem_seed_vm()
for n in [int(v) for v in (sys.argv[2].split(',') if len(sys.argv) > 2 else ['8','16','24'])]:
    t = time.time()
    a1, b1, p1 = (solve(n, gp=False) if n < 32 else (np.nan,)*3); a2, b2, p2 = solve(n, gp=True)
    e = nem.get(n, {'vM dt->0': float('nan'), 'p dt->0': float('nan')})
    print(f"{n:>2}^3  element means: vM {a1:.4e} p {p1:+.4e} | Gauss points: vM(avg comp) {a2:.4e} vM(avg of GP vM) {b2:.4e} p {p2:+.4e}"
          f" | ANSYS dt->0 vM {e['vM dt->0']:.4e} p {e['p dt->0']:+.4e}  ({time.time()-t:.0f} s)", flush=True)
