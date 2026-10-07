#!/usr/bin/env python3
"""Klempt 2024, Table 3: which unit of mu matches the plotted pressure?

Table 2 gives mu = 3.3557 Pa (E = 10 Pa, nu = 0.49); Table 3 plots the
hydrostatic stress p of test case 4.1 in MPa, legend +3e-4 ... -6.5e-4 MPa.
Here the growth field of the fitted 4.1 run (klempt2024_timescale.py: Table 2,
edge-line source, s = 1.5) is put into a small-strain linear-elastic solve:
hex8 on the 1 um grid, E_elem = E max(phi_elem, 1e-6) (no stiffness in the
void, as in the paper), nu = 0.49 with selective reduced integration,
eigenstrain (alpha - 1) I, free cube held against rigid-body motion only.
The result is linear in E, so it is printed as p/E.

    python JAXFEM/klempt2024_pressure_scale.py     (about a minute)

Result (2026-10-03), p/E (tension > 0):
  T* = 0.05   max(alpha-1) 7.5e-5   p/E  -3.4e-5 ... +1.5e-5
  T* = 0.25                3.8e-4        -1.9e-4 ... +7.4e-5
  T* = 0.45                6.8e-4        -3.1e-4 ... +2.3e-4
  T* = 1.00                1.5e-3        -6.8e-4 ... +5.3e-4
  - With E = 10 Pa the pressure is at most 7e-3 Pa = 7e-9 MPa: the whole of
    Table 3 would be one colour on a legend of +3e-4 ... -6.5e-4 MPa.
  - With E = 10 read as MPa (ANSYS's micro-MKS system, um with MPa, where a
    typed 3.3557 is 3.3557 MPa) the pressure at T* = 0.05 is -3.4e-4 ...
    +1.5e-4 MPa, the legend's order, and later times exceed it, which matches
    the saturated interior (blue) and ring (red) in Table 3's later rows.
  - So the plotted stresses fit mu = 3.3557 MPa, not Pa, to within a factor
    of about two, against a factor of 1e6. A question for the authors; it
    matters wherever Klempt 2024's stresses are quoted in Pa.
"""

import sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla
import klempt2024_case1_bc as B, klempt2024_quantitative as K
s = 1.5
K.R, K.BETA, K.K_A = 100*s, 2*s, 1e-3*s
p0,_,_=K.setup("fig4_edge"); solve=B.FirstOrderC(B.nutrient_mask("edge"),1e8)
phi=p0.copy(); al=np.ones_like(phi); c=solve(phi); dt=1e-3; snaps={}
for st in range(1,1001):
    gx,gy,gz=K.grad_c(c); mag=np.sqrt(gx**2+gy**2+gz**2)
    spd=np.where(mag>1e-14,K.R*c/(K.K_M+c)/np.maximum(mag,1e-14),0)
    phi=np.clip(phi+dt*(K.BETA*K.lap(phi)+K.K_A*al+np.abs(K.upwind_dot(phi,(spd*gx,spd*gy,spd*gz)))),0,1)
    al=al+dt*K.K_A*phi; c=solve(phi)
    if st in (50,250,450,1000): snaps[st/1000]=(phi.copy(),al.copy())
# hex8 element on h=1
N=K.N; ne=N-1; h=1.0
g=np.array([-1,1])/np.sqrt(3)
corners=np.array([[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],[-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]],float)
def Bmat(xi):
    dN=np.array([[cx*(1+xi[1]*cy)*(1+xi[2]*cz),cy*(1+xi[0]*cx)*(1+xi[2]*cz),cz*(1+xi[0]*cx)*(1+xi[1]*cy)] for cx,cy,cz in corners])/8*2/h
    Bm=np.zeros((6,24))
    for a in range(8):
        dx,dy,dz=dN[a]; i=3*a
        Bm[0,i]=dx;Bm[1,i+1]=dy;Bm[2,i+2]=dz
        Bm[3,i]=dy;Bm[3,i+1]=dx;Bm[4,i+1]=dz;Bm[4,i+2]=dy;Bm[5,i]=dz;Bm[5,i+2]=dx
    return Bm
nu=0.49
lam=nu/((1+nu)*(1-2*nu)); mu=1/(2*(1+nu))           # per unit E
Dd=np.diag([2*mu]*3+[mu]*3); Dv=np.zeros((6,6)); Dv[:3,:3]=lam
m=np.array([1,1,1,0,0,0.])
Kd=np.zeros((24,24)); 
for a in g:
  for b in g:
    for cc in g:
      Bm=Bmat((a,b,cc)); Kd+=Bm.T@Dd@Bm*(h/2)**3
B0=Bmat((0,0,0)); Kv=B0.T@Dv@B0*h**3
# eigenstrain load per unit (alpha-1): f = int B^T D (e* m)
fd=np.zeros(24)
for a in g:
  for b in g:
    for cc in g:
      Bm=Bmat((a,b,cc)); fd+=Bm.T@Dd@m*(h/2)**3
fv=B0.T@Dv@m*h**3
idx=np.arange(N**3).reshape(N,N,N)
off=[(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)]
conn=np.stack([idx[i:i+ne,j:j+ne,k:k+ne].ravel() for i,j,k in off],1)
dofs=np.concatenate([3*conn[:,[a]]+np.arange(3) for a in range(8)],1)
def solve_p(phi,al,floor=1e-6):
    pe=np.mean([phi[i:i+ne,j:j+ne,k:k+ne] for i,j,k in off],0).ravel()
    ae=np.mean([al[i:i+ne,j:j+ne,k:k+ne] for i,j,k in off],0).ravel()-1
    Ee=np.maximum(pe,floor)
    Ke=(Kd+Kv)[None]*Ee[:,None,None]
    rows=np.repeat(dofs,24,1).ravel(); cols=np.tile(dofs,(1,24)).ravel()
    Kg=sp.csr_matrix((Ke.ravel(),(rows,cols)),shape=(3*N**3,)*2)
    F=np.zeros(3*N**3); np.add.at(F,dofs.ravel(),((fd+fv)[None]*(Ee*ae)[:,None]).ravel())
    fix=[0,1,2, 3*idx[-1,0,0]+1,3*idx[-1,0,0]+2, 3*idx[0,-1,0]+2]
    free=np.setdiff1d(np.arange(3*N**3),fix)
    u=np.zeros(3*N**3); u[free]=spla.spsolve(Kg[free][:,free].tocsc(),F[free])
    ue=u[dofs]
    eps=(B0@ue.T).T - ae[:,None]*m               # elastic strain at centre
    sig=Ee[:,None]*(eps@(Dd+Dv).T)
    p=sig[:,:3].mean(1)                         # hydrostatic stress, tension > 0
    return p.reshape(ne,ne,ne), pe.reshape(ne,ne,ne), ae.max()
for t,(ph,a) in snaps.items():
    t0=time.time(); p,pe,am=solve_p(ph,a)
    inb=pe>0.5
    print(f"T*={t:.2f} max(alpha-1)={am:.2e}  p/E inside biofilm: min {p[inb].min():.2e} max {p[inb].max():.2e}; all: min {p.min():.2e} max {p.max():.2e} ({time.time()-t0:.0f}s)", flush=True)
