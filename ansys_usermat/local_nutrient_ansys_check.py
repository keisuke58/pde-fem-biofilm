#!/usr/bin/env python3
"""Check of the ANSYS step-1 run (5 Oct 2026, 8^3, case 6, consumption 4):
nutrient in the seed 0.435-0.656, share phi_1/(phi_1+phi_2) 0.025-0.142.

Python, same model and coupling step (0.025): the consumption is set so that
the lowest nutrient in the seed is 0.435, once with zero-order consumption as
in the element (d lap c = g phi) and once with first-order consumption as in
composition_local_nutrient.py (d lap c = g phi c).

Result:
  8^3   zero order  c 0.435-0.637, share 0.026-0.139 (corr. with c -0.75)
        first order c 0.435-0.611, share 0.026-0.139
  16^3  zero order  c 0.435-0.729, share 0.024-0.139
        first order c 0.435-0.705, share 0.025-0.139
Python gives the ANSYS range to within 0.003, with either order. The share
depends on the local nutrient only; the order of the consumption matters only
in how low c falls. The wider range of the Python figure (0.03-0.44 at
Lambda = 4) comes from a lower nutrient (down to about 0.2), not from the order.

    python ansys_usermat/local_nutrient_ansys_check.py
"""
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent))
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla
import composition_local_nutrient as L
def nut0(S,D,g):
    N,H=L.N,L.H; n=N**3; ix=np.arange(n).reshape((N,)*3); R=[];C=[];V=[]; b=np.zeros(n)
    for i in range(N):
     for j in range(N):
      for k in range(N):
        r=ix[i,j,k]
        if D[i,j,k]: R+=[r];C+=[r];V+=[1.0]; b[r]=1; continue
        dg=0.0
        for di,dj,dk in ((1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)):
          a,bb,cc=i+di,j+dj,k+dk
          if 0<=a<N and 0<=bb<N and 0<=cc<N: R+=[r];C+=[ix[a,bb,cc]];V+=[-1.0]; dg+=1
        R+=[r];C+=[r];V+=[dg]; b[r]=-g*H*H*S[i,j,k]
    return spla.spsolve(sp.csr_matrix((V,(R,C)),shape=(n,n)).tocsc(),b).reshape((N,)*3)
L.T.DT=0.025
for N in (8,16):
    L.N,L.H=N,2.0/N; S,D=L.grid_masks()
    # find g so that min c in the seed = 0.435 (ANSYS consumption 4)
    lo,hi=0.0,50.0
    for _ in range(50):
        g=(lo+hi)/2; c=nut0(S,D,g)[S]
        lo,hi=(g,hi) if c.min()>0.435 else (lo,g)
    c=nut0(S,D,g)[S]; sh=L.shares("2sp_case6",np.clip(c,0,1))
    print(N,'g',round(g,3),'c',c.min().round(3),c.max().round(3),'chi',sh.min().round(3),sh.max().round(3),'corr',np.corrcoef(c,sh)[0,1].round(2))
    # first order with the same c range
    lo,hi=0.0,10.0
    for _ in range(50):
        lam=(lo+hi)/2; c1=L.nutrient(S,D,lam)[S]
        lo,hi=(lam,hi) if c1.min()>0.435 else (lo,lam)
    sh1=L.shares("2sp_case6",c1); print('  first order Lambda',round(lam,3),'c',c1.min().round(3),c1.max().round(3),'chi',sh1.min().round(3),sh1.max().round(3))
