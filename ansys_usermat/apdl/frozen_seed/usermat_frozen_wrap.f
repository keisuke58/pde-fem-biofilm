C=======================================================================
C  usermat_frozen_wrap.f -- test version for request A-16 (the seed with
C  frozen phase, ANSYS side of the cross-code check with the unified UEL,
C  p3/uel/seed_frozen_reference.py). Not for production runs.
C
C  phi is the temperature: BF,<node>,TEMP,phi on every node, so SOLID185
C  hands each Gauss point the shape-function interpolation of the nodal phi
C  (Temp + dTemp). This wrapper then
C    - advances Eq. 36 of Klempt et al. 2024 with phi frozen:
C        alpha_n+1 = alpha_n + k_alpha phi dTime   (ustatev(10) = alpha - 1)
C    - scales the stiffness by (phi^2 + f): C10 -> (phi^2 + f) C10,
C        D1 -> D1/(phi^2 + f)   (neo-Hooke deviator + (1/D1)(J_e - 1)^2,
C        i.e. K/2 (J_e - 1)^2 with D1 = 2/K)
C  and calls the unchanged thesis core (usermat_biofilm.f, renamed
C  usermat_core in the build directory only) with Fg = (1 + ustatev(10)) I,
C  no viscous branch, no hooks.
C
C  TB,USER constants:  1 C10   2 D1   3 f   4 k_alpha
C  TB,STATE (12):      1-9 Fv (core), 10 alpha - 1, 11 phi at the Gauss
C                      point, 12 phi^2 + f      (11, 12 only as a record)
C=======================================================================
      subroutine usermat(
     &   matId, elemId, kDomIntPt, kLayer, kSectPt,
     &   ldstep, isubst, keycut,
     &   nDirect, nShear, ncomp, nStatev, nProp,
     &   Time, dTime, Temp, dTemp,
     &   stress, ustatev, dsdePl, sedEl, sedPl, epseq,
     &   Strain, dStrain, epsPl, prop, coords,
     &   var0, defGrad_t, defGrad, tsstif, epsZZ, cutFactor,
     &   var1, var2, var3, var4, var5, var6, var7, var8)
      implicit none
      integer          matId, elemId, kDomIntPt, kLayer, kSectPt,
     &                 ldstep, isubst, keycut,
     &                 nDirect, nShear, ncomp, nStatev, nProp
      double precision Time, dTime, Temp, dTemp, sedEl, sedPl, epseq,
     &                 epsZZ, cutFactor
      double precision stress(ncomp), ustatev(nStatev),
     &                 dsdePl(ncomp,ncomp), Strain(ncomp),
     &                 dStrain(ncomp), epsPl(ncomp), prop(nProp),
     &                 coords(3), defGrad_t(3,3), defGrad(3,3),
     &                 tsstif(2)
      double precision var0, var1, var2, var3, var4, var5, var6,
     &                 var7, var8
      double precision PC(8), PHI, FS
      integer          I

      PHI = Temp + dTemp
      FS = PHI*PHI + prop(3)
      do I = 1, 8
        PC(I) = 0.0d0
      end do
      PC(1) = prop(1)*FS
      PC(3) = prop(2)/FS
C     ustatev comes in at the start of the increment on every iteration
      ustatev(10) = ustatev(10) + prop(4)*PHI*dTime
      ustatev(11) = PHI
      ustatev(12) = FS
      call usermat_core(
     &   matId, elemId, kDomIntPt, kLayer, kSectPt,
     &   ldstep, isubst, keycut,
     &   nDirect, nShear, ncomp, 10, 8,
     &   Time, dTime, Temp, dTemp,
     &   stress, ustatev, dsdePl, sedEl, sedPl, epseq,
     &   Strain, dStrain, epsPl, PC, coords,
     &   var0, defGrad_t, defGrad, tsstif, epsZZ, cutFactor,
     &   var1, var2, var3, var4, var5, var6, var7, var8)
      return
      end
