C     usermat_endtoend_driver.f — calls the REAL `usermat` subroutine itself
C     (not BIOFILM_STRESS_CORE), toggling kUsePy via prop(6). This is the
C     proof that the KUSEPY branch in usermat_biofilm.f actually works end
C     to end: Fortran usermat() -> biofilm_py_bridge -> biofilm_py_eval.c ->
C     socket -> material_server.py -> back through the Abaqus->ANSYS Voigt
C     reindex -> stress(ncomp)/dsdePl/ustatev, matching kUsePy=0 (the
C     verified inline core) to numerical precision.
C
C     Every argument usermat() never reads is passed as a sized-but-unused
C     dummy -- this driver exists to exercise the real ABI, not to model an
C     element.
C
C     stdin: F(3x3), Fv_old(3x3), alpha, C10, C01, D1, eta, mtype, dt, kUsePy
C            kStateMat, C10s, C01s, D1s, etas
C              (the trailing line drives prop(7)/ustatev(11:14), the
C               composition-dependent per-IP material path; pass
C               "0 0 0 0 0" to leave it disabled)
C     stdout: stress(6) [ANSYS order], Fv_new(3x3, from ustatev), keycut,
C             cutFactor, dsdePl(6x6) [ANSYS order, row-major]
C
C     Link: gfortran usermat_endtoend_driver.f usermat_py_hook.f
C           ../usermat_biofilm.f biofilm_py_eval.c -o driver
      PROGRAM USERMAT_E2E
      USE biofilm_py_bridge, only: biofilm_py_hook
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INTEGER, PARAMETER :: NCOMP=6, NSTATEV=14, NPROP=7
      DOUBLE PRECISION STRESS(NCOMP), USTATEV(NSTATEV),
     &                 DSDEPL(NCOMP,NCOMP), STRAIN(NCOMP),
     &                 DSTRAIN(NCOMP), EPSPL(NCOMP), PROP(NPROP),
     &                 COORDS(3), DEFGRAD_T(3,3), DEFGRAD(3,3),
     &                 TSSTIF(2), FV_OLD(3,3)
      DOUBLE PRECISION VAR0, VAR1, VAR2, VAR3, VAR4, VAR5, VAR6,
     &                 VAR7, VAR8
      INTEGER MATID, ELEMID, KDOMINTPT, KLAYER, KSECTPT,
     &        LDSTEP, ISUBST, KEYCUT, NDIRECT, NSHEAR
      DOUBLE PRECISION TIME, DTIME, TEMP, DTEMP, SEDEL, SEDPL, EPSEQ,
     &                 EPSZZ, CUTFACTOR
      DOUBLE PRECISION ALPHA, C10, C01, D1, ETA, MTYPE, DT, KUSEPY
      DOUBLE PRECISION KSTMAT, C10S, C01S, D1S, ETAS
      INTEGER I, J, K

      READ(*,*) ((DEFGRAD(I,J),J=1,3),I=1,3)
      READ(*,*) ((FV_OLD(I,J),J=1,3),I=1,3)
      READ(*,*) ALPHA, C10, C01, D1, ETA, MTYPE, DT, KUSEPY
      READ(*,*) KSTMAT, C10S, C01S, D1S, ETAS

      MATID=1; ELEMID=1; KDOMINTPT=1; KLAYER=1; KSECTPT=1
      LDSTEP=1; ISUBST=1; KEYCUT=0
      NDIRECT=3; NSHEAR=3
      TIME=0.0D0; DTIME=DT; TEMP=0.0D0; DTEMP=0.0D0
      SEDEL=-999.0D0; SEDPL=-999.0D0; EPSEQ=0.0D0
      EPSZZ=0.0D0; CUTFACTOR=1.0D0
      VAR0=0.0D0; VAR1=0.0D0; VAR2=0.0D0; VAR3=0.0D0; VAR4=0.0D0
      VAR5=0.0D0; VAR6=0.0D0; VAR7=0.0D0; VAR8=0.0D0
      DO I=1,3
        COORDS(I)=0.0D0
        DO J=1,3
          DEFGRAD_T(I,J)=0.0D0
        END DO
        DEFGRAD_T(I,I)=1.0D0
      END DO
      DO I=1,NCOMP
        STRESS(I)=0.0D0; STRAIN(I)=0.0D0; DSTRAIN(I)=0.0D0
        EPSPL(I)=0.0D0
        DO J=1,NCOMP
          DSDEPL(I,J)=0.0D0
        END DO
      END DO
      TSSTIF(1)=0.0D0; TSSTIF(2)=0.0D0

      PROP(1)=C10; PROP(2)=C01; PROP(3)=D1; PROP(4)=ETA
      PROP(5)=MTYPE; PROP(6)=KUSEPY; PROP(7)=KSTMAT

C     ustatev(1:9) = Fv_old (row-major); ustatev(10) = alpha. usermat()
C     reads these on entry and overwrites 1:9 with Fv_new on return.
C     ustatev(11:14) = per-IP C10,C01,D1,eta, read only when prop(7)>0.5.
      K=0
      DO I=1,3
        DO J=1,3
          K=K+1
          USTATEV(K)=FV_OLD(I,J)
        END DO
      END DO
      USTATEV(10)=ALPHA
      USTATEV(11)=C10S; USTATEV(12)=C01S; USTATEV(13)=D1S; USTATEV(14)=ETAS

      CALL USERMAT(
     &   MATID, ELEMID, KDOMINTPT, KLAYER, KSECTPT,
     &   LDSTEP, ISUBST, KEYCUT,
     &   NDIRECT, NSHEAR, NCOMP, NSTATEV, NPROP,
     &   TIME, DTIME, TEMP, DTEMP,
     &   STRESS, USTATEV, DSDEPL, SEDEL, SEDPL, EPSEQ,
     &   STRAIN, DSTRAIN, EPSPL, PROP, COORDS,
     &   VAR0, DEFGRAD_T, DEFGRAD, TSSTIF, EPSZZ, CUTFACTOR,
     &   VAR1, VAR2, VAR3, VAR4, VAR5, VAR6, VAR7, VAR8)

      WRITE(*,'(6E26.17)') (STRESS(K),K=1,NCOMP)
      WRITE(*,'(9E26.17)') (USTATEV(K),K=1,9)
      WRITE(*,'(I4,1X,E26.17)') KEYCUT, CUTFACTOR
      DO I=1,NCOMP
        WRITE(*,'(6E26.17)') (DSDEPL(I,J),J=1,NCOMP)
      END DO
      END
