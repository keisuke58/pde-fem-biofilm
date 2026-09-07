C     usermat_ecology_e2e_driver.f — calls the REAL `usermat` subroutine to
C     exercise the ecology (0D Hamilton ODE) hook end to end, the same way
C     usermat_endtoend_driver.f exercises kUsePy: Fortran usermat() ->
C     biofilm_py_bridge -> biofilm_ecology_eval (biofilm_py_eval.c) ->
C     socket -> material_server.py's evaluate_ecology -> back through
C     usermat_biofilm.f's own alpha update.
C
C     stdin: F(3x3), Fv_old(3x3), alpha, C10, C01, D1, eta, mtype, dt
C            kUseEcology, k_alpha
C            theta(20)
C            g_old(12)
C     stdout: stress(6) [ANSYS order], ustatev(1:26), keycut, cutFactor
C
C     Link: gfortran usermat_ecology_e2e_driver.f usermat_py_hook.f
C           ../usermat_biofilm.f biofilm_py_eval.c -o eco_driver
      PROGRAM USERMAT_ECO_E2E
      USE biofilm_py_bridge, only: biofilm_py_hook
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INTEGER, PARAMETER :: NCOMP=6, NSTATEV=26, NPROP=29
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
      DOUBLE PRECISION ALPHA, C10, C01, D1, ETA, MTYPE, DT
      DOUBLE PRECISION KUSEECO, KALPHA, THETA(20), GOLD(12)
      INTEGER I, J, K

      READ(*,*) ((DEFGRAD(I,J),J=1,3),I=1,3)
      READ(*,*) ((FV_OLD(I,J),J=1,3),I=1,3)
      READ(*,*) ALPHA, C10, C01, D1, ETA, MTYPE, DT
      READ(*,*) KUSEECO, KALPHA
      READ(*,*) (THETA(I),I=1,20)
      READ(*,*) (GOLD(I),I=1,12)

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
      PROP(5)=MTYPE; PROP(6)=0.0D0; PROP(7)=0.0D0
      PROP(8)=KUSEECO; PROP(9)=KALPHA
      DO I=1,20
        PROP(9+I)=THETA(I)
      END DO

C     ustatev(1:9)=Fv_old, ustatev(10)=alpha, ustatev(11:14) unused here
C     (kStateMat=0), ustatev(15:26)=g_old (the ecology state).
      DO I=1,NSTATEV
        USTATEV(I)=0.0D0
      END DO
      K=0
      DO I=1,3
        DO J=1,3
          K=K+1
          USTATEV(K)=FV_OLD(I,J)
        END DO
      END DO
      USTATEV(10)=ALPHA
      DO I=1,12
        USTATEV(14+I)=GOLD(I)
      END DO

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
      WRITE(*,'(6E26.17)') (USTATEV(K),K=1,NSTATEV)
      WRITE(*,'(I4,1X,E26.17)') KEYCUT, CUTFACTOR
      END
