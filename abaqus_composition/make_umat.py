"""make_umat.py -- one Abaqus user-subroutine file with the same composition
coupling the ANSYS runs use (Keio continuation: the model moves to Abaqus).

    python abaqus_composition/make_umat.py OUT.for

The ANSYS side pastes the call-site fragments (ansys_usermat/apdl/callsite/)
into the partner's usermat. Here the same fragments, unchanged, go into a
small UMAT, so the point model, the composition (prop(28) = 7) and Eq. 36
run the same code in both programs. The file holds, in this order:
  - usermat_py_hook.f (module biofilm_py_bridge, the Fortran side of the
    Python bridge), split_rates.f (module biofilm_split, the per-point cache),
    growth_from_phi.f (Eq. 36), biofilm_material_v01.f (stress from
    F_g = (1+alpha) I with the Klempt 2024 stiffness branch) and, from
    usermat_biofilm.f, BIOFILM_STRESS_CORE and its helpers (the ANSYS
    entry point is left out);
  - UMAT, which checks the sizes and calls BIOFILM_COMP_UMAT;
  - BIOFILM_COMP_UMAT, which sets the names the fragments expect from the
    partner's usermat, includes the two fragments and turns alpha into stress.

Inputs from the Abaqus model:
  field variable 1 = phi, the amount of biofilm (Klempt 2024 phi; the ANSYS
                     runs take it from the partner's field), value at the end
                     of the increment
  field variable 2 = nutrient c (used with prop(33) = c_ref > 0)
  constants 1-42   = the ANSYS prop layout (prop(7) k_alpha, prop(8:27) theta,
                     prop(28) mode, prop(29:42) as in phi_mode_exec.inc:
                     37 number of species, 38-42 growth weights)
  constants 43-46  = E, E_void, nu, nu_void as biofilm_material_v01.f reads
                     them (E_void < 0 selects the Klempt 2024 stiffness)
  *DEPVAR >= 100   = the ANSYS state layout (72-83 point model, 84 alpha - 1),
                     85-93 the viscous state F_v
The C shim (ansys_usermat/coupling/biofilm_py_eval.c) is compiled separately
and linked through the job's abaqus_v6.env (run_comp.ps1 does both).
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
AU = ROOT / "ansys_usermat"
CS = AU / "apdl" / "callsite"
sys.path.insert(0, str(AU / "apdl"))
import paste_fragments as pf  # noqa: E402

PARTS = [AU / "coupling" / "usermat_py_hook.f", CS / "split_rates.f",
         AU / "growth_from_phi.f", AU / "biofilm_material_v01.f",
         AU / "usermat_biofilm.f"]
# usermat_biofilm.f: only from BIOFILM_STRESS_CORE on (the stress core and its
# matrix helpers); its ANSYS entry point usermat would need ANSYS symbols
FROM = {"usermat_biofilm.f": "subroutine BIOFILM_STRESS_CORE"}


def part_lines(p: Path) -> list[str]:
    lines = p.read_text(encoding="latin-1").splitlines()
    key = FROM.get(p.name)
    if key:
        lines = lines[next(i for i, l in enumerate(lines) if l.strip().startswith(key)):]
    return lines

NUT = """\
C=======================================================================
C  biofilm_nut_store: the nutrient c at the integration points, written
C  by the nutrient UEL, read by the UMAT (one iteration later; the
C  partner's element too hands the material the c of the previous
C  sub-step). Indexed by the solid element and its integration point;
C  c = 1 (MY_NUTSTART1) until the UEL has written a value.
C=======================================================================
      MODULE biofilm_nut_store
      IMPLICIT NONE
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: NUT_GP(:,:), NUT_GR(:,:,:)
      CONTAINS
      SUBROUTINE NUT_PUT(IEL, IGP, C, GX, GY, GZ)
      INTEGER IEL, IGP, N
      DOUBLE PRECISION C, GX, GY, GZ
      DOUBLE PRECISION, ALLOCATABLE :: TMP(:,:), TMG(:,:,:)
      IF (IEL .LT. 1) RETURN
      IF (.NOT. ALLOCATED(NUT_GP)) THEN
        ALLOCATE(NUT_GP(8, MAX(IEL, 1000)))
        ALLOCATE(NUT_GR(3, 8, MAX(IEL, 1000)))
        NUT_GP = 1.0D0
        NUT_GR = 0.0D0
      END IF
      N = SIZE(NUT_GP, 2)
      IF (IEL .GT. N) THEN
        ALLOCATE(TMP(8, MAX(IEL, 2 * N)))
        ALLOCATE(TMG(3, 8, MAX(IEL, 2 * N)))
        TMP = 1.0D0
        TMG = 0.0D0
        TMP(:, 1:N) = NUT_GP
        TMG(:, :, 1:N) = NUT_GR
        CALL MOVE_ALLOC(TMP, NUT_GP)
        CALL MOVE_ALLOC(TMG, NUT_GR)
      END IF
      NUT_GP(IGP, IEL) = C
      NUT_GR(1, IGP, IEL) = GX
      NUT_GR(2, IGP, IEL) = GY
      NUT_GR(3, IGP, IEL) = GZ
      END SUBROUTINE
      SUBROUTINE NUT_GRAD(IEL, IGP, G)
      INTEGER IEL, IGP
      DOUBLE PRECISION G(3)
      G = 0.0D0
      IF (.NOT. ALLOCATED(NUT_GR)) RETURN
      IF (IEL .LT. 1 .OR. IEL .GT. SIZE(NUT_GR, 3)) RETURN
      G = NUT_GR(:, IGP, IEL)
      END SUBROUTINE
      DOUBLE PRECISION FUNCTION NUT_GET(IEL, IGP)
      INTEGER IEL, IGP
      NUT_GET = 1.0D0
      IF (.NOT. ALLOCATED(NUT_GP)) RETURN
      IF (IEL .LT. 1 .OR. IEL .GT. SIZE(NUT_GP, 2)) RETURN
      NUT_GET = NUT_GP(IGP, IEL)
      END FUNCTION
      END MODULE biofilm_nut_store

C=======================================================================
C  UEL: the nutrient c (Klempt 2024 Eq. 35 as the partner's element
C  solves it: quasi-static, zero order, d lap(c) = g phi) on the nodes of
C  the C3D8T mesh, overlaid on it. DOF 12 = c; DOF 11 = phi is read only
C  (no residual, the C3D8T elements carry it). 8 nodes, trilinear,
C  2x2x2 Gauss points in Abaqus' C3D8 order. Properties: d, g, offset
C  (UEL element number - offset = the C3D8T element it overlays).
C  *USER ELEMENT, NODES=8, TYPE=U1, PROPERTIES=3, COORDINATES=3,
C  VARIABLES=1, UNSYMM  /  11, 12
C=======================================================================
      SUBROUTINE UEL(RHS, AMATRX, SVARS, ENERGY, NDOFEL, NRHS, NSVARS,
     1 PROPS, NPROPS, COORDS, MCRD, NNODE, U, DU, V, A, JTYPE, TIME,
     2 DTIME, KSTEP, KINC, JELEM, PARAMS, NDLOAD, JDLTYP, ADLMAG,
     3 PREDEF, NPREDF, LFLAGS, MLVARX, DDLMAG, MDLOAD, PNEWDT, JPROPS,
     4 NJPROP, PERIOD)
      USE biofilm_nut_store
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION RHS(MLVARX,*), AMATRX(NDOFEL,NDOFEL), PROPS(*),
     1 SVARS(*), ENERGY(8), COORDS(MCRD,NNODE), U(NDOFEL),
     2 DU(MLVARX,*), V(NDOFEL), A(NDOFEL), TIME(2), PARAMS(*),
     3 JDLTYP(MDLOAD,*), ADLMAG(MDLOAD,*), DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE), LFLAGS(*), JPROPS(*)
      DIMENSION XI(3,8), GN(8), DNX(3,8), DNL(3,8), XJ(3,3), XJI(3,3)
      DATA XI / -1.D0,-1.D0,-1.D0,  1.D0,-1.D0,-1.D0,  1.D0, 1.D0,-1.D0,
     1          -1.D0, 1.D0,-1.D0, -1.D0,-1.D0, 1.D0,  1.D0,-1.D0, 1.D0,
     2           1.D0, 1.D0, 1.D0, -1.D0, 1.D0, 1.D0 /
      D = PROPS(1)
      G = PROPS(2)
      IOFF = NINT(PROPS(3))
      DO I = 1, NDOFEL
        DO J = 1, NRHS
          RHS(I, J) = 0.0D0
        END DO
        DO J = 1, NDOFEL
          AMATRX(I, J) = 0.0D0
        END DO
      END DO
      DO I = 1, 8
        ENERGY(I) = 0.0D0
      END DO
      IF (LFLAGS(3) .NE. 1 .AND. LFLAGS(3) .NE. 2 .AND.
     1    LFLAGS(3) .NE. 5) RETURN
      GP = 1.0D0 / SQRT(3.0D0)
      IGP = 0
      DO KZ = 1, 2
      DO KY = 1, 2
      DO KX = 1, 2
        IGP = IGP + 1
        S1 = GP * (2 * KX - 3)
        S2 = GP * (2 * KY - 3)
        S3 = GP * (2 * KZ - 3)
        DO IA = 1, 8
          GN(IA) = 0.125D0 * (1 + XI(1,IA) * S1) * (1 + XI(2,IA) * S2)
     1             * (1 + XI(3,IA) * S3)
          DNL(1,IA) = 0.125D0 * XI(1,IA) * (1 + XI(2,IA) * S2)
     1                * (1 + XI(3,IA) * S3)
          DNL(2,IA) = 0.125D0 * XI(2,IA) * (1 + XI(1,IA) * S1)
     1                * (1 + XI(3,IA) * S3)
          DNL(3,IA) = 0.125D0 * XI(3,IA) * (1 + XI(1,IA) * S1)
     1                * (1 + XI(2,IA) * S2)
        END DO
        DO I = 1, 3
          DO J = 1, 3
            XJ(I, J) = 0.0D0
            DO IA = 1, 8
              XJ(I, J) = XJ(I, J) + DNL(I, IA) * COORDS(J, IA)
            END DO
          END DO
        END DO
        DET = XJ(1,1) * (XJ(2,2) * XJ(3,3) - XJ(2,3) * XJ(3,2))
     1      - XJ(1,2) * (XJ(2,1) * XJ(3,3) - XJ(2,3) * XJ(3,1))
     2      + XJ(1,3) * (XJ(2,1) * XJ(3,2) - XJ(2,2) * XJ(3,1))
        XJI(1,1) = (XJ(2,2) * XJ(3,3) - XJ(2,3) * XJ(3,2)) / DET
        XJI(1,2) = (XJ(1,3) * XJ(3,2) - XJ(1,2) * XJ(3,3)) / DET
        XJI(1,3) = (XJ(1,2) * XJ(2,3) - XJ(1,3) * XJ(2,2)) / DET
        XJI(2,1) = (XJ(2,3) * XJ(3,1) - XJ(2,1) * XJ(3,3)) / DET
        XJI(2,2) = (XJ(1,1) * XJ(3,3) - XJ(1,3) * XJ(3,1)) / DET
        XJI(2,3) = (XJ(1,3) * XJ(2,1) - XJ(1,1) * XJ(2,3)) / DET
        XJI(3,1) = (XJ(2,1) * XJ(3,2) - XJ(2,2) * XJ(3,1)) / DET
        XJI(3,2) = (XJ(1,2) * XJ(3,1) - XJ(1,1) * XJ(3,2)) / DET
        XJI(3,3) = (XJ(1,1) * XJ(2,2) - XJ(1,2) * XJ(2,1)) / DET
        DO IA = 1, 8
          DO I = 1, 3
            DNX(I, IA) = XJI(I,1) * DNL(1,IA) + XJI(I,2) * DNL(2,IA)
     1                 + XJI(I,3) * DNL(3,IA)
          END DO
        END DO
        PHI = 0.0D0
        C = 0.0D0
        GX = 0.0D0
        GY = 0.0D0
        GZ = 0.0D0
        DO IA = 1, 8
          PHI = PHI + GN(IA) * U(2 * IA - 1)
          C = C + GN(IA) * U(2 * IA)
          GX = GX + DNX(1, IA) * U(2 * IA)
          GY = GY + DNX(2, IA) * U(2 * IA)
          GZ = GZ + DNX(3, IA) * U(2 * IA)
        END DO
        CALL NUT_PUT(JELEM - IOFF, IGP, C, GX, GY, GZ)
        W = DET
        DO IA = 1, 8
          IR = 2 * IA
          RES = D * (DNX(1,IA) * GX + DNX(2,IA) * GY + DNX(3,IA) * GZ)
     1        + G * PHI * GN(IA)
          RHS(IR, 1) = RHS(IR, 1) - RES * W
          DO IB = 1, 8
            AMATRX(IR, 2 * IB) = AMATRX(IR, 2 * IB) + W * D *
     1        (DNX(1,IA) * DNX(1,IB) + DNX(2,IA) * DNX(2,IB)
     2         + DNX(3,IA) * DNX(3,IB))
            AMATRX(IR, 2 * IB - 1) = AMATRX(IR, 2 * IB - 1)
     1        + W * G * GN(IA) * GN(IB)
          END DO
        END DO
      END DO
      END DO
      END DO
      IF (LFLAGS(3) .EQ. 2) THEN
        DO I = 1, NDOFEL
          RHS(I, 1) = 0.0D0
        END DO
      END IF
      RETURN
      END

"""

UMAT = """\
C=======================================================================
C  UMAT: Abaqus entry point (generated by abaqus_composition/make_umat.py)
C=======================================================================
      SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD,
     1 RPL, DDSDDT, DRPLDE, DRPLDT,
     2 STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED,
     3 CMNAME, NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS, COORDS,
     4 DROT, PNEWDT, CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER,
     5 KSPT, JSTEP, KINC)
      USE biofilm_nut_store
      INCLUDE 'ABA_PARAM.INC'
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS), STATEV(NSTATV), DDSDDE(NTENS,NTENS),
     1 DDSDDT(NTENS), DRPLDE(NTENS), STRAN(NTENS), DSTRAN(NTENS),
     2 TIME(2), PREDEF(*), DPRED(*), PROPS(NPROPS), COORDS(3),
     3 DROT(3,3), DFGRD0(3,3), DFGRD1(3,3), JSTEP(4)
      IF (NTENS .NE. 6 .OR. NSTATV .LT. 100 .OR. NPROPS .LT. 46) THEN
        WRITE(7,*) 'BIOFILM UMAT needs NTENS = 6, *DEPVAR >= 100',
     1             ' and 46 constants'
        CALL XIT
      END IF
C     phi: the temperature of a coupled temperature-displacement run
C     (constant 47 = 1, solved by UMATHT below) or field variable 1;
C     c: field variable 2 (not read when phi is the temperature)
      IF (NPROPS .GE. 47 .AND. PROPS(47) .GT. 0.5D0) THEN
        PHI = TEMP + DTEMP
        CN = -1.0D30
      ELSE
        PHI = PREDEF(1) + DPRED(1)
        CN = PREDEF(2) + DPRED(2)
      END IF
C     constant 48 = 1: c from the nutrient UEL overlaid on this element
      IF (NPROPS .GE. 48) THEN
        IF (PROPS(48) .GT. 0.5D0) CN = NUT_GET(NOEL, NPT)
      END IF
      CALL BIOFILM_COMP_UMAT(STRESS, STATEV, DDSDDE, SSE, DTIME,
     1  PHI, CN, NSTATV, PROPS, NPROPS, PNEWDT, DFGRD1,
     2  NOEL, NPT, JSTEP(1), KINC)
C     no thermal coupling terms in the stress (growth enters through
C     alpha, a state variable) and no heat from mechanical work
      DO I = 1, NTENS
        DDSDDT(I) = 0.0D0
        DRPLDE(I) = 0.0D0
      END DO
      RPL = 0.0D0
      DRPLDT = 0.0D0
      RETURN
      END

C=======================================================================
C  UMATHT: phi as the temperature, Klempt et al. 2024 Eq. 34 as the
C  partner's element solves it (its front term is inactive, see
C  FRONT_TERM_FIX.md), with its penalty that keeps phi in [0, 1]:
C      phi_dot = beta lap(phi) + k_alpha alpha_K
C                - P (max(0, phi - 1) + min(0, phi))
C  alpha_K = 1 + alpha (state variable 84, shared with UMAT).
C  Abaqus solves  dU/dt + div(f) = 0  here with f = -beta grad(phi) and
C  U = phi - (integral of the source), so the source sits in U.
C  Constants (*USER MATERIAL, TYPE=THERMAL): beta, k_alpha, P and,
C  optionally, r, k (half-velocity), h (element size): with r > 0 the
C  front term of Eq. 34 as printed is added,
C      - |grad phi| r c/(k+c) n_gradphi . n_gradc = - v . grad(phi),
C      v = r c/(k+c) n_c,  n_c = grad c / |grad c|  (0 where grad c = 0),
C  c and grad c from the nutrient UEL (biofilm_nut_store), and, as the
C  first-order upwind of the finite-difference reproduction, streamline
C  diffusion v h / 2 along n_c (the cell Peclet number is ~10 here).
C  Needs *DENSITY 1.
C=======================================================================
      SUBROUTINE UMATHT(U, DUDT, DUDG, FLUX, DFDT, DFDG,
     1 STATEV, TEMP, DTEMP, DTEMDX, TIME, DTIME, PREDEF, DPRED,
     2 CMNAME, NTGRD, NSTATV, PROPS, NPROPS, COORDS, PNEWDT,
     3 NOEL, NPT, LAYER, KSPT, KSTEP, KINC)
      USE biofilm_nut_store
      INCLUDE 'ABA_PARAM.INC'
      CHARACTER*80 CMNAME
      DIMENSION DUDG(NTGRD), FLUX(NTGRD), DFDT(NTGRD),
     1 DFDG(NTGRD,NTGRD), STATEV(NSTATV), DTEMDX(NTGRD),
     2 TIME(2), PREDEF(1), DPRED(1), PROPS(NPROPS), COORDS(3)
      DIMENSION GC(3), EN(3)
      BETA = PROPS(1)
      AK = PROPS(2)
      PEN = PROPS(3)
      PHI = TEMP + DTEMP
      SRC = AK * (1.0D0 + STATEV(84))
      DSRC = 0.0D0
      IF (PHI .GT. 1.0D0) THEN
        SRC = SRC - PEN * (PHI - 1.0D0)
        DSRC = -PEN
      ELSE IF (PHI .LT. 0.0D0) THEN
        SRC = SRC - PEN * PHI
        DSRC = -PEN
      END IF
      V = 0.0D0
      BS = 0.0D0
      DO I = 1, 3
        EN(I) = 0.0D0
      END DO
      IF (NPROPS .GE. 6) THEN
        IF (PROPS(4) .GT. 0.0D0) THEN
          C = MAX(NUT_GET(NOEL, NPT), 0.0D0)
          CALL NUT_GRAD(NOEL, NPT, GC)
          GM = SQRT(GC(1)**2 + GC(2)**2 + GC(3)**2)
          IF (GM .GT. 1.0D-14) THEN
            V = PROPS(4) * C / (PROPS(5) + C)
            DO I = 1, 3
              EN(I) = GC(I) / GM
            END DO
            BS = 0.5D0 * V * PROPS(6)
          END IF
        END IF
      END IF
      VG = 0.0D0
      DO I = 1, NTGRD
        VG = VG + EN(I) * DTEMDX(I)
      END DO
      SRC = SRC - V * VG
      U = U + DTEMP - SRC * DTIME
      DUDT = 1.0D0 - DSRC * DTIME
      DO I = 1, NTGRD
        DUDG(I) = V * EN(I) * DTIME
        DFDT(I) = 0.0D0
        FLUX(I) = -BETA * DTEMDX(I) - BS * EN(I) * VG
        DO J = 1, NTGRD
          DFDG(I, J) = -BS * EN(I) * EN(J)
        END DO
        DFDG(I, I) = DFDG(I, I) - BETA
      END DO
      RETURN
      END

C=======================================================================
C  BIOFILM_COMP_UMAT: the ANSYS call-site fragments in an Abaqus UMAT
C=======================================================================
      SUBROUTINE BIOFILM_COMP_UMAT(STRESS, ustatev, DDSDDE, SSE,
     1  dTime, AB_PHI, AB_CN, NSTATV, prop, nProp, PNEWDT, DFGRD1,
     2  elemId, kDomIntPt, ldstep, isubst)
      USE biofilm_py_bridge
      USE biofilm_split
      IMPLICIT NONE
      INTEGER NSTATV, nProp, elemId, kDomIntPt, ldstep, isubst
      DOUBLE PRECISION STRESS(6), ustatev(NSTATV), DDSDDE(6,6), SSE
      DOUBLE PRECISION dTime, AB_PHI, AB_CN, prop(nProp)
      DOUBLE PRECISION PNEWDT, DFGRD1(3,3)
C     names the fragments expect from the partner's usermat
      DOUBLE PRECISION Sdp_bio1_n, Sdp_locbio1_n, Sdp_bio2_n
      DOUBLE PRECISION Sdp_locbio2_n, Sdp_sumLocal, Sbio_GrowthConst
      INTEGER keycut
C     Abaqus side: nutrient, mechanics, Voigt map (ANSYS 23/13 -> 13/23)
      DOUBLE PRECISION AB_C, AB_SIG(6), AB_TAN(6,6), AB_W
      DOUBLE PRECISION AB_FV(3,3), AB_FV1(3,3)
      INTEGER AB_KC, AB_I, AB_J, AB_MAP(6)
      DATA AB_MAP /1, 2, 3, 4, 6, 5/
C     trace files: Abaqus runs in a scratch directory that is removed at
C     the end and keeps units below 100 for itself, so the fragments'
C     units 95-99 are moved to 195-199 (make_umat.py) and opened here, once,
C     in the job's output directory
      LOGICAL AB_OPENED
      SAVE AB_OPENED
      DATA AB_OPENED /.FALSE./
      CHARACTER*256 AB_DIR
      INTEGER AB_LD
{decl}
C     --- host values: phi and c as UMAT found them --------------------
      Sdp_bio1_n = AB_PHI
      Sdp_bio2_n = 0.0D0
      Sdp_locbio1_n = 1.0D0
      Sdp_locbio2_n = 1.0D0
      Sdp_sumLocal = 1.0D0
      AB_C = AB_CN
      IF (AB_C .GT. -1.0D29) ustatev(51) = AB_C
      Sbio_GrowthConst = ustatev(84)
      keycut = 0
      IF (.NOT. AB_OPENED) THEN
        CALL GETOUTDIR(AB_DIR, AB_LD)
        OPEN(195, FILE=AB_DIR(1:AB_LD)//'/age_trace.csv',
     1       STATUS='UNKNOWN', POSITION='APPEND')
        OPEN(196, FILE=AB_DIR(1:AB_LD)//'/comp_trace.csv',
     1       STATUS='UNKNOWN', POSITION='APPEND')
        OPEN(197, FILE=AB_DIR(1:AB_LD)//'/phi_trace.csv',
     1       STATUS='UNKNOWN', POSITION='APPEND')
        OPEN(198, FILE=AB_DIR(1:AB_LD)//'/pm_trace.csv',
     1       STATUS='UNKNOWN', POSITION='APPEND')
        OPEN(199, FILE=AB_DIR(1:AB_LD)//'/split_trace.csv',
     1       STATUS='UNKNOWN', POSITION='APPEND')
        AG_TRACE_OPEN = .TRUE.
        CM_TRACE_OPEN = .TRUE.
        TRACE_OPEN = .TRUE.
        PM_TRACE_OPEN = .TRUE.
        SP_TRACE_OPEN = .TRUE.
        AB_OPENED = .TRUE.
      END IF
{exec}
C     --- stress from F_g = (1 + alpha) I, as in the ANSYS runs --------
      DO AB_J = 1, 3
        DO AB_I = 1, 3
          AB_FV(AB_I, AB_J) = ustatev(84 + AB_I + 3 * (AB_J - 1))
        END DO
      END DO
      IF (AB_FV(1,1) .EQ. 0.0D0) THEN
        DO AB_J = 1, 3
          DO AB_I = 1, 3
            AB_FV(AB_I, AB_J) = 0.0D0
          END DO
          AB_FV(AB_J, AB_J) = 1.0D0
        END DO
      END IF
      AB_KC = 0
      CALL BIOFILM_GROWTH_VISCO_V01(DFGRD1, AB_SIG, AB_TAN,
     1  prop(43), prop(44), prop(45), prop(46),
     2  Sdp_bio1_n + Sdp_bio2_n, Sbio_GrowthConst, AB_FV, AB_FV1,
     3  0.0D0, dTime, 0.0D0, 0.0D0, AB_W, AB_KC, elemId)
      DO AB_I = 1, 6
        STRESS(AB_I) = AB_SIG(AB_MAP(AB_I))
        DO AB_J = 1, 6
          DDSDDE(AB_I, AB_J) = AB_TAN(AB_MAP(AB_I), AB_MAP(AB_J))
        END DO
      END DO
      DO AB_J = 1, 3
        DO AB_I = 1, 3
          ustatev(84 + AB_I + 3 * (AB_J - 1)) = AB_FV1(AB_I, AB_J)
        END DO
      END DO
      SSE = AB_W
      IF (keycut .NE. 0 .OR. AB_KC .NE. 0) PNEWDT = 0.5D0
      RETURN
      END
"""


def build() -> str:
    decl = pf.read_lines(CS / "phi_mode_decl.inc")
    exe = pf.set_nut_source(pf.read_lines(CS / "phi_mode_exec.inc"), "AB_C")
    # units 95-99 -> 195-199 (Abaqus keeps units below 100); opened by the wrapper
    exe = [re.sub(r"\b(OPEN|WRITE|FLUSH)\((9[5-9])\b", r"\g<1>(1\2", l) for l in exe]
    assert not any(re.search(r"\b(OPEN|WRITE|FLUSH)\(9[5-9]\b", l) for l in exe)
    out = [f"C  generated by abaqus_composition/make_umat.py from {', '.join(p.name for p in PARTS)},",
           "C  phi_mode_decl.inc, phi_mode_exec.inc -- do not edit, regenerate"]
    for p in PARTS:
        out += [f"C ---- {p.relative_to(ROOT).as_posix()} ----"] + part_lines(p)
    out += NUT.splitlines()
    out += UMAT.format(decl="\n".join(decl), exec="\n".join(exe)).splitlines()
    return "\n".join(out) + "\n"


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    Path(sys.argv[1]).write_text(build(), encoding="latin-1")
    print(f"wrote {sys.argv[1]}")
