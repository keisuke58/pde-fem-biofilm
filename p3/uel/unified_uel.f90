! unified_uel.f90 -- the two-level biofilm model of p3/unified_model.tex as an Abaqus user element
! (UEL), step 3 of P3_STEP3_FEM_PLAN.ja.md.  Free-form Fortran (Abaqus: add '-free' or the
! compiler's free-form flag to compile_fortran in abaqus_v6.env).
!
! Element: 8-node hexahedron, 2x2x2 Gauss points.  Nodal degrees of freedom, in this order per node:
!   1, 2, 3  displacement u            (Abaqus dofs 1-3)
!   4        biofilm phase phi          (Abaqus dof 11)
!   5        nutrient c                 (Abaqus dof 12)
! so the element has 40 degrees of freedom.  State variables per Gauss point (NSV = 14):
!   1-12  composition g = (s_1..5, s_0, psi_1..5, gamma)  (the point model's state)
!   13    growth xi = alpha - 1
!   14    flag: 1 once the point model has started (phi above phi_min)
!
! Equations (weak form, implicit Euler; p3/unified_model.tex, Sections 3-4, beta_s = 0):
!   phase    : int N [eta_phi (phi - phi_n)/dt + 2 K_phi/(1-phi)^3 + lam h'(phi) (Psi_pm - Psi_ref)]
!              + beta grad N . grad phi  dV = 0,   h = phi^2
!   nutrient : int N [(c - c_n)/dt + Y xidot] + D grad N . grad c  dV = 0
!   growth   : xidot = L0 phi sum_i s_i psi_i (Y ln c - mu_b)   (implicit in phi and c)
!   mechanics: small strain, sigma = (phi^2 + f) C : (eps - xi I), int B^T sigma dV = 0
!              (a check of the coupling; the finite-strain law of Klempt et al. 2024 is next)
!   composition: the point model of Klempt et al. 2026 at every Gauss point, advanced once per
!              increment over s*dt with the nutrient of the start of the increment (staggered),
!              by ecology_native.f (biofilm_py_bridge).
! Residual is assembled as Abaqus wants it: RHS = -R.  AMATRX = dR/dU by central differences of R
! with the composition held at its increment value (an exact Newton tangent of the discrete
! equations solved; a closed-form tangent is a later optimisation).
!
! PROPS (NPROPS = 16): E_bio, nu, f_void, beta, eta_phi, K_phi, lam, Psi_ref, L0, Y, mu_b, D,
!                      s (time link), n_sub, phi_min, c_ref (<= 0: constant c*; > 0: local nutrient c/c_ref)
! The point model's theta, c*, alpha* (antibiotic) and eta come from biofilm_py_bridge
! (eco_native_set, or the file named by BIOFILM_ECO_CASE).

module unified_uel_mod
  use biofilm_py_bridge
  implicit none
  integer, parameter :: NN = 8, NDN = 5, NDOF = NN * NDN, NGP = 8, NSV = 14
contains

  subroutine shape(xi, eta, zeta, N, dN)
    double precision, intent(in) :: xi, eta, zeta
    double precision, intent(out) :: N(8), dN(3, 8)
    double precision, parameter :: s(3, 8) = reshape([ -1d0,-1d0,-1d0,  1d0,-1d0,-1d0,  1d0, 1d0,-1d0, -1d0, 1d0,-1d0, &
                                                        -1d0,-1d0, 1d0,  1d0,-1d0, 1d0,  1d0, 1d0, 1d0, -1d0, 1d0, 1d0], [3, 8])
    integer :: a
    do a = 1, 8
      N(a) = 0.125d0 * (1 + s(1,a)*xi) * (1 + s(2,a)*eta) * (1 + s(3,a)*zeta)
      dN(1,a) = 0.125d0 * s(1,a) * (1 + s(2,a)*eta) * (1 + s(3,a)*zeta)
      dN(2,a) = 0.125d0 * s(2,a) * (1 + s(1,a)*xi) * (1 + s(3,a)*zeta)
      dN(3,a) = 0.125d0 * s(3,a) * (1 + s(1,a)*xi) * (1 + s(2,a)*eta)
    end do
  end subroutine shape

  subroutine inv3(J, Ji, det)
    double precision, intent(in) :: J(3,3)
    double precision, intent(out) :: Ji(3,3), det
    det = J(1,1)*(J(2,2)*J(3,3)-J(2,3)*J(3,2)) - J(1,2)*(J(2,1)*J(3,3)-J(2,3)*J(3,1)) &
        + J(1,3)*(J(2,1)*J(3,2)-J(2,2)*J(3,1))
    Ji(1,1) = (J(2,2)*J(3,3)-J(2,3)*J(3,2))/det; Ji(1,2) = (J(1,3)*J(3,2)-J(1,2)*J(3,3))/det
    Ji(1,3) = (J(1,2)*J(2,3)-J(1,3)*J(2,2))/det; Ji(2,1) = (J(2,3)*J(3,1)-J(2,1)*J(3,3))/det
    Ji(2,2) = (J(1,1)*J(3,3)-J(1,3)*J(3,1))/det; Ji(2,3) = (J(1,3)*J(2,1)-J(1,1)*J(2,3))/det
    Ji(3,1) = (J(2,1)*J(3,2)-J(2,2)*J(3,1))/det; Ji(3,2) = (J(1,2)*J(3,1)-J(1,1)*J(3,2))/det
    Ji(3,3) = (J(1,1)*J(2,2)-J(1,2)*J(2,1))/det
  end subroutine inv3

  double precision function psi_pm(g, ast) result(p)
    ! the point model's stored energy (barrier terms included), constants from biofilm_py_bridge
    double precision, intent(in) :: g(12), ast
    double precision :: A(5,5), sb(5)
    integer :: i, j
    call theta_to_a(eco_theta, A)
    p = 0d0
    do i = 1, 5
      p = p + eco_kp1 * (bar(g(i)) + bar(g(6+i)))
      sb(i) = g(i) * g(6+i)
    end do
    p = p + eco_kp1 * bar(g(6))
    do i = 1, 5
      p = p + 0.5d0 * ast * theta_b(i) * g(6+i)**2   ! b_i from theta (4, 5, 9, 10, 16)
    end do
    do i = 1, 5
      do j = 1, 5
        p = p - 0.5d0 * eco_cstar * A(i,j) * sb(i) * sb(j)
      end do
    end do
  contains
    double precision function bar(x)
      double precision, intent(in) :: x
      bar = 1d0 / (x*x * (1d0-x)**2)
    end function bar
  end function psi_pm

  double precision function theta_b(i)
    integer, intent(in) :: i
    integer, parameter :: idx(5) = [4, 5, 9, 10, 16]
    theta_b = eco_theta(idx(i))
  end function theta_b

  subroutine theta_to_a(th, A)
    double precision, intent(in) :: th(20)
    double precision, intent(out) :: A(5,5)
    A = 0d0
    A(1,1)=th(1); A(1,2)=th(2); A(2,1)=th(2); A(2,2)=th(3)
    A(3,3)=th(6); A(3,4)=th(7); A(4,3)=th(7); A(4,4)=th(8)
    A(1,3)=th(11); A(3,1)=th(11); A(1,4)=th(12); A(4,1)=th(12)
    A(2,3)=th(13); A(3,2)=th(13); A(2,4)=th(14); A(4,2)=th(14)
    A(5,5)=th(15); A(1,5)=th(17); A(5,1)=th(17); A(2,5)=th(18); A(5,2)=th(18)
    A(3,5)=th(19); A(5,3)=th(19); A(4,5)=th(20); A(5,4)=th(20)
  end subroutine theta_to_a

  subroutine residual(U, Un, coords, sv, props, dt, R)
    ! element residual R(U) with the Gauss-point composition sv(1:12,:) held fixed (already advanced)
    double precision, intent(in) :: U(NDOF), Un(NDOF), coords(3, NN), sv(NSV, NGP), props(16), dt
    double precision, intent(out) :: R(NDOF)
    double precision :: N(8), dNr(3,8), dN(3,8), Jm(3,3), Ji(3,3), det, w, gpc(2)
    double precision :: ph, phn, c, cn, gph(3), gc(3), eps(6), sig(6), Cm(6,6), lamL, mu, Ef
    double precision :: E, nu, f, beta, etaph, kphi, lam, psiref, L0, Y, mub, D, liv, xid, xin, psi
    integer :: a, i, j, k, q
    E = props(1); nu = props(2); f = props(3); beta = props(4); etaph = props(5); kphi = props(6)
    lam = props(7); psiref = props(8); L0 = props(9); Y = props(10); mub = props(11); D = props(12)
    lamL = E*nu/((1+nu)*(1-2*nu)); mu = E/(2*(1+nu))
    Cm = 0d0
    Cm(1:3,1:3) = lamL
    do i = 1, 3
      Cm(i,i) = lamL + 2*mu; Cm(3+i,3+i) = mu
    end do
    gpc = [-1d0, 1d0] / sqrt(3d0)
    R = 0d0
    q = 0
    do k = 1, 2
      do j = 1, 2
        do i = 1, 2
          q = q + 1
          call shape(gpc(i), gpc(j), gpc(k), N, dNr)
          Jm = matmul(dNr, transpose(coords))
          call inv3(Jm, Ji, det)
          dN = matmul(transpose(Ji), dNr)
          w = det
          ph = 0d0; phn = 0d0; c = 0d0; cn = 0d0; gph = 0d0; gc = 0d0
          do a = 1, NN
            ph = ph + N(a) * U(NDN*(a-1)+4); phn = phn + N(a) * Un(NDN*(a-1)+4)
            c = c + N(a) * U(NDN*(a-1)+5);   cn = cn + N(a) * Un(NDN*(a-1)+5)
            gph = gph + dN(:,a) * U(NDN*(a-1)+4)
            gc = gc + dN(:,a) * U(NDN*(a-1)+5)
          end do
          liv = 0d0
          if (sv(14,q) > 0.5d0) liv = sum(sv(1:5,q) * sv(7:11,q))
          psi = 0d0
          if (sv(14,q) > 0.5d0) psi = psi_pm(sv(1:12,q), eco_alstar)
          xid = L0 * max(ph, 0d0) * liv * max(Y * log(max(c, 1d-12)) - mub, 0d0)
          ! growth: xi_{n+1} = xi_n + dt xidot, with xi_n stored in sv(13,q)
          xin = sv(13,q) + dt * xid
          ! small-strain mechanics with the growth eigenstrain
          eps = 0d0
          do a = 1, NN
            eps(1) = eps(1) + dN(1,a)*U(NDN*(a-1)+1); eps(2) = eps(2) + dN(2,a)*U(NDN*(a-1)+2)
            eps(3) = eps(3) + dN(3,a)*U(NDN*(a-1)+3)
            eps(4) = eps(4) + dN(2,a)*U(NDN*(a-1)+1) + dN(1,a)*U(NDN*(a-1)+2)
            eps(5) = eps(5) + dN(3,a)*U(NDN*(a-1)+2) + dN(2,a)*U(NDN*(a-1)+3)
            eps(6) = eps(6) + dN(3,a)*U(NDN*(a-1)+1) + dN(1,a)*U(NDN*(a-1)+3)
          end do
          eps(1:3) = eps(1:3) - xin
          Ef = max(ph, 0d0)**2 + f
          sig = Ef * matmul(Cm, eps)
          do a = 1, NN
            R(NDN*(a-1)+1) = R(NDN*(a-1)+1) + w * (dN(1,a)*sig(1) + dN(2,a)*sig(4) + dN(3,a)*sig(6))
            R(NDN*(a-1)+2) = R(NDN*(a-1)+2) + w * (dN(2,a)*sig(2) + dN(1,a)*sig(4) + dN(3,a)*sig(5))
            R(NDN*(a-1)+3) = R(NDN*(a-1)+3) + w * (dN(3,a)*sig(3) + dN(2,a)*sig(5) + dN(1,a)*sig(6))
            R(NDN*(a-1)+4) = R(NDN*(a-1)+4) + w * ( N(a) * (etaph*(ph-phn)/dt + 2d0*kphi/(1d0-ph)**3 &
                                                    + lam * 2d0*ph * (psi - psiref)) + beta * dot_product(dN(:,a), gph) )
            R(NDN*(a-1)+5) = R(NDN*(a-1)+5) + w * ( N(a) * ((c-cn)/dt + Y*xid) + D * dot_product(dN(:,a), gc) )
          end do
        end do
      end do
    end do
  end subroutine residual

  subroutine advance_composition(U, Un, coords, sv, props, dt, ok)
    ! once per increment: start the point model where phi has reached phi_min, advance it over
    ! s*dt with the nutrient of the start of the increment, then store xi_n (done after convergence)
    double precision, intent(in) :: U(NDOF), Un(NDOF), coords(3, NN), props(16), dt
    double precision, intent(inout) :: sv(NSV, NGP)
    logical, intent(out) :: ok
    double precision :: N(8), dNr(3,8), gpc(2), phn, cn, gnew(12), pint, crel
    logical :: have
    integer :: a, i, j, k, q, nsub
    gpc = [-1d0, 1d0] / sqrt(3d0)
    nsub = max(1, nint(props(14)))
    ok = .true.
    q = 0
    do k = 1, 2
      do j = 1, 2
        do i = 1, 2
          q = q + 1
          call shape(gpc(i), gpc(j), gpc(k), N, dNr)
          phn = 0d0; cn = 0d0
          do a = 1, NN
            phn = phn + N(a) * Un(NDN*(a-1)+4); cn = cn + N(a) * Un(NDN*(a-1)+5)
          end do
          if (sv(14,q) < 0.5d0 .and. phn > props(15)) then
            call biofilm_ecology_init_state(5, sv(1:12,q), have)
            sv(14,q) = 1d0
          end if
          if (sv(14,q) > 0.5d0) then
            crel = -1d0
            if (props(16) > 0d0) crel = max(cn, 0d0) / props(16)
            call eco_substeps(sv(1:12,q), eco_theta, props(13)*dt, nsub, crel, gnew, pint, ok)
            if (.not. ok) return
            sv(1:12,q) = gnew
          end if
        end do
      end do
    end do
  end subroutine advance_composition

  subroutine store_growth(U, coords, sv, props, dt)
    ! after convergence: xi_n <- xi_{n+1} at every Gauss point
    double precision, intent(in) :: U(NDOF), coords(3, NN), props(16), dt
    double precision, intent(inout) :: sv(NSV, NGP)
    double precision :: N(8), dNr(3,8), gpc(2), ph, c, liv
    integer :: a, i, j, k, q
    gpc = [-1d0, 1d0] / sqrt(3d0)
    q = 0
    do k = 1, 2
      do j = 1, 2
        do i = 1, 2
          q = q + 1
          call shape(gpc(i), gpc(j), gpc(k), N, dNr)
          ph = 0d0; c = 0d0
          do a = 1, NN
            ph = ph + N(a) * U(NDN*(a-1)+4); c = c + N(a) * U(NDN*(a-1)+5)
          end do
          liv = 0d0
          if (sv(14,q) > 0.5d0) liv = sum(sv(1:5,q) * sv(7:11,q))
          sv(13,q) = sv(13,q) + dt * props(9) * max(ph, 0d0) * liv &
                     * max(props(10) * log(max(c, 1d-12)) - props(11), 0d0)
        end do
      end do
    end do
  end subroutine store_growth

  subroutine tangent(U, Un, coords, sv, props, dt, K)
    double precision, intent(in) :: U(NDOF), Un(NDOF), coords(3, NN), sv(NSV, NGP), props(16), dt
    double precision, intent(out) :: K(NDOF, NDOF)
    double precision :: Up(NDOF), Um(NDOF), Rp(NDOF), Rm(NDOF), hstep
    integer :: m
    do m = 1, NDOF
      hstep = 1d-7 * max(1d0, abs(U(m)))
      Up = U; Um = U
      Up(m) = U(m) + hstep; Um(m) = U(m) - hstep
      call residual(Up, Un, coords, sv, props, dt, Rp)
      call residual(Um, Un, coords, sv, props, dt, Rm)
      K(:, m) = (Rp - Rm) / (2d0 * hstep)
    end do
  end subroutine tangent
end module unified_uel_mod


! ---- C-callable entry points for the Python test driver (p3/uel/driver.py) ----------------------
subroutine uel_set_point_model(theta, cstar, alstar, nwmax, nwtol) bind(C, name="uel_set_point_model")
  use biofilm_py_bridge
  use iso_c_binding
  real(c_double), intent(in) :: theta(20), cstar, alstar, nwtol
  integer(c_int), value :: nwmax
  double precision :: eta(5)
  eta = 1d0
  call eco_native_set(5, cstar, alstar, eta, theta, .false.)
  eco_nwmax = nwmax
  eco_nwtol = nwtol
  eco_damp = .true.
end subroutine uel_set_point_model

subroutine uel_set_phi_init(v) bind(C, name="uel_set_phi_init")
  use biofilm_py_bridge
  use iso_c_binding
  real(c_double), intent(in) :: v(5)
  eco_phi_init = v
  eco_haveinit = .true.
end subroutine uel_set_phi_init

subroutine uel_set_alstar(x) bind(C, name="uel_set_alstar")
  use biofilm_py_bridge
  use iso_c_binding
  real(c_double), value :: x
  eco_alstar = x
end subroutine uel_set_alstar

subroutine uel_element(U, Un, coords, sv, props, dt, R, K) bind(C, name="uel_element")
  use unified_uel_mod
  use iso_c_binding
  real(c_double), intent(in) :: U(NDOF), Un(NDOF), coords(3, NN), sv(NSV, NGP), props(16)
  real(c_double), value :: dt
  real(c_double), intent(out) :: R(NDOF), K(NDOF, NDOF)
  call residual(U, Un, coords, sv, props, dt, R)
  call tangent(U, Un, coords, sv, props, dt, K)
end subroutine uel_element

subroutine uel_advance(U, Un, coords, sv, props, dt, ok) bind(C, name="uel_advance")
  use unified_uel_mod
  use iso_c_binding
  real(c_double), intent(in) :: U(NDOF), Un(NDOF), coords(3, NN), props(16)
  real(c_double), intent(inout) :: sv(NSV, NGP)
  real(c_double), value :: dt
  integer(c_int), intent(out) :: ok
  logical :: lok
  call advance_composition(U, Un, coords, sv, props, dt, lok)
  ok = 0
  if (lok) ok = 1
end subroutine uel_advance

subroutine uel_store(U, coords, sv, props, dt) bind(C, name="uel_store")
  use unified_uel_mod
  use iso_c_binding
  real(c_double), intent(in) :: U(NDOF), coords(3, NN), props(16)
  real(c_double), intent(inout) :: sv(NSV, NGP)
  real(c_double), value :: dt
  call store_growth(U, coords, sv, props, dt)
end subroutine uel_store


! ---- Abaqus entry point ---------------------------------------------------------------------------
! *USER ELEMENT, NODES=8, TYPE=U1, PROPERTIES=16, COORDINATES=3, VARIABLES=112
! 1, 2, 3, 11, 12
! SVARS per element: NSV x NGP = 14 x 8 = 112 (composition 1-12, growth 13, started 14).
! Abaqus hands the element the state variables of the start of the increment at every
! iteration; the composition is therefore advanced from them on every call (the same result each
! iteration) and the end-of-increment values are written back.  The nodal dofs in U are ordered
! node by node as listed above (1, 2, 3, 11, 12), which is the order the module uses.
! The point model's constants are read once from BIOFILM_ECO_CASE (ecology_native.f, through
! UEXTERNALDB or at the first call).
subroutine UEL(RHS, AMATRX, SVARS, ENERGY, NDOFEL, NRHS, NSVARS, PROPS, NPROPS, COORDS, MCRD, NNODE, &
               U, DU, V, A, JTYPE, TIME, DTIME, KSTEP, KINC, JELEM, PARAMS, NDLOAD, JDLTYP, ADLMAG, &
               PREDEF, NPREDF, LFLAGS, MLVARX, DDLMAG, MDLOAD, PNEWDT, JPROPS, NJPROP, PERIOD)
  use unified_uel_mod
  implicit none
  integer :: NDOFEL, NRHS, NSVARS, NPROPS, MCRD, NNODE, JTYPE, KSTEP, KINC, JELEM, NDLOAD, NPREDF, &
             MLVARX, MDLOAD, NJPROP
  integer :: JDLTYP(MDLOAD, *), LFLAGS(*), JPROPS(*)
  double precision :: RHS(MLVARX, *), AMATRX(NDOFEL, NDOFEL), SVARS(*), ENERGY(8), PROPS(*), &
                      COORDS(MCRD, NNODE), U(NDOFEL), DU(MLVARX, *), V(NDOFEL), A(NDOFEL), TIME(2), &
                      DTIME, PARAMS(*), ADLMAG(MDLOAD, *), PREDEF(2, NPREDF, NNODE), DDLMAG(MDLOAD, *), &
                      PNEWDT, PERIOD
  double precision :: Un(NDOF), sv(NSV, NGP), R(NDOF), K(NDOF, NDOF)
  logical :: ok
  if (.not. eco_inited) call eco_native_init()
  eco_damp = .true.
  Un = U(1:NDOF) - DU(1:NDOF, 1)
  sv = reshape(SVARS(1:NSV*NGP), [NSV, NGP])
  if (DTIME > 0d0) then
    call advance_composition(U, Un, COORDS, sv, PROPS(1:16), DTIME, ok)
    if (.not. ok) then
      PNEWDT = 0.5d0                       ! ask Abaqus for a smaller increment
      return
    end if
  end if
  call residual(U, Un, COORDS, sv, PROPS(1:16), max(DTIME, 1d-30), R)
  RHS(1:NDOF, 1) = -R
  if (LFLAGS(3) == 1 .or. LFLAGS(3) == 2) then
    call tangent(U, Un, COORDS, sv, PROPS(1:16), max(DTIME, 1d-30), K)
    AMATRX = K
  end if
  if (DTIME > 0d0) call store_growth(U, COORDS, sv, PROPS(1:16), DTIME)
  SVARS(1:NSV*NGP) = reshape(sv, [NSV*NGP])
end subroutine UEL
