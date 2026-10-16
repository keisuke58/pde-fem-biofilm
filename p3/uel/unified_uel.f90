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
! With NSVARS = 8 x 14 + 16 = 128 the element also writes von Mises stress and pressure (Cauchy,
! pressure positive in compression) at each Gauss point to SVARS(113:128) (pairs, Gauss point order).
!
! Equations (weak form, implicit Euler; p3/unified_model.tex, Sections 3-4, beta_s = 0):
!   phase    : int N [eta_phi (phi - phi_n)/dt + 2 K_phi/(1-phi)^3 + lam h'(phi) (Psi_pm - Psi_ref)]
!              + beta grad N . grad phi  dV = 0,   h = phi^2
!   nutrient : int N [(c - c_n)/dt + Y xidot] + D grad N . grad c  dV = 0
!   growth   : xidot = L0 phi sum_i s_i psi_i (Y ln c - mu_b)   (implicit in phi and c)
!   mechanics: finite strain (default), F = F_e F_g, F_g = (1 + xi) I, int P : Grad N dV = 0 with
!              P from W = (1+xi)^3 (phi^2 + f) [mu/2 (J_e^(-2/3) tr C_e - 3) + K/2 (J_e - 1)^2]
!              (nh_stress; Klempt et al. 2024 Eq. 20, volumetric penalty in place of the constraint);
!              or small strain, sigma = (phi^2 + f) C : (eps - xi I) (JPROPS(1) = 0; kin_finite);
!              JPROPS(2) = 1: mean dilatation against volumetric locking (kin_md; bilayer_test.py);
!              JPROPS(3) = 1: nutrient in the current configuration (kin_cur)
!   composition: the point model of Klempt et al. 2026 at every Gauss point, advanced once per
!              increment over s*dt with the nutrient of the start of the increment (staggered),
!              by ecology_native.f (biofilm_py_bridge).
! Residual is assembled as Abaqus wants it: RHS = -R.  AMATRX = dR/dU in closed form (tangent), with
! the composition held at its increment value: the exact Newton tangent of the discrete equations
! solved.  tangent_fd (central differences of R) is kept as its check.
!
! PROPS (NPROPS = 16): E_bio, nu, f_void, beta, eta_phi, K_phi, lam, Psi_ref, L0, Y, mu_b, D,
!                      s (time link), n_sub, phi_min, c_ref (<= 0: constant c*; > 0: local nutrient c/c_ref)
!   optional 17, 18: antibiotic start time and alpha*; 19, 20: homeostatic-pressure law p_ref and form
!   (module header; absent or p_ref <= 0: off)
! The point model's theta, c*, alpha* (antibiotic) and eta come from biofilm_py_bridge
! (eco_native_set, or the file named by BIOFILM_ECO_CASE).

module unified_uel_mod
  use biofilm_py_bridge
  implicit none
  integer, parameter :: NN = 8, NDN = 5, NDOF = NN * NDN, NGP = 8, NSV = 14
  ! kinematics: .true. finite strain (default), .false. small strain; Abaqus: JPROPS(1) = 1 / 0
  logical, save :: kin_finite = .true.
  ! mean dilatation (Q1/P0, Simo, Taylor and Pister 1985) against volumetric locking near
  ! incompressibility (finite strain only): the volumetric energy sees the element mean theta of
  ! J_e = J / alpha^3 instead of J_e at each Gauss point (one pressure per element, condensed, as
  ! Klempt et al. 2024 do with their element pressure); Abaqus: JPROPS(2) = 1
  logical, save :: kin_md = .false.
  ! nutrient in the current configuration (finite strain only): amount J c per reference volume and the
  ! pulled-back flux J D C^-1 Grad c; Abaqus: JPROPS(3) = 1. Tangent by central differences for now.
  logical, save :: kin_cur = .false.
  ! Cauchy stress (11 22 33 12 23 13) at the Gauss points from the last call of residual
  double precision, save :: gp_sig(6, 8) = 0d0
  ! Gauss-point records for the output file (Abaqus cannot print a user element's state variables):
  ! per element and Gauss point x, y, z (reference), phi, c, xi = alpha - 1, living, von Mises, pressure
  integer, parameter :: NREC = 9
  double precision, allocatable, save :: out_rec(:,:,:)
  integer, save :: out_nel = 0
  ! Homeostatic-pressure growth law (keio_wp2/README.ja.md, stages 2B and 4-7; the P3 law of KEIO_PLAN.ja.md;
  ! an assumption of this work, not from a paper): the growth rate xidot is multiplied by
  !   g = min(1, max(0, 1 - p_n / p_h)),
  ! p_n = the Gauss point's Cauchy pressure (positive in compression) at the start of the increment, i.e. the
  ! converged pressure of the previous increment, SVARS(NSV*NGP + 2q) (VARIABLES = 128). The gate is explicit
  ! (staggered, as keio_wp2 and the ANSYS cal5gp build): within an increment g is a constant, so the tangent
  ! keeps its form with xidot, dxidot/dphi and dxidot/dc scaled by g. Because xi never decreases and g is
  ! clipped at 0, the explicit gate cannot oscillate: p overshoots p_h by at most one increment of growth,
  ! dt xidot dp/dxi, a first-order error (homeostatic_test.py shows the order).
  ! PROPS(19) = p_ref (absent or <= 0: off, g = 1 exactly, every result unchanged); PROPS(20) = form:
  !   0  p_h = p_ref                      (keio_wp2 "E": p_ref = P_h E k_alpha T*)
  !   1  p_h = p_ref (phi_n^2 + f)         (keio_wp2 "local stiffness", phi_n of the start of the increment)
  !   2  g = 1 / (1 + max(p_n, 0) / p_ref) (the ANSYS cal5gp build, BIOFILM_HOMEOSTATIC_FACTOR form 2)
  ! gp_gate(q) holds g of the element being processed: set by set_gate (from UEL) or uel_set_gate; 1 otherwise.
  double precision, save :: gp_gate(8) = 1d0
contains

  double precision function growth_gate(pref, iform, pn, phn, f) result(g)
    ! the factor g of the homeostatic-pressure law (see the module header); pref <= 0 gives 1
    double precision, intent(in) :: pref, pn, phn, f
    integer, intent(in) :: iform
    double precision :: ph
    g = 1d0
    if (pref <= 0d0) return
    if (iform == 2) then
      g = 1d0 / (1d0 + max(pn, 0d0) / pref)
      return
    end if
    ph = pref
    if (iform == 1) ph = pref * (min(max(phn, 0d0), 1d0)**2 + f)
    if (ph <= 0d0) return
    g = min(1d0, max(0d0, 1d0 - pn / ph))
  end function growth_gate

  subroutine set_gate(Un, props, nprops, st)
    ! gp_gate from the start-of-increment pressures st(2q) (SVARS(NSV*NGP + 2q)) and phi_n at the Gauss points
    double precision, intent(in) :: Un(NDOF), props(*), st(2*NGP)
    integer, intent(in) :: nprops
    double precision :: N(8), dNr(3,8), gpc(2), phn
    integer :: a, i, j, k, q, iform
    gp_gate = 1d0
    if (nprops < 19) return
    if (props(19) <= 0d0) return
    iform = 0
    if (nprops >= 20) iform = nint(props(20))
    gpc = [-1d0, 1d0] / sqrt(3d0)
    q = 0
    do k = 1, 2
      do j = 1, 2
        do i = 1, 2
          q = q + 1
          call shape(gpc(i), gpc(j), gpc(k), N, dNr)
          phn = 0d0
          do a = 1, NN
            phn = phn + N(a) * Un(NDN*(a-1)+4)
          end do
          gp_gate(q) = growth_gate(props(19), iform, st(2*q), phn, props(3))
        end do
      end do
    end do
  end subroutine set_gate

  subroutine nh_stress(F, sf, ag, mu, Kb, P, AF, dPda)
    ! growth F_g = ag I (ag = alpha = 1 + xi), F_e = F / ag; energy per reference volume
    ! W = ag^3 sf [mu/2 (J_e^(-2/3) tr C_e - 3) + Kb/2 (J_e - 1)^2], sf = phi^2 + f (Klempt et al. 2024
    ! Eq. 20 with a volumetric penalty in place of the incompressibility constraint, as the UMAT core).
    ! P = dW/dF, AF = dP/dF (AF(i,J,k,L)), dPda = dP/d ag.
    double precision, intent(in) :: F(3,3), sf, ag, mu, Kb
    double precision, intent(out) :: P(3,3), AF(3,3,3,3), dPda(3,3)
    double precision :: Fe(3,3), G(3,3), Gt(3,3), J, I1, Jm23, Pe(3,3), Ae(3,3,3,3), d
    integer :: i, jj, k, l
    Fe = F / ag
    call inv3(Fe, Gt, J)                 ! Gt = Fe^-1, J = det Fe
    G = transpose(Gt)                    ! Fe^-T
    I1 = sum(Fe * Fe)
    Jm23 = J**(-2d0/3d0)
    Pe = mu * Jm23 * (Fe - I1 / 3d0 * G) + Kb * (J - 1d0) * J * G
    do l = 1, 3
      do k = 1, 3
        do jj = 1, 3
          do i = 1, 3
            d = 0d0
            if (i == k .and. jj == l) d = 1d0
            Ae(i,jj,k,l) = mu * ( -2d0/3d0 * Jm23 * G(k,l) * (Fe(i,jj) - I1/3d0 * G(i,jj)) &
                                 + Jm23 * (d - 2d0/3d0 * Fe(k,l) * G(i,jj) + I1/3d0 * G(i,l) * G(k,jj)) ) &
                         + Kb * ((2d0*J - 1d0) * J * G(k,l) * G(i,jj) - (J*J - J) * G(i,l) * G(k,jj))
          end do
        end do
      end do
    end do
    P = sf * ag**2 * Pe
    AF = sf * ag * Ae
    do jj = 1, 3
      do i = 1, 3
        dPda(i,jj) = sf * (2d0 * ag * Pe(i,jj) - sum(Ae(i,jj,:,:) * F))
      end do
    end do
  end subroutine nh_stress

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
    double precision :: Fd(3,3), Pk(3,3), AFk(3,3,3,3), dPk(3,3)
    double precision :: Fdi(3,3), Jdet, Fgp(3,3,NGP), dNgp(3,8,NGP), wgp(NGP), agp(NGP), sgp(NGP), Jgp(NGP)
    double precision :: vol, th, pv, Sc(3,3), Cinv(3,3), Fn(3,3), Fni(3,3), Jn
    integer :: a, i, j, k, q, jd
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
    vol = 0d0; th = 0d0; pv = 0d0
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
          xid = L0 * max(ph, 0d0) * liv * max(Y * log(max(c, 1d-12)) - mub, 0d0) * gp_gate(q)
          ! growth: xi_{n+1} = xi_n + dt xidot, with xi_n stored in sv(13,q); gp_gate = 1 unless the
          ! homeostatic-pressure law is on (module header)
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
          if (kin_finite) then
            Fd = 0d0
            do a = 1, NN
              do jd = 1, 3
                Fd(:, jd) = Fd(:, jd) + U(NDN*(a-1)+1:NDN*(a-1)+3) * dN(jd,a)
              end do
            end do
            do jd = 1, 3
              Fd(jd,jd) = Fd(jd,jd) + 1d0
            end do
            if (kin_md) then
              ! deviatoric part here; the volumetric part after the loop with the element mean
              call nh_stress(Fd, Ef, 1d0 + xin, mu, 0d0, Pk, AFk, dPk)
              call inv3(Fd, Fdi, Jdet)
              Fgp(:,:,q) = transpose(Fdi); dNgp(:,:,q) = dN; wgp(q) = w
              agp(q) = 1d0 + xin; sgp(q) = Ef; Jgp(q) = Jdet
              vol = vol + w; th = th + w * Jdet / agp(q)**3
              pv = pv + w * agp(q)**3 * Ef
            else
              call nh_stress(Fd, Ef, 1d0 + xin, mu, lamL + 2d0*mu/3d0, Pk, AFk, dPk)
            end if
          else
            sig = Ef * matmul(Cm, eps)
          end if
          if (kin_finite) then
            call inv3(Fd, Fdi, Jdet)
            Sc = matmul(Pk, transpose(Fd)) / Jdet       ! Cauchy (deviatoric only with mean dilatation)
            if (kin_cur) then
              Cinv = matmul(Fdi, transpose(Fdi))
              Fn = 0d0
              do a = 1, NN
                do jd = 1, 3
                  Fn(:, jd) = Fn(:, jd) + Un(NDN*(a-1)+1:NDN*(a-1)+3) * dN(jd,a)
                end do
              end do
              do jd = 1, 3
                Fn(jd,jd) = Fn(jd,jd) + 1d0
              end do
              call inv3(Fn, Fni, Jn)
            end if
            gp_sig(:, q) = [Sc(1,1), Sc(2,2), Sc(3,3), Sc(1,2), Sc(2,3), Sc(1,3)]
          else
            gp_sig(:, q) = sig
          end if
          do a = 1, NN
            if (kin_finite) then
              R(NDN*(a-1)+1:NDN*(a-1)+3) = R(NDN*(a-1)+1:NDN*(a-1)+3) + w * matmul(Pk, dN(:,a))
            else
            R(NDN*(a-1)+1) = R(NDN*(a-1)+1) + w * (dN(1,a)*sig(1) + dN(2,a)*sig(4) + dN(3,a)*sig(6))
            R(NDN*(a-1)+2) = R(NDN*(a-1)+2) + w * (dN(2,a)*sig(2) + dN(1,a)*sig(4) + dN(3,a)*sig(5))
            R(NDN*(a-1)+3) = R(NDN*(a-1)+3) + w * (dN(3,a)*sig(3) + dN(2,a)*sig(5) + dN(1,a)*sig(6))
            end if
            R(NDN*(a-1)+4) = R(NDN*(a-1)+4) + w * ( N(a) * (etaph*(ph-phn)/dt + 2d0*kphi/(1d0-ph)**3 &
                                                    + lam * 2d0*ph * (psi - psiref)) + beta * dot_product(dN(:,a), gph) )
            if (kin_finite .and. kin_cur) then
              R(NDN*(a-1)+5) = R(NDN*(a-1)+5) + w * ( N(a) * ((Jdet*c - Jn*cn)/dt + Y*xid) &
                                                       + D * Jdet * dot_product(matmul(Cinv, dN(:,a)), gc) )
            else
            R(NDN*(a-1)+5) = R(NDN*(a-1)+5) + w * ( N(a) * ((c-cn)/dt + Y*xid) + D * dot_product(dN(:,a), gc) )
            end if
          end do
        end do
      end do
    end do
    if (kin_finite .and. kin_md) then
      ! U_vol = sum_q w_q alpha_q^3 s_q K/2 (theta - 1)^2, theta = sum_q w_q J_q/alpha_q^3 / V
      th = th / vol
      do q = 1, NGP
        do a = 1, NN
          R(NDN*(a-1)+1:NDN*(a-1)+3) = R(NDN*(a-1)+1:NDN*(a-1)+3) + pv * (lamL + 2d0*mu/3d0) * (th - 1d0) &
               * wgp(q) * Jgp(q) / agp(q)**3 / vol * matmul(Fgp(:,:,q), dNgp(:,a,q))
        end do
        gp_sig(1:3, q) = gp_sig(1:3, q) + pv * (lamL + 2d0*mu/3d0) * (th - 1d0) / vol / agp(q)**3
      end do
    end if
  end subroutine residual

  subroutine vm_p(sg, vm, p)
    ! von Mises stress and pressure (positive in compression) of a Cauchy stress (11 22 33 12 23 13)
    double precision, intent(in) :: sg(6)
    double precision, intent(out) :: vm, p
    p = -(sg(1) + sg(2) + sg(3)) / 3d0
    vm = sqrt(0.5d0 * ((sg(1)-sg(2))**2 + (sg(2)-sg(3))**2 + (sg(3)-sg(1))**2) + 3d0 * (sg(4)**2 + sg(5)**2 + sg(6)**2))
  end subroutine vm_p

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
                     * max(props(10) * log(max(c, 1d-12)) - props(11), 0d0) * gp_gate(q)
        end do
      end do
    end do
  end subroutine store_growth

  subroutine tangent(U, Un, coords, sv, props, dt, K)
    ! closed-form K = dR/dU of residual (composition held fixed); tangent_fd is the check
    double precision, intent(in) :: U(NDOF), Un(NDOF), coords(3, NN), sv(NSV, NGP), props(16), dt
    double precision, intent(out) :: K(NDOF, NDOF)
    double precision :: N(8), dNr(3,8), dN(3,8), Jm(3,3), Ji(3,3), det, w, gpc(2), Bm(6,3,8)
    double precision :: ph, c, eps(6), Cm(6,6), lamL, mu, Ef, dEf, CB(6,3), Cmv(6), sig0(6)
    double precision :: E, nu, f, beta, etaph, kphi, lam, psiref, L0, Y, mub, D, liv, xid, xin, psi
    double precision :: aff, dxp, dxc, kpp, kcc, kcp, m(6), Fd(3,3), Pk(3,3), AFk(3,3,3,3), dPk(3,3)
    integer :: a, b, i, j, l, q, ia, ib, jd, kd
    double precision :: Gq(3,3,NGP), Jq(NGP), aq(NGP), sq(NGP), wq(NGP), Nq(8,NGP), dap(NGP), dac(NGP), dEq(NGP)
    double precision :: gv(3,8,NGP), Bv(3,8), Vv, thv, pvv, Kb, Av, dpv, dth, dB(3), cq, Fi(3,3), Jdet
    if (kin_finite .and. kin_cur) then          ! nutrient in the current configuration: by differences
      call tangent_fd(U, Un, coords, sv, props, dt, K)
      return
    end if
    E = props(1); nu = props(2); f = props(3); beta = props(4); etaph = props(5); kphi = props(6)
    lam = props(7); psiref = props(8); L0 = props(9); Y = props(10); mub = props(11); D = props(12)
    lamL = E*nu/((1+nu)*(1-2*nu)); mu = E/(2*(1+nu))
    Cm = 0d0
    Cm(1:3,1:3) = lamL
    do i = 1, 3
      Cm(i,i) = lamL + 2*mu; Cm(3+i,3+i) = mu
    end do
    m = [1d0, 1d0, 1d0, 0d0, 0d0, 0d0]
    Cmv = matmul(Cm, m)
    gpc = [-1d0, 1d0] / sqrt(3d0)
    K = 0d0
    q = 0
    do l = 1, 2
      do j = 1, 2
        do i = 1, 2
          q = q + 1
          call shape(gpc(i), gpc(j), gpc(l), N, dNr)
          Jm = matmul(dNr, transpose(coords))
          call inv3(Jm, Ji, det)
          dN = matmul(transpose(Ji), dNr)
          w = det
          ph = 0d0; c = 0d0
          do a = 1, NN
            ph = ph + N(a) * U(NDN*(a-1)+4); c = c + N(a) * U(NDN*(a-1)+5)
          end do
          liv = 0d0
          if (sv(14,q) > 0.5d0) liv = sum(sv(1:5,q) * sv(7:11,q))
          psi = 0d0
          if (sv(14,q) > 0.5d0) psi = psi_pm(sv(1:12,q), eco_alstar)
          aff = Y * log(max(c, 1d-12)) - mub
          xid = L0 * max(ph, 0d0) * liv * max(aff, 0d0) * gp_gate(q)
          dxp = 0d0; dxc = 0d0
          if (aff > 0d0 .and. ph > 0d0) dxp = L0 * liv * aff * gp_gate(q)
          if (aff > 0d0 .and. c > 1d-12) dxc = L0 * max(ph, 0d0) * liv * Y / c * gp_gate(q)
          xin = sv(13,q) + dt * xid
          ! B matrices (Voigt: 11 22 33 12 23 13 engineering shear, as in residual)
          Bm = 0d0
          do a = 1, NN
            Bm(1,1,a) = dN(1,a); Bm(2,2,a) = dN(2,a); Bm(3,3,a) = dN(3,a)
            Bm(4,1,a) = dN(2,a); Bm(4,2,a) = dN(1,a)
            Bm(5,2,a) = dN(3,a); Bm(5,3,a) = dN(2,a)
            Bm(6,1,a) = dN(3,a); Bm(6,3,a) = dN(1,a)
          end do
          eps = 0d0
          do a = 1, NN
            eps = eps + matmul(Bm(:,:,a), U(NDN*(a-1)+1:NDN*(a-1)+3))
          end do
          eps(1:3) = eps(1:3) - xin
          Ef = max(ph, 0d0)**2 + f
          dEf = 0d0
          if (ph > 0d0) dEf = 2d0 * ph
          sig0 = matmul(Cm, eps)                       ! sigma / Ef
          kpp = etaph/dt + 6d0*kphi/(1d0-ph)**4 + lam * 2d0 * (psi - psiref)
          kcc = 1d0/dt + Y * dxc
          kcp = Y * dxp
          if (kin_finite) then
            Fd = 0d0
            do a = 1, NN
              do jd = 1, 3
                Fd(:, jd) = Fd(:, jd) + U(NDN*(a-1)+1:NDN*(a-1)+3) * dN(jd,a)
              end do
            end do
            do jd = 1, 3
              Fd(jd,jd) = Fd(jd,jd) + 1d0
            end do
            if (kin_md) then
              ! deviatoric part per Gauss point; the volumetric part with the element mean after the loop
              call nh_stress(Fd, Ef, 1d0 + xin, mu, 0d0, Pk, AFk, dPk)
              call inv3(Fd, Fi, Jdet)
              Gq(:,:,q) = transpose(Fi); Jq(q) = Jdet; aq(q) = 1d0 + xin; sq(q) = Ef; wq(q) = w
              Nq(:,q) = N; dap(q) = dt * dxp; dac(q) = dt * dxc; dEq(q) = dEf
              do a = 1, NN
                gv(:,a,q) = matmul(Gq(:,:,q), dN(:,a))
              end do
            else
              call nh_stress(Fd, Ef, 1d0 + xin, mu, lamL + 2d0*mu/3d0, Pk, AFk, dPk)
            end if
          end if
          do b = 1, NN
            CB = matmul(Cm, Bm(:,:,b))
            do a = 1, NN
              ia = NDN*(a-1); ib = NDN*(b-1)
              ! mechanics: d(R_u,a)/du_b, /dphi_b, /dc_b
              if (kin_finite) then
                do kd = 1, 3
                  do jd = 1, 3
                    K(ia+jd, ib+kd) = K(ia+jd, ib+kd) + w * sum(spread(dN(:,a), 2, 3) * AFk(jd,:,kd,:) &
                                                                * spread(dN(:,b), 1, 3))
                  end do
                end do
                K(ia+1:ia+3, ib+4) = K(ia+1:ia+3, ib+4) + w * N(b) * matmul(dEf / Ef * Pk + dt * dxp * dPk, dN(:,a))
                K(ia+1:ia+3, ib+5) = K(ia+1:ia+3, ib+5) + w * N(b) * dt * dxc * matmul(dPk, dN(:,a))
              else
              K(ia+1:ia+3, ib+1:ib+3) = K(ia+1:ia+3, ib+1:ib+3) + w * Ef * matmul(transpose(Bm(:,:,a)), CB)
              K(ia+1:ia+3, ib+4) = K(ia+1:ia+3, ib+4) + w * N(b) * matmul(transpose(Bm(:,:,a)), &
                                     dEf * sig0 - Ef * dt * dxp * Cmv)
              K(ia+1:ia+3, ib+5) = K(ia+1:ia+3, ib+5) - w * N(b) * Ef * dt * dxc * matmul(transpose(Bm(:,:,a)), Cmv)
              end if
              ! phase
              K(ia+4, ib+4) = K(ia+4, ib+4) + w * (N(a)*N(b)*kpp + beta * dot_product(dN(:,a), dN(:,b)))
              ! nutrient
              K(ia+5, ib+5) = K(ia+5, ib+5) + w * (N(a)*N(b)*kcc + D * dot_product(dN(:,a), dN(:,b)))
              K(ia+5, ib+4) = K(ia+5, ib+4) + w * N(a)*N(b)*kcp
            end do
          end do
        end do
      end do
    end do
    if (kin_finite .and. kin_md) then
      ! volumetric part: R_a += A B_a, A = pv K (theta - 1), B_a = d theta / d u_a
      Kb = lamL + 2d0*mu/3d0
      Vv = sum(wq); thv = sum(wq * Jq / aq**3) / Vv; pvv = sum(wq * aq**3 * sq); Av = pvv * Kb * (thv - 1d0)
      do a = 1, NN
        Bv(:,a) = 0d0
        do q = 1, NGP
          Bv(:,a) = Bv(:,a) + wq(q) * Jq(q) / aq(q)**3 * gv(:,a,q) / Vv
        end do
      end do
      do b = 1, NN
        ib = NDN*(b-1)
        do a = 1, NN
          ia = NDN*(a-1)
          do kd = 1, 3
            do jd = 1, 3
              cq = pvv * Kb * Bv(jd,a) * Bv(kd,b)
              do q = 1, NGP
                cq = cq + Av * wq(q) * Jq(q) / aq(q)**3 / Vv * (gv(jd,a,q)*gv(kd,b,q) - gv(jd,b,q)*gv(kd,a,q))
              end do
              K(ia+jd, ib+kd) = K(ia+jd, ib+kd) + cq
            end do
          end do
          ! phi_b (through s and alpha) and c_b (through alpha)
          dpv = 0d0; dth = 0d0; dB = 0d0
          do q = 1, NGP
            dpv = dpv + wq(q) * Nq(b,q) * (3d0 * aq(q)**2 * sq(q) * dap(q) + aq(q)**3 * dEq(q))
            dth = dth - 3d0 * wq(q) * Nq(b,q) * Jq(q) / aq(q)**4 * dap(q) / Vv
            dB = dB - 3d0 * wq(q) * Nq(b,q) * Jq(q) / aq(q)**4 * dap(q) / Vv * gv(:,a,q)
          end do
          K(ia+1:ia+3, ib+4) = K(ia+1:ia+3, ib+4) + Kb * ((thv - 1d0) * dpv + pvv * dth) * Bv(:,a) + Av * dB
          dpv = 0d0; dth = 0d0; dB = 0d0
          do q = 1, NGP
            dpv = dpv + wq(q) * Nq(b,q) * 3d0 * aq(q)**2 * sq(q) * dac(q)
            dth = dth - 3d0 * wq(q) * Nq(b,q) * Jq(q) / aq(q)**4 * dac(q) / Vv
            dB = dB - 3d0 * wq(q) * Nq(b,q) * Jq(q) / aq(q)**4 * dac(q) / Vv * gv(:,a,q)
          end do
          K(ia+1:ia+3, ib+5) = K(ia+1:ia+3, ib+5) + Kb * ((thv - 1d0) * dpv + pvv * dth) * Bv(:,a) + Av * dB
        end do
      end do
    end if
  end subroutine tangent

  subroutine tangent_fd(U, Un, coords, sv, props, dt, K)
    ! central differences of residual: the check of tangent
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
  end subroutine tangent_fd
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

subroutine uel_set_finite(k, fb) bind(C, name="uel_set_finite")
  use unified_uel_mod
  use iso_c_binding
  integer(c_int), value :: k, fb
  kin_finite = (k /= 0)
  kin_md = (fb /= 0)
end subroutine uel_set_finite

subroutine uel_set_current(k) bind(C, name="uel_set_current")
  use unified_uel_mod
  use iso_c_binding
  integer(c_int), value :: k
  kin_cur = (k /= 0)
end subroutine uel_set_current

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

subroutine uel_set_gate(g) bind(C, name="uel_set_gate")
  ! the growth gate per Gauss point for the C entry points (uel_element, uel_store); 1 = off
  use unified_uel_mod
  use iso_c_binding
  real(c_double), intent(in) :: g(8)
  gp_gate = g
end subroutine uel_set_gate

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
! *USER ELEMENT, NODES=8, TYPE=U1, PROPERTIES=16, I PROPERTIES=2, COORDINATES=3, VARIABLES=128
! 1, 2, 3, 11, 12
! SVARS per element: NSV x NGP = 14 x 8 = 112 (composition 1-12, growth 13, started 14).
! Abaqus hands the element the state variables of the start of the increment at every
! iteration; the composition is therefore advanced from them on every call (the same result each
! iteration) and the end-of-increment values are written back.  The nodal dofs in U are ordered
! node by node as listed above (1, 2, 3, 11, 12), which is the order the module uses.
! The point model's constants are read once from BIOFILM_ECO_CASE (ecology_native.f, through
! UEXTERNALDB or at the first call).  With VARIABLES = 128 the Gauss-point results (phi, c, alpha - 1,
! living biomass, von Mises, pressure) are also written to uel_gp.txt by UEXTERNALDB, every 10th
! increment and at the end of the analysis.
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
  integer :: q
  if (.not. eco_inited) call eco_native_init()
  eco_damp = .true.
  ! kinematics: JPROPS(1) = 0 small strain, otherwise (or no integer property) finite strain
  kin_finite = .true.
  if (NJPROP >= 1) kin_finite = (JPROPS(1) /= 0)
  kin_md = .false.
  if (NJPROP >= 2) kin_md = (JPROPS(2) /= 0)
  kin_cur = .false.
  if (NJPROP >= 3) kin_cur = (JPROPS(3) /= 0)
  ! optional antibiotic: PROPS(17) = total time at which it is switched on, PROPS(18) = alpha*
  if (NPROPS >= 18) then
    if (TIME(2) >= PROPS(17) - 1d-12) then
      eco_alstar = PROPS(18)
    else
      eco_alstar = 0d0
    end if
  end if
  Un = U(1:NDOF) - DU(1:NDOF, 1)
  sv = reshape(SVARS(1:NSV*NGP), [NSV, NGP])
  ! optional homeostatic-pressure law: PROPS(19) = p_ref, PROPS(20) = form (module header); needs the
  ! start-of-increment pressures of VARIABLES = 128
  gp_gate = 1d0
  if (NSVARS >= NSV*NGP + 2*NGP) call set_gate(Un, PROPS, NPROPS, SVARS(NSV*NGP+1:NSV*NGP+2*NGP))
  if (DTIME > 0d0) then
    call advance_composition(U, Un, COORDS, sv, PROPS(1:16), DTIME, ok)
    if (.not. ok) then
      PNEWDT = 0.5d0                       ! ask Abaqus for a smaller increment
      return
    end if
  end if
  call residual(U, Un, COORDS, sv, PROPS(1:16), max(DTIME, 1d-30), R)
  RHS(1:NDOF, 1) = -R
  ! optional output: SVARS(NSV*NGP + 2q - 1 : NSV*NGP + 2q) = von Mises, pressure at Gauss point q
  if (NSVARS >= NSV*NGP + 2*NGP) then
    do q = 1, NGP
      call vm_p(gp_sig(:, q), SVARS(NSV*NGP + 2*q - 1), SVARS(NSV*NGP + 2*q))
    end do
  end if
  if (LFLAGS(3) == 1 .or. LFLAGS(3) == 2) then
    call tangent(U, Un, COORDS, sv, PROPS(1:16), max(DTIME, 1d-30), K)
    AMATRX = K
  end if
  if (DTIME > 0d0) call store_growth(U, COORDS, sv, PROPS(1:16), DTIME)
  ! after store_growth: the record holds the end-of-increment alpha - 1 (it was one increment behind before)
  if (NSVARS >= NSV*NGP + 2*NGP) then
    call keep_record(JELEM, U, COORDS, sv, SVARS(NSV*NGP+1:NSV*NGP+2*NGP), DTIME, PROPS(1:16))
  end if
  SVARS(1:NSV*NGP) = reshape(sv, [NSV*NGP])
end subroutine UEL

subroutine keep_record(jel, U, coords, sv, st, dt, props)
  ! the latest values of an element's Gauss points; the last call of an increment is the converged one
  use unified_uel_mod
  implicit none
  integer, intent(in) :: jel
  double precision, intent(in) :: U(NDOF), coords(3, NN), sv(NSV, NGP), st(2*NGP), dt, props(16)
  double precision, allocatable :: tmp(:,:,:)
  double precision :: N(8), dNr(3,8), gpc(2), ph, c, xin, liv
  integer :: i, j, l, q, a, nnew
  if (jel > out_nel) then
    nnew = max(jel, 2 * out_nel, 64)
    allocate(tmp(NREC, NGP, nnew)); tmp = -1d30
    if (out_nel > 0) tmp(:, :, 1:out_nel) = out_rec
    call move_alloc(tmp, out_rec); out_nel = nnew
  end if
  gpc = [-1d0, 1d0] / sqrt(3d0)
  q = 0
  do l = 1, 2
    do j = 1, 2
      do i = 1, 2
        q = q + 1
        call shape(gpc(i), gpc(j), gpc(l), N, dNr)
        ph = 0d0; c = 0d0
        do a = 1, NN
          ph = ph + N(a) * U(NDN*(a-1)+4); c = c + N(a) * U(NDN*(a-1)+5)
        end do
        liv = 0d0
        if (sv(14,q) > 0.5d0) liv = sum(sv(1:5,q) * sv(7:11,q))
        xin = sv(13,q)
        out_rec(1:3, q, jel) = matmul(coords, N)
        out_rec(4:9, q, jel) = [ph, c, xin, liv, st(2*q-1), st(2*q)]
      end do
    end do
  end do
end subroutine keep_record

subroutine write_records(kstep, kinc, ttot)
  ! appends the Gauss-point records to uel_gp.txt in the job directory
  use unified_uel_mod
  implicit none
  integer, intent(in) :: kstep, kinc
  double precision, intent(in) :: ttot
  integer :: u, e, q, st
  character(len=512) :: path
  if (out_nel == 0) return
  ! the file: BIOFILM_UEL_GP if set (an absolute path), else uel_gp.txt in the solver's directory
  call get_environment_variable('BIOFILM_UEL_GP', path, status=st)
  if (st /= 0) path = 'uel_gp.txt'
  open(newunit=u, file=trim(path), position='append', status='unknown', action='write')
  write(u, '(a,i4,a,i6,a,es16.8)') '# step ', kstep, ' increment ', kinc, ' total time ', ttot
  write(u, '(a)') '# element gp x y z phi c alpha-1 living von_Mises pressure'
  do e = 1, out_nel
    if (out_rec(1, 1, e) < -1d29) cycle
    do q = 1, NGP
      write(u, '(i8,i3,9es16.8)') e, q, out_rec(:, q, e)
    end do
  end do
  close(u)
end subroutine write_records

subroutine UEXTERNALDB(LOP, LRESTART, TIME, DTIME, KSTEP, KINC)
  ! LOP 0: read the point model's constants once, before any element call (one thread);
  ! LOP 2 (end of an increment): every 10th increment, and LOP 3 (end of the analysis), write uel_gp.txt
  use unified_uel_mod
  implicit none
  integer :: LOP, LRESTART, KSTEP, KINC
  double precision :: TIME(2), DTIME
  if (LOP == 0 .and. .not. eco_inited) call eco_native_init()
  if ((LOP == 2 .and. mod(KINC, 10) == 0) .or. LOP == 3) call write_records(KSTEP, KINC, TIME(2))
end subroutine UEXTERNALDB
