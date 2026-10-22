! living_uel.f90 -- the coupling agreed with Oliver (MEETING_2026-10-20.md) as an Abaqus user element (A-46, 22 Oct 2026):
! Klempt et al. 2024 Eqs. 34-36 with the front term, the point model of Klempt et al. 2026 at the Gauss points
! (ecology_native.f, used as published), consumption by the living biomass g phi_3D (sum_i phibar_i) c, and the
! antibiotic field.  The same equations as ansys_usermat/apdl/coupled_living_consumption.py (finite volumes) and the
! ANSYS element (NEM); here standard Galerkin finite elements: the Laplacian is symmetric and conserves phi, and the
! zero-flux walls need no treatment.  Free-form source; make_abaqus_source.py writes the fixed-form file for Abaqus.
!
! Element: 8-node hexahedron, 2x2x2 Gauss points.  Nodal dofs, in this order per node (Abaqus dofs 1, 2, 3, 11, 12, 13):
!   u1, u2, u3, phi (amount of biofilm phi_3D), c (nutrient), a (antibiotic)          -> 48 dofs
! State variables per Gauss point (NSV = 16; VARIABLES = 128):
!   1-12 point-model state g = (phi_1..5, phi_0, psi_1..5, gamma); 13 alpha - 1; 14 started (1) or not (0);
!   15, 16 c and a at the start of the previous increment (the point model is advanced from them to the present ones)
!
! Equations (implicit Euler over an increment dt, Newton by Abaqus):
!   phi : int N [(phi - phi_n)/dt - k_alpha alpha + P(phi) + v . grad phi] + beta grad N . grad phi dV = 0
!         v = r c_n/(K + c_n) grad c_n / sqrt(|grad c_n|^2 + eps)   (the front term of Eq. 34, nutrient of the start
!         of the increment), beta = beta_0 + r h / 4 (the stabilisation of the ANSYS runs, --artdiff), P the penalty
!   alpha: alpha = alpha_n + dt k_alpha phi   (Eq. 36, at the Gauss point; alpha - 1 is stored)
!   c   : int grad N . grad c + N q_c phi_c L c dV = 0   (Eq. 35 quasi-static, d of Klempt 2024 Table 2; q_c = g l^2/d),
!         phi_c = min(max(phi, 0), 1), L = sum_i phi_i psi_i of the point model (start of the increment)
!   a   : the same with q_a
!   u   : int Bbar^T sigma dV = 0, sigma = E_bio (phi_c^2 + f) C(nu) : (eps_bar - (alpha - 1) I), small strain,
!         B-bar (as wp3_stress_reference.py)
!   point model: at the start of an increment, each started Gauss point is advanced over s dt with c and a linearly
!         interpolated from the previous increment's start values to the present ones (biofilm_ecology_hook_lin, as
!         ANSYS mode 11 and the Python reference); a point starts when phi_n >= phi_min, from the start state of its
!         position (separate: species B where |x| and |y| > x_split, A elsewhere; mixed: the mixed start).
! PROPS (NPROPS = 18): E_bio, nu, f_void, beta, k_alpha, penalty, r, K, q_c, q_a, phi_min, s, n_sub, layout (1 separate,
!   2 mixed), x_split, total start amount, mixed start of species 1 and 2

module living_uel_mod
  use biofilm_py_bridge
  implicit none
  integer, parameter :: NN = 8, NDN = 6, NDOF = NN * NDN, NGP = 8, NSV = 16
  integer, parameter :: NREC = 12        ! x y z phi c a alpha-1 living vM p shareA started
  double precision, allocatable, save :: out_rec(:,:,:)
  integer, save :: out_nel = 0
  double precision, parameter :: EPSG = 1d-6
  ! A-50, manufactured solution (mms.py, mms_src.f90): PROPS(19) = 1 switches the point model off (L = 1) and adds the
  ! sources S_phi, S_c at the end of the increment (step time mms_t)
  logical, save :: mms_on = .false.
  double precision, save :: mms_t = 0d0
contains

  subroutine shape(r, s, t, N, dNr)
    double precision, intent(in) :: r, s, t
    double precision, intent(out) :: N(8), dNr(3, 8)
    double precision, parameter :: xr(8) = [-1d0, 1d0, 1d0, -1d0, -1d0, 1d0, 1d0, -1d0]
    double precision, parameter :: xs(8) = [-1d0, -1d0, 1d0, 1d0, -1d0, -1d0, 1d0, 1d0]
    double precision, parameter :: xt(8) = [-1d0, -1d0, -1d0, -1d0, 1d0, 1d0, 1d0, 1d0]
    integer :: a
    do a = 1, 8
      N(a) = 0.125d0 * (1 + xr(a) * r) * (1 + xs(a) * s) * (1 + xt(a) * t)
      dNr(1, a) = 0.125d0 * xr(a) * (1 + xs(a) * s) * (1 + xt(a) * t)
      dNr(2, a) = 0.125d0 * xs(a) * (1 + xr(a) * r) * (1 + xt(a) * t)
      dNr(3, a) = 0.125d0 * xt(a) * (1 + xr(a) * r) * (1 + xs(a) * s)
    end do
  end subroutine shape

  subroutine gauss(coords, q, N, dN, wdj, x)
    ! Gauss point q (1..8): shape functions, their spatial gradients, weight x det J, position
    double precision, intent(in) :: coords(3, NN)
    integer, intent(in) :: q
    double precision, intent(out) :: N(8), dN(3, 8), wdj, x(3)
    double precision :: dNr(3, 8), J(3, 3), Ji(3, 3), det, g
    integer :: i, k, l
    g = 1d0 / sqrt(3d0)
    i = mod(q - 1, 2); k = mod((q - 1) / 2, 2); l = (q - 1) / 4
    call shape(merge(g, -g, i == 1), merge(g, -g, k == 1), merge(g, -g, l == 1), N, dNr)
    J = matmul(dNr, transpose(coords))                  ! J(i, j) = d x_j / d r_i
    det = J(1,1) * (J(2,2) * J(3,3) - J(2,3) * J(3,2)) - J(1,2) * (J(2,1) * J(3,3) - J(2,3) * J(3,1)) &
        + J(1,3) * (J(2,1) * J(3,2) - J(2,2) * J(3,1))
    Ji(1,1) = J(2,2) * J(3,3) - J(2,3) * J(3,2); Ji(1,2) = J(1,3) * J(3,2) - J(1,2) * J(3,3)
    Ji(1,3) = J(1,2) * J(2,3) - J(1,3) * J(2,2); Ji(2,1) = J(2,3) * J(3,1) - J(2,1) * J(3,3)
    Ji(2,2) = J(1,1) * J(3,3) - J(1,3) * J(3,1); Ji(2,3) = J(1,3) * J(2,1) - J(1,1) * J(2,3)
    Ji(3,1) = J(2,1) * J(3,2) - J(2,2) * J(3,1); Ji(3,2) = J(1,2) * J(3,1) - J(1,1) * J(3,2)
    Ji(3,3) = J(1,1) * J(2,2) - J(1,2) * J(2,1)
    Ji = Ji / det
    dN = matmul(Ji, dNr)                                ! dN(j, a) = d N_a / d x_j
    wdj = det
    x = matmul(coords, N)
  end subroutine gauss

  subroutine start_state(props, x, g)
    ! the start state of a newly colonised Gauss point (coupled_living_consumption.py start_state / group_state)
    double precision, intent(in) :: props(*), x(3)
    double precision, intent(out) :: g(12)
    double precision :: tot
    integer :: n
    n = eco_nact
    tot = props(16)
    g = 0d0
    if (nint(props(14)) == 1) then
      g(1:n) = 1d-6
      if (abs(x(1)) > props(15) .and. abs(x(2)) > props(15)) then
        g(2) = tot - 1d-6 * (n - 1)
      else
        g(1) = tot - 1d-6 * (n - 1)
      end if
    else
      g(1) = props(17); g(2) = props(18)
    end if
    g(6) = 1d0 - sum(g(1:5))
    g(7:11) = 0.999d0
    g(12) = 0d0
  end subroutine start_state

  double precision function living(g)
    double precision, intent(in) :: g(12)
    living = sum(g(1:eco_nact) * g(7:6 + eco_nact))
  end function living

  subroutine advance(Un, coords, sv, props, dt, ok)
    ! point model at the start of the increment: advance the started points, start the new ones
    double precision, intent(in) :: Un(NDOF), coords(3, NN), props(*), dt
    double precision, intent(inout) :: sv(NSV, NGP)
    logical, intent(out) :: ok
    double precision :: N(8), dN(3, 8), w, x(3), ph, c, a, gnew(12), pint
    integer :: q, k
    ok = .true.
    do q = 1, NGP
      call gauss(coords, q, N, dN, w, x)
      ph = 0d0; c = 0d0; a = 0d0
      do k = 1, NN
        ph = ph + N(k) * Un(NDN*(k-1)+4); c = c + N(k) * Un(NDN*(k-1)+5); a = a + N(k) * Un(NDN*(k-1)+6)
      end do
      if (mms_on) then
        sv(14, q) = 0d0
      else if (sv(14, q) > 0.5d0 .and. dt > 1d-6) then
        call biofilm_ecology_hook_lin(sv(1:12, q), eco_theta, props(12) * dt, nint(props(13)), &
                                      sv(15, q), c, sv(16, q), a, gnew, pint, ok)
        if (.not. ok) return
        sv(1:12, q) = gnew
      else if (sv(14, q) < 0.5d0 .and. ph >= props(11)) then
        call start_state(props, x, sv(1:12, q))
        sv(14, q) = 1d0
      end if
      sv(15, q) = c; sv(16, q) = a
    end do
  end subroutine advance

  subroutine element(U, Un, coords, sv, props, dt, R, K, want_k, sig)
    ! residual R(U) and tangent K = dR/dU with the point-model state sv held (already advanced)
    double precision, intent(in) :: U(NDOF), Un(NDOF), coords(3, NN), props(*), dt
    double precision, intent(inout) :: sv(NSV, NGP)
    double precision, intent(out) :: R(NDOF), K(NDOF, NDOF), sig(6, NGP)
    logical, intent(in) :: want_k
    double precision :: N(8), dN(3, 8), w, x(3), ph, phn, c, cn, a, gph(3), gc(3), gcn(3), ga(3), v(3)
    double precision :: al, pen, dpen, phc, dphc, L, E, dE, D0(6, 6), lam, mu, eps(6), sm(6), B(6, 24), Bv(3, 8)
    double precision :: vol, m(6), dsig(6), nrm, sp, dt_, sphi, sc, phim, cm
    double precision :: Ebio, nu, f, beta, ka, P, rr, Km, qc, qa
    integer :: q, i, j, ii, jj, d
    Ebio = props(1); nu = props(2); f = props(3); beta = props(4); ka = props(5); P = props(6)
    rr = props(7); Km = props(8); qc = props(9); qa = props(10)
    dt_ = max(dt, 1d-30)
    R = 0d0; K = 0d0; sig = 0d0
    m = [1d0, 1d0, 1d0, 0d0, 0d0, 0d0]
    lam = nu / ((1 + nu) * (1 - 2 * nu)); mu = 0.5d0 / (1 + nu)
    D0 = 0d0
    D0(1:3, 1:3) = lam
    do i = 1, 3
      D0(i, i) = lam + 2 * mu; D0(i + 3, i + 3) = mu
    end do
    ! B-bar: element mean of the shape-function gradients for the volumetric part
    Bv = 0d0; vol = 0d0
    do q = 1, NGP
      call gauss(coords, q, N, dN, w, x)
      Bv = Bv + w * dN; vol = vol + w
    end do
    Bv = Bv / vol
    do q = 1, NGP
      call gauss(coords, q, N, dN, w, x)
      ph = 0d0; phn = 0d0; c = 0d0; cn = 0d0; a = 0d0; gph = 0d0; gc = 0d0; gcn = 0d0; ga = 0d0
      do i = 1, NN
        ph = ph + N(i) * U(NDN*(i-1)+4); phn = phn + N(i) * Un(NDN*(i-1)+4)
        c = c + N(i) * U(NDN*(i-1)+5); cn = cn + N(i) * Un(NDN*(i-1)+5); a = a + N(i) * U(NDN*(i-1)+6)
        gph = gph + dN(:, i) * U(NDN*(i-1)+4); gc = gc + dN(:, i) * U(NDN*(i-1)+5)
        gcn = gcn + dN(:, i) * Un(NDN*(i-1)+5); ga = ga + dN(:, i) * U(NDN*(i-1)+6)
      end do
      nrm = sqrt(sum(gcn**2) + EPSG)
      sp = rr * cn / (Km + cn)
      v = sp * gcn / nrm
      al = 1d0 + sv(13, q) + dt * ka * ph               ! alpha at the end of the increment
      pen = P * (max(ph - 1d0, 0d0) + min(ph, 0d0))
      dpen = 0d0
      if (ph > 1d0 .or. ph < 0d0) dpen = P
      phc = min(max(ph, 0d0), 1d0)
      dphc = 0d0
      if (ph > 0d0 .and. ph < 1d0) dphc = 1d0
      L = 0d0
      if (sv(14, q) > 0.5d0) L = living(sv(1:12, q))
      sphi = 0d0; sc = 0d0
      if (mms_on) then
        L = 1d0
        call mms_eval(x(1), x(2), x(3), mms_t, beta, rr, Km, qc, sphi, sc, phim, cm)
      end if
      ! phase
      do i = 1, NN
        ii = NDN*(i-1)+4
        R(ii) = R(ii) + w * (N(i) * ((ph - phn) / dt_ - ka * al + pen + dot_product(v, gph) - sphi) &
                             + beta * dot_product(dN(:, i), gph))
      end do
      ! nutrient and antibiotic
      do i = 1, NN
        ii = NDN*(i-1)+5
        R(ii) = R(ii) + w * (dot_product(dN(:, i), gc) + N(i) * qc * phc * L * c - N(i) * sc)
        ii = NDN*(i-1)+6
        R(ii) = R(ii) + w * (dot_product(dN(:, i), ga) + N(i) * qa * phc * L * a)
      end do
      ! mechanics
      B = 0d0
      do i = 1, NN
        do d = 1, 3
          j = 3*(i-1)+d
          B(d, j) = dN(d, i)
          B(1:3, j) = B(1:3, j) + (Bv(d, i) - dN(d, i)) / 3d0
        end do
        B(4, 3*(i-1)+1) = dN(2, i); B(4, 3*(i-1)+2) = dN(1, i)
        B(5, 3*(i-1)+2) = dN(3, i); B(5, 3*(i-1)+3) = dN(2, i)
        B(6, 3*(i-1)+1) = dN(3, i); B(6, 3*(i-1)+3) = dN(1, i)
      end do
      eps = 0d0
      do i = 1, NN
        do d = 1, 3
          eps = eps + B(:, 3*(i-1)+d) * U(NDN*(i-1)+d)
        end do
      end do
      E = Ebio * (phc**2 + f)
      dE = Ebio * 2d0 * phc * dphc
      sm = matmul(D0, eps - (al - 1d0) * m)
      sig(:, q) = E * sm
      do i = 1, NN
        do d = 1, 3
          ii = NDN*(i-1)+d
          R(ii) = R(ii) + w * dot_product(B(:, 3*(i-1)+d), sig(:, q))
        end do
      end do
      if (.not. want_k) cycle
      dsig = dE * sm - E * dt * ka * matmul(D0, m)       ! d sigma / d phi
      do i = 1, NN
        do j = 1, NN
          ! phase rows
          ii = NDN*(i-1)+4; jj = NDN*(j-1)+4
          K(ii, jj) = K(ii, jj) + w * (N(i) * N(j) * (1d0 / dt_ - ka * ka * dt + dpen) &
                                       + N(i) * dot_product(v, dN(:, j)) + beta * dot_product(dN(:, i), dN(:, j)))
          ! nutrient rows
          ii = NDN*(i-1)+5
          K(ii, NDN*(j-1)+5) = K(ii, NDN*(j-1)+5) + w * (dot_product(dN(:, i), dN(:, j)) + N(i) * N(j) * qc * phc * L)
          K(ii, jj) = K(ii, jj) + w * N(i) * N(j) * qc * dphc * L * c
          ii = NDN*(i-1)+6
          K(ii, NDN*(j-1)+6) = K(ii, NDN*(j-1)+6) + w * (dot_product(dN(:, i), dN(:, j)) + N(i) * N(j) * qa * phc * L)
          K(ii, jj) = K(ii, jj) + w * N(i) * N(j) * qa * dphc * L * a
          ! mechanics rows
          do d = 1, 3
            ii = NDN*(i-1)+d
            K(ii, jj) = K(ii, jj) + w * dot_product(B(:, 3*(i-1)+d), dsig) * N(j)
            do jj = 1, 3
              K(ii, NDN*(j-1)+jj) = K(ii, NDN*(j-1)+jj) &
                + w * E * dot_product(B(:, 3*(i-1)+d), matmul(D0, B(:, 3*(j-1)+jj)))
            end do
            jj = NDN*(j-1)+4
          end do
        end do
      end do
    end do
    ! alpha at the end of the increment (stored by the caller after the last, converged call)
  end subroutine element

  subroutine alpha_store(U, coords, sv, props, dt)
    double precision, intent(in) :: U(NDOF), coords(3, NN), props(*), dt
    double precision, intent(inout) :: sv(NSV, NGP)
    double precision :: N(8), dN(3, 8), w, x(3), ph
    integer :: q, i
    do q = 1, NGP
      call gauss(coords, q, N, dN, w, x)
      ph = 0d0
      do i = 1, NN
        ph = ph + N(i) * U(NDN*(i-1)+4)
      end do
      sv(13, q) = sv(13, q) + dt * props(5) * ph
    end do
  end subroutine alpha_store

  subroutine keep_record(jel, U, coords, sv, sig)
    integer, intent(in) :: jel
    double precision, intent(in) :: U(NDOF), coords(3, NN), sv(NSV, NGP), sig(6, NGP)
    double precision, allocatable :: tmp(:,:,:)
    double precision :: N(8), dN(3, 8), w, x(3), ph, c, a, p, s(6), vm, sh
    integer :: q, i, nnew
    if (jel > out_nel) then
      nnew = max(jel, 2 * out_nel, 64)
      allocate(tmp(NREC, NGP, nnew)); tmp = -1d30
      if (out_nel > 0) tmp(:, :, 1:out_nel) = out_rec
      call move_alloc(tmp, out_rec); out_nel = nnew
    end if
    do q = 1, NGP
      call gauss(coords, q, N, dN, w, x)
      ph = 0d0; c = 0d0; a = 0d0
      do i = 1, NN
        ph = ph + N(i) * U(NDN*(i-1)+4); c = c + N(i) * U(NDN*(i-1)+5); a = a + N(i) * U(NDN*(i-1)+6)
      end do
      p = sum(sig(1:3, q)) / 3d0
      s = sig(:, q); s(1:3) = s(1:3) - p
      vm = sqrt(1.5d0 * (sum(s(1:3)**2) + 2d0 * sum(s(4:6)**2)))
      sh = 0d0
      if (sum(sv(1:eco_nact, q)) > 0d0) sh = sv(1, q) / sum(sv(1:eco_nact, q))
      out_rec(:, q, jel) = [x(1), x(2), x(3), ph, c, a, sv(13, q), living(sv(1:12, q)) * sv(14, q), vm, p, sh, sv(14, q)]
    end do
  end subroutine keep_record

  subroutine write_records(kstep, kinc, ttot)
    integer, intent(in) :: kstep, kinc
    double precision, intent(in) :: ttot
    integer :: u, e, q, st
    character(len=512) :: path
    if (out_nel == 0) return
    call get_environment_variable('BIOFILM_UEL_GP', path, status=st)
    if (st /= 0) path = 'living_gp.txt'
    open(newunit=u, file=trim(path), position='append', status='unknown', action='write')
    write(u, '(a,i4,a,i6,a,es16.8)') '# step ', kstep, ' increment ', kinc, ' total time ', ttot
    write(u, '(a)') '# element gp x y z phi c a alpha-1 living von_Mises pressure shareA started'
    do e = 1, out_nel
      if (out_rec(1, 1, e) < -1d29) cycle
      do q = 1, NGP
        write(u, '(i8,i3,12es16.8)') e, q, out_rec(:, q, e)
      end do
    end do
    close(u)
  end subroutine write_records
end module living_uel_mod

! ---- Abaqus entry point ---------------------------------------------------------------------------
! *USER ELEMENT, NODES=8, TYPE=U1, UNSYMM, PROPERTIES=18, COORDINATES=3, VARIABLES=128
! 1, 2, 3, 11, 12, 13
subroutine UEL(RHS, AMATRX, SVARS, ENERGY, NDOFEL, NRHS, NSVARS, PROPS, NPROPS, COORDS, MCRD, NNODE, &
               U, DU, V, A, JTYPE, TIME, DTIME, KSTEP, KINC, JELEM, PARAMS, NDLOAD, JDLTYP, ADLMAG, &
               PREDEF, NPREDF, LFLAGS, MLVARX, DDLMAG, MDLOAD, PNEWDT, JPROPS, NJPROP, PERIOD)
  use living_uel_mod
  implicit none
  integer :: NDOFEL, NRHS, NSVARS, NPROPS, MCRD, NNODE, JTYPE, KSTEP, KINC, JELEM, NDLOAD, NPREDF, &
             MLVARX, MDLOAD, NJPROP
  integer :: JDLTYP(MDLOAD, *), LFLAGS(*), JPROPS(*)
  double precision :: RHS(MLVARX, *), AMATRX(NDOFEL, NDOFEL), SVARS(*), ENERGY(8), PROPS(*), &
                      COORDS(MCRD, NNODE), U(NDOFEL), DU(MLVARX, *), V(NDOFEL), A(NDOFEL), TIME(2), &
                      DTIME, PARAMS(*), ADLMAG(MDLOAD, *), PREDEF(2, NPREDF, NNODE), DDLMAG(MDLOAD, *), &
                      PNEWDT, PERIOD
  double precision :: Un(NDOF), sv(NSV, NGP), R(NDOF), K(NDOF, NDOF), sig(6, NGP)
  logical :: ok, wk
  if (.not. eco_inited) call eco_native_init()
  mms_on = .false.
  if (NPROPS >= 19) mms_on = (PROPS(19) > 0.5d0)
  mms_t = TIME(1) + DTIME
  Un = U(1:NDOF) - DU(1:NDOF, 1)
  sv = reshape(SVARS(1:NSV*NGP), [NSV, NGP])
  ! the state variables are those of the start of the increment at every iteration: the point model is advanced
  ! from them on every call (the same result each iteration) and the end-of-increment values are written back
  if (DTIME > 0d0) then
    call advance(Un, COORDS, sv, PROPS, DTIME, ok)
    if (.not. ok) then
      PNEWDT = 0.5d0
      return
    end if
  end if
  wk = (LFLAGS(3) == 1 .or. LFLAGS(3) == 2)
  call element(U, Un, COORDS, sv, PROPS, DTIME, R, K, wk, sig)
  RHS(1:NDOF, 1) = -R
  if (wk) AMATRX = K
  if (DTIME > 0d0) call alpha_store(U, COORDS, sv, PROPS, DTIME)
  call keep_record(JELEM, U, COORDS, sv, sig)
  SVARS(1:NSV*NGP) = reshape(sv, [NSV*NGP])
end subroutine UEL

subroutine UEXTERNALDB(LOP, LRESTART, TIME, DTIME, KSTEP, KINC)
  ! LOP 0: the point model's constants, once, before any element call; LOP 2 (end of an increment): the Gauss-point
  ! records every BIOFILM_UEL_EVERY-th increment (default 40) of step 2; LOP 3 (end of the analysis): the records
  use living_uel_mod
  implicit none
  integer :: LOP, LRESTART, KSTEP, KINC, st, every
  double precision :: TIME(2), DTIME
  character(len=32) :: s
  if (LOP == 0) then
    if (.not. eco_inited) call eco_native_init()
  else if (LOP == 2) then
    every = 40
    call get_environment_variable('BIOFILM_UEL_EVERY', s, status=st)
    if (st == 0) read(s, *) every
    if (KSTEP >= 2 .and. mod(KINC, every) == 0) call write_records(KSTEP, KINC, TIME(2))
  else if (LOP == 3) then
    call write_records(KSTEP, KINC, TIME(2))
  end if
end subroutine UEXTERNALDB
