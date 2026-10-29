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
!   c   : int N (c - c_n)/(d dt) + grad N . grad c + N q_c phi_c L c dV = 0   (Eq. 35 with its time derivative, backward
!         Euler; d of Klempt 2024 Table 2 in mm^2/T*, PROPS(25); PROPS(25) <= 0 or absent: the limit d -> infinity;
!         q_c = g l^2/d),
!         phi_c = min(max(phi, 0), 1), L = sum_i phi_i psi_i of the point model (start of the increment)
!   a   : the same with q_a
!   u   : int Bbar^T sigma dV = 0, sigma = E_bio (phi_c^2 + f) C(nu) : (eps_bar - (alpha - 1) I), small strain,
!         B-bar (as wp3_stress_reference.py)
!   point model: at the start of an increment, each started Gauss point is advanced over s dt with c and a linearly
!         interpolated from the previous increment's start values to the present ones (biofilm_ecology_hook_lin, as
!         ANSYS mode 11 and the Python reference); a point starts when phi_n >= phi_min, from the start state of its
!         position (separate: species B where |x| and |y| > x_split, A elsewhere; mixed: the mixed start).
! PROPS (NPROPS = 18): E_bio, nu, f_void, beta, k_alpha, penalty, r, K, q_c, q_a, phi_min, s, n_sub, layout (1 separate,
!   2 mixed, 3 mixed start, then inherited from the neighbours: A-74; 4 as 3 with the living composition and a fresh
!   psi: A-82), x_split, total start amount,
!   mixed start of species 1 and 2
! Optional PROPS (absent or 0: off): 19 manufactured solution (mms), 20 SUPG, 21 rule for newly reached points,
!   22, 23 homeostatic-pressure growth law p_ref and form (below; a Keio extension, not part of the thesis model),
!   24 front speed times the living amount (below; a Keio extension, not part of the thesis model)
!
! Front speed times the living amount (28 Oct 2026; coupled_living_consumption.py --front-living, KEIO_PLAN.ja.md
! section 10): PROPS(24) = 1 multiplies the front speed r c/(K + c) of Eq. 34 by sum_i phibar_i, taken as the largest
! value of the Gauss point and its six nearest Gauss points on the regular grid (0 where a point has not started), so
! that a front of dead biofilm does not advance.  The value is that of the end of the last increment (UEXTERNALDB,
! LOP 2, front_living_fill, with the lattice of inherit_fill), i.e. the start of this one, as the Python reference
! (living fraction of the start of the outer step, the cell and its six neighbours).  It is constant within the
! increment, so the tangent keeps its form.  Difference to the Python stencil: Gauss points sit h/sqrt(3) and
! h (1 - 1/sqrt(3)) apart (the cells of the Python code h), so the max reaches 0.42 h to 0.58 h instead of h.
! Before the first record (the 1e-8 start step) the factor is 1.
!
! Homeostatic-pressure growth law (27 Oct 2026; the law of p3/uel/unified_uel.f90 PROPS(19:20), same convention,
! here PROPS(22:23) because 19-21 are taken): Eq. 36 becomes
!   alpha_dot = k_alpha phi g,   g = min(1, max(0, 1 - p / p_h)),   p = -tr(sigma)/3 (positive in compression),
! PROPS(22) = p_ref (<= 0 or absent: off, g = 1 exactly, every result unchanged to the last bit), PROPS(23) = form:
!   0  p_h = p_ref;  1  p_h = p_ref (phi_n^2 + f) (phi_n of the start of the increment);  2  g = 1/(1 + max(p, 0)/p_ref)
!   form + 10 (10, 11, 12): the implicit gate.
! Explicit gate (forms 0-2): p of the start of the increment, recomputed at the Gauss point from U_n, phi_n and
!   alpha_n (small strain: sigma is a function of these, so this is the converged pressure of the previous increment;
!   no extra state variable).  g is a constant within the increment and the tangent keeps its form with dt k_alpha
!   replaced by dt k_alpha g.  p overshoots p_h by at most one increment of growth (first order).
! Implicit gate (forms 10-12): p of the end of the increment, p = p_0 + 3 E' Kb xi (linear in xi = alpha - 1 for
!   small strain, E' = E_bio (phi_c^2 + f), Kb = 1/(3 (1 - 2 nu)), p_0 = -E' Kb tr eps with the B-bar strain), so
!   xi = xi_n + dt k_alpha phi g(p_0 + 3 E' Kb xi) is a monotone scalar equation per Gauss point (solve_xi, safeguarded
!   Newton to round-off; no coupling between Gauss points); the tangent takes d xi / d phi and d xi / d u from the
!   implicit function theorem (one-sided dg/dp at the kinks).  No overshoot: p <= p_h where g reached 0.

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
  ! A-54, SUPG (Brooks and Hughes 1982) for the front term instead of the added diffusion r h/4: PROPS(20) = 1 adds
  ! tau (v . grad N) times the strong residual of the phi equation (consistent, so no O(h) term); beta = beta_0 then
  logical, save :: supg_on = .false.
  ! A-74, layout 3 (PROPS(14) = 3): a Gauss point reached in step 2 takes the phi-weighted mean state of its started
  ! neighbours (the six nearest Gauss points on the regular grid), as the Python codes (coupled_living_consumption.py,
  ! maze_coupled_2d.py); while no point has started yet (the colony), the mixed start state.  The element does not know its neighbours, so the
  ! end of each increment (UEXTERNALDB, LOP 2) collects every Gauss point (keep_record), fills layer by layer until
  ! every reached point connected to started ones has a state (inherit), and the next increment's advance takes it.
  double precision, allocatable, save :: gp_g(:,:,:), inh_g(:,:,:)
  logical, allocatable, save :: inh_ok(:,:)
  integer, save :: layout_mode = 0
  ! rule for newly reached points (A-82): 0 start state (layout 1 by position, 2 mixed), 1 the neighbours' state,
  ! 2 the neighbours' living composition with a fresh psi; PROPS(21), or layout 3 / 4 (= mixed start with rule 1 / 2)
  integer, save :: rule_mode = 0
  ! homeostatic-pressure growth law (header): PROPS(22) p_ref, PROPS(23) form (>= 10: implicit); set in UEL or by a driver
  double precision, save :: gate_pref = 0d0
  integer, save :: gate_form = 0
  logical, save :: gate_implicit = .false.
  ! front speed times the living amount (header): PROPS(24); liv_front(q, e) is written only in UEXTERNALDB (serial)
  logical, save :: front_liv_on = .false.
  ! Eq. 35 with its time derivative (29 Oct 2026): PROPS(25) = d in mm^2/T* adds N (c - c_n)/(d dt) to the nutrient and
  ! antibiotic rows (backward Euler, as coupled_living_consumption.py --transient-c); <= 0 or absent: the limit
  double precision, save :: cdot_d = 0d0
  double precision, allocatable, save :: liv_front(:,:)
  double precision, save :: phi_min_s = 1d-3, h_el = 0d0
  logical, save :: any_started = .false.               ! some Gauss point has started (end of the last increment)
  ! Thread-parallel Abaqus (cpus > 1, mp_mode threads; 25 Oct 2026): the elements run UEL at the same time, and
  ! keep_record grows and writes the shared arrays out_rec, gp_g (and h_el), so it is under Abaqus mutex 1 (initialised
  ! in UEXTERNALDB, LOP 0).  inh_g, inh_ok, any_started are written only in UEXTERNALDB (LOP 2, serial) and read in UEL.
  ! The gfortran builds (test_tangent.f90) have no Abaqus library: without ABQ_FORTRAN (defined by Abaqus's
  ! compile_fortran) the three wrappers are empty.
#ifdef ABQ_FORTRAN
#include <SMAAspUserSubroutines.hdr>
#endif
contains

  subroutine shared_init()
#ifdef ABQ_FORTRAN
    call MutexInit(1)
#endif
  end subroutine shared_init

  subroutine shared_lock()
#ifdef ABQ_FORTRAN
    call MutexLock(1)
#endif
  end subroutine shared_lock

  subroutine shared_unlock()
#ifdef ABQ_FORTRAN
    call MutexUnlock(1)
#endif
  end subroutine shared_unlock

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
      if (n == 4) then
        ! four species (A-79, as the Python GROUP_A): species 1, 2 play A (seed), 3, 4 play B (corners)
        if (abs(x(1)) > props(15) .and. abs(x(2)) > props(15)) then
          g(3:4) = (tot - 2d-6) / 2d0
        else
          g(1:2) = (tot - 2d-6) / 2d0
        end if
      else if (abs(x(1)) > props(15) .and. abs(x(2)) > props(15)) then
        g(2) = tot - 1d-6 * (n - 1)
      else
        g(1) = tot - 1d-6 * (n - 1)
      end if
    else if (n == 4) then
      g(1:4) = tot / 4d0                        ! the start of Klempt 2026 Table 2: 0.02 each
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

  subroutine bmatrix(dN, Bv, B)
    ! the B-bar strain-displacement matrix of a Gauss point (volumetric part from the element mean Bv)
    double precision, intent(in) :: dN(3, 8), Bv(3, 8)
    double precision, intent(out) :: B(6, 24)
    integer :: i, d, j
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
  end subroutine bmatrix

  subroutine gate_eval(p, phh, g, dg)
    ! g(p) and the one-sided dg/dp of the homeostatic law (gate_form 0, 1: p_h = phh; 2: hyperbolic with gate_pref)
    double precision, intent(in) :: p, phh
    double precision, intent(out) :: g, dg
    g = 1d0; dg = 0d0
    if (gate_pref <= 0d0) return
    if (gate_form == 2) then
      g = 1d0 / (1d0 + max(p, 0d0) / gate_pref)
      if (p > 0d0) dg = -g * g / gate_pref
      return
    end if
    if (phh <= 0d0) return
    g = min(1d0, max(0d0, 1d0 - p / phh))
    if (p > 0d0 .and. p < phh) dg = -1d0 / phh
  end subroutine gate_eval

  double precision function gate_ph(phn, f)
    ! p_h of forms 0 and 1 from phi_n of the start of the increment
    double precision, intent(in) :: phn, f
    gate_ph = gate_pref
    if (gate_form == 1) gate_ph = gate_pref * (min(max(phn, 0d0), 1d0)**2 + f)
  end function gate_ph

  subroutine solve_xi(xin, A, p0, B3, phh, xi, g, dg)
    ! implicit gate: the root of h(xi) = xi - xin - A g(p0 + B3 xi) (A = dt k_alpha phi, B3 = 3 E' Kb >= 0); h increases
    ! with xi, the root lies between xin and xin + A; safeguarded Newton (bisection when the step leaves the bracket)
    double precision, intent(in) :: xin, A, p0, B3, phh
    double precision, intent(out) :: xi, g, dg
    double precision :: lo, hi, h, dh, xn, tol
    integer :: it
    lo = min(xin, xin + A); hi = max(xin, xin + A)
    tol = 1d-15 * (abs(xin) + abs(A)) + 1d-300
    call gate_eval(p0 + B3 * xin, phh, g, dg)
    xi = xin + A * g
    do it = 1, 100
      call gate_eval(p0 + B3 * xi, phh, g, dg)
      h = xi - xin - A * g
      if (abs(h) <= tol) return
      if (h > 0d0) then
        hi = xi
      else
        lo = xi
      end if
      if (hi - lo <= tol) exit
      dh = 1d0 - A * dg * B3
      xn = 0.5d0 * (lo + hi)
      if (dh > 0d0) then
        xn = xi - h / dh
        if (xn <= lo .or. xn >= hi) xn = 0.5d0 * (lo + hi)
      end if
      xi = xn
    end do
    call gate_eval(p0 + B3 * xi, phh, g, dg)
  end subroutine solve_xi

  subroutine advance(Un, coords, sv, props, dt, ok, jel, kstep)
    ! point model at the start of the increment: advance the started points, start the new ones
    double precision, intent(in) :: Un(NDOF), coords(3, NN), props(*), dt
    double precision, intent(inout) :: sv(NSV, NGP)
    logical, intent(out) :: ok
    integer, intent(in) :: jel, kstep
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
        if (rule_mode >= 1 .and. kstep >= 2 .and. any_started) then
          ! layout 3: the state inherited at the end of the last increment (none yet: wait, as the Python fill)
          if (allocated(inh_ok)) then
            if (jel <= size(inh_ok, 2)) then
              if (inh_ok(q, jel)) then
                sv(1:12, q) = inh_g(:, q, jel)
                sv(14, q) = 1d0
              end if
            end if
          end if
        else
          call start_state(props, x, sv(1:12, q))
          sv(14, q) = 1d0
        end if
      end if
      sv(15, q) = c; sv(16, q) = a
    end do
  end subroutine advance

  subroutine element(U, Un, coords, sv, props, dt, R, K, want_k, sig, fl)
    ! residual R(U) and tangent K = dR/dU with the point-model state sv held (already advanced); fl (optional): the
    ! factor of the front speed per Gauss point (PROPS(24), header)
    double precision, intent(in) :: U(NDOF), Un(NDOF), coords(3, NN), props(*), dt
    double precision, intent(inout) :: sv(NSV, NGP)
    double precision, intent(out) :: R(NDOF), K(NDOF, NDOF), sig(6, NGP)
    logical, intent(in) :: want_k
    double precision, intent(in), optional :: fl(NGP)
    double precision :: N(8), dN(3, 8), w, x(3), ph, phn, c, cn, a, an, cdt, gph(3), gc(3), gcn(3), ga(3), v(3)
    double precision :: al, pen, dpen, phc, dphc, L, E, dE, D0(6, 6), lam, mu, eps(6), sm(6), B(6, 24), Bv(3, 8)
    double precision :: vol, m(6), dsig(6), nrm, sp, dt_, sphi, sc, phim, cm, tau, vm, hs, pe, res, vdn(8)
    double precision :: Ebio, nu, f, beta, ka, P, rr, Km, qc, qa
    double precision :: g, dg, kal, dalp, Kb, epsn(6), phcn, En, pn, Ag, B3, p0, xi, dhdxi, dalu(3, 8), D0m(6)
    integer :: q, i, j, ii, jj, d
    Ebio = props(1); nu = props(2); f = props(3); beta = props(4); ka = props(5); P = props(6)
    rr = props(7); Km = props(8); qc = props(9); qa = props(10)
    dt_ = max(dt, 1d-30)
    cdt = 0d0
    if (cdot_d > 0d0 .and. dt > 0d0 .and. .not. mms_on) cdt = 1d0 / (cdot_d * dt)
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
      ph = 0d0; phn = 0d0; c = 0d0; cn = 0d0; a = 0d0; an = 0d0; gph = 0d0; gc = 0d0; gcn = 0d0; ga = 0d0
      do i = 1, NN
        ph = ph + N(i) * U(NDN*(i-1)+4); phn = phn + N(i) * Un(NDN*(i-1)+4)
        c = c + N(i) * U(NDN*(i-1)+5); cn = cn + N(i) * Un(NDN*(i-1)+5); a = a + N(i) * U(NDN*(i-1)+6)
        an = an + N(i) * Un(NDN*(i-1)+6)
        gph = gph + dN(:, i) * U(NDN*(i-1)+4); gc = gc + dN(:, i) * U(NDN*(i-1)+5)
        gcn = gcn + dN(:, i) * Un(NDN*(i-1)+5); ga = ga + dN(:, i) * U(NDN*(i-1)+6)
      end do
      nrm = sqrt(sum(gcn**2) + EPSG)
      sp = rr * cn / (Km + cn)
      if (present(fl)) sp = sp * fl(q)                 ! PROPS(24): times the living amount (header)
      v = sp * gcn / nrm
      pen = P * (max(ph - 1d0, 0d0) + min(ph, 0d0))
      dpen = 0d0
      if (ph > 1d0 .or. ph < 0d0) dpen = P
      phc = min(max(ph, 0d0), 1d0)
      dphc = 0d0
      if (ph > 0d0 .and. ph < 1d0) dphc = 1d0
      call bmatrix(dN, Bv, B)
      eps = 0d0
      do i = 1, NN
        do d = 1, 3
          eps = eps + B(:, 3*(i-1)+d) * U(NDN*(i-1)+d)
        end do
      end do
      E = Ebio * (phc**2 + f)
      dE = Ebio * 2d0 * phc * dphc
      ! growth over the increment (Eq. 36; the homeostatic law of the header when gate_pref > 0)
      g = 1d0; dalp = 0d0; dalu = 0d0; kal = ka * ka * dt
      if (gate_pref > 0d0) then
        Kb = lam + 2d0 * mu / 3d0
        if (.not. gate_implicit) then
          ! explicit gate: p of the start of the increment from U_n, phi_n, alpha_n
          phcn = min(max(phn, 0d0), 1d0)
          En = Ebio * (phcn**2 + f)
          epsn = 0d0
          do i = 1, NN
            do d = 1, 3
              epsn = epsn + B(:, 3*(i-1)+d) * Un(NDN*(i-1)+d)
            end do
          end do
          pn = -En * Kb * (sum(epsn(1:3)) - 3d0 * sv(13, q))
          call gate_eval(pn, gate_ph(phn, f), g, dg)
          kal = ka * ka * dt * g
        else
          Ag = dt * ka * ph
          B3 = 3d0 * E * Kb
          p0 = -E * Kb * sum(eps(1:3))
          call solve_xi(sv(13, q), Ag, p0, B3, gate_ph(phn, f), xi, g, dg)
          dhdxi = 1d0 - Ag * dg * B3
          dalp = (dt * ka * g + Ag * dg * (-dE * Kb * sum(eps(1:3)) + xi * 3d0 * dE * Kb)) / dhdxi
          dalu = -Ag * dg * E * Kb * Bv / dhdxi
          kal = ka * dalp
        end if
      end if
      if (gate_implicit .and. gate_pref > 0d0) then
        al = 1d0 + xi
      else
        al = 1d0 + sv(13, q) + dt * ka * ph * g           ! alpha at the end of the increment
      end if
      L = 0d0
      if (sv(14, q) > 0.5d0) L = living(sv(1:12, q))
      sphi = 0d0; sc = 0d0
      if (mms_on) then
        L = 1d0
        call mms_eval(x(1), x(2), x(3), mms_t, beta, rr, Km, qc, sphi, sc, phim, cm)
      end if
      ! phase
      tau = 0d0
      do i = 1, NN
        vdn(i) = dot_product(v, dN(:, i))
      end do
      vm = sqrt(sum(v**2))
      if (supg_on .and. vm > 1d-12) then
        hs = 2d0 * vm / sum(abs(vdn))                    ! element length along v (Tezduyar)
        pe = vm * hs / (2d0 * beta)
        if (pe > 1d-3) then
          tau = hs / (2d0 * vm) * (1d0 / tanh(pe) - 1d0 / pe)
        else
          tau = hs / (2d0 * vm) * pe / 3d0
        end if
      end if
      ! strong residual of the phi equation (the Laplacian of the trilinear field is dropped, as usual)
      res = (ph - phn) / dt_ - ka * al + pen + dot_product(v, gph) - sphi
      do i = 1, NN
        ii = NDN*(i-1)+4
        R(ii) = R(ii) + w * (N(i) * res + beta * dot_product(dN(:, i), gph) + tau * vdn(i) * res)
      end do
      ! nutrient and antibiotic
      do i = 1, NN
        ii = NDN*(i-1)+5
        R(ii) = R(ii) + w * (N(i) * cdt * (c - cn) + dot_product(dN(:, i), gc) + N(i) * qc * phc * L * c - N(i) * sc)
        ii = NDN*(i-1)+6
        R(ii) = R(ii) + w * (N(i) * cdt * (a - an) + dot_product(dN(:, i), ga) + N(i) * qa * phc * L * a)
      end do
      ! mechanics
      sm = matmul(D0, eps - (al - 1d0) * m)
      sig(:, q) = E * sm
      do i = 1, NN
        do d = 1, 3
          ii = NDN*(i-1)+d
          R(ii) = R(ii) + w * dot_product(B(:, 3*(i-1)+d), sig(:, q))
        end do
      end do
      if (.not. want_k) cycle
      if (gate_implicit .and. gate_pref > 0d0) then
        dsig = dE * sm - E * dalp * matmul(D0, m)        ! d sigma / d phi, d alpha / d phi of the implicit gate
      else
        dsig = dE * sm - E * dt * ka * g * matmul(D0, m)   ! d sigma / d phi
      end if
      D0m = matmul(D0, m)
      do i = 1, NN
        do j = 1, NN
          ! phase rows
          ii = NDN*(i-1)+4; jj = NDN*(j-1)+4
          K(ii, jj) = K(ii, jj) + w * (N(i) * N(j) * (1d0 / dt_ - kal + dpen) &
                                       + N(i) * dot_product(v, dN(:, j)) + beta * dot_product(dN(:, i), dN(:, j)) &
                                       + tau * vdn(i) * (N(j) * (1d0 / dt_ - kal + dpen) + vdn(j)))
          if (gate_implicit .and. gate_pref > 0d0) then
            ! implicit gate: alpha depends on the displacements through p (phase and mechanics rows)
            do d = 1, 3
              K(ii, NDN*(j-1)+d) = K(ii, NDN*(j-1)+d) - w * ka * dalu(d, j) * (N(i) + tau * vdn(i))
              do jj = 1, 3
                K(NDN*(i-1)+jj, NDN*(j-1)+d) = K(NDN*(i-1)+jj, NDN*(j-1)+d) &
                  - w * E * dot_product(B(:, 3*(i-1)+jj), D0m) * dalu(d, j)
              end do
            end do
            jj = NDN*(j-1)+4
          end if
          ! nutrient rows
          ii = NDN*(i-1)+5
          K(ii, NDN*(j-1)+5) = K(ii, NDN*(j-1)+5) + w * (dot_product(dN(:, i), dN(:, j)) + N(i) * N(j) * (qc * phc * L + cdt))
          K(ii, jj) = K(ii, jj) + w * N(i) * N(j) * qc * dphc * L * c
          ii = NDN*(i-1)+6
          K(ii, NDN*(j-1)+6) = K(ii, NDN*(j-1)+6) + w * (dot_product(dN(:, i), dN(:, j)) + N(i) * N(j) * (qa * phc * L + cdt))
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

  subroutine alpha_store(U, Un, coords, sv, props, dt)
    ! alpha - 1 at the end of the increment (after the last, converged call); the same growth as element()
    double precision, intent(in) :: U(NDOF), Un(NDOF), coords(3, NN), props(*), dt
    double precision, intent(inout) :: sv(NSV, NGP)
    double precision :: N(8), dN(3, 8), w, x(3), ph, phn, g, dg, Bv(3, 8), B(6, 24), eps(6), vol, Kb, lam, mu
    double precision :: phc, E, A, B3, p0, xi
    integer :: q, i, d
    if (gate_pref > 0d0) then
      lam = props(2) / ((1 + props(2)) * (1 - 2 * props(2))); mu = 0.5d0 / (1 + props(2))
      Kb = lam + 2d0 * mu / 3d0
      Bv = 0d0; vol = 0d0
      do q = 1, NGP
        call gauss(coords, q, N, dN, w, x)
        Bv = Bv + w * dN; vol = vol + w
      end do
      Bv = Bv / vol
    end if
    do q = 1, NGP
      call gauss(coords, q, N, dN, w, x)
      ph = 0d0; phn = 0d0
      do i = 1, NN
        ph = ph + N(i) * U(NDN*(i-1)+4); phn = phn + N(i) * Un(NDN*(i-1)+4)
      end do
      g = 1d0
      if (gate_pref > 0d0) then
        call bmatrix(dN, Bv, B)
        if (.not. gate_implicit) then
          phc = min(max(phn, 0d0), 1d0)
          E = props(1) * (phc**2 + props(3))
          eps = 0d0
          do i = 1, NN
            do d = 1, 3
              eps = eps + B(:, 3*(i-1)+d) * Un(NDN*(i-1)+d)
            end do
          end do
          call gate_eval(-E * Kb * (sum(eps(1:3)) - 3d0 * sv(13, q)), gate_ph(phn, props(3)), g, dg)
        else
          phc = min(max(ph, 0d0), 1d0)
          E = props(1) * (phc**2 + props(3))
          eps = 0d0
          do i = 1, NN
            do d = 1, 3
              eps = eps + B(:, 3*(i-1)+d) * U(NDN*(i-1)+d)
            end do
          end do
          A = dt * props(5) * ph
          B3 = 3d0 * E * Kb
          p0 = -E * Kb * sum(eps(1:3))
          call solve_xi(sv(13, q), A, p0, B3, gate_ph(phn, props(3)), xi, g, dg)
          sv(13, q) = xi
          cycle
        end if
      end if
      sv(13, q) = sv(13, q) + dt * props(5) * ph * g
    end do
  end subroutine alpha_store

  subroutine keep_record(jel, U, coords, sv, sig)
    integer, intent(in) :: jel
    double precision, intent(in) :: U(NDOF), coords(3, NN), sv(NSV, NGP), sig(6, NGP)
    double precision, allocatable :: tmp(:,:,:)
    double precision :: N(8), dN(3, 8), w, x(3), ph, c, a, p, s(6), vm, sh
    integer :: q, i, nnew
    ! the whole record under the mutex: another element may be growing (move_alloc) the arrays this one writes
    call shared_lock()
    if (jel > out_nel) then
      nnew = max(jel, 2 * out_nel, 64)
      allocate(tmp(NREC, NGP, nnew)); tmp = -1d30
      if (out_nel > 0) tmp(:, :, 1:out_nel) = out_rec
      call move_alloc(tmp, out_rec)
      allocate(tmp(12, NGP, nnew)); tmp = 0d0
      if (out_nel > 0) tmp(:, :, 1:out_nel) = gp_g
      call move_alloc(tmp, gp_g); out_nel = nnew
    end if
    if (h_el <= 0d0) h_el = maxval(coords(1, :)) - minval(coords(1, :))
    gp_g(:, :, jel) = sv(1:12, :)
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
      if (sum(sv(1:eco_nact, q)) > 0d0) then
        if (eco_nact == 4) then                 ! share of group A (species 1, 2)
          sh = sum(sv(1:2, q)) / sum(sv(1:eco_nact, q))
        else
          sh = sv(1, q) / sum(sv(1:eco_nact, q))
        end if
      end if
      out_rec(:, q, jel) = [x(1), x(2), x(3), ph, c, a, sv(13, q), living(sv(1:12, q)) * sv(14, q), vm, p, sh, sv(14, q)]
    end do
    call shared_unlock()
  end subroutine keep_record

  subroutine living_share(g)
    ! layout 4 (25 Oct): the inherited amount shared among the species as their living biomass phi_i psi_i (growth
    ! spreads the biofilm, what arrives is alive), psi of the start state and gamma = 0; without living biomass the
    ! shares of phi_i (as coupled_living_consumption.py living_start)
    double precision, intent(inout) :: g(12)
    double precision :: tot, lv(5), sl
    integer :: m
    m = eco_nact
    tot = sum(g(1:m))
    lv = 0d0
    lv(1:m) = g(1:m) * g(7:6 + m)
    sl = sum(lv(1:m))
    if (sl > 1d-12) then
      g(1:m) = tot * lv(1:m) / sl
    end if
    g(7:11) = 0.999d0
    g(12) = 0d0
  end subroutine living_share

  subroutine inherit_fill()
    ! layout 3: Gauss points with phi >= phi_min that have not started take the phi-weighted mean state of their
    ! started (or already filled) neighbours; repeated until nothing changes.  Neighbours: the six nearest Gauss points
    ! on the regular grid, found by the lattice index of each point (element index times 2 plus the side of the
    ! element centre, per axis; the Gauss points sit at +-h/(2 sqrt 3) from the centre)
    integer, allocatable :: lat(:,:,:), ke(:,:,:)
    double precision :: wsum, g(12), w
    integer :: e, q, d, n(3), kk(3), s, changed, m, ne, qn
    logical, allocatable :: done(:,:), done0(:,:)
    if (out_nel == 0 .or. h_el <= 0d0) return
    any_started = any(out_rec(12, :, 1:out_nel) > 0.5d0 .and. out_rec(1, :, 1:out_nel) > -1d29)
    call gp_lattice(ke, lat, n)
    if (allocated(inh_ok)) deallocate(inh_ok, inh_g)
    allocate(inh_ok(NGP, out_nel), inh_g(12, NGP, out_nel), done(NGP, out_nel), done0(NGP, out_nel))
    inh_ok = .false.; inh_g = 0d0
    do e = 1, out_nel
      do q = 1, NGP
        done(q, e) = (out_rec(1, 1, e) > -1d29) .and. (out_rec(12, q, e) > 0.5d0)
      end do
    end do
    changed = 1
    do while (changed > 0)                                ! one layer per sweep (from the points done before it)
      changed = 0
      done0 = done
      do e = 1, out_nel
        if (out_rec(1, 1, e) < -1d29) cycle
        do q = 1, NGP
          if (done(q, e) .or. out_rec(4, q, e) < phi_min_s) cycle
          wsum = 0d0; g = 0d0
          do d = 1, 3
            do s = -1, 1, 2
              kk = ke(:, q, e); kk(d) = kk(d) + s
              if (any(kk < 0) .or. any(kk >= n)) cycle
              m = lat(kk(1), kk(2), kk(3))
              if (m == 0) cycle
              ne = (m - 1) / NGP + 1; qn = m - NGP * (ne - 1)
              if (.not. done0(qn, ne)) cycle
              w = max(out_rec(4, qn, ne), 0d0)
              if (out_rec(12, qn, ne) > 0.5d0) then
                g = g + w * gp_g(:, qn, ne)
              else
                g = g + w * inh_g(:, qn, ne)
              end if
              wsum = wsum + w
            end do
          end do
          if (wsum <= 0d0) cycle
          g = g / wsum
          g(1:5) = max(g(1:5), 0d0)
          if (rule_mode == 2) call living_share(g)
          g(1:eco_nact) = max(g(1:eco_nact), 1d-10)
          g(6) = 1d0 - sum(g(1:5))                       ! the point model's own constraint
          inh_g(:, q, e) = g; inh_ok(q, e) = .true.
          done(q, e) = .true.; changed = changed + 1
        end do
      end do
    end do
  end subroutine inherit_fill

  subroutine gp_lattice(ke, lat, n)
    ! the lattice index of every recorded Gauss point (element index times 2 plus the side of the element centre, per
    ! axis; the Gauss points sit at +-h/(2 sqrt 3) from the centre): ke(:, q, e), and lat(k) = NGP (e - 1) + q (0: none)
    integer, allocatable, intent(out) :: ke(:,:,:), lat(:,:,:)
    integer, intent(out) :: n(3)
    double precision :: lo(3)
    integer :: e, q, d, k(3), m
    lo = 1d30
    do e = 1, out_nel
      if (out_rec(1, 1, e) < -1d29) cycle
      do q = 1, NGP
        lo = min(lo, out_rec(1:3, q, e))
      end do
    end do
    lo = lo - (0.5d0 - 0.5d0 / sqrt(3d0)) * h_el            ! the lowest node
    allocate(ke(3, NGP, out_nel)); ke = -1
    n = 0
    do e = 1, out_nel
      if (out_rec(1, 1, e) < -1d29) cycle
      do q = 1, NGP
        do d = 1, 3
          m = floor((out_rec(d, q, e) - lo(d)) / h_el)
          k(d) = 2 * m + merge(1, 0, out_rec(d, q, e) - lo(d) > (m + 0.5d0) * h_el)
        end do
        ke(:, q, e) = k
        n = max(n, k + 1)
      end do
    end do
    allocate(lat(0:n(1)-1, 0:n(2)-1, 0:n(3)-1)); lat = 0
    do e = 1, out_nel
      if (out_rec(1, 1, e) < -1d29) cycle
      do q = 1, NGP
        k = ke(:, q, e)
        lat(k(1), k(2), k(3)) = NGP * (e - 1) + q
      end do
    end do
  end subroutine gp_lattice

  subroutine front_living_fill()
    ! PROPS(24) (header): liv_front(q, e) = the largest living amount (sum_i phibar_i, 0 where not started) of the Gauss
    ! point and its six nearest Gauss points on the regular grid, from the records of the increment that just ended
    integer, allocatable :: lat(:,:,:), ke(:,:,:)
    double precision :: lv
    integer :: e, q, d, s, n(3), kk(3), m, ne, qn
    if (out_nel == 0 .or. h_el <= 0d0) return
    call gp_lattice(ke, lat, n)
    if (allocated(liv_front)) deallocate(liv_front)
    allocate(liv_front(NGP, out_nel)); liv_front = 1d0
    do e = 1, out_nel
      if (out_rec(1, 1, e) < -1d29) cycle
      do q = 1, NGP
        lv = out_rec(8, q, e)
        do d = 1, 3
          do s = -1, 1, 2
            kk = ke(:, q, e); kk(d) = kk(d) + s
            if (any(kk < 0) .or. any(kk >= n)) cycle
            m = lat(kk(1), kk(2), kk(3))
            if (m == 0) cycle
            ne = (m - 1) / NGP + 1; qn = m - NGP * (ne - 1)
            lv = max(lv, out_rec(8, qn, ne))
          end do
        end do
        liv_front(q, e) = lv
      end do
    end do
  end subroutine front_living_fill

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
! *USER ELEMENT, NODES=8, TYPE=U1, UNSYMM, PROPERTIES=18, COORDINATES=3, VARIABLES=128   (PROPERTIES up to 24, header)
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
  logical :: ok, wk, lmms, lsupg, limp, lfl
  integer :: lay, rule, iform
  double precision :: pref, dcd
  if (.not. eco_inited) call eco_native_init()     ! never reached in Abaqus: UEXTERNALDB LOP 0 runs first, serial
  ! the module switches from PROPS: every element finds the same values (one property set), so with threads they are
  ! stored only when they differ (the first call) and read afterwards; mms_t changes each increment and is written
  ! by every element (the same value, before its own use; an aligned 8-byte store)
  lmms = .false.
  if (NPROPS >= 19) lmms = (PROPS(19) > 0.5d0)
  lsupg = .false.
  if (NPROPS >= 20) lsupg = (PROPS(20) > 0.5d0)
  lay = nint(PROPS(14))
  rule = 0
  if (lay == 3) rule = 1
  if (lay == 4) rule = 2
  if (NPROPS >= 21) then
    if (PROPS(21) > 0.5d0) rule = nint(PROPS(21))
  end if
  ! homeostatic-pressure growth law (header): PROPS(22) p_ref, PROPS(23) form (>= 10: implicit); absent: off
  pref = 0d0; iform = 0; limp = .false.
  if (NPROPS >= 22) pref = PROPS(22)
  if (NPROPS >= 23) iform = nint(PROPS(23))
  if (iform >= 10) then
    limp = .true.; iform = iform - 10
  end if
  if (gate_pref /= pref) gate_pref = pref
  if (gate_form /= iform) gate_form = iform
  if (gate_implicit .neqv. limp) gate_implicit = limp
  lfl = .false.
  if (NPROPS >= 24) lfl = (PROPS(24) > 0.5d0)
  if (front_liv_on .neqv. lfl) front_liv_on = lfl
  dcd = 0d0
  if (NPROPS >= 25) dcd = PROPS(25)
  if (cdot_d /= dcd) cdot_d = dcd
  if (mms_on .neqv. lmms) mms_on = lmms
  if (supg_on .neqv. lsupg) supg_on = lsupg
  if (lmms) mms_t = TIME(1) + DTIME
  if (layout_mode /= lay) layout_mode = lay
  if (phi_min_s /= PROPS(11)) phi_min_s = PROPS(11)
  if (rule_mode /= rule) rule_mode = rule
  Un = U(1:NDOF) - DU(1:NDOF, 1)
  sv = reshape(SVARS(1:NSV*NGP), [NSV, NGP])
  ! the state variables are those of the start of the increment at every iteration: the point model is advanced
  ! from them on every call (the same result each iteration) and the end-of-increment values are written back
  if (DTIME > 0d0) then
    call advance(Un, COORDS, sv, PROPS, DTIME, ok, JELEM, KSTEP)
    if (.not. ok) then
      PNEWDT = 0.5d0
      return
    end if
  end if
  wk = (LFLAGS(3) == 1 .or. LFLAGS(3) == 2)
  ! PROPS(24): the factor of the end of the last increment (none yet, in the 1e-8 start step: factor 1)
  lfl = .false.
  if (front_liv_on .and. allocated(liv_front)) lfl = (JELEM <= size(liv_front, 2))
  if (lfl) then
    call element(U, Un, COORDS, sv, PROPS, DTIME, R, K, wk, sig, liv_front(:, JELEM))
  else
    call element(U, Un, COORDS, sv, PROPS, DTIME, R, K, wk, sig)
  end if
  RHS(1:NDOF, 1) = -R
  if (wk) AMATRX = K
  if (DTIME > 0d0) call alpha_store(U, Un, COORDS, sv, PROPS, DTIME)
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
    call shared_init()                                  ! mutex 1 of keep_record, before any element call
    if (.not. eco_inited) call eco_native_init()
  else if (LOP == 2) then
    every = 40
    call get_environment_variable('BIOFILM_UEL_EVERY', s, status=st)
    if (st == 0) read(s, *) every
    if (KSTEP >= 2 .and. mod(KINC, every) == 0) call write_records(KSTEP, KINC, TIME(2))
    if (rule_mode >= 1) call inherit_fill()
    if (front_liv_on) call front_living_fill()
  else if (LOP == 3) then
    call write_records(KSTEP, KINC, TIME(2))
  end if
end subroutine UEXTERNALDB
