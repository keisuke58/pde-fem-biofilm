! test_tangent.f90 -- the closed-form tangent of living_uel.f90 against central differences of its residual, for one
! distorted element with random fields, started Gauss points and phi inside (0, 1) (A-46).
!   gfortran -cpp -ffixed-form -ffixed-line-length-132 -c living_uel_abaqus.for
!   gfortran -ffree-line-length-none test_tangent.f90 living_uel_abaqus.o -o test_tangent && ./test_tangent [supg] [gate|gate1|gate2|igate|igate1|igate2]
! Arguments (any order): supg = SUPG (A-54); h=<step of the displacement dofs> (default 1e-9, implicit gate 1e-7); gate, gate1, gate2 = the homeostatic-pressure law with the explicit gate,
! forms 0, 1, 2, p_ref chosen so that 0 < g < 1 at the Gauss points; igate, igate1, igate2 = the implicit gate (27 Oct 2026).
! With the law on, the effective gate g of every Gauss point ((alpha - alpha_n) / (dt k_alpha phi)) and the relative
! change of the tangent against the law switched off are printed, so a check sees that the gate was active.
program test_tangent
  use living_uel_mod
  implicit none
  double precision :: U(NDOF), Un(NDOF), coords(3, NN), sv(NSV, NGP), sv2(NSV, NGP), props(18), R(NDOF), K(NDOF, NDOF)
  double precision :: Rp(NDOF), Rm(NDOF), Kfd(NDOF, NDOF), Kd(NDOF, NDOF), K0(NDOF, NDOF), sig(6, NGP), h, err, sc
  double precision :: N(8), dN(3, 8), w, x(3), ph, geff(NGP), pref0, dt, hu
  integer :: i, j, q, na
  character(len=16) :: arg
  double precision, parameter :: xr(8) = [-1d0, 1d0, 1d0, -1d0, -1d0, 1d0, 1d0, -1d0]
  double precision, parameter :: xs(8) = [-1d0, -1d0, 1d0, 1d0, -1d0, -1d0, 1d0, 1d0]
  double precision, parameter :: xt(8) = [-1d0, -1d0, -1d0, -1d0, 1d0, 1d0, 1d0, 1d0]
  eco_nact = 2
  dt = 0.005d0
  call random_seed()
  do i = 1, NN
    call random_number(coords(:, i))
    coords(:, i) = 0.125d0 * [xr(i), xs(i), xt(i)] + 0.02d0 * (coords(:, i) - 0.5d0)
  end do
  call random_number(U); call random_number(Un)
  do i = 1, NN
    U(NDN*(i-1)+1:NDN*(i-1)+3) = 1d-4 * (U(NDN*(i-1)+1:NDN*(i-1)+3) - 0.5d0)
    Un(NDN*(i-1)+1:NDN*(i-1)+3) = 1d-5 * (Un(NDN*(i-1)+1:NDN*(i-1)+3) - 0.5d0)
    U(NDN*(i-1)+4) = 0.2d0 + 0.6d0 * U(NDN*(i-1)+4)
    Un(NDN*(i-1)+4) = 0.2d0 + 0.6d0 * Un(NDN*(i-1)+4)
    U(NDN*(i-1)+5) = 0.8d0 + 0.2d0 * U(NDN*(i-1)+5); Un(NDN*(i-1)+5) = 0.8d0 + 0.2d0 * Un(NDN*(i-1)+5)
  end do
  call random_number(sv)
  sv(14, :) = 1d0
  sv(13, :) = 1d-4 * (0.3d0 + 0.4d0 * sv(13, :))       ! p_n of about 0.005-0.012 Pa at the Gauss points (E = 10 Pa)
  props = [10d0, 0.49d0, 1d-3, 0.645d0, 1d-3, 5d0, 10d0, 1d0, 1d0, 1d0, 1d-3, 0.1d0, 5d0, 1d0, 0.5625d0, 0.55d0, &
           0.275d0, 0.275d0]
  supg_on = .false.; gate_pref = 0d0; gate_form = 0; gate_implicit = .false.
  hu = 1d-9                                             ! step of the displacement dofs
  do na = 1, command_argument_count()
    call get_command_argument(na, arg)
    if (arg == 'supg') supg_on = .true.
    if (arg(1:4) == 'gate') gate_pref = 0.03d0
    if (arg(1:5) == 'igate') then
      gate_pref = 0.02d0; gate_implicit = .true.
      hu = 1d-7      ! the residual is linear in u between the kinks, so a larger step only lowers the round-off of the local solve
    end if
    if (arg(1:2) == 'h=') read(arg(3:), *) hu
    if (arg == 'gate1' .or. arg == 'igate1') gate_form = 1
    if (arg == 'gate2' .or. arg == 'igate2') gate_form = 2
  end do
  if (gate_form == 1) gate_pref = gate_pref / 0.3d0    ! p_h = p_ref (phi_n^2 + f), phi_n about 0.2-0.8
  call element(U, Un, coords, sv, props, dt, R, K, .true., sig)
  do j = 1, NDOF
    h = 1d-6 * max(1d0, abs(U(j)))
    if (mod(j - 1, NDN) < 3) h = hu
    U(j) = U(j) + h; call element(U, Un, coords, sv, props, dt, Rp, Kd, .false., sig)
    U(j) = U(j) - 2 * h; call element(U, Un, coords, sv, props, dt, Rm, Kd, .false., sig)
    U(j) = U(j) + h
    Kfd(:, j) = (Rp - Rm) / (2 * h)
  end do
  sc = maxval(abs(K))
  err = maxval(abs(K - Kfd)) / sc
  print '(a,es10.3,a,es10.3)', 'max |K - K_fd| / max |K| = ', err, ',  max |K| = ', sc
  do i = 1, NDN
    print '(a,i2,a,es10.3)', '  dof type ', i, ': max |K - K_fd| in its rows ', &
      maxval(abs(K(i:NDOF:NDN, :) - Kfd(i:NDOF:NDN, :))) / max(maxval(abs(K(i:NDOF:NDN, :))), 1d-300)
  end do
  if (gate_pref > 0d0) then
    ! the effective gate of every Gauss point, from the stored growth, and the tangent against the law off
    sv2 = sv
    call alpha_store(U, Un, coords, sv2, props, dt)
    do q = 1, NGP
      call gauss(coords, q, N, dN, w, x)
      ph = 0d0
      do i = 1, NN
        ph = ph + N(i) * U(NDN*(i-1)+4)
      end do
      geff(q) = (sv2(13, q) - sv(13, q)) / (dt * props(5) * ph)
    end do
    pref0 = gate_pref; gate_pref = 0d0
    call element(U, Un, coords, sv, props, dt, Rp, K0, .true., sig)
    gate_pref = pref0
    print '(a,8f8.4)', 'gate g at the Gauss points:', geff
    print '(a,es10.3)', 'tangent changed by the law (relative): ', maxval(abs(K - K0)) / sc
  end if
  if (err > 1d-6) stop 1
end program test_tangent
