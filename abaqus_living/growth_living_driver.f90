! growth_living_driver.f90 -- one element of living_uel.f90 under the growth of the living biomass (PROPS(26), 1 Nov
! 2026; header of living_uel.f90), without Abaqus, in the style of homeostatic_driver.f90: phi = 1, c = 1, a = 0 at
! every node (phi, c, a held), every Gauss point started with a fixed point-model state (phi_1 = <phi_1>, psi_1 = <psi_1>,
! so L = sum_i phibar_i = phi_1 psi_1 is known and constant: the point model is not advanced, consumption is off), no
! front term, no diffusion.  The growth is then
!   alpha_dot = k_alpha phi L = k_alpha L,   alpha(0) = 1,   so alpha - 1 = k_alpha L t exactly (Euler is exact for a
! constant rate), and the source of Eq. 34 is k_alpha alpha L: with phi = 1 held, the phase residual summed over the
! nodes is -k_alpha alpha L times the element volume (1 mm^3).  Newton on the displacement dofs (3-2-1 support, p = 0).
!   gfortran -cpp -ffixed-form -ffixed-line-length-132 -c living_uel_abaqus.for
!   gfortran -ffree-line-length-none growth_living_driver.f90 living_uel_abaqus.o -o growth_living_driver
!   ./growth_living_driver <on 0|1> <phi_1> <psi_1> <steps> <T> [k_alpha]
! Prints one line per increment: k, t, alpha - 1 (Gauss-point mean), its expected value k_alpha L t (L = phi_1 psi_1
! with the law on, 1 off), the summed phase residual, its expected value -k_alpha alpha L, Newton iterations; the last
! line gives the largest relative deviations.
program growth_living_driver
  use living_uel_mod
  implicit none
  double precision :: U(NDOF), Un(NDOF), coords(3, NN), sv(NSV, NGP), props(18), R(NDOF), K(NDOF, NDOF), sig(6, NGP)
  double precision :: T, dt, ka, rn, Kf(NDOF, NDOF), Rf(NDOF), dU(NDOF), ph1, ps1, L, xi, xie, rph, rphe, e1, e2
  integer :: steps, i, kk, it, nf, idx(NDOF), ion
  logical :: free(NDOF)
  character(len=16) :: arg
  double precision, parameter :: xr(8) = [0d0, 1d0, 1d0, 0d0, 0d0, 1d0, 1d0, 0d0]
  double precision, parameter :: xs(8) = [0d0, 0d0, 1d0, 1d0, 0d0, 0d0, 1d0, 1d0]
  double precision, parameter :: xt(8) = [0d0, 0d0, 0d0, 0d0, 1d0, 1d0, 1d0, 1d0]
  if (command_argument_count() < 5) stop 'usage: growth_living_driver on phi_1 psi_1 steps T [k_alpha]'
  call get_command_argument(1, arg); read(arg, *) ion
  call get_command_argument(2, arg); read(arg, *) ph1
  call get_command_argument(3, arg); read(arg, *) ps1
  call get_command_argument(4, arg); read(arg, *) steps
  call get_command_argument(5, arg); read(arg, *) T
  ka = 1d-3
  if (command_argument_count() >= 6) then
    call get_command_argument(6, arg); read(arg, *) ka
  end if
  eco_nact = 2
  growth_liv_on = (ion /= 0)
  gate_pref = 0d0; gate_form = 0; gate_implicit = .false.; supg_on = .false.; mms_on = .false.; cdot_d = 0d0
  ! E = 10 Pa, nu = 0.49, f = 1e-3 (Klempt et al. 2024 Table 2); no front term, no diffusion, no consumption
  props = [10d0, 0.49d0, 1d-3, 0d0, ka, 5d0, 0d0, 1d0, 0d0, 0d0, 1d-3, 0.1d0, 1d0, 2d0, 0.5625d0, 0.55d0, 0.275d0, 0.275d0]
  do i = 1, NN
    coords(:, i) = [xr(i), xs(i), xt(i)]
  end do
  U = 0d0
  do i = 1, NN
    U(NDN*(i-1)+4) = 1d0; U(NDN*(i-1)+5) = 1d0
  end do
  ! every Gauss point started, with the fixed state (phi_1, phi_0 = 1 - phi_1, psi_1; species 2 empty): L = phi_1 psi_1
  sv = 0d0
  sv(1, :) = ph1; sv(6, :) = 1d0 - ph1; sv(7, :) = ps1; sv(8, :) = 0.999d0; sv(14, :) = 1d0
  sv(15, :) = 1d0
  L = living(sv(1:12, 1))
  if (.not. growth_liv_on) L = 1d0
  free = .false.
  do i = 1, NN
    free(NDN*(i-1)+1:NDN*(i-1)+3) = .true.
  end do
  free(1:3) = .false.                       ! node 1 (0,0,0): u1 = u2 = u3 = 0
  free(NDN+2:NDN+3) = .false.               ! node 2 (1,0,0): u2 = u3 = 0
  free(3*NDN+3) = .false.                   ! node 4 (0,1,0): u3 = 0
  nf = 0
  do i = 1, NDOF
    if (free(i)) then
      nf = nf + 1; idx(nf) = i
    end if
  end do
  dt = T / steps
  e1 = 0d0; e2 = 0d0
  print '(a,f10.6)', '# L = sum_i phibar_i at the Gauss points = ', L
  print '(a)', '#    k         t              alpha-1 (mean)        k_alpha L t         sum R_phi            -k_alpha alpha L  it'
  do kk = 1, steps
    Un = U
    ! the point model is held (no advance): L is the given, constant field of this check
    do it = 1, 30
      call element(U, Un, coords, sv, props, dt, R, K, .true., sig)
      rn = maxval(abs(R(idx(1:nf))))
      if (rn < 1d-13) exit
      Kf(1:nf, 1:nf) = K(idx(1:nf), idx(1:nf)); Rf(1:nf) = -R(idx(1:nf))
      call solve(nf, Kf(1:nf, 1:nf), Rf(1:nf))
      dU = 0d0; dU(idx(1:nf)) = Rf(1:nf)
      U = U + dU
    end do
    rph = 0d0
    do i = 1, NN
      rph = rph + R(NDN*(i-1)+4)
    end do
    call alpha_store(U, Un, coords, sv, props, dt)
    xi = sum(sv(13, :)) / NGP
    xie = ka * L * kk * dt
    rphe = -ka * (1d0 + xie) * L                        ! alpha of the end of the increment, element volume 1
    print '(i6,es16.8,4es22.14,i4)', kk, kk * dt, xi, xie, rph, rphe, it - 1
    e1 = max(e1, abs(xi - xie) / max(abs(xie), 1d-300))
    e2 = max(e2, abs(rph - rphe) / max(abs(rphe), 1d-300))
  end do
  print '(a,es10.3,a,es10.3)', '# max relative deviation: alpha - 1 ', e1, ',  phase source ', e2
contains
  subroutine solve(n, A, b)
    ! A x = b by Gaussian elimination with partial pivoting; x returned in b
    integer, intent(in) :: n
    double precision, intent(inout) :: A(n, n), b(n)
    double precision :: t, row(n)
    integer :: i, j, piv
    do i = 1, n
      piv = i - 1 + maxloc(abs(A(i:n, i)), 1)
      if (piv /= i) then
        row = A(i, :); A(i, :) = A(piv, :); A(piv, :) = row
        t = b(i); b(i) = b(piv); b(piv) = t
      end if
      do j = i + 1, n
        t = A(j, i) / A(i, i)
        A(j, i:n) = A(j, i:n) - t * A(i, i:n); b(j) = b(j) - t * b(i)
      end do
    end do
    do i = n, 1, -1
      b(i) = (b(i) - dot_product(A(i, i+1:n), b(i+1:n))) / A(i, i)
    end do
  end subroutine solve
end program growth_living_driver
