! homeostatic_driver.f90 -- one element of living_uel.f90 under the homeostatic-pressure growth law (27 Oct 2026), the
! check of p3/uel/homeostatic_test.py for this element, without Abaqus: phi = 1, c = 1, a = 0 at every node, no point
! model (phi_min = 2, so no Gauss point starts and nothing is consumed), so the growth is Eq. 36 with the gate,
! alpha_dot = k_alpha g.  Newton on the displacement dofs with the element's tangent (dense solve).
!   gfortran -cpp -ffixed-form -ffixed-line-length-132 -c living_uel_abaqus.for
!   gfortran -ffree-line-length-none homeostatic_driver.f90 living_uel_abaqus.o -o homeostatic_driver
!   ./homeostatic_driver <support> <p_ref> <form> <steps> <T> [k_alpha]
!     support: free (3-2-1 support: p = 0, Eq. 36) or all (every u = 0: p = E (1 + f) xi / (1 - 2 nu))
!     p_ref <= 0: the law off; form 0, 1, 2 explicit, 10, 11, 12 implicit (living_uel.f90 header)
! Prints one line per increment: k, t, alpha - 1 (Gauss-point mean), p = -tr(sigma)/3 (mean), Newton iterations.
program homeostatic_driver
  use living_uel_mod
  implicit none
  double precision :: U(NDOF), Un(NDOF), coords(3, NN), sv(NSV, NGP), props(18), R(NDOF), K(NDOF, NDOF), sig(6, NGP)
  double precision :: pref, T, dt, ka, p, rn, Kf(NDOF, NDOF), Rf(NDOF), dU(NDOF)
  integer :: steps, iform, i, kk, it, nf, idx(NDOF), q
  logical :: free(NDOF), ok
  character(len=16) :: arg
  double precision, parameter :: xr(8) = [0d0, 1d0, 1d0, 0d0, 0d0, 1d0, 1d0, 0d0]
  double precision, parameter :: xs(8) = [0d0, 0d0, 1d0, 1d0, 0d0, 0d0, 1d0, 1d0]
  double precision, parameter :: xt(8) = [0d0, 0d0, 0d0, 0d0, 1d0, 1d0, 1d0, 1d0]
  if (command_argument_count() < 5) stop 'usage: homeostatic_driver free|all p_ref form steps T [k_alpha]'
  call get_command_argument(1, arg)
  call get_command_argument(2, arg); read(arg, *) pref
  call get_command_argument(3, arg); read(arg, *) iform
  call get_command_argument(4, arg); read(arg, *) steps
  call get_command_argument(5, arg); read(arg, *) T
  ka = 1d-3
  if (command_argument_count() >= 6) then
    call get_command_argument(6, arg); read(arg, *) ka
  end if
  call get_command_argument(1, arg)
  eco_nact = 1
  gate_pref = pref; gate_form = iform; gate_implicit = .false.
  if (iform >= 10) then
    gate_implicit = .true.; gate_form = iform - 10
  end if
  ! E = 10 Pa, nu = 0.49, f = 1e-3 (Klempt et al. 2024 Table 2); no front term, no diffusion, no consumption
  props = [10d0, 0.49d0, 1d-3, 0d0, ka, 5d0, 0d0, 1d0, 0d0, 0d0, 2d0, 0.1d0, 1d0, 2d0, 0.5625d0, 0.55d0, 0.275d0, 0.275d0]
  do i = 1, NN
    coords(:, i) = [xr(i), xs(i), xt(i)]
  end do
  U = 0d0
  do i = 1, NN
    U(NDN*(i-1)+4) = 1d0; U(NDN*(i-1)+5) = 1d0
  end do
  sv = 0d0
  free = .false.
  if (arg == 'free') then
    do i = 1, NN
      free(NDN*(i-1)+1:NDN*(i-1)+3) = .true.
    end do
    free(1:3) = .false.                       ! node 1 (0,0,0): u1 = u2 = u3 = 0
    free(NDN+2:NDN+3) = .false.               ! node 2 (1,0,0): u2 = u3 = 0
    free(3*NDN+3) = .false.                   ! node 4 (0,1,0): u3 = 0
  end if
  nf = 0
  do i = 1, NDOF
    if (free(i)) then
      nf = nf + 1; idx(nf) = i
    end if
  end do
  dt = T / steps
  do kk = 1, steps
    Un = U
    call advance(Un, coords, sv, props, dt, ok, 1, 2)
    do it = 1, 30
      call element(U, Un, coords, sv, props, dt, R, K, .true., sig)
      rn = 0d0
      if (nf > 0) rn = maxval(abs(R(idx(1:nf))))
      if (rn < 1d-13 .or. nf == 0) exit
      Kf(1:nf, 1:nf) = K(idx(1:nf), idx(1:nf)); Rf(1:nf) = -R(idx(1:nf))
      call solve(nf, Kf(1:nf, 1:nf), Rf(1:nf))
      dU = 0d0; dU(idx(1:nf)) = Rf(1:nf)
      U = U + dU
    end do
    call alpha_store(U, Un, coords, sv, props, dt)
    p = 0d0
    do q = 1, NGP
      p = p - sum(sig(1:3, q)) / 3d0
    end do
    p = p / NGP
    print '(i6,es24.16,es24.16,es24.16,i4)', kk, kk * dt, sum(sv(13, :)) / NGP, p, it - 1
  end do
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
end program homeostatic_driver
