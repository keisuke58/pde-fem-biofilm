! test_tangent.f90 -- the closed-form tangent of living_uel.f90 against central differences of its residual, for one
! distorted element with random fields, started Gauss points and phi inside (0, 1) (A-46).
!   gfortran -cpp -ffixed-form -ffixed-line-length-132 -c living_uel_abaqus.for
!   gfortran -ffree-line-length-none test_tangent.f90 living_uel_abaqus.o -o test_tangent && ./test_tangent [supg]
program test_tangent
  use living_uel_mod
  implicit none
  double precision :: U(NDOF), Un(NDOF), coords(3, NN), sv(NSV, NGP), props(18), R(NDOF), K(NDOF, NDOF)
  double precision :: Rp(NDOF), Rm(NDOF), Kfd(NDOF, NDOF), Kd(NDOF, NDOF), sig(6, NGP), h, err, sc
  integer :: i, j
  double precision, parameter :: xr(8) = [-1d0, 1d0, 1d0, -1d0, -1d0, 1d0, 1d0, -1d0]
  double precision, parameter :: xs(8) = [-1d0, -1d0, 1d0, 1d0, -1d0, -1d0, 1d0, 1d0]
  double precision, parameter :: xt(8) = [-1d0, -1d0, -1d0, -1d0, 1d0, 1d0, 1d0, 1d0]
  eco_nact = 2
  call random_seed()
  do i = 1, NN
    call random_number(coords(:, i))
    coords(:, i) = 0.125d0 * [xr(i), xs(i), xt(i)] + 0.02d0 * (coords(:, i) - 0.5d0)
  end do
  call random_number(U); call random_number(Un)
  do i = 1, NN
    U(NDN*(i-1)+1:NDN*(i-1)+3) = 1d-4 * (U(NDN*(i-1)+1:NDN*(i-1)+3) - 0.5d0)
    U(NDN*(i-1)+4) = 0.2d0 + 0.6d0 * U(NDN*(i-1)+4)
    Un(NDN*(i-1)+4) = 0.2d0 + 0.6d0 * Un(NDN*(i-1)+4)
    U(NDN*(i-1)+5) = 0.8d0 + 0.2d0 * U(NDN*(i-1)+5); Un(NDN*(i-1)+5) = 0.8d0 + 0.2d0 * Un(NDN*(i-1)+5)
  end do
  call random_number(sv)
  sv(14, :) = 1d0
  sv(13, :) = 1d-4 * sv(13, :)
  props = [10d0, 0.49d0, 1d-3, 0.645d0, 1d-3, 5d0, 10d0, 1d0, 1d0, 1d0, 1d-3, 0.1d0, 5d0, 1d0, 0.5625d0, 0.55d0, &
           0.275d0, 0.275d0]
  supg_on = (command_argument_count() > 0)           ! any argument: with SUPG (A-54)
  call element(U, Un, coords, sv, props, 0.005d0, R, K, .true., sig)
  do j = 1, NDOF
    h = 1d-6 * max(1d0, abs(U(j)))
    if (mod(j - 1, NDN) < 3) h = 1d-9
    U(j) = U(j) + h; call element(U, Un, coords, sv, props, 0.005d0, Rp, Kd, .false., sig)
    U(j) = U(j) - 2 * h; call element(U, Un, coords, sv, props, 0.005d0, Rm, Kd, .false., sig)
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
  if (err > 1d-6) stop 1
end program test_tangent
