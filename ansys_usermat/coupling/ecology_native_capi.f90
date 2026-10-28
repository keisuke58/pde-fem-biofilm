! ecology_native_capi.f90 -- C entry points to the Fortran point model of the ANSYS runs (ecology_native.f,
! module biofilm_py_bridge), so that the cloud can run exactly the code ANSYS runs (ctypes, see
! ecology_native_ctypes.py). Nothing of ecology_native.f is changed; this file only calls it.
!   gfortran -O2 -fPIC -shared ecology_native.f ecology_native_capi.f90 -o libecology_native.so
module ecology_native_capi
  use, intrinsic :: iso_c_binding
  use biofilm_py_bridge
  implicit none
contains
  subroutine eco_capi_init(nact, cstar, alstar, nwmax, nwtol, theta) bind(c, name="eco_capi_init")
    ! reads the file named by BIOFILM_ECO_CASE (as in ANSYS) and returns its constants and theta
    integer(c_int), intent(out) :: nact, nwmax
    real(c_double), intent(out) :: cstar, alstar, nwtol, theta(20)
    if (.not. eco_inited) call eco_native_init()
    nact = eco_nact; cstar = eco_cstar; alstar = eco_alstar
    nwmax = eco_nwmax; nwtol = eco_nwtol; theta = eco_theta
  end subroutine eco_capi_init

  subroutine eco_capi_lin(n, g, theta, dt, n_sub, c0, c1, a0, a1, gnew, ok) bind(c, name="eco_capi_lin")
    ! n points, each one call of biofilm_ecology_hook_lin, the routine the ANSYS material routine calls
    integer(c_int), value, intent(in) :: n, n_sub
    real(c_double), value, intent(in) :: dt
    real(c_double), intent(in) :: g(12, n), theta(20), c0(n), c1(n), a0(n), a1(n)
    real(c_double), intent(out) :: gnew(12, n)
    integer(c_int), intent(out) :: ok(n)
    real(c_double) :: phi_int
    logical :: lok
    integer :: k
    do k = 1, n
      call biofilm_ecology_hook_lin(g(:, k), theta, dt, int(n_sub), c0(k), c1(k), a0(k), a1(k), &
                                    gnew(:, k), phi_int, lok)
      ok(k) = merge(1, 0, lok)
    end do
  end subroutine eco_capi_lin
end module ecology_native_capi
