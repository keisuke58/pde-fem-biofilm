! test_uel_wrapper.f90 -- calls the Abaqus entry point UEL twice from the same start-of-increment
! state (as Abaqus does in two iterations) and checks: identical RHS, AMATRX and end-of-increment
! SVARS both times; RHS = -R of the module's residual with the advanced composition; growth written.
program test_uel_wrapper
  use unified_uel_mod
  use iso_c_binding
  implicit none
  interface
    subroutine uel_set_point_model(theta, cstar, alstar, nwmax, nwtol) bind(C, name="uel_set_point_model")
      import :: c_double, c_int
      real(c_double), intent(in) :: theta(20), cstar, alstar, nwtol
      integer(c_int), value :: nwmax
    end subroutine
    subroutine uel_set_phi_init(v) bind(C, name="uel_set_phi_init")
      import :: c_double
      real(c_double), intent(in) :: v(5)
    end subroutine
  end interface
  integer, parameter :: NDOFEL = 40, NSVARS = 112, NPROPS = 16
  double precision :: RHS(NDOFEL,1), AM(NDOFEL,NDOFEL), SV0(NSVARS), SV1(NSVARS), SV2(NSVARS), ENERGY(8), &
       PROPS(NPROPS), COORDS(3,8), U(NDOFEL), DU(NDOFEL,1), V(NDOFEL), A(NDOFEL), TIME(2), DTIME, PARAMS(3), &
       ADLMAG(1,1), PREDEF(2,1,8), DDLMAG(1,1), PNEWDT, PERIOD, RHS1(NDOFEL,1), AM1(NDOFEL,NDOFEL), &
       th(20), R(NDOF), sv(NSV,NGP), Un(NDOF), g(5)
  integer :: JDLTYP(1,1), LFLAGS(5), JPROPS(1), ia, q
  logical :: okc
  double precision, parameter :: cs = 25d0, al = 0d0, tol = 1d-26
  th = [ -0.6929614183633133d0, -0.3903645866494718d0, -0.4610407137003207d0, 0d0, 0d0, -0.7769305773355909d0, &
         -0.9298104492087518d0, -0.8590640026947014d0, 0d0, 0d0, 0.5189105952206899d0, -4.449056646982786d0, &
          1.3483384494191148d0, 0.8332096197960628d0, -0.9920945055134962d0, 0d0, -4.545884027121279d0, &
          0.5733977664159978d0, -0.5955397035054063d0, 0.14858871650880234d0 ]
  call uel_set_point_model(th, cs, al, 60, tol)
  g = [0.75d0, 0.05d0, 0.06d0, 0.03d0, 0.01d0]
  call uel_set_phi_init(g)
  COORDS = reshape([0d0,0d0,0d0, 1d0,0d0,0d0, 1d0,1d0,0d0, 0d0,1d0,0d0, &
                    0d0,0d0,1d0, 1d0,0d0,1d0, 1d0,1d0,1d0, 0d0,1d0,1d0], [3,8]) * 0.1d0
  PROPS = [1d0, 0.3d0, 1d-3, 0.02d0, 1d0, 1d-5, 1d0, 0d0, 0.7d0, 1d0, log(0.1d0), 1d0, 1d0, 1d0, 0.05d0, 0d0]
  U = 0d0
  do ia = 1, 8
    U(5*(ia-1)+4) = 0.6d0 + 0.05d0*ia/8d0
    U(5*(ia-1)+5) = 0.9d0
  end do
  DU = 0d0; DU(4,1) = 0.01d0
  SV0 = 0d0
  DTIME = 1d-3; LFLAGS = [1, 1, 1, 0, 0]; PNEWDT = 1d0
  SV1 = SV0
  call UEL(RHS, AM, SV1, ENERGY, NDOFEL, 1, NSVARS, PROPS, NPROPS, COORDS, 3, 8, U, DU, V, A, 1, TIME, DTIME, &
           1, 1, 1, PARAMS, 0, JDLTYP, ADLMAG, PREDEF, 1, LFLAGS, NDOFEL, DDLMAG, 1, PNEWDT, JPROPS, 0, PERIOD)
  RHS1 = RHS; AM1 = AM
  SV2 = SV0
  call UEL(RHS, AM, SV2, ENERGY, NDOFEL, 1, NSVARS, PROPS, NPROPS, COORDS, 3, 8, U, DU, V, A, 1, TIME, DTIME, &
           1, 1, 1, PARAMS, 0, JDLTYP, ADLMAG, PREDEF, 1, LFLAGS, NDOFEL, DDLMAG, 1, PNEWDT, JPROPS, 0, PERIOD)
  print '(a,es10.2,a,es10.2,a,es10.2)', 'second call vs first: RHS ', maxval(abs(RHS-RHS1)), ', AMATRX ', &
        maxval(abs(AM-AM1)), ', SVARS ', maxval(abs(SV2-SV1))
  ! residual of the module with the composition advanced from SV0
  sv = reshape(SV0, [NSV, NGP]); Un = U - DU(:,1)
  call advance_composition(U, Un, COORDS, sv, PROPS, DTIME, okc)
  call residual(U, Un, COORDS, sv, PROPS, DTIME, R)
  print '(a,es10.2)', 'RHS + R (module residual) ', maxval(abs(RHS(:,1) + R))
  q = 0
  do ia = 1, NGP
    if (SV1(14*(ia-1)+14) > 0.5d0) q = q + 1
  end do
  print '(a,i2,a,es10.3,a,f8.4)', 'Gauss points started ', q, ', growth xi at Gauss point 1 ', SV1(13), &
        ', sum of shares ', sum(SV1(1:6))
end program test_uel_wrapper
