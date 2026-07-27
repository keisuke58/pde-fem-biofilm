C     usermat_py_hook.f  —  SKELETON of the Fortran side of the Gauss-point
C     Python material bridge (the `PYTHON MATERIAL HOOK` in usermat_biofilm.f).
C
C     Two supported mechanisms (choose one at integration time):
C       (A) ISO_C_BINDING  — link a tiny C shim `biofilm_py_eval` that talks to
C           the Python material server (in-process or over a local socket).
C           This module declares that interface and does the array marshalling.
C       (B) pure-socket    — the C shim opens a TCP client to material_server.py
C           on 127.0.0.1:8765 and exchanges one newline-delimited JSON frame per
C           call (see coupling/protocol.py).
C
C     Call site in the USERMAT (pseudocode):
C         if (kUsePy .gt. 0.5d0) then
C             call biofilm_py_hook(F, Fv, alpha, C10, C01, D1, eta, mtype, dt,
C        &                         stress, Fvnew, dsde, ok)
C             if (.not. ok) <fall back to the inline BIOFILM_STRESS_CORE>
C         else
C             <inline BIOFILM_STRESS_CORE>            ! verified reference/fallback
C         end if
C
C     Voigt order: Abaqus 11,22,33,12,13,23 (matches material_server.py).
C     Syntax check:  gfortran -c -fsyntax-only -ffixed-line-length-132 usermat_py_hook.f

      module biofilm_py_bridge
        use, intrinsic :: iso_c_binding
        implicit none
        interface
C         C shim (to be provided): returns 0 on success, nonzero on failure.
C         int biofilm_py_eval(const double* F9, const double* Fv9,
C                             const double* params7, double* stress6,
C                             double* Fvnew9, double* dsde36);
          function biofilm_py_eval(F9, Fv9, params7, stress6, Fvnew9, dsde36)
     &             bind(C, name="biofilm_py_eval") result(ierr)
            import :: c_int, c_double
            real(c_double), intent(in)  :: F9(9), Fv9(9), params7(7)
            real(c_double), intent(out) :: stress6(6), Fvnew9(9), dsde36(36)
            integer(c_int) :: ierr
          end function biofilm_py_eval
        end interface
      contains
        subroutine biofilm_py_hook(F, Fv, alpha, C10, C01, D1, eta, mtype,
     &                             dt, stress, Fvnew, dsde, ok)
          real(c_double), intent(in)  :: F(3,3), Fv(3,3)
          real(c_double), intent(in)  :: alpha, C10, C01, D1, eta, mtype, dt
          real(c_double), intent(out) :: stress(6), Fvnew(9), dsde(6,6)
          logical, intent(out)        :: ok
          real(c_double) :: F9(9), Fv9(9), params7(7), d36(36)
          integer(c_int) :: ierr
          integer :: i, j, k
          k = 0
          do i = 1, 3
            do j = 1, 3
              k = k + 1
              F9(k)  = F(i,j)
              Fv9(k) = Fv(i,j)
            end do
          end do
          params7(1) = alpha
          params7(2) = C10
          params7(3) = C01
          params7(4) = D1
          params7(5) = eta
          params7(6) = mtype
          params7(7) = dt
          ierr  = biofilm_py_eval(F9, Fv9, params7, stress, Fvnew, d36)
          ok    = (ierr .eq. 0)
          dsde  = reshape(d36, [6, 6])
        end subroutine biofilm_py_hook
      end module biofilm_py_bridge
