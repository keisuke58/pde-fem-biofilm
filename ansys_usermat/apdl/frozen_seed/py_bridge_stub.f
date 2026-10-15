C  py_bridge_stub.f -- request A-16 test build only: the module
C  biofilm_py_bridge with the two hooks usermat_biofilm.f uses, both
C  returning ok = .false. (never reached: prop(6) = prop(8) = 0 in the
C  wrapper), so the build needs neither the C shim nor a material server.
      module biofilm_py_bridge
      contains
        subroutine biofilm_py_hook(F, Fv, alpha, C10, C01, D1, eta, mtype,
     &                             dt, stress, Fvnew, dsde, ok)
          double precision, intent(in)  :: F(3,3), Fv(3,3)
          double precision, intent(in)  :: alpha, C10, C01, D1, eta,
     &                                     mtype, dt
          double precision, intent(out) :: stress(6), Fvnew(9),
     &                                     dsde(6,6)
          logical, intent(out)          :: ok
          stress = 0.0d0
          Fvnew = 0.0d0
          dsde = 0.0d0
          ok = .false.
        end subroutine biofilm_py_hook
        subroutine biofilm_ecology_hook(g, theta, dt, n_sub, g_new,
     &                                   phi_int, ok)
          double precision, intent(in)  :: g(12), theta(20), dt
          integer, intent(in)           :: n_sub
          double precision, intent(out) :: g_new(12), phi_int
          logical, intent(out)          :: ok
          g_new = g
          phi_int = 0.0d0
          ok = .false.
        end subroutine biofilm_ecology_hook
      end module biofilm_py_bridge
