C     ecology_native.f -- the point model (0D Hamilton ODE, Klempt et al.
C     2026; ecology_jax.ecology_substeps) in Fortran, as a drop-in for
C     usermat_py_hook.f: the same module name biofilm_py_bridge and the
C     same hooks biofilm_ecology_hook / biofilm_ecology_hook_c, so the
C     call-site fragments stay unchanged, but no material server, no
C     socket and no C shim (6 Oct 2026: the server call per integration
C     point was ~100x the rest of an Abaqus run).
C
C     Same scheme as jax_hamilton_0d_5species_demo.newton_step: per
C     sub-step 6 Newton iterations on the 12 residuals in
C     g = (phi_1..5, phi_0, psi_1..5, gamma), each preceded and followed
C     by clip_state; the Jacobian is the exact derivative (jax.jacfwd
C     there, written out here); LU with partial pivoting. Differences to
C     the server are round-off (tests/test_ecology_native.py: <= 1e-12).
C     The Hill gate is off (ecology_constants.K_HILL = 0); a case with
C     K_hill > 0 is not supported here.
C
C     Constants: as material_server.py --case (c*, alpha*, eta_i of Klempt
C     et al. 2026 Table 1, n species, theta checked against the deck) or,
C     without a case, ecology_constants.py (c* = 25, alpha* = 0, eta = 1,
C     5 species). Set once, single-threaded: eco_native_init reads the
C     file named by the environment variable BIOFILM_ECO_CASE
C     (abaqus_composition/write_eco_cfg.py; Abaqus calls it in
C     UEXTERNALDB, ANSYS at the first hook call), or eco_native_set sets
C     them directly.

      module biofilm_py_bridge
        use, intrinsic :: iso_c_binding
        implicit none
        integer, save :: eco_nact = 5
        logical, save :: eco_hascase = .false.
        double precision, save :: eco_cstar = 25.0d0
        double precision, save :: eco_alstar = 0.0d0
        double precision, save :: eco_eta(5) = 1.0d0
        double precision, save :: eco_theta(20) = 0.0d0
        double precision, parameter :: eco_kp1 = 1.0d-4
        logical, save :: eco_inited = .false.
      contains

        subroutine eco_native_set(nact, cstar, alstar, eta, theta,
     &                            hascase)
          integer, intent(in) :: nact
          double precision, intent(in) :: cstar, alstar, eta(5)
          double precision, intent(in) :: theta(20)
          logical, intent(in) :: hascase
          eco_inited = .true.
          eco_nact = nact
          eco_cstar = cstar
          eco_alstar = alstar
          eco_eta = eta
          eco_theta = theta
          eco_hascase = hascase
        end subroutine eco_native_set

        subroutine eco_native_init()
C         file: n_active / c* alpha* / eta(5) / theta(20), free format
          character(len=512) :: path
          integer :: st, u, nact
          double precision :: cs, als, eta(5), th(20)
          eco_inited = .true.
          call get_environment_variable('BIOFILM_ECO_CASE', path,
     &                                  status=st)
          if (st .ne. 0) return
          open(newunit=u, file=trim(path), status='old', iostat=st)
          if (st .ne. 0) return
          read(u, *, iostat=st) nact
          if (st .eq. 0) read(u, *, iostat=st) cs, als
          if (st .eq. 0) read(u, *, iostat=st) eta
          if (st .eq. 0) read(u, *, iostat=st) th
          close(u)
          if (st .eq. 0) call eco_native_set(nact, cs, als, eta, th,
     &                                       .true.)
        end subroutine eco_native_init

        subroutine eco_clip(g, mask)
          double precision, intent(inout) :: g(12)
          double precision, intent(in) :: mask(5)
          double precision, parameter :: e = 1.0d-10
          integer :: i
          do i = 1, 5
            g(i) = mask(i) * min(max(g(i), e), 1.0d0 - e)
            g(6+i) = mask(i) * min(max(g(6+i), e), 1.0d0 - e)
          end do
          g(6) = min(max(g(6), e), 1.0d0 - e)
          g(12) = min(max(g(12), -1.0d6), 1.0d6)
        end subroutine eco_clip

        subroutine eco_resjac(g, gp, dt, a, b, c, al, eta, mask, q, jac)
C         residual Q(g; g_prev) of the demo's residual() and dQ/dg
          double precision, intent(in) :: g(12), gp(12), dt, a(5,5)
          double precision, intent(in) :: b(5), c, al, eta(5), mask(5)
          double precision, intent(out) :: q(12), jac(12,12)
          double precision :: p(5), s(5), pd(5), sd(5), ia(5), w(5)
          double precision :: x, d, dd, t1, dt1, u1, du1, k
          integer :: i, j
          k = eco_kp1
          q = 0.0d0
          jac = 0.0d0
          do i = 1, 5
            p(i) = g(i)
            s(i) = g(6+i)
            pd(i) = (g(i) - gp(i)) / dt
            sd(i) = (g(6+i) - gp(6+i)) / dt
            w(i) = p(i) * s(i)
          end do
          do i = 1, 5
            ia(i) = 0.0d0
            do j = 1, 5
              ia(i) = ia(i) + a(i,j) * w(j)
            end do
          end do
          do i = 1, 5
            if (mask(i) .gt. 0.5d0) then
              x = p(i)
              d = (x - 1.0d0)**3 * x**3
              t1 = k * (2.0d0 - 4.0d0 * x) / d
              dd = 3.0d0 * (x - 1.0d0)**2 * x**2 * (2.0d0 * x - 1.0d0)
              dt1 = k * (-4.0d0 * d - (2.0d0 - 4.0d0 * x) * dd) / d**2
              q(i) = t1 + (1.0d0 / eta(i)) * (g(12)
     &             + (eta(i) + eta(i) * s(i)**2) * pd(i)
     &             + eta(i) * p(i) * s(i) * sd(i))
     &             - (c / eta(i)) * s(i) * ia(i)
              jac(i,i) = dt1 + (1.0d0 / eta(i)) * ((eta(i)
     &                 + eta(i) * s(i)**2) / dt + eta(i) * s(i) * sd(i))
              jac(i,6+i) = (1.0d0 / eta(i)) * (2.0d0 * eta(i) * s(i)
     &                   * pd(i) + eta(i) * p(i) * sd(i) + eta(i) * p(i)
     &                   * s(i) / dt) - (c / eta(i)) * ia(i)
              jac(i,12) = 1.0d0 / eta(i)
              do j = 1, 5
                jac(i,j) = jac(i,j) - (c / eta(i)) * s(i) * a(i,j)
     &                   * s(j)
                jac(i,6+j) = jac(i,6+j) - (c / eta(i)) * s(i)
     &                     * a(i,j) * p(j)
              end do
            else
              q(i) = p(i)
              jac(i,i) = 1.0d0
            end if
          end do
          x = g(6)
          d = (x - 1.0d0)**3 * x**3
          dd = 3.0d0 * (x - 1.0d0)**2 * x**2 * (2.0d0 * x - 1.0d0)
          q(6) = g(12) + k * (2.0d0 - 4.0d0 * x) / d
     &         + (x - gp(6)) / dt
          jac(6,6) = k * (-4.0d0 * d - (2.0d0 - 4.0d0 * x) * dd) / d**2
     &             + 1.0d0 / dt
          jac(6,12) = 1.0d0
          do i = 1, 5
            if (mask(i) .gt. 0.5d0) then
              x = s(i)
              u1 = (-2.0d0 * k) / ((x - 1.0d0)**2 * x**3)
     &           - (2.0d0 * k) / ((x - 1.0d0)**3 * x**2)
              du1 = -2.0d0 * k * (-2.0d0 / ((x - 1.0d0)**3 * x**3)
     &            - 3.0d0 / ((x - 1.0d0)**2 * x**4))
     &            - 2.0d0 * k * (-3.0d0 / ((x - 1.0d0)**4 * x**2)
     &            - 2.0d0 / ((x - 1.0d0)**3 * x**3))
              q(6+i) = u1 + (b(i) * al / eta(i)) * s(i)
     &               + p(i) * s(i) * pd(i) + p(i)**2 * sd(i)
     &               - (c / eta(i)) * p(i) * ia(i)
              jac(6+i,6+i) = du1 + b(i) * al / eta(i) + p(i) * pd(i)
     &                     + p(i)**2 / dt
              jac(6+i,i) = s(i) * pd(i) + p(i) * s(i) / dt
     &                   + 2.0d0 * p(i) * sd(i) - (c / eta(i)) * ia(i)
              do j = 1, 5
                jac(6+i,j) = jac(6+i,j) - (c / eta(i)) * p(i)
     &                     * a(i,j) * s(j)
                jac(6+i,6+j) = jac(6+i,6+j) - (c / eta(i)) * p(i)
     &                       * a(i,j) * p(j)
              end do
            else
              q(6+i) = s(i)
              jac(6+i,6+i) = 1.0d0
            end if
          end do
          q(12) = p(1) + p(2) + p(3) + p(4) + p(5) + g(6) - 1.0d0
          do j = 1, 6
            jac(12,j) = 1.0d0
          end do
        end subroutine eco_resjac

        subroutine eco_solve(a, r, x, ok)
C         x = a^{-1} r, Gaussian elimination with partial pivoting
          double precision, intent(inout) :: a(12,12), r(12)
          double precision, intent(out) :: x(12)
          logical, intent(out) :: ok
          double precision :: f, t
          integer :: i, j, kk, piv
          ok = .true.
          do kk = 1, 12
            piv = kk
            do i = kk + 1, 12
              if (abs(a(i,kk)) .gt. abs(a(piv,kk))) piv = i
            end do
            if (a(piv,kk) .eq. 0.0d0) then
              ok = .false.
              return
            end if
            if (piv .ne. kk) then
              do j = 1, 12
                t = a(kk,j)
                a(kk,j) = a(piv,j)
                a(piv,j) = t
              end do
              t = r(kk)
              r(kk) = r(piv)
              r(piv) = t
            end if
            do i = kk + 1, 12
              f = a(i,kk) / a(kk,kk)
              a(i,kk) = 0.0d0
              do j = kk + 1, 12
                a(i,j) = a(i,j) - f * a(kk,j)
              end do
              r(i) = r(i) - f * r(kk)
            end do
          end do
          do i = 12, 1, -1
            t = r(i)
            do j = i + 1, 12
              t = t - a(i,j) * x(j)
            end do
            x(i) = t / a(i,i)
          end do
        end subroutine eco_solve

        subroutine eco_substeps(g, theta, dt_h, n_sub, crel, g_new,
     &                          phi_int, ok)
C         ecology_jax.ecology_substeps: n_sub equal steps over dt_h;
C         phi_int = sum_k dt_sub sum_i phi_i psi_i (state after step k).
C         crel < 0: no local nutrient (c = c*).
          double precision, intent(in) :: g(12), theta(20), dt_h, crel
          integer, intent(in) :: n_sub
          double precision, intent(out) :: g_new(12), phi_int
          logical, intent(out) :: ok
          double precision :: a(5,5), b(5), mask(5), gp(12), gg(12)
          double precision :: q(12), jac(12,12), dl(12), dts, c, s
          integer :: i, it, ks, n
          ok = .false.
          g_new = g
          phi_int = 0.0d0
C         ANSYS has no UEXTERNALDB: read the constants at the first call
C         (single-threaded with -np 1, as the wired decks run)
          if (.not. eco_inited) call eco_native_init()
          if (n_sub .lt. 1) return
          if (eco_hascase) then
            do i = 1, 20
              if (theta(i) .ne. eco_theta(i)) return
            end do
          end if
          a = 0.0d0
          b = 0.0d0
          a(1,1) = theta(1)
          a(1,2) = theta(2)
          a(2,1) = theta(2)
          a(2,2) = theta(3)
          b(1) = theta(4)
          b(2) = theta(5)
          a(3,3) = theta(6)
          a(3,4) = theta(7)
          a(4,3) = theta(7)
          a(4,4) = theta(8)
          b(3) = theta(9)
          b(4) = theta(10)
          a(1,3) = theta(11)
          a(3,1) = theta(11)
          a(1,4) = theta(12)
          a(4,1) = theta(12)
          a(2,3) = theta(13)
          a(3,2) = theta(13)
          a(2,4) = theta(14)
          a(4,2) = theta(14)
          a(5,5) = theta(15)
          b(5) = theta(16)
          a(1,5) = theta(17)
          a(5,1) = theta(17)
          a(2,5) = theta(18)
          a(5,2) = theta(18)
          a(3,5) = theta(19)
          a(5,3) = theta(19)
          a(4,5) = theta(20)
          a(5,4) = theta(20)
          n = eco_nact
          mask = 0.0d0
          do i = 1, n
            mask(i) = 1.0d0
          end do
          c = eco_cstar
          if (crel .ge. 0.0d0) c = eco_cstar * crel
          gp = g
C         seed_inactive: zero the switched-off species once
          if (n .lt. 5) then
            s = 0.0d0
            do i = n + 1, 5
              s = s + abs(gp(i)) + abs(gp(6+i))
            end do
            if (s .ne. 0.0d0) then
              do i = n + 1, 5
                gp(i) = 0.0d0
                gp(6+i) = 0.0d0
              end do
              s = 0.0d0
              do i = 1, n
                s = s + gp(i)
              end do
              gp(6) = 1.0d0 - s
            end if
          end if
          dts = dt_h / dble(n_sub)
          do ks = 1, n_sub
            gg = gp
            call eco_clip(gg, mask)
            do it = 1, 6
              call eco_clip(gg, mask)
              call eco_resjac(gg, gp, dts, a, b, c, eco_alstar, eco_eta,
     &                        mask, q, jac)
              q = -q
              call eco_solve(jac, q, dl, ok)
              if (.not. ok) return
              gg = gg + dl
              call eco_clip(gg, mask)
            end do
            gp = gg
            s = 0.0d0
            do i = 1, 5
              s = s + gg(i) * gg(6+i)
            end do
            phi_int = phi_int + dts * s
          end do
          g_new = gp
          ok = .true.
        end subroutine eco_substeps

        subroutine biofilm_py_hook(F, Fv, alpha, C10, C01, D1, eta,
     &                             mtype, dt, stress, Fvnew, dsde, ok)
C         the material bridge is not part of the native build: ok =
C         .false. sends the caller to its inline BIOFILM_STRESS_CORE
          real(c_double), intent(in)  :: F(3,3), Fv(3,3)
          real(c_double), intent(in)  :: alpha, C10, C01, D1, eta
          real(c_double), intent(in)  :: mtype, dt
          real(c_double), intent(out) :: stress(6), Fvnew(9), dsde(6,6)
          logical, intent(out)        :: ok
          stress = 0.0d0
          Fvnew = 0.0d0
          dsde = 0.0d0
          ok = .false.
        end subroutine biofilm_py_hook

        subroutine biofilm_ecology_hook(g, theta, dt, n_sub, g_new,
     &                                   phi_int, ok)
          real(c_double), intent(in)  :: g(12), theta(20), dt
          integer, intent(in)         :: n_sub
          real(c_double), intent(out) :: g_new(12)
          real(c_double), intent(out) :: phi_int
          logical, intent(out)        :: ok
          call eco_substeps(g, theta, dt, n_sub, -1.0d0, g_new,
     &                      phi_int, ok)
        end subroutine biofilm_ecology_hook

        subroutine biofilm_ecology_hook_c(g, theta, dt, n_sub, c_rel,
     &                                     g_new, phi_int, ok)
          real(c_double), intent(in)  :: g(12), theta(20), dt, c_rel
          integer, intent(in)         :: n_sub
          real(c_double), intent(out) :: g_new(12)
          real(c_double), intent(out) :: phi_int
          logical, intent(out)        :: ok
          if (c_rel .lt. 0.0d0) then
            call eco_substeps(g, theta, dt, n_sub, -1.0d0, g_new,
     &                        phi_int, ok)
          else
            call eco_substeps(g, theta, dt, n_sub, c_rel, g_new,
     &                        phi_int, ok)
          end if
        end subroutine biofilm_ecology_hook_c
      end module biofilm_py_bridge
