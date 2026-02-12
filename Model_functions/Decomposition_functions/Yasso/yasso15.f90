! =============================================================================
! Yasso15 Fortran Implementation for R Interface
! =============================================================================
!
! Based on:
!   Järvenpää (2015) Yasso15 core code, Finnish Meteorological Institute
!   GPL-3.0 licence
!
! Key differences from Yasso07:
!   - Pool-group-specific climate responses:
!       AWE pools: beta1, beta2, gamma        (params 22-24 in original, 17-19 here)
!       N pool:    betaN1, betaN2, gammaN      (params 25-27 here)
!       H pool:    betaH1, betaH2, gammaH      (params 28-30 here)
!   - Leaching term on AWEN diagonal (requires leac + precip)
!   - 35-parameter vector (original indexing preserved in comments)
!
! Architecture mirrors Yasso07:
!   - xi computation in R, passed as pre-computed arrays (xi_awe, xi_n, xi_h)
!   - Two R entry points: yasso15_steady_state_r, yasso15_run_r
!
! Compile with:
!   R CMD SHLIB yasso15.f90
!
! =============================================================================

MODULE yasso15_mod
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 307)
  INTEGER, PARAMETER :: NPOOLS = 5
  REAL(dp), PARAMETER :: TOL = 1.0E-12_dp
CONTAINS

  ! =========================================================================
  ! BUILD DECOMPOSITION MATRIX A
  ! =========================================================================

  SUBROUTINE build_matrix_A(params, xi_awe, xi_n, xi_h, woody_diam, leac, precip, A)
    ! Constructs the 5x5 decomposition matrix.
    !
    ! Parameter vector (35 elements, standardised):
    !   1-4:   alpha_A, alpha_W, alpha_E, alpha_N   (base decomp rates)
    !   5-16:  p_WA, p_EA, p_NA, p_AW, p_EW, p_NW,
    !          p_AE, p_WE, p_NE, p_AN, p_WN, p_EN  (transfer fractions, same order as Yasso07)
    !   17-19: beta1, beta2, gamma                   (AWE climate response)
    !   20-22: betaN1, betaN2, gammaN                (N climate response)
    !   23-25: betaH1, betaH2, gammaH                (H climate response)
    !   26-27: p_H, alpha_H                          (humus formation & decomp)
    !   28-30: delta1, delta2, r                     (size modifier)
    !   31-35: w1, w2, w3, w4, w5                   (leaching; w-params, stored but leac scalar used)
    !
    ! xi_awe, xi_n, xi_h: pre-computed climate modifiers (from R)
    ! leac:    leaching scalar (site-level, passed from R)
    ! precip:  annual precipitation (mm), needed for leaching term

    REAL(dp), INTENT(IN)  :: params(35)
    REAL(dp), INTENT(IN)  :: xi_awe, xi_n, xi_h
    REAL(dp), INTENT(IN)  :: woody_diam, leac, precip
    REAL(dp), INTENT(OUT) :: A(NPOOLS, NPOOLS)

    REAL(dp) :: alpha(NPOOLS), k(NPOOLS), size_mod
    REAL(dp) :: p12, p13, p14
    REAL(dp) :: p21, p23, p24
    REAL(dp) :: p31, p32, p34
    REAL(dp) :: p41, p42, p43
    REAL(dp) :: p_H
    INTEGER  :: i

    ! Base decomposition rates
    alpha(1) = params(1)   ! A
    alpha(2) = params(2)   ! W
    alpha(3) = params(3)   ! E
    alpha(4) = params(4)   ! N
    alpha(5) = params(27)  ! alpha_H

    ! Transfer fractions (same positional order as Yasso07 params 5-16)
    p12 = params(5)    ! W -> A
    p13 = params(6)    ! E -> A
    p14 = params(7)    ! N -> A
    p21 = params(8)    ! A -> W
    p23 = params(9)    ! E -> W
    p24 = params(10)   ! N -> W
    p31 = params(11)   ! A -> E
    p32 = params(12)   ! W -> E
    p34 = params(13)   ! N -> E
    p41 = params(14)   ! A -> N
    p42 = params(15)   ! W -> N
    p43 = params(16)   ! E -> N
    p_H = params(26)

    ! Size modifier: only affects AWEN pools
    IF (woody_diam > 0.0_dp) THEN
      size_mod = MIN(1.0_dp, (1.0_dp + params(28) * woody_diam + &
                               params(29) * woody_diam**2)**(-ABS(params(30))))
    ELSE
      size_mod = 1.0_dp
    END IF

    ! Effective decomposition rates (pool-group-specific xi)
    k(1) = -ABS(alpha(1)) * xi_awe * size_mod   ! A
    k(2) = -ABS(alpha(2)) * xi_awe * size_mod   ! W
    k(3) = -ABS(alpha(3)) * xi_awe * size_mod   ! E
    k(4) = -ABS(alpha(4)) * xi_n   * size_mod   ! N
    k(5) = -ABS(alpha(5)) * xi_h                ! H (no size effect)

    ! Build A: diagonal = pool loss rates
    A = 0.0_dp
    DO i = 1, NPOOLS
      A(i,i) = k(i)
    END DO

    ! Off-diagonal: transfers between AWEN pools
    A(1,2) = p12 * ABS(k(2))   ! W -> A
    A(1,3) = p13 * ABS(k(3))   ! E -> A
    A(1,4) = p14 * ABS(k(4))   ! N -> A

    A(2,1) = p21 * ABS(k(1))   ! A -> W
    A(2,3) = p23 * ABS(k(3))   ! E -> W
    A(2,4) = p24 * ABS(k(4))   ! N -> W

    A(3,1) = p31 * ABS(k(1))   ! A -> E
    A(3,2) = p32 * ABS(k(2))   ! W -> E
    A(3,4) = p34 * ABS(k(4))   ! N -> E

    A(4,1) = p41 * ABS(k(1))   ! A -> N
    A(4,2) = p42 * ABS(k(2))   ! W -> N
    A(4,3) = p43 * ABS(k(3))   ! E -> N

    ! Humus formation: AWEN -> H (size effect present via k)
    A(5,1) = p_H * ABS(k(1))
    A(5,2) = p_H * ABS(k(2))
    A(5,3) = p_H * ABS(k(3))
    A(5,4) = p_H * ABS(k(4))

    ! No flows from H back to AWEN
    A(1,5) = 0.0_dp
    A(2,5) = 0.0_dp
    A(3,5) = 0.0_dp
    A(4,5) = 0.0_dp

    ! Leaching: additive term on AWEN diagonal only (no leaching for H)
    DO i = 1, 4
      A(i,i) = A(i,i) + leac * precip / 1000.0_dp
    END DO

  END SUBROUTINE build_matrix_A


  ! =========================================================================
  ! SINGLE TIMESTEP: C(t+dt) from C(t)
  ! =========================================================================

  SUBROUTINE model_step(A, C_prev, litter_in, dt, C_next)
    REAL(dp), INTENT(IN)  :: A(NPOOLS, NPOOLS)
    REAL(dp), INTENT(IN)  :: C_prev(NPOOLS), litter_in(NPOOLS), dt
    REAL(dp), INTENT(OUT) :: C_next(NPOOLS)

    REAL(dp) :: At(NPOOLS, NPOOLS), mexpAt(NPOOLS, NPOOLS)
    REAL(dp) :: z1(NPOOLS), z2(NPOOLS)

    IF (ABS(A(1,1)) < TOL .AND. ABS(A(2,2)) < TOL) THEN
      C_next = C_prev + litter_in * dt
      RETURN
    END IF

    z1 = MATMUL(A, C_prev) + litter_in
    At = A * dt
    CALL matrixexp(At, mexpAt)
    z2 = MATMUL(mexpAt, z1) - litter_in
    CALL solve5(A, z2, C_next)

  END SUBROUTINE model_step


  ! =========================================================================
  ! STEADY STATE: C* = -A^-1 * b
  ! =========================================================================

  SUBROUTINE model_steady_state(A, litter_in, C_ss)
    REAL(dp), INTENT(IN)  :: A(NPOOLS, NPOOLS), litter_in(NPOOLS)
    REAL(dp), INTENT(OUT) :: C_ss(NPOOLS)

    REAL(dp) :: neg_A(NPOOLS, NPOOLS)

    neg_A = -A
    CALL solve5(neg_A, litter_in, C_ss)
    C_ss = MAX(0.0_dp, C_ss)

  END SUBROUTINE model_steady_state


  ! =========================================================================
  ! STEADY STATE — three litter types, called from R entry point
  ! =========================================================================

  SUBROUTINE yasso15_steady_state(params, &
                                   nwl_mean, fwl_mean, cwl_mean, &
                                   xi_awe_mean, xi_n_mean, xi_h_mean, &
                                   leac, precip_mean, &
                                   diam_fwl, diam_cwl, &
                                   C_init)
    ! C_init(15): [C_nwl(1:5) | C_fwl(6:10) | C_cwl(11:15)]

    REAL(dp), INTENT(IN)  :: params(35)
    REAL(dp), INTENT(IN)  :: nwl_mean(4), fwl_mean(4), cwl_mean(4)
    REAL(dp), INTENT(IN)  :: xi_awe_mean, xi_n_mean, xi_h_mean
    REAL(dp), INTENT(IN)  :: leac, precip_mean
    REAL(dp), INTENT(IN)  :: diam_fwl, diam_cwl
    REAL(dp), INTENT(OUT) :: C_init(3 * NPOOLS)

    REAL(dp) :: A(NPOOLS, NPOOLS), litter(NPOOLS), C_ss(NPOOLS)

    ! --- Non-woody litter ---
    litter(1:4) = nwl_mean
    litter(5)   = 0.0_dp
    CALL build_matrix_A(params, xi_awe_mean, xi_n_mean, xi_h_mean, &
                         0.0_dp, leac, precip_mean, A)
    CALL model_steady_state(A, litter, C_ss)
    C_init(1:5) = C_ss

    ! --- Fine woody litter ---
    litter(1:4) = fwl_mean
    litter(5)   = 0.0_dp
    CALL build_matrix_A(params, xi_awe_mean, xi_n_mean, xi_h_mean, &
                         diam_fwl, leac, precip_mean, A)
    CALL model_steady_state(A, litter, C_ss)
    C_init(6:10) = C_ss

    ! --- Coarse woody litter ---
    litter(1:4) = cwl_mean
    litter(5)   = 0.0_dp
    CALL build_matrix_A(params, xi_awe_mean, xi_n_mean, xi_h_mean, &
                         diam_cwl, leac, precip_mean, A)
    CALL model_steady_state(A, litter, C_ss)
    C_init(11:15) = C_ss

  END SUBROUTINE yasso15_steady_state


  ! =========================================================================
  ! TRANSIENT FORWARD SIMULATION — called from R entry point
  ! =========================================================================

  SUBROUTINE yasso15_run(n_years, params, &
                          nwl_awen, fwl_awen, cwl_awen, &
                          xi_awe_array, xi_n_array, xi_h_array, &
                          leac, precip_array, &
                          diam_fwl, diam_cwl, &
                          C_init, &
                          C_out, resp_out)
    ! Receives pre-computed xi arrays (one per pool group per year) and
    ! precip_array (needed for leaching term). No climate logic here.

    INTEGER,  INTENT(IN)    :: n_years
    REAL(dp), INTENT(IN)    :: params(35)
    REAL(dp), INTENT(IN)    :: nwl_awen(n_years, 4)
    REAL(dp), INTENT(IN)    :: fwl_awen(n_years, 4)
    REAL(dp), INTENT(IN)    :: cwl_awen(n_years, 4)
    REAL(dp), INTENT(IN)    :: xi_awe_array(n_years)
    REAL(dp), INTENT(IN)    :: xi_n_array(n_years)
    REAL(dp), INTENT(IN)    :: xi_h_array(n_years)
    REAL(dp), INTENT(IN)    :: leac               ! leaching scalar (site-level)
    REAL(dp), INTENT(IN)    :: precip_array(n_years)
    REAL(dp), INTENT(IN)    :: diam_fwl, diam_cwl
    REAL(dp), INTENT(IN)    :: C_init(3 * NPOOLS)
    REAL(dp), INTENT(INOUT) :: C_out(n_years, NPOOLS)
    REAL(dp), INTENT(INOUT) :: resp_out(n_years)

    REAL(dp) :: dt
    REAL(dp) :: A_nwl(NPOOLS,NPOOLS), A_fwl(NPOOLS,NPOOLS), A_cwl(NPOOLS,NPOOLS)
    REAL(dp) :: C_nwl(NPOOLS), C_fwl(NPOOLS), C_cwl(NPOOLS)
    REAL(dp) :: C_nwl_next(NPOOLS), C_fwl_next(NPOOLS), C_cwl_next(NPOOLS)
    REAL(dp) :: litter_nwl(NPOOLS), litter_fwl(NPOOLS), litter_cwl(NPOOLS)
    REAL(dp) :: C_total_prev, C_total_next, input_total
    REAL(dp) :: xi_awe, xi_n, xi_h, precip
    INTEGER  :: yr

    dt = 1.0_dp

    ! Unpack initial states
    C_nwl = C_init(1:5)
    C_fwl = C_init(6:10)
    C_cwl = C_init(11:15)

    DO yr = 1, n_years

      xi_awe = xi_awe_array(yr)
      xi_n   = xi_n_array(yr)
      xi_h   = xi_h_array(yr)
      precip = precip_array(yr)

      litter_nwl(1:4) = nwl_awen(yr, :)
      litter_nwl(5)   = 0.0_dp
      litter_fwl(1:4) = fwl_awen(yr, :)
      litter_fwl(5)   = 0.0_dp
      litter_cwl(1:4) = cwl_awen(yr, :)
      litter_cwl(5)   = 0.0_dp

      C_total_prev = SUM(C_nwl) + SUM(C_fwl) + SUM(C_cwl)
      input_total  = SUM(litter_nwl) + SUM(litter_fwl) + SUM(litter_cwl)

      CALL build_matrix_A(params, xi_awe, xi_n, xi_h, 0.0_dp,   leac, precip, A_nwl)
      CALL build_matrix_A(params, xi_awe, xi_n, xi_h, diam_fwl, leac, precip, A_fwl)
      CALL build_matrix_A(params, xi_awe, xi_n, xi_h, diam_cwl, leac, precip, A_cwl)

      CALL model_step(A_nwl, C_nwl, litter_nwl, dt, C_nwl_next)
      CALL model_step(A_fwl, C_fwl, litter_fwl, dt, C_fwl_next)
      CALL model_step(A_cwl, C_cwl, litter_cwl, dt, C_cwl_next)

      C_nwl = C_nwl_next
      C_fwl = C_fwl_next
      C_cwl = C_cwl_next

      C_out(yr, :) = C_nwl + C_fwl + C_cwl

      C_total_next = SUM(C_out(yr, :))
      resp_out(yr) = C_total_prev + input_total - C_total_next

    END DO

  END SUBROUTINE yasso15_run


  ! =========================================================================
  ! MATRIX EXPONENTIAL (Taylor series with scaling & squaring)
  ! =========================================================================

  SUBROUTINE matrixexp(A, B)
    REAL(dp), INTENT(IN)  :: A(NPOOLS, NPOOLS)
    REAL(dp), INTENT(OUT) :: B(NPOOLS, NPOOLS)

    REAL(dp) :: C(NPOOLS, NPOOLS), D(NPOOLS, NPOOLS)
    REAL(dp) :: p, normiter
    INTEGER  :: i, j, q

    q = 10

    B = 0.0_dp
    DO i = 1, NPOOLS
      B(i,i) = 1.0_dp
    END DO

    normiter = 2.0_dp
    j = 1
    CALL matrixnorm(A, p)
    DO
      IF (p < normiter) EXIT
      normiter = normiter * 2.0_dp
      j = j + 1
    END DO

    C = A / normiter
    B = B + C
    D = C

    DO i = 2, q
      D = MATMUL(C, D) / REAL(i, dp)
      B = B + D
    END DO

    DO i = 1, j
      B = MATMUL(B, B)
    END DO

  END SUBROUTINE matrixexp


  ! =========================================================================
  ! FROBENIUS NORM
  ! =========================================================================

  SUBROUTINE matrixnorm(A, nrm)
    REAL(dp), INTENT(IN)  :: A(NPOOLS, NPOOLS)
    REAL(dp), INTENT(OUT) :: nrm

    INTEGER :: i

    nrm = 0.0_dp
    DO i = 1, NPOOLS
      nrm = nrm + SUM(A(:,i)**2)
    END DO
    nrm = SQRT(nrm)

  END SUBROUTINE matrixnorm


  ! =========================================================================
  ! 5x5 LINEAR SOLVER (Gaussian elimination with partial pivoting)
  ! =========================================================================

  SUBROUTINE solve5(A, b, x)
    REAL(dp), INTENT(IN)  :: A(NPOOLS, NPOOLS), b(NPOOLS)
    REAL(dp), INTENT(OUT) :: x(NPOOLS)

    REAL(dp) :: U(NPOOLS, NPOOLS), c(NPOOLS)
    INTEGER  :: i

    CALL pgauss5(A, b, U, c)

    x(NPOOLS) = c(NPOOLS) / U(NPOOLS, NPOOLS)
    DO i = NPOOLS - 1, 1, -1
      x(i) = (c(i) - DOT_PRODUCT(U(i, i+1:NPOOLS), x(i+1:NPOOLS))) / U(i,i)
    END DO

  END SUBROUTINE solve5


  SUBROUTINE pgauss5(A, b, U, c)
    REAL(dp), INTENT(IN)  :: A(NPOOLS, NPOOLS), b(NPOOLS)
    REAL(dp), INTENT(OUT) :: U(NPOOLS, NPOOLS), c(NPOOLS)

    INTEGER  :: k, j, pk, q
    REAL(dp) :: tmp_row(NPOOLS), tmp_val

    U = A
    c = b

    DO k = 1, NPOOLS - 1
      q = MAXLOC(ABS(U(k:NPOOLS, k)), 1)
      IF (q > 1) THEN
        pk = k - 1 + q
        tmp_row = U(k, :)
        U(k, :) = U(pk, :)
        U(pk, :) = tmp_row
        tmp_val = c(k)
        c(k) = c(pk)
        c(pk) = tmp_val
      END IF

      IF (ABS(U(k,k)) <= TOL) CYCLE

      U(k+1:NPOOLS, k) = U(k+1:NPOOLS, k) / U(k, k)
      DO j = k + 1, NPOOLS
        U(j, k+1:NPOOLS) = U(j, k+1:NPOOLS) - U(j, k) * U(k, k+1:NPOOLS)
      END DO
      c(k+1:NPOOLS) = c(k+1:NPOOLS) - c(k) * U(k+1:NPOOLS, k)
    END DO

  END SUBROUTINE pgauss5

END MODULE yasso15_mod


! =============================================================================
! R-CALLABLE ENTRY POINTS
! =============================================================================

SUBROUTINE yasso15_steady_state_r(params, &
                                   nwl_mean, fwl_mean, cwl_mean, &
                                   xi_awe_mean, xi_n_mean, xi_h_mean, &
                                   leac, precip_mean, &
                                   diam_fwl, diam_cwl, &
                                   C_init)
  USE yasso15_mod
  IMPLICIT NONE

  REAL(dp), INTENT(IN)  :: params(35)
  REAL(dp), INTENT(IN)  :: nwl_mean(4), fwl_mean(4), cwl_mean(4)
  REAL(dp), INTENT(IN)  :: xi_awe_mean, xi_n_mean, xi_h_mean
  REAL(dp), INTENT(IN)  :: leac, precip_mean
  REAL(dp), INTENT(IN)  :: diam_fwl, diam_cwl
  REAL(dp), INTENT(OUT) :: C_init(15)

  CALL yasso15_steady_state(params, &
                             nwl_mean, fwl_mean, cwl_mean, &
                             xi_awe_mean, xi_n_mean, xi_h_mean, &
                             leac, precip_mean, &
                             diam_fwl, diam_cwl, &
                             C_init)

END SUBROUTINE yasso15_steady_state_r


SUBROUTINE yasso15_run_r(n_years, params, &
                          nwl_awen, fwl_awen, cwl_awen, &
                          xi_awe_array, xi_n_array, xi_h_array, &
                          leac, precip_array, &
                          diam_fwl, diam_cwl, &
                          C_init, &
                          C_out, resp_out)
  USE yasso15_mod
  IMPLICIT NONE

  INTEGER,  INTENT(IN)    :: n_years
  REAL(dp), INTENT(IN)    :: params(35)
  REAL(dp), INTENT(IN)    :: nwl_awen(n_years, 4)
  REAL(dp), INTENT(IN)    :: fwl_awen(n_years, 4)
  REAL(dp), INTENT(IN)    :: cwl_awen(n_years, 4)
  REAL(dp), INTENT(IN)    :: xi_awe_array(n_years)
  REAL(dp), INTENT(IN)    :: xi_n_array(n_years)
  REAL(dp), INTENT(IN)    :: xi_h_array(n_years)
  REAL(dp), INTENT(IN)    :: leac
  REAL(dp), INTENT(IN)    :: precip_array(n_years)
  REAL(dp), INTENT(IN)    :: diam_fwl, diam_cwl
  REAL(dp), INTENT(IN)    :: C_init(15)
  REAL(dp), INTENT(INOUT) :: C_out(n_years, 5)
  REAL(dp), INTENT(INOUT) :: resp_out(n_years)

  CALL yasso15_run(n_years, params, &
                   nwl_awen, fwl_awen, cwl_awen, &
                   xi_awe_array, xi_n_array, xi_h_array, &
                   leac, precip_array, &
                   diam_fwl, diam_cwl, &
                   C_init, &
                   C_out, resp_out)

END SUBROUTINE yasso15_run_r
