! =============================================================================
! Yasso07 Fortran Implementation for R Interface
! =============================================================================
!
! Based on:
!   Tuomi et al. (2009) Eco. Mod. 220: 3362-3371
!   Original R implementation by Taru Palosuo (FMI, December 2011)
!   Matrix exponential routines adapted from Yasso15 Fortran (Järvenpää, 2015)
!
! Architecture:
!   xi computation is handled entirely in R and passed as a pre-computed array.
!   Two entry points are exposed to R:
!     yasso07_steady_state_r  -- computes C* from mean inputs and mean xi
!     yasso07_run_r           -- transient forward simulation from C_init
!
! Compile with:
!   R CMD SHLIB yasso07.f90
!
! =============================================================================

MODULE yasso07_mod
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 307)  ! double precision
  INTEGER, PARAMETER :: NPOOLS = 5
  REAL(dp), PARAMETER :: TOL = 1.0E-12_dp
CONTAINS

  ! =========================================================================
  ! BUILD DECOMPOSITION MATRIX A
  ! =========================================================================

  SUBROUTINE build_matrix_A(params, xi, woody_diam, A)
    ! Constructs the 5x5 decomposition matrix.
    !
    ! params(1:4)   = alpha_A, alpha_W, alpha_E, alpha_N (base decomp rates)
    ! params(5:16)  = p_WA, p_EA, p_NA, p_AW, p_EW, p_NW,
    !                 p_AE, p_WE, p_NE, p_AN, p_WN, p_EN  (transfer fractions)
    ! params(17:18) = beta1, beta2  (used in R for xi; not used here)
    ! params(19)    = gamma          (used in R for xi; not used here)
    ! params(20:22) = delta1, delta2, r (size modifier)
    ! params(23)    = p_H (humus formation fraction)
    ! params(24)    = alpha_H (humus decomposition rate)

    REAL(dp), INTENT(IN)  :: params(24)
    REAL(dp), INTENT(IN)  :: xi, woody_diam
    REAL(dp), INTENT(OUT) :: A(NPOOLS, NPOOLS)

    REAL(dp) :: alpha(NPOOLS), k(NPOOLS), size_mod
    REAL(dp) :: p12, p13, p14
    REAL(dp) :: p21, p23, p24
    REAL(dp) :: p31, p32, p34
    REAL(dp) :: p41, p42, p43
    REAL(dp) :: p_H

    ! Base decomposition rates
    alpha(1) = params(1)   ! A
    alpha(2) = params(2)   ! W
    alpha(3) = params(3)   ! E
    alpha(4) = params(4)   ! N
    alpha(5) = params(24)  ! H

    ! Transfer fractions (positional, matching R original PA[5:16])
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
    p_H = params(23)

    ! Size modifier (Eq. 3.1): only affects AWEN pools, not H
    IF (woody_diam > 0.0_dp) THEN
      size_mod = (1.0_dp + params(20) * woody_diam + &
                  params(21) * woody_diam**2)**params(22)
    ELSE
      size_mod = 1.0_dp
    END IF

    ! Effective decomposition rates: k_i = -alpha_i * xi * size_mod
    k(1:4) = -alpha(1:4) * xi * size_mod
    k(5)   = -alpha(5) * xi   ! no size effect on H

    ! Build A: diagonal = pool loss rates
    A = 0.0_dp
    A(1,1) = k(1)
    A(2,2) = k(2)
    A(3,3) = k(3)
    A(4,4) = k(4)
    A(5,5) = k(5)

    ! Off-diagonal: transfers between AWEN pools
    ! Convention: A(dest, source) = p_{source->dest} * |k(source)|
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

    ! Humus formation: all AWEN pools -> H
    A(5,1) = p_H * ABS(k(1))
    A(5,2) = p_H * ABS(k(2))
    A(5,3) = p_H * ABS(k(3))
    A(5,4) = p_H * ABS(k(4))

    ! No flows from H back to AWEN
    A(1,5) = 0.0_dp
    A(2,5) = 0.0_dp
    A(3,5) = 0.0_dp
    A(4,5) = 0.0_dp

  END SUBROUTINE build_matrix_A


  ! =========================================================================
  ! SINGLE TIMESTEP: C(t+dt) from C(t)
  ! =========================================================================

  SUBROUTINE model_step(A, C_prev, litter_in, dt, C_next)
    ! Analytical solution for one timestep.

    REAL(dp), INTENT(IN)  :: A(NPOOLS, NPOOLS)
    REAL(dp), INTENT(IN)  :: C_prev(NPOOLS), litter_in(NPOOLS), dt
    REAL(dp), INTENT(OUT) :: C_next(NPOOLS)

    REAL(dp) :: At(NPOOLS, NPOOLS), mexpAt(NPOOLS, NPOOLS)
    REAL(dp) :: z1(NPOOLS), z2(NPOOLS)

    ! Degenerate case: no decomposition
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

  SUBROUTINE yasso07_steady_state(params, &
                                   nwl_mean, fwl_mean, cwl_mean, &
                                   xi_mean, &
                                   diam_fwl, diam_cwl, &
                                   C_init)
    ! Computes steady-state carbon pools for three litter types.
    ! xi_mean is pre-computed in R from mean climate and passed in.
    ! C_init(15): first 5 = nwl pools, next 5 = fwl pools, last 5 = cwl pools.

    REAL(dp), INTENT(IN)  :: params(24)
    REAL(dp), INTENT(IN)  :: nwl_mean(4), fwl_mean(4), cwl_mean(4)  ! mean AWEN inputs
    REAL(dp), INTENT(IN)  :: xi_mean
    REAL(dp), INTENT(IN)  :: diam_fwl, diam_cwl
    REAL(dp), INTENT(OUT) :: C_init(3 * NPOOLS)   ! [C_nwl | C_fwl | C_cwl]

    REAL(dp) :: A(NPOOLS, NPOOLS)
    REAL(dp) :: litter(NPOOLS)
    REAL(dp) :: C_ss(NPOOLS)

    ! --- Non-woody litter ---
    litter(1:4) = nwl_mean
    litter(5)   = 0.0_dp
    CALL build_matrix_A(params, xi_mean, 0.0_dp, A)
    CALL model_steady_state(A, litter, C_ss)
    C_init(1:5) = C_ss

    ! --- Fine woody litter ---
    litter(1:4) = fwl_mean
    litter(5)   = 0.0_dp
    CALL build_matrix_A(params, xi_mean, diam_fwl, A)
    CALL model_steady_state(A, litter, C_ss)
    C_init(6:10) = C_ss

    ! --- Coarse woody litter ---
    litter(1:4) = cwl_mean
    litter(5)   = 0.0_dp
    CALL build_matrix_A(params, xi_mean, diam_cwl, A)
    CALL model_steady_state(A, litter, C_ss)
    C_init(11:15) = C_ss

  END SUBROUTINE yasso07_steady_state


  ! =========================================================================
  ! TRANSIENT FORWARD SIMULATION — called from R entry point
  ! =========================================================================

  SUBROUTINE yasso07_run(n_years, params, &
                          nwl_awen, fwl_awen, cwl_awen, &
                          xi_array, &
                          diam_fwl, diam_cwl, &
                          C_init, &
                          C_out, resp_out)
    ! Transient forward simulation. Receives pre-computed xi_array and
    ! initial state C_init — no climate data or steady-state logic here.
    !
    ! C_init(15): [C_nwl(1:5) | C_fwl(6:10) | C_cwl(11:15)]
    ! C_out(n_years, 5): summed pool states across litter types per year

    INTEGER,  INTENT(IN)    :: n_years
    REAL(dp), INTENT(IN)    :: params(24)
    REAL(dp), INTENT(IN)    :: nwl_awen(n_years, 4)
    REAL(dp), INTENT(IN)    :: fwl_awen(n_years, 4)
    REAL(dp), INTENT(IN)    :: cwl_awen(n_years, 4)
    REAL(dp), INTENT(IN)    :: xi_array(n_years)
    REAL(dp), INTENT(IN)    :: diam_fwl, diam_cwl
    REAL(dp), INTENT(IN)    :: C_init(3 * NPOOLS)
    REAL(dp), INTENT(INOUT) :: C_out(n_years, NPOOLS)
    REAL(dp), INTENT(INOUT) :: resp_out(n_years)

    REAL(dp) :: xi, dt
    REAL(dp) :: A_nwl(NPOOLS,NPOOLS), A_fwl(NPOOLS,NPOOLS), A_cwl(NPOOLS,NPOOLS)
    REAL(dp) :: C_nwl(NPOOLS), C_fwl(NPOOLS), C_cwl(NPOOLS)
    REAL(dp) :: C_nwl_next(NPOOLS), C_fwl_next(NPOOLS), C_cwl_next(NPOOLS)
    REAL(dp) :: litter_nwl(NPOOLS), litter_fwl(NPOOLS), litter_cwl(NPOOLS)
    REAL(dp) :: C_total_prev, C_total_next, input_total
    INTEGER  :: yr

    dt = 1.0_dp  ! annual timestep

    ! Unpack initial states
    C_nwl = C_init(1:5)
    C_fwl = C_init(6:10)
    C_cwl = C_init(11:15)

    DO yr = 1, n_years

      ! xi pre-computed in R — just index into array
      xi = xi_array(yr)

      ! Current year litter inputs (H pool input = 0)
      litter_nwl(1:4) = nwl_awen(yr, :)
      litter_nwl(5)   = 0.0_dp
      litter_fwl(1:4) = fwl_awen(yr, :)
      litter_fwl(5)   = 0.0_dp
      litter_cwl(1:4) = cwl_awen(yr, :)
      litter_cwl(5)   = 0.0_dp

      ! Total C before step (for respiration by mass balance)
      C_total_prev = SUM(C_nwl) + SUM(C_fwl) + SUM(C_cwl)
      input_total  = SUM(litter_nwl) + SUM(litter_fwl) + SUM(litter_cwl)

      ! Build matrices with current xi
      CALL build_matrix_A(params, xi, 0.0_dp,   A_nwl)
      CALL build_matrix_A(params, xi, diam_fwl, A_fwl)
      CALL build_matrix_A(params, xi, diam_cwl, A_cwl)

      ! Step each litter type
      CALL model_step(A_nwl, C_nwl, litter_nwl, dt, C_nwl_next)
      CALL model_step(A_fwl, C_fwl, litter_fwl, dt, C_fwl_next)
      CALL model_step(A_cwl, C_cwl, litter_cwl, dt, C_cwl_next)

      ! Update pools
      C_nwl = C_nwl_next
      C_fwl = C_fwl_next
      C_cwl = C_cwl_next

      ! Store output: sum across litter types
      C_out(yr, :) = C_nwl + C_fwl + C_cwl

      ! Respiration by mass balance
      C_total_next = SUM(C_out(yr, :))
      resp_out(yr) = C_total_prev + input_total - C_total_next

    END DO

  END SUBROUTINE yasso07_run


  ! =========================================================================
  ! MATRIX EXPONENTIAL (Taylor series with scaling & squaring)
  ! Adapted from Yasso15 Fortran (Järvenpää, 2015)
  ! =========================================================================

  SUBROUTINE matrixexp(A, B)
    REAL(dp), INTENT(IN)  :: A(NPOOLS, NPOOLS)
    REAL(dp), INTENT(OUT) :: B(NPOOLS, NPOOLS)

    REAL(dp) :: C(NPOOLS, NPOOLS), D(NPOOLS, NPOOLS)
    REAL(dp) :: p, normiter
    INTEGER  :: i, j, q

    q = 10  ! terms in Taylor expansion

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
  ! Adapted from Yasso15 Fortran (Järvenpää, 2015)
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

END MODULE yasso07_mod


! =============================================================================
! R-CALLABLE ENTRY POINTS (outside module, C-compatible interface)
! =============================================================================

! -----------------------------------------------------------------------------
! Steady-state initialisation
! -----------------------------------------------------------------------------
SUBROUTINE yasso07_steady_state_r(params, &
                                   nwl_mean, fwl_mean, cwl_mean, &
                                   xi_mean, &
                                   diam_fwl, diam_cwl, &
                                   C_init)
  USE yasso07_mod
  IMPLICIT NONE

  REAL(dp), INTENT(IN)  :: params(24)
  REAL(dp), INTENT(IN)  :: nwl_mean(4), fwl_mean(4), cwl_mean(4)
  REAL(dp), INTENT(IN)  :: xi_mean
  REAL(dp), INTENT(IN)  :: diam_fwl, diam_cwl
  REAL(dp), INTENT(OUT) :: C_init(15)

  CALL yasso07_steady_state(params, &
                             nwl_mean, fwl_mean, cwl_mean, &
                             xi_mean, &
                             diam_fwl, diam_cwl, &
                             C_init)

END SUBROUTINE yasso07_steady_state_r


! -----------------------------------------------------------------------------
! Transient forward simulation
! -----------------------------------------------------------------------------
SUBROUTINE yasso07_run_r(n_years, params, &
                          nwl_awen, fwl_awen, cwl_awen, &
                          xi_array, &
                          diam_fwl, diam_cwl, &
                          C_init, &
                          C_out, resp_out)
  USE yasso07_mod
  IMPLICIT NONE

  INTEGER,  INTENT(IN)    :: n_years
  REAL(dp), INTENT(IN)    :: params(24)
  REAL(dp), INTENT(IN)    :: nwl_awen(n_years, 4)
  REAL(dp), INTENT(IN)    :: fwl_awen(n_years, 4)
  REAL(dp), INTENT(IN)    :: cwl_awen(n_years, 4)
  REAL(dp), INTENT(IN)    :: xi_array(n_years)
  REAL(dp), INTENT(IN)    :: diam_fwl, diam_cwl
  REAL(dp), INTENT(IN)    :: C_init(15)
  REAL(dp), INTENT(INOUT) :: C_out(n_years, 5)
  REAL(dp), INTENT(INOUT) :: resp_out(n_years)

  CALL yasso07_run(n_years, params, &
                   nwl_awen, fwl_awen, cwl_awen, &
                   xi_array, &
                   diam_fwl, diam_cwl, &
                   C_init, &
                   C_out, resp_out)

END SUBROUTINE yasso07_run_r
