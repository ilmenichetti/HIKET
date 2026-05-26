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
!   xi (climate modifier) is computed entirely in R and passed as a pre-computed
!   array. This avoids duplicating the climate response logic in Fortran.
!
!   Two entry points are exposed to R:
!     yasso07_steady_state_r  -- computes equilibrium C pools from mean inputs
!     yasso07_run_r           -- transient forward simulation from C_init
!
! Changes from original (2025):
!   - matrixexp: capped scaling steps at MAX_SCALE_STEPS = 64 to prevent
!     infinite loop when matrix norm is very large (near-singular A matrices)
!   - matrixexp: increased Taylor terms from 10 to 20 for better accuracy
!     when scaling is heavy
!   - model_step: added diagonal dominance check before matrixexp call;
!     falls back to Euler step for near-singular matrices (pathological params)
!   - yasso07_run / yasso07_run_r: added C_final(15) output argument that
!     returns the terminal per-cohort pool state [C_nwl | C_fwl | C_cwl].
!     Without this, the R caller only has the 5 summed pools and cannot
!     correctly initialise a chained projection run (the 15-element C_init
!     required by the next yasso07_run call cannot be reconstructed losslessly
!     from 5 totals). C_final exposes the exact terminal state at zero cost.
!
! Compile with:
!   R CMD SHLIB yasso07.f90
!
! =============================================================================

MODULE yasso07_mod
  IMPLICIT NONE

  ! Double precision kind -- 15 significant digits, exponent up to 307
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 307)

  INTEGER, PARAMETER :: NPOOLS = 5      ! A, W, E, N, H

  ! Tolerance for treating a value as zero (used in solver and step guard)
  REAL(dp), PARAMETER :: TOL = 1.0E-12_dp

  ! Maximum number of scaling-and-squaring steps in matrixexp.
  ! Without this cap, near-singular matrices cause j to grow without bound,
  ! making the squaring loop run thousands of times and effectively hanging.
  ! 64 steps = scaling by 2^64 ~ 1.8e19, more than enough for any physical A.
  INTEGER, PARAMETER :: MAX_SCALE_STEPS = 64

CONTAINS

  ! =========================================================================
  ! BUILD DECOMPOSITION MATRIX A
  ! =========================================================================
  !
  ! Constructs the 5x5 matrix A for the ODE: dC/dt = A*C + b
  !
  ! The matrix encodes:
  !   - Diagonal:     net loss rate of each pool  A(i,i) = -alpha_i * xi * size_mod
  !   - Off-diagonal: carbon transfer from pool j to pool i
  !                   A(i,j) = p_{j->i} * |k_j|
  !                   where p_{j->i} is the fraction of decomposed carbon
  !                   from pool j that moves to pool i (not respired)
  !   - Row 5:        humus formation from all AWEN pools
  !                   A(5,j) = p_H * |k_j|  for j = 1..4
  !
  ! This convention matches the original R code: A = p_matrix %*% diag(k)
  ! where p_matrix has -1 on diagonal and transfer fractions off-diagonal.
  !
  ! Parameters layout (matches R parameter vector order):
  !   params(1:4)   = alpha_A, alpha_W, alpha_E, alpha_N  (base decomp rates)
  !   params(5:16)  = p_WA, p_EA, p_NA,   <- fractions going INTO pool A
  !                   p_AW, p_EW, p_NW,   <- fractions going INTO pool W
  !                   p_AE, p_WE, p_NE,   <- fractions going INTO pool E
  !                   p_AN, p_WN, p_EN    <- fractions going INTO pool N
  !   params(17:19) = beta1, beta2, gamma  (used in R for xi; ignored here)
  !   params(20:22) = delta1, delta2, r    (woody size modifier)
  !   params(23)    = p_H                  (humus formation fraction)
  !   params(24)    = alpha_H              (humus decomp rate)
  !
  ! =========================================================================

  SUBROUTINE build_matrix_A(params, xi, woody_diam, A)

    REAL(dp), INTENT(IN)  :: params(24)
    REAL(dp), INTENT(IN)  :: xi           ! climate modifier (pre-computed in R)
    REAL(dp), INTENT(IN)  :: woody_diam   ! diameter of woody litter (0 = non-woody)
    REAL(dp), INTENT(OUT) :: A(NPOOLS, NPOOLS)

    REAL(dp) :: alpha(NPOOLS)  ! base decomposition rates
    REAL(dp) :: k(NPOOLS)      ! effective rates = alpha * xi * size_mod
    REAL(dp) :: size_mod       ! woody size modifier (1.0 for non-woody)

    ! Transfer fractions: p_ij means "fraction of decomposed carbon
    ! from pool j that moves to pool i"
    ! Named as p{dest}{source}: p12 = fraction from pool 2 going to pool 1
    REAL(dp) :: p12, p13, p14   ! into A from W, E, N
    REAL(dp) :: p21, p23, p24   ! into W from A, E, N
    REAL(dp) :: p31, p32, p34   ! into E from A, W, N
    REAL(dp) :: p41, p42, p43   ! into N from A, W, E
    REAL(dp) :: p_H             ! fraction going to humus from each AWEN pool

    ! -- Unpack decomposition rates --
    alpha(1) = params(1)   ! pool A (acid-hydrolysable)
    alpha(2) = params(2)   ! pool W (water-soluble)
    alpha(3) = params(3)   ! pool E (ethanol-soluble)
    alpha(4) = params(4)   ! pool N (non-soluble)
    alpha(5) = params(24)  ! pool H (humus)

    ! -- Unpack transfer fractions (positional, matching R PA[5:16]) --
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
    p_H = params(23)   ! all AWEN -> H (same fraction for all pools)

    ! -- Woody size modifier (Eq. 3.1 in model description) --
    ! Reduces decomposition rate for larger woody pieces.
    ! Only applied to AWEN pools, not H.
    ! No effect when woody_diam = 0 (non-woody litter).
    IF (woody_diam > 0.0_dp) THEN
      size_mod = (1.0_dp + params(20) * woody_diam + &
                  params(21) * woody_diam**2)**params(22)
    ELSE
      size_mod = 1.0_dp
    END IF

    ! -- Effective decomposition rates --
    ! k_i = alpha_i * xi * size_mod  (negative = loss from pool)
    k(1:4) = -alpha(1:4) * xi * size_mod
    k(5)   = -alpha(5)   * xi          ! humus: no size effect

    ! -- Build matrix A --
    ! Start with zeros
    A = 0.0_dp

    ! Diagonal: net loss rate of each pool
    ! A(i,i) = k_i < 0  (carbon leaves pool i through decomposition)
    A(1,1) = k(1)
    A(2,2) = k(2)
    A(3,3) = k(3)
    A(4,4) = k(4)
    A(5,5) = k(5)

    ! Off-diagonal AWEN transfers: A(dest, source) = p * |k(source)|
    ! Fraction p of the carbon decomposed from source goes to dest.
    ! The remainder (1 - sum_of_p - p_H) is respired as CO2.
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

    ! Humus formation: all AWEN pools contribute to H at rate p_H
    A(5,1) = p_H * ABS(k(1))   ! A -> H
    A(5,2) = p_H * ABS(k(2))   ! W -> H
    A(5,3) = p_H * ABS(k(3))   ! E -> H
    A(5,4) = p_H * ABS(k(4))   ! N -> H

    ! No flows from H back to AWEN (humus is a terminal pool)
    A(1,5) = 0.0_dp
    A(2,5) = 0.0_dp
    A(3,5) = 0.0_dp
    A(4,5) = 0.0_dp

  END SUBROUTINE build_matrix_A


  ! =========================================================================
  ! SINGLE TIMESTEP: advance C by one year
  ! =========================================================================
  !
  ! Solves the ODE: dC/dt = A*C + b  analytically over interval [0, dt].
  !
  ! Analytical solution (Eq. 1.3 in model description):
  !   C(t+dt) = A^{-1} * (exp(A*dt) * (A*C(t) + b) - b)
  !
  ! Two special cases handled before the main path:
  !
  !   1. No decomposition (diagonal near zero): simple accumulation
  !      C(t+dt) = C(t) + b*dt
  !
  !   2. Near-singular matrix (column sums approach diagonal):
  !      This happens when transfer fractions nearly exhaust the pool budget,
  !      leaving almost no carbon respired. The matrix exponential then
  !      requires excessive scaling steps and effectively hangs.
  !      Fallback: Euler step C(t+dt) = C(t) + (A*C(t) + b)*dt
  !      Less accurate but always terminates. These are physically unrealistic
  !      parameter combinations (near-zero respiration) that the likelihood
  !      will strongly penalise anyway.
  !
  ! =========================================================================

  SUBROUTINE model_step(A, C_prev, litter_in, dt, C_next)

    REAL(dp), INTENT(IN)  :: A(NPOOLS, NPOOLS)
    REAL(dp), INTENT(IN)  :: C_prev(NPOOLS)   ! C pools at start of timestep
    REAL(dp), INTENT(IN)  :: litter_in(NPOOLS) ! litter input rate (b vector)
    REAL(dp), INTENT(IN)  :: dt               ! timestep length (1 year)
    REAL(dp), INTENT(OUT) :: C_next(NPOOLS)   ! C pools at end of timestep

    REAL(dp) :: At(NPOOLS, NPOOLS)    ! A * dt
    REAL(dp) :: mexpAt(NPOOLS, NPOOLS) ! matrix exponential exp(A*dt)
    REAL(dp) :: z1(NPOOLS), z2(NPOOLS) ! intermediate vectors
    REAL(dp) :: diag_abs, offdiag_sum  ! for diagonal dominance check
    INTEGER  :: i
    LOGICAL  :: well_conditioned

    ! -- Special case 1: no decomposition --
    ! If the first two diagonal entries are essentially zero, xi ~ 0
    ! (no precipitation or extreme conditions). Just accumulate litter.
    IF (ABS(A(1,1)) < TOL .AND. ABS(A(2,2)) < TOL) THEN
      C_next = C_prev + litter_in * dt
      RETURN
    END IF

    ! -- Special case 2: near-singular matrix check --
    ! The A matrix is diagonally dominant when |A(i,i)| > sum of |A(j,i)| for j/=i.
    ! This holds when the respiration fraction per pool > 0.
    ! When transfer fractions nearly exhaust the budget (respiration -> 0),
    ! diagonal dominance is lost and matrixexp requires huge scaling steps.
    ! We check: offdiag column sum < 0.9999 * |diagonal|
    ! (0.9999 rather than 1.0 gives a small safety margin)
    well_conditioned = .TRUE.
    DO i = 1, NPOOLS
      diag_abs    = ABS(A(i,i))
      offdiag_sum = SUM(ABS(A(:,i))) - diag_abs  ! sum of off-diagonal in column i
      IF (diag_abs < TOL .OR. offdiag_sum > diag_abs * 0.9999_dp) THEN
        well_conditioned = .FALSE.
        EXIT
      END IF
    END DO

    IF (.NOT. well_conditioned) THEN
      ! Euler fallback: first-order approximation, always terminates
      ! C(t+dt) ~ C(t) + dC/dt * dt  where dC/dt = A*C + b
      C_next = MAX(0.0_dp, C_prev + (MATMUL(A, C_prev) + litter_in) * dt)
      RETURN
    END IF

    ! -- Main path: analytical solution --
    ! z1 = A*C(t) + b
    z1 = MATMUL(A, C_prev) + litter_in

    ! Compute exp(A*dt)
    At = A * dt
    CALL matrixexp(At, mexpAt)

    ! z2 = exp(A*dt) * z1 - b
    z2 = MATMUL(mexpAt, z1) - litter_in

    ! C(t+dt) = A^{-1} * z2
    CALL solve5(A, z2, C_next)

  END SUBROUTINE model_step


  ! =========================================================================
  ! STEADY STATE: C* = -A^{-1} * b
  ! =========================================================================
  !
  ! At equilibrium dC/dt = 0, so A*C* + b = 0 => C* = -A^{-1}*b.
  ! Negative pool values are clipped to zero (can occur with singular A).
  !
  ! =========================================================================

  SUBROUTINE model_steady_state(A, litter_in, C_ss)

    REAL(dp), INTENT(IN)  :: A(NPOOLS, NPOOLS)
    REAL(dp), INTENT(IN)  :: litter_in(NPOOLS)
    REAL(dp), INTENT(OUT) :: C_ss(NPOOLS)

    REAL(dp) :: neg_A(NPOOLS, NPOOLS)

    ! Solve -A * C_ss = b  =>  C_ss = (-A)^{-1} * b
    neg_A = -A
    CALL solve5(neg_A, litter_in, C_ss)

    ! Clip negatives (numerical artefact from near-singular A)
    C_ss = MAX(0.0_dp, C_ss)

  END SUBROUTINE model_steady_state


  ! =========================================================================
  ! STEADY STATE INITIALISATION — three litter types
  ! =========================================================================
  !
  ! Computes equilibrium C pools separately for non-woody, fine woody, and
  ! coarse woody litter, then concatenates into C_init(15).
  !
  ! Each litter type has its own size modifier (via woody_diam), so three
  ! separate A matrices are built. xi_mean is the long-run average climate
  ! modifier, pre-computed in R.
  !
  ! Output layout: C_init = [C_nwl(1:5) | C_fwl(6:10) | C_cwl(11:15)]
  !
  ! =========================================================================

  SUBROUTINE yasso07_steady_state(params, &
                                   nwl_mean, fwl_mean, cwl_mean, &
                                   xi_mean, &
                                   diam_fwl, diam_cwl, &
                                   C_init)

    REAL(dp), INTENT(IN)  :: params(24)
    REAL(dp), INTENT(IN)  :: nwl_mean(4)   ! mean annual AWEN litter input, non-woody
    REAL(dp), INTENT(IN)  :: fwl_mean(4)   ! mean annual AWEN litter input, fine woody
    REAL(dp), INTENT(IN)  :: cwl_mean(4)   ! mean annual AWEN litter input, coarse woody
    REAL(dp), INTENT(IN)  :: xi_mean       ! mean climate modifier
    REAL(dp), INTENT(IN)  :: diam_fwl      ! diameter of fine woody litter [cm]
    REAL(dp), INTENT(IN)  :: diam_cwl      ! diameter of coarse woody litter [cm]
    REAL(dp), INTENT(OUT) :: C_init(3 * NPOOLS)

    REAL(dp) :: A(NPOOLS, NPOOLS)
    REAL(dp) :: litter(NPOOLS)
    REAL(dp) :: C_ss(NPOOLS)

    ! -- Non-woody litter (no size effect) --
    litter(1:4) = nwl_mean
    litter(5)   = 0.0_dp           ! no direct litter input to H pool
    CALL build_matrix_A(params, xi_mean, 0.0_dp, A)
    CALL model_steady_state(A, litter, C_ss)
    C_init(1:5) = C_ss

    ! -- Fine woody litter --
    litter(1:4) = fwl_mean
    litter(5)   = 0.0_dp
    CALL build_matrix_A(params, xi_mean, diam_fwl, A)
    CALL model_steady_state(A, litter, C_ss)
    C_init(6:10) = C_ss

    ! -- Coarse woody litter --
    litter(1:4) = cwl_mean
    litter(5)   = 0.0_dp
    CALL build_matrix_A(params, xi_mean, diam_cwl, A)
    CALL model_steady_state(A, litter, C_ss)
    C_init(11:15) = C_ss

  END SUBROUTINE yasso07_steady_state


  ! =========================================================================
  ! TRANSIENT FORWARD SIMULATION
  ! =========================================================================
  !
  ! Advances C pools forward in time, one year at a time.
  ! Each year uses the current xi from xi_array (pre-computed in R).
  !
  ! Three litter types (nwl, fwl, cwl) are tracked separately because they
  ! have different size modifiers, then summed for output.
  !
  ! Respiration is computed by mass balance:
  !   resp = C_total_before + litter_input - C_total_after
  !
  ! Input layout:
  !   C_init(15) = [C_nwl(1:5) | C_fwl(6:10) | C_cwl(11:15)]
  !   nwl_awen, fwl_awen, cwl_awen: annual litter inputs, rows = years
  !   xi_array: annual climate modifier, length = n_years
  !
  ! Output layout:
  !   C_out(n_years, 5): total C per pool (summed across litter types)
  !   resp_out(n_years): total annual respiration
  !   C_final(15): terminal per-cohort pool state [C_nwl | C_fwl | C_cwl]
  !     at the end of the last simulated year. Returned so the R caller can
  !     initialise a chained run (e.g. projection) with the exact 15-element
  !     state vector, avoiding the lossless-reconstruction problem that arises
  !     when only the 5 summed pools in C_out are available.
  !
  ! =========================================================================

  SUBROUTINE yasso07_run(n_years, params, &
                          nwl_awen, fwl_awen, cwl_awen, &
                          xi_array, &
                          diam_fwl, diam_cwl, &
                          C_init, &
                          C_out, resp_out, C_final)

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
    REAL(dp), INTENT(OUT)   :: C_final(3 * NPOOLS)

    REAL(dp) :: xi, dt
    REAL(dp) :: A_nwl(NPOOLS,NPOOLS), A_fwl(NPOOLS,NPOOLS), A_cwl(NPOOLS,NPOOLS)
    REAL(dp) :: C_nwl(NPOOLS), C_fwl(NPOOLS), C_cwl(NPOOLS)
    REAL(dp) :: C_nwl_next(NPOOLS), C_fwl_next(NPOOLS), C_cwl_next(NPOOLS)
    REAL(dp) :: litter_nwl(NPOOLS), litter_fwl(NPOOLS), litter_cwl(NPOOLS)
    REAL(dp) :: C_total_prev, C_total_next, input_total
    INTEGER  :: yr

    dt = 1.0_dp   ! annual timestep

    ! Unpack initial state into three litter-type pools
    C_nwl = C_init(1:5)
    C_fwl = C_init(6:10)
    C_cwl = C_init(11:15)

    DO yr = 1, n_years

      ! Climate modifier for this year (pre-computed in R)
      xi = xi_array(yr)

      ! Litter inputs for this year (H pool receives no direct input)
      litter_nwl(1:4) = nwl_awen(yr, :)
      litter_nwl(5)   = 0.0_dp
      litter_fwl(1:4) = fwl_awen(yr, :)
      litter_fwl(5)   = 0.0_dp
      litter_cwl(1:4) = cwl_awen(yr, :)
      litter_cwl(5)   = 0.0_dp

      ! Total C before step (for mass-balance respiration)
      C_total_prev = SUM(C_nwl) + SUM(C_fwl) + SUM(C_cwl)
      input_total  = SUM(litter_nwl) + SUM(litter_fwl) + SUM(litter_cwl)

      ! Build decomposition matrices for this year's xi
      CALL build_matrix_A(params, xi, 0.0_dp,   A_nwl)
      CALL build_matrix_A(params, xi, diam_fwl, A_fwl)
      CALL build_matrix_A(params, xi, diam_cwl, A_cwl)

      ! Advance each litter type by one year
      CALL model_step(A_nwl, C_nwl, litter_nwl, dt, C_nwl_next)
      CALL model_step(A_fwl, C_fwl, litter_fwl, dt, C_fwl_next)
      CALL model_step(A_cwl, C_cwl, litter_cwl, dt, C_cwl_next)

      ! Update state
      C_nwl = C_nwl_next
      C_fwl = C_fwl_next
      C_cwl = C_cwl_next

      ! Sum across litter types for output
      C_out(yr, :) = C_nwl + C_fwl + C_cwl

      ! Respiration by mass balance (avoids accumulating numerical error)
      C_total_next = SUM(C_out(yr, :))
      resp_out(yr) = C_total_prev + input_total - C_total_next

    END DO

    ! -- Write terminal per-cohort state --
    ! C_nwl, C_fwl, C_cwl hold the pool values at the end of the last year.
    ! Returned as C_final so the R caller can chain a projection run without
    ! losing the per-cohort breakdown (which C_out does not preserve).
    C_final(1:5)   = C_nwl
    C_final(6:10)  = C_fwl
    C_final(11:15) = C_cwl

  END SUBROUTINE yasso07_run


  ! =========================================================================
  ! MATRIX EXPONENTIAL
  ! =========================================================================
  !
  ! Computes B = exp(A) using Taylor series with scaling and squaring.
  ! Algorithm:
  !   1. Find integer j such that ||A|| < 2^j  (scaling step)
  !   2. Compute Taylor expansion of exp(A/2^j) using q terms
  !   3. Recover exp(A) by squaring j times: exp(A) = exp(A/2^j)^{2^j}
  !
  ! FIX (2025): The original code had no cap on j. For near-singular matrices
  ! (which arise when transfer fractions nearly exhaust the pool budget),
  ! ||A|| can be very large, making j grow without bound and causing the
  ! squaring loop to run thousands of times -- effectively hanging.
  ! We now cap j at MAX_SCALE_STEPS = 64, which covers any physically
  ! meaningful matrix (||A|| up to 2^64 ~ 1.8e19).
  !
  ! FIX (2025): Taylor terms increased from 10 to 20 for better accuracy
  ! when heavy scaling is used (large j means exp(A/2^j) is nearly I,
  ! so more terms are needed to capture the deviation accurately).
  !
  ! =========================================================================

  SUBROUTINE matrixexp(A, B)

    REAL(dp), INTENT(IN)  :: A(NPOOLS, NPOOLS)
    REAL(dp), INTENT(OUT) :: B(NPOOLS, NPOOLS)

    REAL(dp) :: C(NPOOLS, NPOOLS)   ! scaled matrix A/2^j
    REAL(dp) :: D(NPOOLS, NPOOLS)   ! current Taylor term
    REAL(dp) :: p                   ! Frobenius norm of A
    REAL(dp) :: normiter            ! current threshold (powers of 2)
    INTEGER  :: i, j, q

    q = 20   ! number of Taylor terms (increased from 10 for accuracy)

    ! Initialise B = identity matrix (zeroth Taylor term)
    B = 0.0_dp
    DO i = 1, NPOOLS
      B(i,i) = 1.0_dp
    END DO

    ! -- Step 1: find scaling factor 2^j such that ||A|| < 2^j --
    ! j counts how many times we need to halve A.
    ! normiter doubles each iteration: 2, 4, 8, 16, ...
    normiter = 2.0_dp
    j = 1
    CALL matrixnorm(A, p)
    DO
      IF (p < normiter) EXIT          ! found sufficient scaling
      normiter = normiter * 2.0_dp
      j = j + 1
      IF (j >= MAX_SCALE_STEPS) EXIT  ! cap: never exceed MAX_SCALE_STEPS
    END DO

    ! -- Step 2: Taylor expansion of exp(A/2^j) --
    ! C = A / 2^j  (normiter = 2^j at this point)
    C = A / normiter

    ! Taylor: exp(C) = I + C + C^2/2! + C^3/3! + ...
    ! B already = I, add first-order term
    B = B + C

    ! D tracks the current power C^i / i!
    D = C
    DO i = 2, q
      D = MATMUL(C, D) / REAL(i, dp)   ! D = C^i / i!
      B = B + D
    END DO

    ! -- Step 3: recover exp(A) by repeated squaring --
    ! exp(A) = exp(A/2^j)^{2^j} = B^{2^j}
    ! Each squaring doubles the exponent: B -> B^2 -> B^4 -> ... -> B^{2^j}
    DO i = 1, j
      B = MATMUL(B, B)
    END DO

  END SUBROUTINE matrixexp


  ! =========================================================================
  ! FROBENIUS NORM
  ! =========================================================================
  !
  ! ||A||_F = sqrt(sum of squared elements)
  ! Used by matrixexp to determine the scaling factor.
  !
  ! =========================================================================

  SUBROUTINE matrixnorm(A, nrm)

    REAL(dp), INTENT(IN)  :: A(NPOOLS, NPOOLS)
    REAL(dp), INTENT(OUT) :: nrm

    INTEGER :: i

    nrm = 0.0_dp
    DO i = 1, NPOOLS
      nrm = nrm + SUM(A(:,i)**2)   ! sum squared elements of column i
    END DO
    nrm = SQRT(nrm)

  END SUBROUTINE matrixnorm


  ! =========================================================================
  ! 5x5 LINEAR SOLVER
  ! =========================================================================
  !
  ! Solves A*x = b using Gaussian elimination with partial pivoting.
  ! Used in two places:
  !   - model_step:        A*x = z2  to recover C(t+dt)
  !   - model_steady_state: (-A)*x = b  to compute C* = -A^{-1}*b
  !
  ! Partial pivoting (row swapping) improves numerical stability when
  ! diagonal elements are small.
  !
  ! =========================================================================

  SUBROUTINE solve5(A, b, x)

    REAL(dp), INTENT(IN)  :: A(NPOOLS, NPOOLS)
    REAL(dp), INTENT(IN)  :: b(NPOOLS)
    REAL(dp), INTENT(OUT) :: x(NPOOLS)

    REAL(dp) :: U(NPOOLS, NPOOLS)   ! upper triangular form of A
    REAL(dp) :: c(NPOOLS)           ! transformed right-hand side
    INTEGER  :: i

    ! Reduce to upper triangular form
    CALL pgauss5(A, b, U, c)

    ! Back substitution: solve U*x = c from bottom up
    x(NPOOLS) = c(NPOOLS) / U(NPOOLS, NPOOLS)
    DO i = NPOOLS - 1, 1, -1
      x(i) = (c(i) - DOT_PRODUCT(U(i, i+1:NPOOLS), x(i+1:NPOOLS))) / U(i,i)
    END DO

  END SUBROUTINE solve5


  ! =========================================================================
  ! GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
  ! =========================================================================
  !
  ! Transforms A*x = b into upper triangular form U*x = c.
  ! At each step k, swaps rows to bring the largest element in column k
  ! to the pivot position -- this is partial pivoting, which avoids
  ! dividing by very small numbers that would amplify rounding errors.
  !
  ! =========================================================================

  SUBROUTINE pgauss5(A, b, U, c)

    REAL(dp), INTENT(IN)  :: A(NPOOLS, NPOOLS)
    REAL(dp), INTENT(IN)  :: b(NPOOLS)
    REAL(dp), INTENT(OUT) :: U(NPOOLS, NPOOLS)
    REAL(dp), INTENT(OUT) :: c(NPOOLS)

    INTEGER  :: k, j, pk, q
    REAL(dp) :: tmp_row(NPOOLS), tmp_val

    U = A
    c = b

    DO k = 1, NPOOLS - 1

      ! Find row with largest absolute value in column k (partial pivot)
      q = MAXLOC(ABS(U(k:NPOOLS, k)), 1)

      ! If pivot is not already in row k, swap rows k and k+q-1
      IF (q > 1) THEN
        pk      = k - 1 + q
        tmp_row = U(k, :)
        U(k, :) = U(pk, :)
        U(pk,:) = tmp_row
        tmp_val = c(k)
        c(k)    = c(pk)
        c(pk)   = tmp_val
      END IF

      ! Skip elimination if pivot is effectively zero (singular column)
      IF (ABS(U(k,k)) <= TOL) CYCLE

      ! Eliminate entries below pivot in column k
      U(k+1:NPOOLS, k) = U(k+1:NPOOLS, k) / U(k, k)   ! multipliers
      DO j = k + 1, NPOOLS
        U(j, k+1:NPOOLS) = U(j, k+1:NPOOLS) - U(j, k) * U(k, k+1:NPOOLS)
      END DO
      c(k+1:NPOOLS) = c(k+1:NPOOLS) - c(k) * U(k+1:NPOOLS, k)

    END DO

  END SUBROUTINE pgauss5

END MODULE yasso07_mod


! =============================================================================
! R-CALLABLE ENTRY POINTS
! =============================================================================
!
! These subroutines are outside the module so they have a C-compatible
! interface and can be called from R via .Fortran().
! They are thin wrappers that simply delegate to the module subroutines.
!
! =============================================================================

! -----------------------------------------------------------------------------
! Steady-state initialisation
! Called from R: yasso07_steady_state(params, nwl_mean, fwl_mean, cwl_mean,
!                                      xi_mean, diam_fwl, diam_cwl)
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
! Called from R: yasso07_run(n_years, params, nwl_awen, fwl_awen, cwl_awen,
!                             xi_array, diam_fwl, diam_cwl, C_init,
!                             C_out, resp_out, C_final)
! C_final(15): terminal per-cohort state [C_nwl | C_fwl | C_cwl].
!   Use as C_init for a chained projection run.
! -----------------------------------------------------------------------------
SUBROUTINE yasso07_run_r(n_years, params, &
                          nwl_awen, fwl_awen, cwl_awen, &
                          xi_array, &
                          diam_fwl, diam_cwl, &
                          C_init, &
                          C_out, resp_out, C_final)
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
  REAL(dp), INTENT(OUT)   :: C_final(15)

  CALL yasso07_run(n_years, params, &
                   nwl_awen, fwl_awen, cwl_awen, &
                   xi_array, &
                   diam_fwl, diam_cwl, &
                   C_init, &
                   C_out, resp_out, C_final)

END SUBROUTINE yasso07_run_r