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
!       AWE pools: beta1, beta2, gamma        (params 17-19)
!       N pool:    betaN1, betaN2, gammaN      (params 20-22)
!       H pool:    betaH1, betaH2, gammaH      (params 23-25)
!   - Leaching term on AWEN diagonal (requires leac scalar + precip)
!   - 35-parameter vector (original Järvenpää indexing preserved in comments)
!
! Architecture mirrors Yasso07:
!   - xi computation entirely in R, passed as pre-computed arrays
!     (xi_awe, xi_n, xi_h -- one per pool group per year)
!   - Two R entry points: yasso15_steady_state_r, yasso15_run_r
!
! Changes from original Järvenpää (2025):
!   - matrixexp: capped scaling steps at MAX_SCALE_STEPS = 64 to prevent
!     effective hang when matrix norm is very large (near-singular A matrices)
!   - matrixexp: increased Taylor terms from 10 to 20 for better accuracy
!     when scaling is heavy
!   - model_step: added diagonal dominance check before matrixexp call;
!     falls back to Euler step for near-singular matrices (pathological params)
!
! Compile with:
!   R CMD SHLIB yasso15.f90
!
! =============================================================================

MODULE yasso15_mod
  IMPLICIT NONE

  ! Double precision kind -- 15 significant digits, exponent up to 307
  !INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 307)
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(6, 37)   ! single precision
  
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
  ! Identical structure to Yasso07 except:
  !   - Three separate xi values (xi_awe, xi_n, xi_h) replace the single xi,
  !     reflecting pool-group-specific climate responses
  !   - Size modifier uses MIN(1, ...) with negative exponent convention
  !     (Yasso15 Eq. vs Yasso07 Eq. 3.1 -- different sign convention for r)
  !   - Leaching adds an extra term to the AWEN diagonal:
  !     A(i,i) += leac * precip / 1000  for i = 1..4
  !     This represents dissolved organic carbon loss via water percolation.
  !     No leaching for H pool.
  !
  ! Parameter vector (35 elements, standardised for this wrapper):
  !   1-4:   alpha_A, alpha_W, alpha_E, alpha_N   (base decomp rates)
  !   5-16:  p_WA, p_EA, p_NA, p_AW, p_EW, p_NW,
  !          p_AE, p_WE, p_NE, p_AN, p_WN, p_EN  (transfer fractions, same order as Yasso07)
  !   17-19: beta1, beta2, gamma                   (AWE climate response; orig params 22-24, 28)
  !   20-22: betaN1, betaN2, gammaN                (N climate response;   orig params 24-25, 29)
  !   23-25: betaH1, betaH2, gammaH                (H climate response;   orig params 26-27, 30)
  !   26-27: p_H, alpha_H                          (humus formation & decomp; orig 31-32)
  !   28-30: delta1, delta2, r                     (size modifier; orig 33-35)
  !   31-35: w1, w2, w3, w4, w5                   (leaching weights; stored but leac scalar used)
  !
  ! =========================================================================

  SUBROUTINE build_matrix_A(params, xi_awe, xi_n, xi_h, woody_diam, leac, precip, A)

    REAL(dp), INTENT(IN)  :: params(35)
    REAL(dp), INTENT(IN)  :: xi_awe      ! climate modifier for AWE pools (pre-computed in R)
    REAL(dp), INTENT(IN)  :: xi_n        ! climate modifier for N pool
    REAL(dp), INTENT(IN)  :: xi_h        ! climate modifier for H pool
    REAL(dp), INTENT(IN)  :: woody_diam  ! diameter of woody litter (0 = non-woody)
    REAL(dp), INTENT(IN)  :: leac        ! leaching scalar (site-level)
    REAL(dp), INTENT(IN)  :: precip      ! annual precipitation [mm]
    REAL(dp), INTENT(OUT) :: A(NPOOLS, NPOOLS)

    REAL(dp) :: alpha(NPOOLS)  ! base decomposition rates
    REAL(dp) :: k(NPOOLS)      ! effective rates (alpha * xi * size_mod)
    REAL(dp) :: size_mod       ! woody size modifier (1.0 for non-woody)

    ! Transfer fractions: p_ij = fraction of decomposed carbon
    ! from pool j that moves to pool i (not respired, not to humus)
    REAL(dp) :: p12, p13, p14   ! into A from W, E, N
    REAL(dp) :: p21, p23, p24   ! into W from A, E, N
    REAL(dp) :: p31, p32, p34   ! into E from A, W, N
    REAL(dp) :: p41, p42, p43   ! into N from A, W, E
    REAL(dp) :: p_H             ! fraction going to humus from each AWEN pool

    INTEGER  :: i

    ! -- Unpack decomposition rates --
    alpha(1) = params(1)   ! pool A (acid-hydrolysable)
    alpha(2) = params(2)   ! pool W (water-soluble)
    alpha(3) = params(3)   ! pool E (ethanol-soluble)
    alpha(4) = params(4)   ! pool N (non-soluble)
    alpha(5) = params(27)  ! pool H (humus)

    ! -- Unpack transfer fractions (same positional order as Yasso07 params 5-16) --
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
    p_H = params(26)   ! all AWEN -> H (same fraction for all pools)

    ! -- Woody size modifier --
    ! Yasso15 convention: size_mod = MIN(1, (1 + d1*d + d2*d^2)^(-|r|))
    ! Note the negative exponent and MIN(1,...) cap -- different from Yasso07.
    ! Larger pieces decompose more slowly (size_mod <= 1).
    ! No effect when woody_diam = 0 (non-woody litter).
    IF (woody_diam > 0.0_dp) THEN
      size_mod = MIN(1.0_dp, (1.0_dp + params(28) * woody_diam + &
                               params(29) * woody_diam**2)**(-ABS(params(30))))
    ELSE
      size_mod = 1.0_dp
    END IF

    ! -- Effective decomposition rates --
    ! Pool-group-specific xi: AWE share one xi, N and H have their own.
    ! ABS() on alpha ensures sign does not matter in parameter vector.
    k(1) = -ABS(alpha(1)) * xi_awe * size_mod   ! A
    k(2) = -ABS(alpha(2)) * xi_awe * size_mod   ! W
    k(3) = -ABS(alpha(3)) * xi_awe * size_mod   ! E
    k(4) = -ABS(alpha(4)) * xi_n   * size_mod   ! N (different climate response)
    k(5) = -ABS(alpha(5)) * xi_h                ! H (no size effect)

    ! -- Build matrix A --
    ! Start with zeros
    A = 0.0_dp

    ! Diagonal: net loss rate of each pool
    DO i = 1, NPOOLS
      A(i,i) = k(i)
    END DO

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

    ! -- Leaching: additive loss from AWEN diagonal --
    ! Represents dissolved organic carbon lost via water percolation.
    ! Proportional to precipitation (converted mm -> m).
    ! H pool is not subject to leaching.
    DO i = 1, 4
      A(i,i) = A(i,i) + leac * precip / 1000.0_dp
    END DO

  END SUBROUTINE build_matrix_A


  ! =========================================================================
  ! SINGLE TIMESTEP: advance C by one year
  ! =========================================================================
  !
  ! Identical to Yasso07 model_step except for the additional variables
  ! needed for the diagonal dominance check (diag_abs, offdiag_sum,
  ! well_conditioned). The logic and fallback are the same.
  !
  ! Solves the ODE: dC/dt = A*C + b  analytically over interval [0, dt].
  !
  ! Analytical solution:
  !   C(t+dt) = A^{-1} * (exp(A*dt) * (A*C(t) + b) - b)
  !
  ! Special cases:
  !   1. No decomposition (diagonal ~ 0): simple accumulation
  !   2. Near-singular matrix (diagonal dominance lost): Euler fallback
  !      These correspond to near-zero respiration -- physically unrealistic,
  !      strongly penalised by the likelihood during MCMC.
  !
  ! =========================================================================

  SUBROUTINE model_step(A, C_prev, litter_in, dt, C_next)

    REAL(dp), INTENT(IN)  :: A(NPOOLS, NPOOLS)
    REAL(dp), INTENT(IN)  :: C_prev(NPOOLS)    ! C pools at start of timestep
    REAL(dp), INTENT(IN)  :: litter_in(NPOOLS) ! litter input rate (b vector)
    REAL(dp), INTENT(IN)  :: dt                ! timestep length (1 year)
    REAL(dp), INTENT(OUT) :: C_next(NPOOLS)    ! C pools at end of timestep

    REAL(dp) :: At(NPOOLS, NPOOLS)     ! A * dt
    REAL(dp) :: mexpAt(NPOOLS, NPOOLS) ! matrix exponential exp(A*dt)
    REAL(dp) :: z1(NPOOLS), z2(NPOOLS) ! intermediate vectors
    REAL(dp) :: diag_abs, offdiag_sum  ! for diagonal dominance check
    INTEGER  :: i
    LOGICAL  :: well_conditioned

    ! -- Special case 1: no decomposition --
    ! If the first two diagonal entries are essentially zero, xi ~ 0.
    ! Just accumulate litter without decomposition.
    IF (ABS(A(1,1)) < TOL .AND. ABS(A(2,2)) < TOL) THEN
      C_next = C_prev + litter_in * dt
      RETURN
    END IF

    ! -- Special case 2: near-singular matrix check --
    ! Diagonal dominance: |A(i,i)| > sum of |A(j,i)| for j/=i.
    ! Lost when transfer fractions nearly exhaust the pool budget (respiration -> 0).
    ! Near-singular matrices cause matrixexp to require excessive scaling steps.
    ! Threshold 0.9999 gives a small safety margin below exact singularity.
    well_conditioned = .TRUE.
    DO i = 1, NPOOLS
      diag_abs    = ABS(A(i,i))
      offdiag_sum = SUM(ABS(A(:,i))) - diag_abs
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
  ! separate A matrices are built. xi_*_mean are long-run average climate
  ! modifiers, pre-computed in R.
  !
  ! Output layout: C_init = [C_nwl(1:5) | C_fwl(6:10) | C_cwl(11:15)]
  !
  ! =========================================================================

  SUBROUTINE yasso15_steady_state(params, &
                                   nwl_mean, fwl_mean, cwl_mean, &
                                   xi_awe_mean, xi_n_mean, xi_h_mean, &
                                   leac, precip_mean, &
                                   diam_fwl, diam_cwl, &
                                   C_init)

    REAL(dp), INTENT(IN)  :: params(35)
    REAL(dp), INTENT(IN)  :: nwl_mean(4)    ! mean annual AWEN litter input, non-woody
    REAL(dp), INTENT(IN)  :: fwl_mean(4)    ! mean annual AWEN litter input, fine woody
    REAL(dp), INTENT(IN)  :: cwl_mean(4)    ! mean annual AWEN litter input, coarse woody
    REAL(dp), INTENT(IN)  :: xi_awe_mean    ! mean climate modifier for AWE pools
    REAL(dp), INTENT(IN)  :: xi_n_mean      ! mean climate modifier for N pool
    REAL(dp), INTENT(IN)  :: xi_h_mean      ! mean climate modifier for H pool
    REAL(dp), INTENT(IN)  :: leac           ! leaching scalar (site-level)
    REAL(dp), INTENT(IN)  :: precip_mean    ! mean annual precipitation [mm]
    REAL(dp), INTENT(IN)  :: diam_fwl       ! diameter of fine woody litter [cm]
    REAL(dp), INTENT(IN)  :: diam_cwl       ! diameter of coarse woody litter [cm]
    REAL(dp), INTENT(OUT) :: C_init(3 * NPOOLS)

    REAL(dp) :: A(NPOOLS, NPOOLS)
    REAL(dp) :: litter(NPOOLS)
    REAL(dp) :: C_ss(NPOOLS)

    ! -- Non-woody litter (no size effect) --
    litter(1:4) = nwl_mean
    litter(5)   = 0.0_dp           ! no direct litter input to H pool
    CALL build_matrix_A(params, xi_awe_mean, xi_n_mean, xi_h_mean, &
                         0.0_dp, leac, precip_mean, A)
    CALL model_steady_state(A, litter, C_ss)
    C_init(1:5) = C_ss

    ! -- Fine woody litter --
    litter(1:4) = fwl_mean
    litter(5)   = 0.0_dp
    CALL build_matrix_A(params, xi_awe_mean, xi_n_mean, xi_h_mean, &
                         diam_fwl, leac, precip_mean, A)
    CALL model_steady_state(A, litter, C_ss)
    C_init(6:10) = C_ss

    ! -- Coarse woody litter --
    litter(1:4) = cwl_mean
    litter(5)   = 0.0_dp
    CALL build_matrix_A(params, xi_awe_mean, xi_n_mean, xi_h_mean, &
                         diam_cwl, leac, precip_mean, A)
    CALL model_steady_state(A, litter, C_ss)
    C_init(11:15) = C_ss

  END SUBROUTINE yasso15_steady_state


  ! =========================================================================
  ! TRANSIENT FORWARD SIMULATION
  ! =========================================================================
  !
  ! Advances C pools forward in time, one year at a time.
  ! Receives pre-computed xi arrays (one value per pool group per year)
  ! and precip_array (needed for the leaching term in build_matrix_A).
  !
  ! Three litter types (nwl, fwl, cwl) are tracked separately because they
  ! have different size modifiers, then summed for output.
  !
  ! Respiration is computed by mass balance:
  !   resp = C_total_before + litter_input - C_total_after
  !
  ! Input layout:
  !   C_init(15) = [C_nwl(1:5) | C_fwl(6:10) | C_cwl(11:15)]
  !   *_awen: annual litter inputs, rows = years, cols = AWEN
  !   xi_*_array: annual climate modifiers per pool group
  !   precip_array: annual precipitation [mm]
  !
  ! Output layout:
  !   C_out(n_years, 5): total C per pool (summed across litter types)
  !   resp_out(n_years): total annual respiration
  !
  ! =========================================================================

  SUBROUTINE yasso15_run(n_years, params, &
                          nwl_awen, fwl_awen, cwl_awen, &
                          xi_awe_array, xi_n_array, xi_h_array, &
                          leac, precip_array, &
                          diam_fwl, diam_cwl, &
                          C_init, &
                          C_out, resp_out)

    INTEGER,  INTENT(IN)    :: n_years
    REAL(dp), INTENT(IN)    :: params(35)
    REAL(dp), INTENT(IN)    :: nwl_awen(n_years, 4)
    REAL(dp), INTENT(IN)    :: fwl_awen(n_years, 4)
    REAL(dp), INTENT(IN)    :: cwl_awen(n_years, 4)
    REAL(dp), INTENT(IN)    :: xi_awe_array(n_years)  ! climate modifier for AWE
    REAL(dp), INTENT(IN)    :: xi_n_array(n_years)    ! climate modifier for N
    REAL(dp), INTENT(IN)    :: xi_h_array(n_years)    ! climate modifier for H
    REAL(dp), INTENT(IN)    :: leac                   ! leaching scalar (site-level)
    REAL(dp), INTENT(IN)    :: precip_array(n_years)  ! annual precipitation [mm]
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

    dt = 1.0_dp   ! annual timestep

    ! Unpack initial state into three litter-type pools
    C_nwl = C_init(1:5)
    C_fwl = C_init(6:10)
    C_cwl = C_init(11:15)

    DO yr = 1, n_years

      ! Climate modifiers and precipitation for this year (pre-computed in R)
      xi_awe = xi_awe_array(yr)
      xi_n   = xi_n_array(yr)
      xi_h   = xi_h_array(yr)
      precip = precip_array(yr)

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

      ! Build decomposition matrices for this year's climate
      CALL build_matrix_A(params, xi_awe, xi_n, xi_h, 0.0_dp,   leac, precip, A_nwl)
      CALL build_matrix_A(params, xi_awe, xi_n, xi_h, diam_fwl, leac, precip, A_fwl)
      CALL build_matrix_A(params, xi_awe, xi_n, xi_h, diam_cwl, leac, precip, A_cwl)

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

  END SUBROUTINE yasso15_run


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
  ! to the pivot position -- avoids dividing by very small numbers.
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

END MODULE yasso15_mod


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
! Called from R: yasso15_steady_state(params, nwl_mean, fwl_mean, cwl_mean,
!                                      xi_awe_mean, xi_n_mean, xi_h_mean,
!                                      leac, precip_mean, diam_fwl, diam_cwl)
! -----------------------------------------------------------------------------
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


! -----------------------------------------------------------------------------
! Transient forward simulation
! Called from R: yasso15_run(n_years, params, nwl_awen, fwl_awen, cwl_awen,
!                             xi_awe_array, xi_n_array, xi_h_array,
!                             leac, precip_array, diam_fwl, diam_cwl, C_init)
! -----------------------------------------------------------------------------
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