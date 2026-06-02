# =============================================================================
# yasso15_wrapper_transient.R
#
# Standalone transient wrapper for YASSO15.
# Complete copy of yasso15_wrapper.R with yasso*_transient_init appended.
# Does NOT source the original wrapper -- self-contained.
# =============================================================================

# =============================================================================
# Yasso15 R Wrapper for Fortran Implementation
# =============================================================================
#
# Compile:  system("R CMD SHLIB yasso15.f90")
# Load:     dyn.load("yasso15.so")
#
# Three-stage call pattern:
#   1. compute_xi_yasso15()    -- three climate modifier time series (R)
#   2. yasso15_steady_state()  -- initial pool state C* (Fortran)
#   3. yasso15_run()           -- transient forward simulation (Fortran)
#
# Key difference from Yasso07: pool-group-specific climate responses.
# AWE, N, and H pools each have their own beta1/beta2/gamma parameters,
# producing three separate xi arrays passed to Fortran.
#
# Precision note:
#   Despite REAL(dp) declarations throughout the source, dp is defined as
#   SELECTED_REAL_KIND(6, 37) = REAL(4) (single precision) in yasso15_mod.
#   All .Fortran() calls use as.single() / single() accordingly.
#   Outputs are cast back to double for R-side use.
#
# Steady-state xi note:
#   The reference (mod5c with steady=TRUE) computes xi internally from the
#   mean climate vector c(mean_T, mean_precip, mean_Tamp). This is NOT the
#   same as mean(xi_annual) because xi is nonlinear in temperature and
#   precipitation. yasso15_steady_state() therefore accepts raw mean climate
#   scalars and computes xi internally, matching the reference behaviour.
#
# =============================================================================

# Parameter vector: 35-element standardised format
# Indices:
#   1-4:   alpha_A, alpha_W, alpha_E, alpha_N   (base decomp rates)
#   5-16:  p_WA, p_EA, p_NA, p_AW, p_EW, p_NW,
#          p_AE, p_WE, p_NE, p_AN, p_WN, p_EN  (transfer fractions)
#   17-19: beta1, beta2, gamma                   (AWE climate response)
#   20-22: betaN1, betaN2, gammaN                (N climate response)
#   23-25: betaH1, betaH2, gammaH                (H climate response)
#   26-27: p_H, alpha_H                          (humus)
#   28-30: delta1, delta2, r                     (size modifier)
#   31-35: w1, w2, w3, w4, w5                   (leaching weights; leac scalar used in practice)

YASSO15_PARAM_NAMES <- c(
  "alpha_A", "alpha_W", "alpha_E", "alpha_N",
  "p_WA", "p_EA", "p_NA",
  "p_AW", "p_EW", "p_NW",
  "p_AE", "p_WE", "p_NE",
  "p_AN", "p_WN", "p_EN",
  "beta1",  "beta2",  "gamma",
  "betaN1", "betaN2", "gammaN",
  "betaH1", "betaH2", "gammaH",
  "p_H", "alpha_H",
  "delta1", "delta2", "r",
  "w1", "w2", "w3", "w4", "w5"
)

YASSO15_DEFAULT_PARAMS <- c(
  0.48971473, 4.9138734, 0.24197346, 0.094876416,          # alpha_A/W/E/N (Fortran takes ABS)
  0.43628932, 0.24997402, 0.91512685,                      # p_WA, p_EA, p_NA
  0.99258227, 0.083853738, 0.011476783,                    # p_AW, p_EW, p_NW
  6.08e-04, 4.76e-04, 0.066037729,                        # p_AE, p_WE, p_NE
  7.71e-04, 0.10401742, 0.64880756,                       # p_AN, p_WN, p_EN
  0.090598047, -2.14e-04, -1.8089202,                     # beta1,  beta2,  gamma  (AWE)
  0.048772465, -7.91e-05, -1.1725473,                     # betaN1, betaN2, gammaN (N)
  0.035185492, -2.09e-04, -12.535951,                     # betaH1, betaH2, gammaH (H)
  0.004596472, 0.001302583,                               # p_H, alpha_H
  -0.43892271, 1.2674668, 0.25691424,                     # delta1, delta2, r
  -0.15487177, -0.019568024, -0.9171713, -4.04e-04, -1.67e-04  # w1-w5
)
names(YASSO15_DEFAULT_PARAMS) <- YASSO15_PARAM_NAMES


# -----------------------------------------------------------------------------
# Stage 1: Xi computation (R, vectorised, three pool groups)
# -----------------------------------------------------------------------------

#' Compute Yasso15 climate modifiers
#'
#' Returns three xi vectors (AWE, N, H), one per pool group, using
#' pool-group-specific beta and gamma parameters.
#' Computed in double precision in R -- matches reference behaviour.
#'
#' Can be called with vectors (one value per year) or scalars (mean climate).
#'
#' @param temp_mean   Numeric vector. Mean annual temperature (°C).
#' @param temp_amp    Numeric vector. Annual temperature amplitude (°C).
#' @param precip      Numeric vector. Annual precipitation (mm).
#' @param params      Numeric vector (length 35).
#' @return List with elements xi_awe, xi_n, xi_h (each same length as inputs).
compute_xi_yasso15 <- function(temp_mean, temp_amp, precip, params) {
  
  .xi_one_group <- function(temp_mean, temp_amp, precip, beta1, beta2, gamma) {
    sq2 <- sqrt(2)
    T1 <- temp_mean + 4 * temp_amp / pi * (1 / sq2 - 1)
    T2 <- temp_mean - 4 * temp_amp / (sq2 * pi)
    T3 <- temp_mean + 4 * temp_amp / pi * (1 - 1 / sq2)
    T4 <- temp_mean + 4 * temp_amp / (sq2 * pi)
    
    temp_mod <- (exp(beta1 * T1 + beta2 * T1^2) +
                   exp(beta1 * T2 + beta2 * T2^2) +
                   exp(beta1 * T3 + beta2 * T3^2) +
                   exp(beta1 * T4 + beta2 * T4^2)) / 4
    
    temp_mod * (1 - exp(gamma * precip / 1000))
  }
  
  list(
    xi_awe = .xi_one_group(temp_mean, temp_amp, precip,
                           params[17], params[18], params[19]),
    xi_n   = .xi_one_group(temp_mean, temp_amp, precip,
                           params[20], params[21], params[22]),
    xi_h   = .xi_one_group(temp_mean, temp_amp, precip,
                           params[23], params[24], params[25])
  )
}


# -----------------------------------------------------------------------------
# Stage 1b: Mean-climate xi for steady-state initialisation (R)
# -----------------------------------------------------------------------------
#
# IMPORTANT: steady-state xi must be xi(mean_climate), NOT mean(xi_annual).
# Because xi is nonlinear in T and P (Jensen's inequality), averaging xi over
# years gives a biased estimate of xi at mean climate. This function takes the
# steady-state climate window (clim_ss), averages the three climate scalars,
# then calls compute_xi_yasso15 once with scalar inputs, returning a list of
# three scalars (xi_awe, xi_n, xi_h).
#
# Mirrors compute_xi_mean_yasso07 in yasso07_wrapper.R. Called from the
# engine via the compute_xi_mean argument to make_likelihood().
# -----------------------------------------------------------------------------

compute_xi_mean_yasso15 <- function(clim_ss, params) {
  compute_xi_yasso15(
    temp_mean = mean(clim_ss$temp_mean),
    temp_amp  = mean(clim_ss$temp_amplitude),
    precip    = mean(clim_ss$precip),
    params    = params
  )
}


# -----------------------------------------------------------------------------
# Stage 2: Steady-state initialisation (Fortran)
# -----------------------------------------------------------------------------

#' Compute steady-state initial carbon pools for Yasso15
#'
#' Accepts pre-computed xi values (from compute_xi_mean_yasso15) rather than
#' raw climate scalars. The Fortran entry point has always taken pre-computed
#' xi; this refactor makes the R interface match that contract and keeps the
#' engine's compute_xi_mean -> steady_state call pattern uniform across models.
#'
#' precip_mean is still required: the Fortran uses it separately for the
#' leaching term and it cannot be derived from xi alone.
#'
#' @param params      Numeric vector (length 35).
#' @param nwl_mean    Numeric vector (length 4). Mean AWEN inputs, non-woody.
#' @param fwl_mean    Numeric vector (length 4). Mean AWEN inputs, fine woody.
#' @param cwl_mean    Numeric vector (length 4). Mean AWEN inputs, coarse woody.
#' @param xi_ss       Named list with xi_awe, xi_n, xi_h (scalars) from
#'                    compute_xi_mean_yasso15().
#' @param precip_mean Scalar. Mean annual precipitation (mm) for leaching term.
#' @param leac        Scalar. Leaching parameter (site-level).
#' @param diam_fwl    Fine woody litter diameter (cm).
#' @param diam_cwl    Coarse woody litter diameter (cm).
#' @return Numeric vector of length 15 (initial pool states, double precision).
yasso15_steady_state <- function(params,
                                 nwl_mean, fwl_mean, cwl_mean,
                                 xi_ss,
                                 precip_mean,
                                 leac     = 0.0,
                                 diam_fwl = 2.0, diam_cwl = 15.0) {
  
  result <- .Fortran("yasso15_steady_state_r",
                     params      = as.single(params),
                     nwl_mean    = as.single(nwl_mean),
                     fwl_mean    = as.single(fwl_mean),
                     cwl_mean    = as.single(cwl_mean),
                     xi_awe_mean = as.single(xi_ss$xi_awe),
                     xi_n_mean   = as.single(xi_ss$xi_n),
                     xi_h_mean   = as.single(xi_ss$xi_h),
                     leac        = as.single(leac),
                     precip_mean = as.single(precip_mean),
                     diam_fwl    = as.single(diam_fwl),
                     diam_cwl    = as.single(diam_cwl),
                     C_init      = single(15)
  )
  
  as.double(result$C_init)
}


# -----------------------------------------------------------------------------
# Stage 3: Transient forward simulation (Fortran)
# -----------------------------------------------------------------------------

#' Run Yasso15 transient forward simulation via Fortran
#'
#' All inputs are cast to single precision before the Fortran call to match
#' the original Järvenpää reference implementation. Output is cast back to
#' double for R-side use.
#'
#' @param input_df   Data frame with columns: year, nwl_A/W/E/N,
#'                   fwl_A/W/E/N, cwl_A/W/E/N.
#' @param params     Numeric vector (length 35).
#' @param C_init     Numeric vector (length 15) from yasso15_steady_state().
#' @param xi_arrays  List with xi_awe, xi_n, xi_h vectors (from compute_xi_yasso15).
#' @param precip     Numeric vector. Annual precipitation (mm).
#' @param leac       Scalar. Leaching parameter (site-level).
#' @param diam_fwl   Fine woody litter diameter (cm).
#' @param diam_cwl   Coarse woody litter diameter (cm).
#' @return Data frame: year, A, W, E, N, H, total_soc, respiration.
yasso15_run <- function(input_df, params, C_init, xi_arrays,
                        precip,
                        leac = 0.0,
                        diam_fwl = 2.0, diam_cwl = 15.0) {
  
  n_years      <- nrow(input_df)
  precip_array <- precip
  nwl_awen <- as.matrix(input_df[, c("nwl_A", "nwl_W", "nwl_E", "nwl_N")])
  fwl_awen <- as.matrix(input_df[, c("fwl_A", "fwl_W", "fwl_E", "fwl_N")])
  cwl_awen <- as.matrix(input_df[, c("cwl_A", "cwl_W", "cwl_E", "cwl_N")])
  
  result <- .Fortran("yasso15_run_r",
                     n_years       = as.integer(n_years),
                     params        = as.single(params),
                     nwl_awen      = as.single(nwl_awen),
                     fwl_awen      = as.single(fwl_awen),
                     cwl_awen      = as.single(cwl_awen),
                     xi_awe_array  = as.single(xi_arrays$xi_awe),
                     xi_n_array    = as.single(xi_arrays$xi_n),
                     xi_h_array    = as.single(xi_arrays$xi_h),
                     leac          = as.single(leac),
                     precip_array  = as.single(precip_array),
                     diam_fwl      = as.single(diam_fwl),
                     diam_cwl      = as.single(diam_cwl),
                     C_init        = as.single(C_init),
                     C_out         = single(n_years * 5),
                     resp_out      = single(n_years),
                     C_final       = single(15)        # terminal per-cohort state
  )
  
  C_mat <- matrix(as.double(result$C_out), nrow = n_years, ncol = 5)
  colnames(C_mat) <- c("A", "W", "E", "N", "H")
  
  out <- data.frame(
    year        = input_df$year,
    C_mat,
    total_soc   = rowSums(C_mat),
    respiration = as.double(result$resp_out)
  )
  attr(out, "C_final") <- as.double(result$C_final)
  out
}

# =============================================================================
# yasso15_transient_init
#
# 68-year pre-run (1917 -> 1985), same structure as yasso07_transient_init.
# xi_mean is the list (xi_awe, xi_n, xi_h) from compute_xi_mean_yasso15.
# precip held constant at lm$precip_mean (inactive leaching term, leac = 0).
# =============================================================================

yasso15_transient_init <- function(model_params, lm, xi_mean, ...) {
  sigma_init  <- unname(model_params["sigma_init"])
  sigma_input <- unname(model_params["sigma_input"])
  params      <- model_params[YASSO15_PARAM_NAMES]

  # Fully-scaled AWEN litter endpoints
  nwl_1917 <- lm$nwl_full_mean * sigma_init * sigma_input
  fwl_1917 <- lm$fwl_full_mean * sigma_init * sigma_input
  cwl_1917 <- lm$cwl_full_mean * sigma_init * sigma_input
  nwl_1985 <- lm$nwl_t0_mean   * sigma_input
  fwl_1985 <- lm$fwl_t0_mean   * sigma_input
  cwl_1985 <- lm$cwl_t0_mean   * sigma_input

  # Start at steady state under 1917 litter (15-element per-cohort state)
  C_init <- yasso15_steady_state(params      = params,
                                 nwl_mean    = nwl_1917,
                                 fwl_mean    = fwl_1917,
                                 cwl_mean    = cwl_1917,
                                 xi_ss       = xi_mean,
                                 precip_mean = lm$precip_mean)

  # Build 68-row input_df with linearly interpolated AWEN columns
  n_pre <- 68L
  fracs <- (seq_len(n_pre) - 1L) / (n_pre - 1L)
  nwl_mat <- outer(1 - fracs, nwl_1917) + outer(fracs, nwl_1985)
  fwl_mat <- outer(1 - fracs, fwl_1917) + outer(fracs, fwl_1985)
  cwl_mat <- outer(1 - fracs, cwl_1917) + outer(fracs, cwl_1985)
  input_df <- data.frame(year = seq_len(n_pre), nwl_mat, fwl_mat, cwl_mat)
  names(input_df) <- c("year",
                       "nwl_A","nwl_W","nwl_E","nwl_N",
                       "fwl_A","fwl_W","fwl_E","fwl_N",
                       "cwl_A","cwl_W","cwl_E","cwl_N")

  # Constant xi and precip arrays (no pre-1985 climate observations)
  xi_const     <- list(xi_awe = rep(xi_mean$xi_awe, n_pre),
                       xi_n   = rep(xi_mean$xi_n,   n_pre),
                       xi_h   = rep(xi_mean$xi_h,   n_pre))
  precip_const <- rep(lm$precip_mean, n_pre)

  # Run pre-run; C_final carries the full 15-element terminal state
  pre_out <- yasso15_run(input_df  = input_df,
                         params    = params,
                         C_init    = C_init,
                         xi_arrays = xi_const,
                         precip    = precip_const)
  attr(pre_out, "C_final")
}