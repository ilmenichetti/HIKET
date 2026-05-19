# =============================================================================
# yasso20_wrapper_transient.R
#
# Standalone transient wrapper for YASSO20.
# Complete copy of yasso20_wrapper.R with yasso*_transient_init appended.
# Does NOT source the original wrapper -- self-contained.
# =============================================================================

# =============================================================================
# Yasso20 R Wrapper
# =============================================================================
#
# Reuses yasso15.f90 compiled as yasso15.so -- no separate Fortran needed.
# The only difference from Yasso15 is compute_xi_yasso20(), which averages
# over all 12 monthly temperatures directly instead of the four-point Gaussian
# approximation used in Yasso07/15.
#
# Three-stage call pattern:
#   1. compute_xi_yasso20()    -- climate modifier time series (R)
#   2. yasso20_steady_state()  -- initial pool state C* (Fortran, via yasso15)
#   3. yasso20_run()           -- transient forward simulation (Fortran, via yasso15)
#
# Input: monthly_climate from map_climate_monthly(), which provides
#   columns: plot_id, year, month, temp_air, precip
#
# Inputs: Yasso20_inputs from map_inputs_yasso20(), which provides
#   columns: year, nwl_A/W/E/N, fwl_A/W/E/N, cwl_A/W/E/N
# name change
# =============================================================================

# Parameter vector: identical 35-element format as Yasso15
# (see yasso15_wrapper.R for full index documentation)
YASSO20_PARAM_NAMES   <- YASSO15_PARAM_NAMES
YASSO20_DEFAULT_PARAMS <- YASSO15_DEFAULT_PARAMS


# -----------------------------------------------------------------------------
# Stage 1: Xi computation (R, vectorised, three pool groups)
# -----------------------------------------------------------------------------

#' Compute Yasso20 climate modifiers
#'
#' Averages the temperature response over all 12 monthly temperatures
#' rather than the four-point Gaussian approximation used in Yasso15.
#' Precipitation modifier uses the annual total (sum of monthly precip).
#'
#' @param climate_df  Data frame with columns: year, month, temp_air, precip.
#'                    Long monthly format (12 rows per year).
#' @param params      Numeric vector (length 35).
#' @return List with xi_awe, xi_n, xi_h (each length n_years).
compute_xi_yasso20 <- function(climate_df, params) {
  
  .xi_one_group <- function(temp_monthly, precip_annual, beta1, beta2, gamma) {
    # Temperature modifier: mean of exp(beta1*T + beta2*T^2) over 12 months
    temp_mod <- mean(exp(beta1 * temp_monthly + beta2 * temp_monthly^2))
    # Precipitation modifier: uses annual total (mm -> m)
    temp_mod * (1 - exp(gamma * precip_annual / 1000))
  }
  
  years        <- sort(unique(climate_df$year))
  precip_annual <- as.numeric(tapply(climate_df$precip, climate_df$year, sum, na.rm = TRUE))
  
  xi_awe <- xi_n <- xi_h <- numeric(length(years))
  
  for (i in seq_along(years)) {
    yr      <- years[i]
    T_month <- climate_df$temp_air[climate_df$year == yr]
    PR      <- precip_annual[i]
    
    xi_awe[i] <- .xi_one_group(T_month, PR, params[17], params[18], params[19])
    xi_n[i]   <- .xi_one_group(T_month, PR, params[20], params[21], params[22])
    xi_h[i]   <- .xi_one_group(T_month, PR, params[23], params[24], params[25])
  }
  
  list(xi_awe = xi_awe, xi_n = xi_n, xi_h = xi_h)
}

# -----------------------------------------------------------------------------
# Stage 1b: Mean-climate xi for steady-state initialisation (R)
# -----------------------------------------------------------------------------
#
# Builds a synthetic single mean-year from the monthly climate window clim_ss,
# then calls compute_xi_yasso20 once to get scalar xi values for each pool
# group. The mean year has 12 rows (one per calendar month), with:
#   temp_air : per-month mean across all years in clim_ss
#   precip   : mean monthly value across all rows (summed to annual total
#              inside compute_xi_yasso20 via tapply)
#
# IMPORTANT: xi(mean_climate) != mean(xi_annual). By averaging the monthly
# temperatures first and calling compute_xi once, we correctly compute xi at
# mean climate rather than averaging xi across years.
#
# Mirrors compute_xi_mean_yasso15 in yasso15_wrapper.R.
# Called from the engine via the compute_xi_mean argument to make_likelihood().
# -----------------------------------------------------------------------------

compute_xi_mean_yasso20 <- function(clim_ss, params) {
  mean_clim <- data.frame(
    year     = 9999L,
    month    = 1:12,
    temp_air = as.numeric(tapply(clim_ss$temp_air, clim_ss$month,
                                 mean, na.rm = TRUE)),
    precip   = mean(clim_ss$precip)   # monthly mean; summed to annual inside compute_xi_yasso20
  )
  compute_xi_yasso20(mean_clim, params)
}

# -----------------------------------------------------------------------------
# Stage 2: Steady-state initialisation (Fortran, via yasso15)
# -----------------------------------------------------------------------------

#' Compute steady-state initial carbon pools for Yasso20
#'
#' Accepts pre-computed xi values (from compute_xi_mean_yasso20) and
#' precip_mean, then delegates to the refactored yasso15_steady_state().
#' The internal xi and precip computation that previously lived here has
#' moved to compute_xi_mean_yasso20() and the run script's litter_means
#' construction, keeping the engine's compute_xi_mean -> steady_state
#' call pattern uniform across all three models.
#'
#' @param params      Numeric vector (length 35).
#' @param nwl_mean    Numeric vector (length 4). Mean AWEN inputs, non-woody.
#' @param fwl_mean    Numeric vector (length 4). Mean AWEN inputs, fine woody.
#' @param cwl_mean    Numeric vector (length 4). Mean AWEN inputs, coarse woody.
#' @param xi_ss       Named list with xi_awe, xi_n, xi_h (scalars) from
#'                    compute_xi_mean_yasso20().
#' @param precip_mean Scalar. Mean annual precipitation (mm) for leaching term.
#'                    For Yasso20: sum of mean monthly values over clim_ss.
#' @param leac        Scalar. Leaching parameter (site-level).
#' @param diam_fwl    Fine woody litter diameter (cm).
#' @param diam_cwl    Coarse woody litter diameter (cm).
#' @return Numeric vector of length 15 (initial pool states, double precision).
yasso20_steady_state <- function(params,
                                 nwl_mean, fwl_mean, cwl_mean,
                                 xi_ss,
                                 precip_mean,
                                 leac     = 0.0,
                                 diam_fwl = 2.0, diam_cwl = 15.0) {
  
  yasso15_steady_state(
    params      = params,
    nwl_mean    = nwl_mean,
    fwl_mean    = fwl_mean,
    cwl_mean    = cwl_mean,
    xi_ss       = xi_ss,
    precip_mean = precip_mean,
    leac        = leac,
    diam_fwl    = diam_fwl,
    diam_cwl    = diam_cwl
  )
}


# -----------------------------------------------------------------------------
# Stage 3: Transient forward simulation (Fortran, via yasso15)
# -----------------------------------------------------------------------------

#' Run Yasso20 transient forward simulation
#'
#' Computes xi and annual precip from monthly climate internally,
#' then delegates to the Yasso15 Fortran run routine.
#'
#' @param climate_df  Data frame with columns: year, month, temp_air, precip.
#' @param inputs_df   Data frame from map_inputs_yasso20().
#'                    Columns: year, nwl_A/W/E/N, fwl_A/W/E/N, cwl_A/W/E/N.
#' @param params      Numeric vector (length 35).
#' @param C_init      Numeric vector (length 15) from yasso20_steady_state().
#' @param leac        Scalar. Leaching parameter.
#' @param diam_fwl    Fine woody litter diameter (cm).
#' @param diam_cwl    Coarse woody litter diameter (cm).
#' @return Data frame: year, A, W, E, N, H, total_soc, respiration.
yasso20_run <- function(climate_df, inputs_df, params, C_init,
                        leac = 0.0,
                        diam_fwl = 2.0, diam_cwl = 15.0) {
  
  xi_arrays <- compute_xi_yasso20(climate_df, params)
  precip    <- as.numeric(tapply(climate_df$precip, climate_df$year, sum, na.rm = TRUE))
  
  yasso15_run(
    input_df  = inputs_df,
    params    = params,
    C_init    = C_init,
    xi_arrays = xi_arrays,
    precip    = precip,
    leac      = leac,
    diam_fwl  = diam_fwl,
    diam_cwl  = diam_cwl
  )
}
# =============================================================================
# yasso20_transient_init
# See yasso07_transient_init for full rationale.
# Yasso20 shares the same 15-element per-cohort C_init structure.
# =============================================================================

yasso20_transient_init <- function(model_params, lm, xi_mean, ...) {
  sigma_init  <- unname(model_params["sigma_init"])
  sigma_input <- unname(model_params["sigma_input"])
  
  yasso20_steady_state(
    params      = model_params[YASSO20_PARAM_NAMES],
    nwl_mean    = lm$nwl_mean * sigma_input * sigma_init,
    fwl_mean    = lm$fwl_mean * sigma_input * sigma_init,
    cwl_mean    = lm$cwl_mean * sigma_input * sigma_init,
    xi_ss       = xi_mean,
    precip_mean = lm$precip_mean
  )
}
