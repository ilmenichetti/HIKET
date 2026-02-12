# =============================================================================
# Yasso20 R Wrapper
# =============================================================================
#
# Reuses yasso15.f90 compiled as yasso15.so — no separate Fortran needed.
# The only difference from Yasso15 is compute_xi_yasso20(), which averages
# over 12 monthly temperatures directly instead of the four-point Gaussian
# approximation.
#
# Three-stage call pattern:
#   1. compute_xi_yasso20()    -- climate modifier time series (R)
#   2. yasso20_steady_state()  -- initial pool state C* (Fortran, via yasso15)
#   3. yasso20_run()           -- transient forward simulation (Fortran, via yasso15)
#
# Climate input: monthly_climate from map_climate_monthly()
#   columns: plot_id, year, month, temp_air, precip, evap
#
# Litter input: Yasso20_inputs from map_inputs_yasso20()
#   columns: plot_id, year, nwl_A/W/E/N, fwl_A/W/E/N, cwl_A/W/E/N
#
# =============================================================================

YASSO20_PARAM_NAMES    <- YASSO15_PARAM_NAMES
YASSO20_DEFAULT_PARAMS <- YASSO15_DEFAULT_PARAMS


# -----------------------------------------------------------------------------
# Stage 1: Xi computation
# -----------------------------------------------------------------------------

#' Compute Yasso20 climate modifiers
#'
#' Sums the Gaussian temperature response over all 12 monthly temperatures
#' per year, then applies the precipitation modifier. Matches the original
#' Yasso20 Fortran exactly.
#'
#' @param climate_df  Data frame from map_climate_monthly().
#'                    Columns: year, month, temp_air, precip.
#' @param params      Numeric vector (length 35).
#' @return List with xi_awe, xi_n, xi_h (each a numeric vector, length n_years).
compute_xi_yasso20 <- function(climate_df, params) {
  
  if (!"month" %in% names(climate_df))
    stop("climate_df must be long monthly format (year, month, temp_air, precip)", call. = FALSE)
  
  years <- sort(unique(climate_df$year))
  n_yrs <- length(years)
  
  # Build n_years x 12 temperature matrix
  Tmat <- matrix(NA_real_, nrow = n_yrs, ncol = 12,
                 dimnames = list(as.character(years), paste0("T_", 1:12)))
  for (y in years)
    Tmat[as.character(y), climate_df$month[climate_df$year == y]] <-
    climate_df$temp_air[climate_df$year == y]
  if (anyNA(Tmat)) stop("Missing monthly temperatures in climate_df", call. = FALSE)
  
  # Annual precip (sum over 12 months)
  precip <- as.numeric(tapply(climate_df$precip, climate_df$year,
                              sum, na.rm = TRUE)[as.character(years)])
  
  .xi_one <- function(beta1, beta2, gamma)
    rowSums(exp(beta1 * Tmat + beta2 * Tmat^2)) * (1 - exp(gamma * precip / 1000)) / 12
  
  list(
    xi_awe = .xi_one(params[17], params[18], params[19]),
    xi_n   = .xi_one(params[20], params[21], params[22]),
    xi_h   = .xi_one(params[23], params[24], params[25])
  )
}


# -----------------------------------------------------------------------------
# Stage 2: Steady-state initialisation
# -----------------------------------------------------------------------------

#' Compute steady-state initial carbon pools for Yasso20
#'
#' @param climate_df   Data frame from map_climate_monthly(). Used to derive
#'                     mean xi and mean annual precip for steady-state.
#' @param params       Numeric vector (length 35).
#' @param nwl_mean     Numeric vector (length 4). Mean AWEN inputs, non-woody.
#' @param fwl_mean     Numeric vector (length 4). Mean AWEN inputs, fine woody.
#' @param cwl_mean     Numeric vector (length 4). Mean AWEN inputs, coarse woody.
#' @param leac         Scalar. Leaching parameter.
#' @param diam_fwl     Fine woody litter diameter (cm).
#' @param diam_cwl     Coarse woody litter diameter (cm).
#' @return Numeric vector of length 15.
yasso20_steady_state <- function(climate_df, params,
                                 nwl_mean, fwl_mean, cwl_mean,
                                 leac = 0.0,
                                 diam_fwl = 2.0, diam_cwl = 15.0) {
  
  xi_arrays   <- compute_xi_yasso20(climate_df, params)
  xi_means    <- list(
    xi_awe = mean(xi_arrays$xi_awe),
    xi_n   = mean(xi_arrays$xi_n),
    xi_h   = mean(xi_arrays$xi_h)
  )
  precip_mean <- mean(tapply(climate_df$precip, climate_df$year, sum, na.rm = TRUE))
  
  yasso15_steady_state(
    params      = params,
    nwl_mean    = nwl_mean,
    fwl_mean    = fwl_mean,
    cwl_mean    = cwl_mean,
    xi_means    = xi_means,
    leac        = leac,
    precip_mean = precip_mean,
    diam_fwl    = diam_fwl,
    diam_cwl    = diam_cwl
  )
}


# -----------------------------------------------------------------------------
# Stage 3: Transient forward simulation
# -----------------------------------------------------------------------------

#' Run Yasso20 transient forward simulation
#'
#' @param climate_df  Data frame from map_climate_monthly().
#'                    Used to derive xi_arrays and annual precip.
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
  precip    <- as.numeric(tapply(climate_df$precip, climate_df$year,
                                 sum, na.rm = TRUE))
  
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