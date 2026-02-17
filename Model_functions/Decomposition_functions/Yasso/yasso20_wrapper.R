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
#
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
# Stage 2: Steady-state initialisation (Fortran, via yasso15)
# -----------------------------------------------------------------------------

#' Compute steady-state initial carbon pools for Yasso20
#'
#' Builds a mean-climate year from monthly_climate and computes xi from it,
#' then delegates to the Yasso15 Fortran steady-state routine.
#'
#' @param climate_df  Data frame with columns: year, month, temp_air, precip.
#' @param params      Numeric vector (length 35).
#' @param nwl_mean    Numeric vector (length 4). Mean AWEN inputs, non-woody.
#' @param fwl_mean    Numeric vector (length 4). Mean AWEN inputs, fine woody.
#' @param cwl_mean    Numeric vector (length 4). Mean AWEN inputs, coarse woody.
#' @param leac        Scalar. Leaching parameter (site-level).
#' @param diam_fwl    Fine woody litter diameter (cm).
#' @param diam_cwl    Coarse woody litter diameter (cm).
#' @return Numeric vector of length 15 (initial pool states).
yasso20_steady_state <- function(climate_df, params,
                                 nwl_mean, fwl_mean, cwl_mean,
                                 leac = 0.0,
                                 diam_fwl = 2.0, diam_cwl = 15.0) {
  
  # Build synthetic mean year: monthly mean temperature, mean annual precip
  # This is the Yasso20 equivalent of (mean_temp, mean_temp_amp, mean_precip)
  mean_clim <- data.frame(
    year     = 9999L,
    month    = 1:12,
    temp_air = tapply(climate_df$temp_air, climate_df$month, mean, na.rm = TRUE),
    precip   = mean(climate_df$precip)   # monthly mean, summed to annual inside compute_xi
  )
  
  xi_ss        <- compute_xi_yasso20(mean_clim, params)
  precip_mean  <- sum(mean_clim$precip)   # annual total from mean monthly values
  
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