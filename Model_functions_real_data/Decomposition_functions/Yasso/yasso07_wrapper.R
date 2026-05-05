# =============================================================================
# Yasso07 R Wrapper for Fortran Implementation
# =============================================================================
#
# Compile:  R CMD SHLIB yasso07.f90
# Load:     dyn.load("yasso07.so")
#
# Three-stage call pattern:
#   1. compute_xi_yasso07()    -- climate modifier time series (R)
#   2. yasso07_steady_state()  -- initial pool state C* (Fortran)
#   3. yasso07_run()           -- transient forward simulation (Fortran)
#
# =============================================================================

YASSO07_PARAM_NAMES <- c(
  "alpha_A", "alpha_W", "alpha_E", "alpha_N",
  "p_WA", "p_EA", "p_NA",
  "p_AW", "p_EW", "p_NW",
  "p_AE", "p_WE", "p_NE",
  "p_AN", "p_WN", "p_EN",
  "beta1", "beta2",
  "gamma",
  "delta1", "delta2", "r",
  "p_H",
  "alpha_H"
)

YASSO07_DEFAULT_PARAMS <- c(
  0.7035942674, 5.681055546, 0.2613542378, 0.02810959704,  # alpha_A/W/E/N (positive; Fortran negates)
  0.4888527989, 0.01905768365, 0.9696374536,               # p_WA, p_EA, p_NA
  0.9872559905, 0.002843263559, 0.003396461252,            # p_AW, p_EW, p_NW
  1.399370376e-05, 1.796692413e-05, 0.01218125969,         # p_AE, p_WE, p_NE
  0.002777846763, 0.01269555371, 0.9713827968,             # p_AN, p_WN, p_EN
  0.09873183817, -0.001571640489, -1.271691799,            # beta1, beta2, gamma
  -1.708411336, 0.8585553765, -0.3068014085,               # delta1, delta2, r
  0.004270370584, 0.001496617449                           # p_H, alpha_H (positive; Fortran negates)
)
names(YASSO07_DEFAULT_PARAMS) <- YASSO07_PARAM_NAMES


# -----------------------------------------------------------------------------
# Stage 1: Xi computation (R, vectorised)
# -----------------------------------------------------------------------------

compute_xi_yasso07 <- function(temp_mean, temp_amp, precip,
                               beta1, beta2, gamma) {
  sq2 <- sqrt(2)
  T1 <- temp_mean + 4 * temp_amp / pi * (1 / sq2 - 1)
  T2 <- temp_mean - 4 * temp_amp / (sq2 * pi)
  T3 <- temp_mean + 4 * temp_amp / pi * (1 - 1 / sq2)
  T4 <- temp_mean + 4 * temp_amp / (sq2 * pi)
  
  temp_mod <- (exp(beta1 * T1 + beta2 * T1^2) +
                 exp(beta1 * T2 + beta2 * T2^2) +
                 exp(beta1 * T3 + beta2 * T3^2) +
                 exp(beta1 * T4 + beta2 * T4^2)) / 4
  
  precip_mod <- 1 - exp(gamma * precip / 1000)
  
  temp_mod * precip_mod
}


# -----------------------------------------------------------------------------
# Stage 2: Steady-state initialisation (Fortran)
# -----------------------------------------------------------------------------

yasso07_steady_state <- function(params,
                                 nwl_mean, fwl_mean, cwl_mean,
                                 xi_mean,
                                 diam_fwl = 2.0, diam_cwl = 15.0) {
  
  result <- .Fortran("yasso07_steady_state_r",
                     params   = as.double(params),
                     nwl_mean = as.double(nwl_mean),
                     fwl_mean = as.double(fwl_mean),
                     cwl_mean = as.double(cwl_mean),
                     xi_mean  = as.double(xi_mean),
                     diam_fwl = as.double(diam_fwl),
                     diam_cwl = as.double(diam_cwl),
                     C_init   = double(15)
  )
  
  result$C_init
}


# -----------------------------------------------------------------------------
# Stage 3: Transient forward simulation (Fortran)
# -----------------------------------------------------------------------------

yasso07_run <- function(input_df, params, C_init, xi_array,
                        diam_fwl = 2.0, diam_cwl = 15.0) {
  
  n_years  <- nrow(input_df)
  nwl_awen <- as.matrix(input_df[, c("nwl_A", "nwl_W", "nwl_E", "nwl_N")])
  fwl_awen <- as.matrix(input_df[, c("fwl_A", "fwl_W", "fwl_E", "fwl_N")])
  cwl_awen <- as.matrix(input_df[, c("cwl_A", "cwl_W", "cwl_E", "cwl_N")])
  
  result <- .Fortran("yasso07_run_r",
                     n_years  = as.integer(n_years),
                     params   = as.double(params),
                     nwl_awen = as.double(nwl_awen),
                     fwl_awen = as.double(fwl_awen),
                     cwl_awen = as.double(cwl_awen),
                     xi_array = as.double(xi_array),
                     diam_fwl = as.double(diam_fwl),
                     diam_cwl = as.double(diam_cwl),
                     C_init   = as.double(C_init),
                     C_out    = double(n_years * 5),
                     resp_out = double(n_years)
  )
  
  C_mat <- matrix(result$C_out, nrow = n_years, ncol = 5)
  colnames(C_mat) <- c("A", "W", "E", "N", "H")
  
  data.frame(
    year        = input_df$year,
    C_mat,
    total_soc   = rowSums(C_mat),
    respiration = result$resp_out
  )
}



# -----------------------------------------------------------------------------
# Stage 1b: Mean-climate xi for steady-state initialisation (R)
# -----------------------------------------------------------------------------
#
# IMPORTANT (Apr 2026 fix):
#   The steady-state climate modifier must be xi(mean_climate), NOT
#   mean(xi_annual). Because xi is convex in T (via exp) and concave in P
#   (via 1 - exp(-gamma*P/1000)), Jensen's inequality gives:
#       mean(xi_annual) >= xi(mean_climate)
#   for typical Finnish inter-annual variability, with bias 5-15%.
#
#   This function takes a subset of climate data (the steady-state window)
#   and returns the SCALAR xi computed from the mean of each climate variable
#   over that window. compute_xi_yasso07 already handles the within-year
#   sinusoidal temperature decomposition (T1-T4) correctly, so calling it
#   once with mean inputs is exactly what we want.
#
#   Called from inside the calibration likelihood (Stage 2, steady_state
#   initialisation). NOT called from the transient simulation, which uses
#   the full xi_array from compute_xi_yasso07 over all years.
# -----------------------------------------------------------------------------

compute_xi_mean_yasso07 <- function(clim_ss, beta1, beta2, gamma) {
  compute_xi_yasso07(
    temp_mean = mean(clim_ss$temp_mean),
    temp_amp  = mean(clim_ss$temp_amplitude),
    precip    = mean(clim_ss$precip),
    beta1     = beta1,
    beta2     = beta2,
    gamma     = gamma
  )
}
