# =============================================================================
# sp1_wrapper_transient.R
#
# Standalone transient wrapper for SP1.
# Complete copy of sp1_wrapper.R with compute_xi functions (from Yasso07)
# and sp1_transient_init appended. Does NOT source any other file.
# =============================================================================

# -----------------------------------------------------------------------------
# Climate modifier functions (identical to Yasso07; included here to avoid
# sourcing yasso07_wrapper.R -- kept standalone)
# -----------------------------------------------------------------------------



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

# =============================================================================
# sp1_wrapper.R
#
# Single-Pool SOC decomposition model (SP1) for the HIKET intercomparison.
#
# STRUCTURE:
#   One carbon pool C receiving all litter inputs (NWL + FWL + CWL, all AWEN
#   fractions collapsed to a single scalar flux J). Decomposition follows a
#   first-order rate equation modified by the climate scalar xi:
#
#     dC/dt = J(t) - alpha * xi(t) * C(t)
#
#   This is the simplest possible linear SOC model. It serves as the
#   low-complexity anchor of the HIKET intercomparison (SP1 -> TP2 -> Yasso07).
#
# CLIMATE MODIFIER:
#   Identical formula and parameters (beta1, beta2, gamma) to Yasso07.
#   xi depends on the seasonal temperature decomposition from annual mean
#   temperature and amplitude:
#
#     T1 = T_mean + T_amp   (warmest season, ~July)
#     T2 = T_mean           (autumn)
#     T3 = T_mean - T_amp   (coldest season, ~January)
#     T4 = T_mean           (spring)
#
#     xi = exp( beta1*(T1+T2+T3+T4) + beta2*(T1^2+T2^2+T3^2+T4^2) )
#          * (1 - exp(gamma * P))
#
#   Substituting the T definitions:
#     T1+T2+T3+T4 = 4*T_mean
#     T1^2+T2^2+T3^2+T4^2 = 4*T_mean^2 + 2*T_amp^2
#
#   So: xi = exp( 4*beta1*T_mean + beta2*(4*T_mean^2 + 2*T_amp^2) )
#            * (1 - exp(gamma * P))
#
#   This matches the Yasso07 formulation (Tuomi et al. 2009, Eq. 4), preserving
#   the seasonal amplitude's contribution through the beta2 quadratic term.
#   Using the same functional form enables fair comparison of model STRUCTURE
#   independent of climate response -- xi parameters are independently
#   calibrated per model.
#
# TRANSIENT SIMULATION:
#   Each annual timestep solved exactly (no approximation beyond floating point):
#     C(t+1) = C_ss(t) + (C(t) - C_ss(t)) * exp(-alpha * xi(t))
#   where C_ss(t) = J(t) / (alpha * xi(t)) is the instantaneous steady state.
#
# PARAMETERS (all free in calibration):
#   alpha       : decomposition rate (yr^-1). Strictly positive (log transform).
#   beta1       : temperature response, linear coefficient. Strictly positive.
#   beta2       : temperature response, quadratic coefficient. Unconstrained.
#   gamma       : precipitation response. Unconstrained (sign enforced by prior).
#   sigma_init  : initial-condition scaling (nuisance). Log-transformed.
#   sigma_input : global litter input multiplier (nuisance). Log-transformed.
#
# ENGINE INTERFACE:
#   sp1_steady_state(model_params, lm, xi_mean)  -> named vector c(C = C_ss)
#   sp1_run(inputs, model_params, C_init, xi_array) -> data.frame with year,
#                                                       C, total_soc
#   compute_xi_sp1(temp_mean, temp_amp, precip, beta1, beta2, gamma) -> numeric
#   compute_xi_mean_sp1(clim_ss, model_params) -> scalar xi
#
# DEPENDENCIES: base R only (no Fortran, no external packages).
# =============================================================================


# Parameter name vector in the order expected by assemble_model_params() in
# the calibration script. Included here for reference and cross-validation.
SP1_PARAM_NAMES <- c("alpha", "beta1", "beta2", "gamma",
                     "sigma_init", "sigma_input")


# =============================================================================
# Steady-state initialisation
# =============================================================================

# Compute the analytical steady-state pool size C* for a given set of
# parameters, mean litter input, and mean climate modifier.
#
# From dC/dt = J - alpha*xi*C = 0:
#   C* = J / (alpha * xi_mean)
#
# sigma_input is applied to J here so that the calibrated litter multiplier
# is consistent between steady-state initialisation and the transient run.
#
# ARGUMENTS:
#   model_params : named numeric -- must contain alpha, sigma_input
#   lm           : named list -- must contain $J_total_mean (scalar, tC/ha/yr)
#   xi_mean      : scalar -- climate modifier at mean climate (from compute_xi_mean_sp1)
#
# RETURNS:
#   named numeric vector c(C = C_ss) for consumption by sp1_run()
sp1_steady_state <- function(model_params, lm, xi_mean) {
  alpha <- model_params["alpha"]
  si    <- model_params["sigma_input"]
  J     <- lm$J_total_mean * si
  C_ss  <- J / (alpha * xi_mean)
  c(C = unname(C_ss))
}


# =============================================================================
# Transient simulation
# =============================================================================

# Run the SP1 model from C_init over the full input time series.
#
# Each annual step uses the exact per-step solution of the linear ODE:
#   C(t+1) = C_ss(t) + (C(t) - C_ss(t)) * exp(-alpha * xi(t))
# where C_ss(t) = J(t)/(alpha*xi(t)) is the (time-varying) steady state.
# This is exact for piecewise-constant J and xi within each year (dt = 1).
#
# sigma_input scales litter inputs, consistent with steady_state above.
#
# ARGUMENTS:
#   inputs      : data.frame -- annual inputs, must have columns year, J_total
#   model_params : named numeric -- alpha, sigma_input, and xi parameters
#   C_init      : named numeric -- initial pool sizes, from sp1_steady_state()
#   xi_array    : numeric vector -- xi(t) for each row of inputs, from compute_xi_sp1()
#
# RETURNS:
#   data.frame with columns: year, C, total_soc
#   (total_soc = C for consistency with the multi-pool engine interface)
sp1_run <- function(inputs, model_params, C_init, xi_array) {
  alpha <- model_params["alpha"]
  si    <- model_params["sigma_input"]
  n     <- nrow(inputs)

  C    <- numeric(n + 1L)    # C[1] = C_init; C[2..n+1] = post-step values
  C[1] <- C_init["C"]

  for (t in seq_len(n)) {
    xi    <- xi_array[t]
    J     <- inputs$J_total[t] * si
    k     <- alpha * xi          # effective decay rate this year
    C_ss  <- J / k               # instantaneous steady state under this year's conditions
    C[t + 1L] <- C_ss + (C[t] - C_ss) * exp(-k)
  }

  data.frame(
    year      = inputs$year,
    C         = C[-1L],
    total_soc = C[-1L]    # single pool: total = C
  )
}

# =============================================================================
# sp1_transient_init
#
# Transient-initialization replacement for sp1_steady_state.
# sigma_init scales the litter level used for the steady-state computation,
# replacing the original likelihood patch (sigma_init in sd_vec).
#   sigma_init = 1.0 : initial SOC at steady state under contemporary litter
#   sigma_init > 1   : historically more productive -> larger initial C stock
#   sigma_init < 1   : historically less productive -> smaller initial C stock
# =============================================================================

sp1_transient_init <- function(model_params, lm, xi_mean, ...) {
  # Prevent name propagation from compute_xi_yasso07() arithmetic chain
  xi_mean <- unname(xi_mean)

  alpha       <- unname(model_params["alpha"])
  sigma_init  <- unname(model_params["sigma_init"])
  sigma_input <- unname(model_params["sigma_input"])

  J    <- lm$J_total_mean * sigma_input * sigma_init
  C_ss <- J / (alpha * xi_mean)
  c(C = unname(C_ss))
}
