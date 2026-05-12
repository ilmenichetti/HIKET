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
