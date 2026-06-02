# =============================================================================
# tp2_wrapper_transient.R
#
# Standalone transient wrapper for TP2.
# Complete copy of tp2_wrapper.R with compute_xi functions (from Yasso07)
# and tp2_transient_init appended. Does NOT source any other file.
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
# tp2_wrapper.R
#
# Two-Pool SOC decomposition model (TP2) for the HIKET intercomparison.
#
# STRUCTURE (ICBM-style):
#   Two pools: A (active/fast) and H (humus/slow). All litter inputs enter
#   the A pool; the H pool receives only transferred mass from A. This is
#   the ICBM topology (Andren & Katterer 1997), adapted to use Yasso07
#   parameter naming and the same xi functional form.
#
#     dA/dt = J(t) - alpha_A * xi(t) * A(t)
#     dH/dt = p_H * alpha_A * xi(t) * A(t) - alpha_H * xi(t) * H(t)
#
#   Notation follows Yasso07 directly: alpha_A, alpha_H are pool decay rates
#   and p_H is the humification fraction (fraction of A decomposition routed
#   to H). The TP2 parameter set is therefore a strict SUBSET of Yasso07's
#   parameters, which makes the structural comparison clear.
#
# WHY ICBM TOPOLOGY (all inputs to A only):
#   Yasso15/20 and Q-model also follow this convention for the fast-to-slow
#   transfer. Routing inputs only to A avoids a free input-partitioning
#   parameter that would be weakly constrained and would confound the
#   structural comparison.
#
# CLIMATE MODIFIER:
#   Identical formula to SP1 and Yasso07 (same functional form, independently
#   calibrated beta1/beta2/gamma). xi multiplies BOTH pool decay rates --
#   this is a key structural difference from ICBM's separate re/h factors.
#   See sp1_wrapper.R for xi derivation.
#
# TRANSIENT SIMULATION:
#   A is solved exactly each step (independent of H). H uses the analytical
#   integration of A's trajectory as its time-varying input term -- no Euler
#   approximation. The combined formula per step (for k_A = alpha_A*xi,
#   k_H = alpha_H*xi, dt = 1 yr):
#
#     A_ss     = J / k_A          (A's instantaneous steady state)
#     delta_A  = A(t) - A_ss
#     A(t+1)   = A_ss + delta_A * exp(-k_A)           [exact]
#
#     H(t+1)   = H(t)*exp(-k_H)                        [decay]
#              + p_H * J * (1 - exp(-k_H)) / k_H       [from A_ss component]
#              + p_H * k_A * delta_A * (exp(-k_A) - exp(-k_H)) / (k_H - k_A)
#                                                        [from transient A term]
#
#   Derivation: solve dH/dt = p_H*k_A*A(t') - k_H*H by integrating the
#   variation-of-constants formula with A(t') expanded analytically.
#   The degenerate case k_A == k_H (unlikely in practice) uses the limiting
#   form: last H term becomes p_H * k_A * delta_A * exp(-k_A).
#
# STEADY STATE:
#   A* = J / (k_A)                           ... trivially
#   H* = p_H * alpha_A * xi * A* / (alpha_H * xi) = p_H * J / (alpha_H * xi)
#   Note xi cancels in H*: H* = p_H * A_ss * alpha_A / alpha_H.
#   Total SOC: C* = J/xi * (1/alpha_A + p_H/alpha_H)
#
# PARAMETERS (all free in calibration):
#   alpha_A    : A pool decay rate (yr^-1). Positive (log transform).
#   alpha_H    : H pool decay rate (yr^-1). Positive (log transform).
#   p_H        : humification fraction, A -> H. In (0,1) (logit transform).
#   beta1      : temperature response, linear. Positive (log transform).
#   beta2      : temperature response, quadratic. Unconstrained.
#   gamma      : precipitation response. Unconstrained.
#   sigma_init : initial-condition scaling. Log-transformed.
#   sigma_input: global litter input multiplier. Log-transformed.
#
# DEPENDENCIES: base R only (no Fortran, no external packages).
# =============================================================================


TP2_PARAM_NAMES <- c("alpha_A", "alpha_H", "p_H",
                     "beta1", "beta2", "gamma",
                     "sigma_init", "sigma_input")


# =============================================================================
# Steady-state initialisation
# =============================================================================

# Analytical steady-state pool distribution for TP2.
#
# At steady state (dA/dt = 0, dH/dt = 0) with inputs sigma_input * J_mean
# and climate scalar xi_mean:
#
#   A* = J / (alpha_A * xi_mean)
#   H* = p_H * J / (alpha_H * xi_mean)
#
# Note that xi appears in both denominators, so H* has the same xi dependence
# as A*. This is because p_H * alpha_A * xi * A* = p_H * J at steady state
# (all A production balanced by input), and then H* = p_H*J / (alpha_H*xi).
#
# ARGUMENTS:
#   model_params : named numeric -- alpha_A, alpha_H, p_H, sigma_input
#   lm           : named list -- $J_total_mean (scalar tC/ha/yr)
#   xi_mean      : scalar
#
# RETURNS:
#   named numeric vector c(A = A_ss, H = H_ss)
tp2_steady_state <- function(model_params, lm, xi_mean) {
  alpha_A <- unname(model_params["alpha_A"])
  alpha_H <- unname(model_params["alpha_H"])
  p_H     <- unname(model_params["p_H"])
  si      <- unname(model_params["sigma_input"])
  J       <- lm$J_total_mean * si

  A_ss <- J / (alpha_A * xi_mean)
  H_ss <- p_H * J / (alpha_H * xi_mean)

  c(A = unname(A_ss), H = unname(H_ss))
}


# =============================================================================
# Per-step exact solver (internal helper)
# =============================================================================

# Advance A and H by one year (dt = 1) under constant J and xi.
# See file header for full derivation.
#
# ARGUMENTS:
#   A, H    : current pool sizes (scalars)
#   J       : litter input this year (already sigma_input-scaled)
#   k_A, k_H: effective decay rates (alpha_A/H * xi)
#   p_H     : humification fraction
#
# RETURNS:
#   named numeric vector c(A_new, H_new)
tp2_step <- function(A, H, J, k_A, k_H, p_H) {
  e_A     <- exp(-k_A)
  e_H     <- exp(-k_H)
  A_ss    <- J / k_A
  delta_A <- A - A_ss

  A_new <- A_ss + delta_A * e_A

  # H integrates two input components:
  #   (i)  p_H * J * (1 - e_H) / k_H     : from A's steady-state component
  #   (ii) p_H * k_A * delta_A * (e_A - e_H) / (k_H - k_A) : from transient A
  # Degenerate case k_A == k_H (limiting form via L'Hopital):
  #   term (ii) = p_H * k_A * delta_A * e_A   [limit as k_A -> k_H]
  if (abs(k_H - k_A) < 1e-8) {
    H_new <- H * e_H +
             p_H * J       * (1 - e_H) / k_H +
             p_H * k_A * delta_A * e_A
  } else {
    H_new <- H * e_H +
             p_H * J       * (1 - e_H) / k_H +
             p_H * k_A * delta_A * (e_A - e_H) / (k_H - k_A)
  }

  c(A = A_new, H = H_new)
}


# =============================================================================
# Transient simulation
# =============================================================================

# Run the TP2 model from C_init over the full input time series.
#
# ARGUMENTS:
#   inputs       : data.frame -- annual inputs; must have columns year, J_total
#   model_params : named numeric -- TP2_PARAM_NAMES subset needed here
#   C_init       : named numeric c(A=..., H=...) from tp2_steady_state()
#   xi_array     : numeric vector -- one xi per row of inputs
#
# RETURNS:
#   data.frame with columns: year, A, H, total_soc
tp2_run <- function(inputs, model_params, C_init, xi_array) {
  alpha_A <- unname(model_params["alpha_A"])
  alpha_H <- unname(model_params["alpha_H"])
  p_H     <- unname(model_params["p_H"])
  si      <- unname(model_params["sigma_input"])
  n       <- nrow(inputs)

  A    <- numeric(n + 1L)
  H    <- numeric(n + 1L)
  A[1] <- C_init["A"]
  H[1] <- C_init["H"]

  for (t in seq_len(n)) {
    xi  <- xi_array[t]
    J   <- inputs$J_total[t] * si
    k_A <- alpha_A * xi
    k_H <- alpha_H * xi

    new_pools <- tp2_step(A[t], H[t], J, k_A, k_H, p_H)
    A[t + 1L] <- new_pools["A"]
    H[t + 1L] <- new_pools["H"]
  }

  A_out <- A[-1L]
  H_out <- H[-1L]

  data.frame(
    year      = inputs$year,
    A         = A_out,
    H         = H_out,
    total_soc = A_out + H_out
  )
}

# =============================================================================
# tp2_transient_init
#
# 68-year pre-run (1917 -> 1985) identical in structure to sp1_transient_init.
# Uses tp2_step for the exact per-year A/H analytical update.
# =============================================================================

tp2_transient_init <- function(model_params, lm, xi_mean, ...) {
  xi_mean     <- unname(xi_mean)
  alpha_A     <- unname(model_params["alpha_A"])
  alpha_H     <- unname(model_params["alpha_H"])
  p_H         <- unname(model_params["p_H"])
  sigma_init  <- unname(model_params["sigma_init"])
  sigma_input <- unname(model_params["sigma_input"])

  # Fully-scaled litter endpoints
  J_1917 <- lm$J_full_mean * sigma_init * sigma_input
  J_1985 <- lm$J_t0_mean   * sigma_input

  # Effective decay rates (constant xi throughout pre-run)
  k_A <- alpha_A * xi_mean
  k_H <- alpha_H * xi_mean

  # Start at analytical steady state under 1917 litter
  A <- unname(J_1917 / k_A)
  H <- unname(p_H * J_1917 / k_H)

  # 68-year pre-run: linearly interpolated J, constant xi
  n_pre <- 68L
  for (i in seq_len(n_pre)) {
    frac   <- (i - 1L) / (n_pre - 1L)
    J      <- J_1917 + (J_1985 - J_1917) * frac
    new_AH <- tp2_step(A, H, J, k_A, k_H, p_H)
    A      <- unname(new_AH["A"])
    H      <- unname(new_AH["H"])
  }

  c(A = A, H = H)
}
