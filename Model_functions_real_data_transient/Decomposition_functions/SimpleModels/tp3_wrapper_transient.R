# =============================================================================
# tp3_wrapper_transient.R
#
# Three-pool sequential cascade (ASH) model for the HIKET pipeline.
#
# Pool topology:
#   Active (A) → Slow (S) → Humus (H)
#
# Annual discrete dynamics:
#   A(t+1) = A(t) * (1 - alpha_A * xi(t)) + J(t)
#   S(t+1) = S(t) * (1 - alpha_S * xi(t)) + p_S * alpha_A * xi(t) * A(t)
#   H(t+1) = H(t) * (1 - alpha_H)         + p_H * alpha_S * xi(t) * S(t)
#
# Notes:
#   - H has no climate modifier (xi), consistent with TP2 and Yasso conventions.
#   - p_S: fraction of A decomposition flux routed to S; remainder is respired.
#   - p_H: fraction of S decomposition flux routed to H; remainder is respired.
#   - xi is computed externally via compute_xi_yasso07 (same function as TP2/SP1).
#   - sigma_input: global litter multiplier (calibrated); applied to J(t) at each step.
#   - sigma_init:  scales pre-run (1917) litter relative to the contemporary mean;
#                  same interpretation as TP2 (see below).
#
# Analytical steady state:
#   A_ss = J / (alpha_A * xi)
#   S_ss = p_S * J / (alpha_S * xi)
#   H_ss = p_H * p_S * J / alpha_H
#
# Transient initialisation (tp3_transient_init):
#   68-year pre-run (1917 → 1985) with linearly interpolated litter:
#     J_1917 = J_full_mean * sigma_init * sigma_input  (historical anchor)
#     J_1985 = J_t0_mean   * sigma_input               (observed-period start)
#   Starts at analytical steady state under J_1917; returns terminal pool state.
#
#   sigma_init = 1.0  →  1917 litter equal to the contemporary series mean.
#   sigma_init > 1.0  →  historically more productive (old-growth / unmanaged).
#   sigma_init < 1.0  →  historically less productive (post-disturbance recovery).
#
# Exported functions:
#   tp3_steady_state   (model_params, lm, xi_mean)       → named numeric(3) [A,S,H]
#   tp3_transient_init (model_params, lm, xi_mean)       → named numeric(3) [A,S,H]
#   tp3_run            (inputs, model_params, C_init, xi_array) → data.frame
# =============================================================================


# -----------------------------------------------------------------------------
# tp3_steady_state
#
# Analytical steady-state pool sizes at a given mean litter flux J and mean
# climate modifier xi_mean.  Used as the starting point for the pre-run inside
# tp3_transient_init and as a fallback when transient initialisation is not used.
# -----------------------------------------------------------------------------
tp3_steady_state <- function(model_params, lm, xi_mean) {
  J    <- lm$J_total_mean * model_params["sigma_input"]
  A_ss <- J / (model_params["alpha_A"] * xi_mean)
  S_ss <- model_params["p_S"] * J / (model_params["alpha_S"] * xi_mean)
  H_ss <- model_params["p_H"] * model_params["p_S"] * J / model_params["alpha_H"]
  c(A = unname(A_ss), S = unname(S_ss), H = unname(H_ss))
}


# -----------------------------------------------------------------------------
# tp3_transient_init
#
# 68-year pre-run (1917 → 1985): litter linearly interpolated from
# J_1917 (J_full_mean * sigma_init * sigma_input) to J_1985 (J_t0_mean *
# sigma_input). Climate constant at xi_mean (no pre-1985 observations).
# Starts at analytical steady state under J_1917; returns terminal state.
# -----------------------------------------------------------------------------
tp3_transient_init <- function(model_params, lm, xi_mean) {
  xi_mean     <- unname(xi_mean)
  alpha_A     <- unname(model_params["alpha_A"])
  alpha_S     <- unname(model_params["alpha_S"])
  alpha_H     <- unname(model_params["alpha_H"])
  p_S         <- unname(model_params["p_S"])
  p_H         <- unname(model_params["p_H"])
  sigma_init  <- unname(model_params["sigma_init"])
  sigma_input <- unname(model_params["sigma_input"])

  # Fully-scaled litter endpoints
  J_1917 <- lm$J_full_mean * sigma_init * sigma_input
  J_1985 <- lm$J_t0_mean   * sigma_input

  # Start at analytical steady state under 1917 litter
  C <- c(A = unname(J_1917 / (alpha_A * xi_mean)),
         S = unname(p_S * J_1917 / (alpha_S * xi_mean)),
         H = unname(p_H * p_S * J_1917 / alpha_H))

  # 68-year pre-run: linearly interpolated J, constant xi
  n_pre <- 68L
  for (i in seq_len(n_pre)) {
    frac   <- (i - 1L) / (n_pre - 1L)
    J      <- J_1917 + (J_1985 - J_1917) * frac
    flux_A <- alpha_A * xi_mean * C["A"]
    flux_S <- alpha_S * xi_mean * C["S"]
    C["A"] <- C["A"] - flux_A + J
    C["S"] <- C["S"] + p_S * flux_A - flux_S
    C["H"] <- C["H"] + p_H * flux_S - alpha_H * C["H"]
  }
  C
}


# -----------------------------------------------------------------------------
# tp3_run
#
# Runs the ASH model forward over the observed period given initial pool state
# C_init, annual litter inputs (J_total per year), and the pre-computed annual
# xi array.
#
# Arguments:
#   inputs      data.frame with columns year, J_total (annual, tC/ha/yr)
#   model_params named vector including alpha_A, alpha_S, alpha_H, p_S, p_H,
#               sigma_input (and any other params; extra names are ignored)
#   C_init      named numeric(3): A, S, H pool sizes at t=0 (tC/ha)
#   xi_array    numeric vector length nrow(inputs): annual climate modifier
#
# Returns:
#   data.frame with columns year, A, S, H, total_soc (all in tC/ha)
# -----------------------------------------------------------------------------
tp3_run <- function(inputs, model_params, C_init, xi_array) {
  n  <- nrow(inputs)
  A  <- numeric(n)
  S  <- numeric(n)
  H  <- numeric(n)

  alpha_A <- unname(model_params["alpha_A"])
  alpha_S <- unname(model_params["alpha_S"])
  alpha_H <- unname(model_params["alpha_H"])
  p_S     <- unname(model_params["p_S"])
  p_H     <- unname(model_params["p_H"])
  sig_inp <- unname(model_params["sigma_input"])

  C_A <- unname(C_init["A"])
  C_S <- unname(C_init["S"])
  C_H <- unname(C_init["H"])

  for (t in seq_len(n)) {
    J      <- inputs$J_total[t] * sig_inp
    xi_t   <- xi_array[t]
    flux_A <- alpha_A * xi_t * C_A
    flux_S <- alpha_S * xi_t * C_S
    C_A    <- C_A - flux_A + J
    C_S    <- C_S + p_S * flux_A - flux_S
    C_H    <- C_H + p_H * flux_S - alpha_H * C_H
    A[t]   <- C_A
    S[t]   <- C_S
    H[t]   <- C_H
  }

  data.frame(
    year      = inputs$year,
    A         = A,
    S         = S,
    H         = H,
    total_soc = A + S + H
  )
}
