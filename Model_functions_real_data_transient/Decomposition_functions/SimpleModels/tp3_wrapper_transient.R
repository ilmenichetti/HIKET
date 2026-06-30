# =============================================================================
# tp3_wrapper_transient.R
#
# Three-pool sequential cascade (ASH) model for the HIKET pipeline.
#
# Pool topology:
#   Active (A) -> Slow (S) -> Humus (H)
#
# Continuous dynamics (constant J and xi within each annual step):
#   dA/dt = J            - k_A A
#   dS/dt = p_S k_A A    - k_S S
#   dH/dt = p_H k_S S    - k_H H
#   with k_A = alpha_A * xi,  k_S = alpha_S * xi,  k_H = alpha_H  (H has no xi).
#
# INTEGRATION (changed 2026-06: was explicit forward Euler).
#   Each annual step is now solved EXACTLY, by the matrix exponential of the
#   linear generator, evaluated in closed form. For the lower-triangular
#   cascade generator M the exponential is itself lower-triangular and its
#   entries are divided differences of exp() over the eigenvalues
#   (-k_A, -k_S, -k_H) -- so the step costs a handful of scalars, not an
#   iterative expm. This makes TP3 analytically exact like SP1, TP2 (closed
#   form) and the Yasso family (matrix exponential); the explicit Euler step
#   previously used here was the sole approximate integrator in the ensemble
#   and produced a spurious interannual oscillation when the calibrated active
#   pool turned over faster than the 1-year step (alpha_A * xi > 1). Exact
#   integration removes that artefact and leaves the inter-model differences
#   purely structural. See HIKET_calibration.Rmd section 13.8.
#
# Notes:
#   - H has no climate modifier (xi), consistent with TP2 and Yasso conventions.
#   - p_S: fraction of A decomposition flux routed to S; remainder is respired.
#   - p_H: fraction of S decomposition flux routed to H; remainder is respired.
#   - xi is computed externally via compute_xi_yasso07 (same function as TP2/SP1).
#   - sigma_input: global litter multiplier (calibrated); applied to J(t) at each step.
#   - sigma_init:  scales pre-run (1917) litter relative to the contemporary mean.
#
# Analytical steady state (within-year equilibrium of the dynamics above):
#   A_ss = J / (alpha_A * xi)
#   S_ss = p_S * J / (alpha_S * xi)
#   H_ss = p_H * p_S * J / alpha_H
#
# Exported functions:
#   tp3_steady_state   (model_params, lm, xi_mean)            -> named numeric(3) [A,S,H]
#   tp3_transient_init (model_params, lm, xi_mean)            -> named numeric(3) [A,S,H]
#   tp3_run            (inputs, model_params, C_init, xi_array) -> data.frame
# =============================================================================


# -----------------------------------------------------------------------------
# .tp3_step  (internal)
#
# Exact one-year update of the ASH cascade with J and xi held constant over the
# step. Solves C(t+1) = C_ss + exp(M) (C(t) - C_ss), where C_ss is the
# within-year equilibrium and exp(M) is the (closed-form) exponential of the
# lower-triangular generator. The divided-difference formula is undefined when
# two eigenvalues coincide; since alpha_A >> alpha_S >> alpha_H here that never
# occurs in practice, but a 1e-6 nudge guards the degenerate case.
# -----------------------------------------------------------------------------
.tp3_step <- function(cA, cS, cH, kA, kS, kH, pS, pH, J) {
  # strip any names off the scalar inputs: if cA/cS/cH arrive named (e.g. via
  # C["A"]), the named arithmetic propagates into c(A=...) as compound names
  # ("A.A", ...) and the caller's next C["A"] then returns NA. Unname guards
  # both call sites (tp3_transient_init passes C["A"]; tp3_run passes scalars).
  cA <- unname(cA); cS <- unname(cS); cH <- unname(cH)

  # keep eigenvalues distinct for the divided-difference exponential
  if (abs(kA - kS) < 1e-6) kS <- kS + 1e-6
  if (abs(kA - kH) < 1e-6) kH <- kH + 1e-6
  if (abs(kS - kH) < 1e-6) kH <- kH + 2e-6

  # within-year equilibrium (litter J enters A only)
  Ass <- J / kA
  Sss <- pS * J / kS
  Hss <- pH * pS * J / kH

  l1 <- -kA; l2 <- -kS; l3 <- -kH
  e1 <- exp(l1); e2 <- exp(l2); e3 <- exp(l3)

  # exp(M) entries: diagonal e^li; off-diagonals = (path product) x divided diff
  d21  <- (e2 - e1) / (l2 - l1)                       # 1st divided diff, nodes l1,l2
  d32  <- (e3 - e2) / (l3 - l2)                       # 1st divided diff, nodes l2,l3
  dd31 <- e1 / ((l1 - l2) * (l1 - l3)) +              # 2nd divided diff, nodes l1,l2,l3
          e2 / ((l2 - l1) * (l2 - l3)) +
          e3 / ((l3 - l1) * (l3 - l2))
  a <- pS * kA            # M[2,1]
  cc <- pH * kS           # M[3,2]

  dA <- cA - Ass; dS <- cS - Sss; dH <- cH - Hss
  c(A = Ass + e1 * dA,
    S = Sss + a * d21 * dA + e2 * dS,
    H = Hss + a * cc * dd31 * dA + cc * d32 * dS + e3 * dH)
}


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
# 68-year pre-run (1917 -> 1985): litter linearly interpolated from
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

  # 68-year pre-run: linearly interpolated J, constant xi, exact annual step
  kA <- alpha_A * xi_mean
  kS <- alpha_S * xi_mean
  kH <- alpha_H
  n_pre <- 68L
  for (i in seq_len(n_pre)) {
    frac <- (i - 1L) / (n_pre - 1L)
    J    <- J_1917 + (J_1985 - J_1917) * frac
    C    <- .tp3_step(C["A"], C["S"], C["H"], kA, kS, kH, p_S, p_H, J)
  }
  C
}


# -----------------------------------------------------------------------------
# tp3_run
#
# Runs the ASH model forward over the observed period given initial pool state
# C_init, annual litter inputs (J_total per year), and the pre-computed annual
# xi array. Each annual step is integrated exactly (see .tp3_step).
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
    J    <- inputs$J_total[t] * sig_inp
    xi_t <- xi_array[t]
    C    <- .tp3_step(C_A, C_S, C_H,
                      alpha_A * xi_t, alpha_S * xi_t, alpha_H,
                      p_S, p_H, J)
    C_A <- C[["A"]]; C_S <- C[["S"]]; C_H <- C[["H"]]
    A[t] <- C_A; S[t] <- C_S; H[t] <- C_H
  }

  data.frame(
    year      = inputs$year,
    A         = A,
    S         = S,
    H         = H,
    total_soc = A + S + H
  )
}
