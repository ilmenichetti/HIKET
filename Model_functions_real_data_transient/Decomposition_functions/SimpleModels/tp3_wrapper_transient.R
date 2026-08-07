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
#   with k_A = alpha_A * xi,  k_S = alpha_S * xi,  k_H = alpha_H * xi.
#   (C2, 2026-07-16: xi now multiplies ALL pool rates including humus, so TP3's
#    climate treatment matches TP2 and the Yasso family; previously H had no xi.)
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
#   - H IS climate-modified (k_H = alpha_H * xi), matching TP2 (which applies xi
#     to both its pool rates) and the Yasso family. (Corrected C2 2026-07-16 — the
#     earlier "H has no xi, consistent with TP2" note was wrong: TP2 has xi on H.)
#   - p_S: fraction of A decomposition flux routed to S; remainder is respired.
#   - p_H: fraction of S decomposition flux routed to H; remainder is respired.
#   - xi is computed externally via compute_xi_yasso07 (same function as TP2/SP1).
#   - sigma_input: global litter multiplier (calibrated); applied to J(t) at each step.
#   - sigma_init:  scales pre-run (1917) litter relative to the contemporary mean.
#
# Analytical steady state (within-year equilibrium of the dynamics above):
#   A_ss = J / (alpha_A * xi)
#   S_ss = p_S * J / (alpha_S * xi)
#   H_ss = p_H * p_S * J / (alpha_H * xi)
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
# lower-triangular generator. The divided-difference exponential is computed in
# a confluence-safe recursive form (expm1), so it stays exact when two eigenvalues
# coincide. This matters after C1+C2: alpha_S and alpha_H are now BOTH at k2-scale
# (kS ~ kH is a regular regime, not the "alpha_S >> alpha_H, never coincides" case
# the earlier crude 1e-6 nudge assumed — that nudge gave ~1e-4 error near kS==kH).
# The one truly catastrophic coincidence (kA == kH -> the l3-l1 denominator below)
# is physically impossible (fast pool vs slow humus differ by ~100x) and guarded.
# -----------------------------------------------------------------------------
.tp3_step <- function(cA, cS, cH, kA, kS, kH, pS, pH, J) {
  # strip any names off the scalar inputs: if cA/cS/cH arrive named (e.g. via
  # C["A"]), the named arithmetic propagates into c(A=...) as compound names
  # ("A.A", ...) and the caller's next C["A"] then returns NA. Unname guards
  # both call sites (tp3_transient_init passes C["A"]; tp3_run passes scalars).
  cA <- unname(cA); cS <- unname(cS); cH <- unname(cH)

  # Guard ONLY the fast-vs-slow coincidences (kA==kS, kA==kH): never physical
  # (fast pool >> slow pair), but they would divide-by-zero in dd31's outer
  # denominator (l3-l1 = kA-kH). kS~kH needs NO nudge — the recursion below is
  # exact through it.
  if (abs(kA - kS) < 1e-9) kS <- kS + 1e-9
  if (abs(kA - kH) < 1e-9) kH <- kH + 1e-9

  # within-year equilibrium (litter J enters A only)
  Ass <- J / kA
  Sss <- pS * J / kS
  Hss <- pH * pS * J / kH

  l1 <- -kA; l2 <- -kS; l3 <- -kH
  e1 <- exp(l1); e2 <- exp(l2); e3 <- exp(l3)

  # Divided differences of exp, confluence-safe:
  #   f[li,lj] = (e_j - e_i)/(l_j - l_i) = e_i * expm1(l_j-l_i)/(l_j-l_i) -> e_i as l_j->l_i
  #   f[l1,l2,l3] = (f[l2,l3] - f[l1,l2])/(l3 - l1)   (recursive 2nd divided diff)
  ddx  <- function(li, lj, ei) { d <- lj - li; if (abs(d) < 1e-12) ei else ei * expm1(d) / d }
  d21  <- ddx(l1, l2, e1)                 # 1st divided diff, nodes l1,l2
  d32  <- ddx(l2, l3, e2)                 # 1st divided diff, nodes l2,l3 (exact when kS~kH)
  dd31 <- (d32 - d21) / (l3 - l1)         # 2nd divided diff; l3-l1 = kA-kH always large
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
  H_ss <- model_params["p_H"] * model_params["p_S"] * J / (model_params["alpha_H"] * xi_mean)  # C2: xi on H
  c(A = unname(A_ss), S = unname(S_ss), H = unname(H_ss))
}


# -----------------------------------------------------------------------------
# tp3_transient_init
#
# 68-year pre-run (1917 -> 1985): litter linearly interpolated from
# J_1917 (J_t0_mean * sigma_init * sigma_input; P1 2026-08-07, was J_full_mean) to J_1985 (J_t0_mean *
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
# P1 (2026-08-07): the 1917 anchor now uses J_t0_mean, NOT J_full_mean.
# WHY: J_full_mean is the mean over the WHOLE 1986-2024 litter record and J_t0_mean
# the mean over its first five years. Anchoring the two ends of the pre-run to
# DIFFERENT aggregates meant sigma_init was not the 1917/1985 flux ratio but that
# ratio times J_t0/J_full (median 0.818) -- so sigma_init = 1 silently asserted a
# 1917 flux 22% ABOVE 1985, contradicting the growing-stock record the ramp shape
# (C3) is built from. With a common anchor, sigma_init IS the ratio and the
# pre-run inverts at sigma_init > 1 rather than at 0.818.
  J_1917 <- lm$J_t0_mean * sigma_init * sigma_input
  J_1985 <- lm$J_t0_mean   * sigma_input

  # Start at analytical steady state under 1917 litter
  C <- c(A = unname(J_1917 / (alpha_A * xi_mean)),
         S = unname(p_S * J_1917 / (alpha_S * xi_mean)),
         H = unname(p_H * p_S * J_1917 / (alpha_H * xi_mean)))  # C2: xi on H

  # 68-year pre-run: linearly interpolated J, constant xi, exact annual step
  kA <- alpha_A * xi_mean
  kS <- alpha_S * xi_mean
  kH <- alpha_H * xi_mean   # C2: xi on H (was alpha_H)
  # C3: J follows the growing-stock shape (linear fallback if the bundle lacks it)
  n_pre <- 68L
  shape <- if (!is.null(lm$preinit_shape) && length(lm$preinit_shape) == n_pre)
             lm$preinit_shape else (seq_len(n_pre) - 1L) / (n_pre - 1L)
  for (i in seq_len(n_pre)) {
    J    <- J_1917 + (J_1985 - J_1917) * shape[i]
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
                      alpha_A * xi_t, alpha_S * xi_t, alpha_H * xi_t,   # C2: xi on H
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
