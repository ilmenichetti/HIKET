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
