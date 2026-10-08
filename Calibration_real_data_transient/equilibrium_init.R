# =============================================================================
# equilibrium_init.R -- the equilibrium-init counterfactual arm (2026-10-08)
# -----------------------------------------------------------------------------
# HIKET_EQUILIBRIUM_INIT=1 turns every run_*_transient_{calibration,predictive}.R
# into the equilibrium arm (on Roihu: SINGULARITYENV_HIKET_EQUILIBRIUM_INIT=1).
# Unset => production, unchanged. Design: NEXT_RUN_equilibrium_init.md.
#
# WHAT THE ARM IS: the production calibration with sigma_init FIXED AT 1 and
# removed from the free set. The pre-run flux is
#   J(i) = J_1917 + (J_1985 - J_1917) * shape[i],  J_1917 = J_t0_mean*sigma_init*sigma_input
# so sigma_init = 1 makes the ramp flat WHATEVER the shape, and the initialiser
# returns the steady state at the 1985-89 litter and the 1985-2004 climate.
# No initialiser is re-implemented; the one factor that changes is whether the
# soil carries a deficit into 1985.
#
# WHAT STAYS: make_likelihood(..., transient_init = TRUE). The FALSE branch is
# the old engine, which adds sigma_init to the error of the first observation:
# a second factor.
#
# PRIOR: flux_pair -> flux_now, whose Jacobian keeps the -log(F_now) term so the
# sigma_input prior is the same as in the transient arm (calibration_engine.R,
# flux_now header; doublechecks/eqinit_prior_check.R).
#
# NAMING: outputs carry MODEL_NAME = "<model>_eqinit", so the figure layer's
# "^<model>_posterior_<RUN_ID>.rds$" auto-selection never picks them up.
# =============================================================================

.hiket_eq_init <- identical(Sys.getenv("HIKET_EQUILIBRIUM_INIT", "0"), "1")
if (.hiket_eq_init)
  message("[EQUILIBRIUM INIT] ON -- sigma_init fixed at 1, removed from the free set; ",
          "outputs tagged _eqinit")

# Output tag: "TP2" -> "TP2_eqinit" in the equilibrium arm.
hiket_eq_tag <- function(model) if (.hiket_eq_init) paste0(model, "_eqinit") else model

# flux_pair (sigma_input, sigma_init) -> flux_now (sigma_input), same window/J_bar.
hiket_eq_param_spec <- function(param_spec) {
  if (!.hiket_eq_init) return(param_spec)
  out <- lapply(param_spec, function(g) {
    if (!identical(g$type, "flux_pair")) return(g)
    stopifnot(identical(g$names, c("sigma_input", "sigma_init")))
    list(names = "sigma_input", type = "flux_now", window = g$window, J_bar = g$J_bar)
  })
  n_old <- length(unlist(lapply(param_spec, `[[`, "names")))
  n_new <- length(unlist(lapply(out,        `[[`, "names")))
  stopifnot(n_new == n_old - 1L)
  message(sprintf("[EQUILIBRIUM INIT] free parameters: %d (production %d)", n_new, n_old))
  out
}

# Wrap a model's assemble_model_params so sigma_init = 1 always (also overrides a
# sigma_init carried in free_defaults, as used by the pre-MCMC sanity run).
hiket_eq_assemble <- function(f) {
  if (!.hiket_eq_init) return(f)
  force(f)
  function(p_free) { p <- f(p_free); p["sigma_init"] <- 1; p }
}

# Drop sigma_init from diagnostic plot lists.
hiket_eq_names <- function(nms) if (.hiket_eq_init && !is.null(nms)) setdiff(nms, "sigma_init") else nms
