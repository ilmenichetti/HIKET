# =============================================================================
# forward_scenarios.R -- helpers for the two side analyses of the projection
# -----------------------------------------------------------------------------
# Sourced by all six run_*_transient_predictive.R (Section 4b). Adds, on the SAME
# posterior draws, plots and projection forcing as the production projection:
#
#   (1) INPUT SCENARIO: the projection repeated with the litter input multiplied
#       by INPUT_SCEN_MULT from the first projected year (default 1.2 = +20%).
#   (2) EQUILIBRIUM STOCK: the stock each plot would reach if the litter input of
#       the last observed year and the recycled recent climate were held forever,
#       from the model's OWN initialiser: sigma_init <- 1 and preinit_shape <- 0
#       hold the pre-run input flat at the chosen level, so the initialiser
#       returns its own steady state (no re-implementation of any model).
#
# Uses no random numbers: every existing output of the predictive stage is
# unchanged. Added 2026-09-23 (manuscript round 5).
# =============================================================================

INPUT_SCEN_MULT <- suppressWarnings(as.numeric(Sys.getenv("HIKET_INPUT_SCEN_MULT", "1.2")))
if (!is.finite(INPUT_SCEN_MULT) || INPUT_SCEN_MULT <= 0)
  stop("HIKET_INPUT_SCEN_MULT must be a positive number")

# Multiply the litter columns of a projection input table (not year, plot, precip).
scale_litter <- function(inputs_proj, mult) {
  lit <- setdiff(names(inputs_proj), c("plot_id", "year", "precip"))
  inputs_proj[lit] <- inputs_proj[lit] * mult
  inputs_proj
}

# Litter summary that makes the initialiser return the equilibrium at the input of
# the last observed year: flat pre-run at that input (the caller sets sigma_init = 1).
equilibrium_lm <- function(lm, inputs, recycle_years) {
  l <- lm
  l$preinit_shape <- rep(0, length(lm$preinit_shape))
  last <- inputs[nrow(inputs), , drop = FALSE]
  if (!is.null(l$J_t0_mean)) l$J_t0_mean <- last$J_total
  for (s in c("nwl", "fwl", "cwl")) {
    nm <- paste0(s, "_t0_mean")
    if (!is.null(l[[nm]])) {
      cols <- paste0(s, "_", c("A", "W", "E", "N"))
      l[[nm]] <- setNames(as.numeric(unlist(last[cols])), cols)
    }
  }
  if (!is.null(l$precip_mean) && "precip" %in% names(inputs))
    l$precip_mean <- mean(tail(inputs$precip, recycle_years))
  l
}
