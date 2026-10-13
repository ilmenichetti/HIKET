# =============================================================================
# scenario_switches.R -- the scenario arms of the recalibrated decomposition (2026-10-13)
# -----------------------------------------------------------------------------
# Each scenario of the sink decomposition becomes its OWN calibration, so the parts are
# differences between calibrated models -- a comparison between assumptions only
# (Lorenzo, 2026-10-13). The scenario definitions are those of
# doublechecks/attribution_decomposition.R, unchanged; only the calibration changes.
#
#   switch                       scenario property                       definition
#   HIKET_EQUILIBRIUM_INIT=1     balanced start in 1985                  equilibrium_init.R
#   HIKET_CLIMATE_HELD=1         climate held at the reference           each plot's forward climate modifier
#                                                                        = its pre-run modifier, every year
#   HIKET_LITTER_HELD=1          litter held at its 1985-89 level        the litter the pre-run ends on
#
# All inert unless set (on Roihu: the SINGULARITYENV_ prefix). They compose:
#   production                   scenario 1   (estimated start, observed litter and climate)
#   EQUILIBRIUM_INIT             scenario 2   (balanced start)
#   EQUILIBRIUM_INIT+CLIMATE     scenario 3   (balanced start, reference climate)  -> three parts
#   EQUILIBRIUM_INIT+LITTER      scenario 4   (balanced start, litter held)        -> five terms
#   CLIMATE                      scenario 1'  (estimated start, reference climate) -> five terms
#
# Climate enters every model only through the climate modifier (Yasso15/20 also take
# precipitation for a leaching term, but leac = 0), so holding the modifier holds the
# climate completely. The reference is the run's pre-run climate (HIKET_PRERUN_CLIM_YEARS,
# 1961-1984 in the current runs), exactly as in the decomposition.
#
# WHERE THE HOOKS ARE: the likelihood (calibration_engine_transient.R, one per-plot path
# shared by the independent and the correlated likelihood). The predictive scripts REFUSE
# to run under a held switch: the recalibrated decomposition runs each arm's posterior
# directly (attribution_decomposition.R, HIKET_ATTR_MODE=recalibrated).
#
# NAMING: hiket_eq_tag() is redefined here to compose every switch, e.g.
# "Yasso15" -> "Yasso15_eqinit_climheld"; the figure layer's run selection keys on the
# exact tag, so arms never mix.
# Tests: doublechecks/scenario_switch_tests.R (no-op identity; hooks vs an independent
# construction of the same scenario).
# =============================================================================

.hiket_climate_held <- identical(Sys.getenv("HIKET_CLIMATE_HELD", "0"), "1")
.hiket_litter_held  <- identical(Sys.getenv("HIKET_LITTER_HELD",  "0"), "1")
if (.hiket_climate_held)
  message("[CLIMATE HELD] ON -- forward climate modifier held at the pre-run (reference) value; ",
          "outputs tagged _climheld")
if (.hiket_litter_held)
  message("[LITTER HELD] ON -- forward litter held at its 1985-89 level; outputs tagged _litheld")

# Output tag, composing every switch (overrides the one in equilibrium_init.R).
hiket_arm_suffix <- function()
  paste0(if (.hiket_eq_init) "_eqinit" else "",
         if (.hiket_climate_held) "_climheld" else "",
         if (.hiket_litter_held) "_litheld" else "")
hiket_eq_tag <- function(model) paste0(model, hiket_arm_suffix())
hiket_arm_label <- function() { s <- hiket_arm_suffix(); if (nzchar(s)) sub("^_", "", s) else "production" }

# Climate held: the forward modifier becomes the pre-run modifier, every year.
# Same construction as xi_rep() in attribution_decomposition.R.
hiket_xi_held <- function(xi_array, xi_ss) {
  if (!.hiket_climate_held) return(xi_array)
  n <- length(if (is.list(xi_array)) xi_array[[1]] else xi_array)
  if (is.list(xi_ss)) lapply(xi_ss, function(v) rep(unname(v)[1], n)) else rep(unname(xi_ss)[1], n)
}

# Litter held: every forward year at the 1985-89 level the pre-run ends on
# (litter_rows(f = 1) in attribution_decomposition.R). Unscaled; the run function
# applies sigma_input, as always.
hiket_litter_held <- function(inputs, lm) {
  if (!.hiket_litter_held) return(inputs)
  r <- inputs
  if ("J_total" %in% names(r)) r$J_total <- lm$J_t0_mean
  else for (s in c("nwl", "fwl", "cwl")) {
    v <- lm[[paste0(s, "_t0_mean")]]
    for (k in c("A", "W", "E", "N")) r[[paste0(s, "_", k)]] <- unname(v[[paste0(s, "_", k)]])
  }
  r
}

# The predictive stage has no held hooks: refuse rather than mislabel an arm.
hiket_arm_predictive_guard <- function() {
  if (.hiket_climate_held || .hiket_litter_held)
    stop("[SCENARIO] the predictive stage does not implement HIKET_CLIMATE_HELD / ",
         "HIKET_LITTER_HELD. The held arms are used only by the recalibrated decomposition: ",
         "Rscript doublechecks/attribution_decomposition.R with HIKET_ATTR_MODE=recalibrated.",
         call. = FALSE)
}
