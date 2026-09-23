# =============================================================================
# build_forward_scenarios.R -- numbers for the two side analyses of the projection
# -----------------------------------------------------------------------------
# Reads the columns added to projection_predictions by the predictive stage
# (Calibration_real_data_transient/forward_scenarios.R, 2026-09-23):
#   total_soc      baseline projection: litter held at the last observed year,
#                  last 20 years of climate recycled, 60 years
#   total_soc_up   the same with the litter input x HIKET_INPUT_SCEN_MULT (1.2)
#   C_last         stock in the last observed year (2024)
#   C_eq           equilibrium stock at the last observed input and recycled climate
# Plot set: the balanced set (observed in all three campaigns), as every observed
# rate in the paper. Summaries are median and 5-95% range over the posterior draws
# of the national (unweighted cross-plot) mean.
#
# Run from repo root:  Rscript manuscript/figures/build_forward_scenarios.R
# =============================================================================
suppressMessages({ source("manuscript/figures/run_ids.R"); source("manuscript/figures/obs_basis.R") })
cat("RUN_IDs:", paste(names(RID), RID, collapse = " | "), "\n")
MULT <- suppressWarnings(as.numeric(Sys.getenv("HIKET_INPUT_SCEN_MULT", "1.2")))
q <- function(x) sprintf("%.3g [%.3g, %.3g]", median(x), quantile(x, .05), quantile(x, .95))

out <- list()
for (m in FIG_MODELS) {
  pp <- readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_predictive_%s.rds", m, RID[[m]]))$projection_predictions
  if (!all(c("total_soc_up", "C_last", "C_eq") %in% names(pp)))
    stop(m, ": projection_predictions lacks the side-analysis columns -- rerun the predictive stage")
  om <- readRDS(sprintf("Data/model_inputs/%s_inputs_%s.rds", m, RID[[m]]))$obs_meta
  pp <- pp[as.integer(pp$plot_id) %in% balanced_plots(om), ]
  y0 <- min(pp$year) - 1L; y1 <- max(pp$year)
  # national mean per draw and year
  S <- tapply(pp$total_soc,    list(pp$draw, pp$year), mean)
  U <- tapply(pp$total_soc_up, list(pp$draw, pp$year), mean)
  first <- !duplicated(pp[c("draw", "plot_id")])
  C0 <- tapply(pp$C_last[first], pp$draw[first], mean)
  CE <- tapply(pp$C_eq[first],   pp$draw[first], mean)
  yr <- as.integer(colnames(S)); at <- function(M, y) M[, match(y, yr)]
  d <- data.frame(
    model        = m, draw = as.integer(rownames(S)), n_plots = length(unique(pp$plot_id)),
    C_last       = C0, C_eq = CE,
    headroom     = CE / C0 - 1,
    # rates over WHOLE 20-year climate cycles: decades alternate between the replayed
    # 2005-14 and 2015-24 climates (doublechecks/yasso20_late_source_check.R)
    sink_first20 = (at(S, y0 + 20L) - C0) / 20,
    sink_last20  = (at(S, y1) - at(S, y1 - 20L)) / 20,
    realised_84  = (at(S, y1) - C0) / (CE - C0),
    up_gain_eq   = (MULT - 1) * CE,                       # eventual extra stock (linear models)
    up_extra_50  = at(U, 2050L) - at(S, 2050L),
    up_extra_84  = at(U, y1) - at(S, y1))
  d$up_frac_50 <- d$up_extra_50 / d$up_gain_eq
  d$up_frac_84 <- d$up_extra_84 / d$up_gain_eq
  out[[m]] <- d
  cat(sprintf("\n%s  (%d plots, %d draws, %d-%d)\n", m, d$n_plots[1], nrow(d), y0, y1))
  cat("  stock", y0, ":", q(d$C_last), "  equilibrium:", q(d$C_eq), "  headroom:", q(d$headroom), "\n")
  cat("  baseline sink, first 20-yr cycle:", q(d$sink_first20), "  last 20-yr cycle:", q(d$sink_last20),
      "  share of headroom realised by", y1, ":", q(d$realised_84), "\n")
  cat(sprintf("  input x%.2f: extra stock 2050 %s, %d %s; eventual %s; share realised 2050 %s, %d %s\n",
              MULT, q(d$up_extra_50), y1, q(d$up_extra_84), q(d$up_gain_eq), q(d$up_frac_50), y1, q(d$up_frac_84)))
}
res <- do.call(rbind, out); attr(res, "run_ids") <- RID
saveRDS(res, "manuscript/figures/forward_scenarios.rds")
cat("\nsaved manuscript/figures/forward_scenarios.rds\n")
