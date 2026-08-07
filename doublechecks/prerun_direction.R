# =============================================================================
# prerun_direction.R   (2026-08-07)
#
# WHY THIS EXISTS.
# The transient initialisation starts each plot at steady state under the 1917
# litter flux and ramps it (C3 growing-stock shape) to the 1985 flux:
#
#     J_1917 = J_full_mean * sigma_init * sigma_input
#     J_1985 = J_t0_mean   *              sigma_input
#
# Because the pre-run STARTS at equilibrium for J_1917 and the flux then moves
# monotonically to J_1985, SOC moves monotonically in the same direction. So the
# reconstructed 1917->1985 pre-run DECLINES exactly when
#
#     sigma_init > J_t0_mean / J_full_mean            (sigma_input cancels)
#
# A declining pre-run contradicts the growing-stock history that C3 encodes
# (Finnish forest stocks rose over this period). This script reports, per model,
# how much of the sigma_init posterior sits on the declining side -- using the
# per-plot threshold rather than a single national number.
#
# Usage:  Rscript doublechecks/prerun_direction.R
# =============================================================================

DIR_RUNS <- "Calibration_real_data_transient/runs"
DIR_INP  <- "Data/model_inputs"
MODELS   <- c("SP1", "TP2", "TP3", "Yasso07", "Yasso15", "Yasso20")

latest_run_id <- function(model) {
  fs <- list.files(DIR_RUNS,
                   pattern = sprintf("^%s_posterior_[0-9]{8}_[0-9]{6}\\.rds$", model))
  if (!length(fs)) return(NA_character_)
  sub(sprintf("^%s_posterior_(.+)\\.rds$", model), "\\1", sort(fs, decreasing = TRUE)[1])
}

cat("\n=== Pre-run direction diagnostic ===\n")
cat("declining  <=>  sigma_init > J_t0_mean / J_full_mean\n\n")

thr_all <- NULL
res <- lapply(MODELS, function(m) {
  rid <- latest_run_id(m)
  if (is.na(rid)) return(NULL)
  post <- readRDS(file.path(DIR_RUNS, sprintf("%s_posterior_%s.rds", m, rid)))
  inp  <- readRDS(file.path(DIR_INP,  sprintf("%s_inputs_%s.rds",    m, rid)))

  # Per-plot threshold ratio. Simple models carry scalar J_*; the Yasso models
  # carry AWEN vectors per litter class (nwl/fwl/cwl) -- both sigmas scale every
  # component identically, so the totals give the same ratio.
  thr <- vapply(inp$litter_means, function(x) {
    if (!is.null(x$J_t0_mean) && !is.null(x$J_full_mean))
      return(x$J_t0_mean / x$J_full_mean)
    t0   <- sum(x$nwl_t0_mean,   x$fwl_t0_mean,   x$cwl_t0_mean)
    full <- sum(x$nwl_full_mean, x$fwl_full_mean, x$cwl_full_mean)
    if (!is.finite(t0) || !is.finite(full) || full <= 0) return(NA_real_)
    t0 / full
  }, numeric(1))
  thr <- thr[is.finite(thr)]
  if (is.null(thr_all)) thr_all <<- thr

  si <- unname(post[, "sigma_init"])
  thr <- unname(thr)

  # Fraction of (draw x plot) combinations that decline. Evaluated on a posterior
  # subsample against every plot threshold.
  set.seed(2025)
  si_s <- sample(si, min(4000L, length(si)))
  frac_declining <- mean(outer(si_s, thr, ">"))

  data.frame(
    model          = m,
    run_id         = rid,
    sigma_init_med = median(si),
    sigma_init_q025 = unname(quantile(si, 0.025)),
    sigma_init_q975 = unname(quantile(si, 0.975)),
    thr_median     = median(thr),
    # share of PLOTS declining at the posterior-median sigma_init
    plots_declining_at_med = mean(median(si) > thr),
    # share of draws declining for the median plot
    draws_declining_med_plot = mean(si > median(thr)),
    frac_declining = frac_declining,
    stringsAsFactors = FALSE)
})

out <- do.call(rbind, res)

cat(sprintf("Per-plot threshold J_t0/J_full: median %.3f  [%.3f .. %.3f]  (n = %d plots)\n\n",
            median(thr_all), min(thr_all), max(thr_all), length(thr_all)))

print(within(out, {
  sigma_init_med           <- round(sigma_init_med, 3)
  sigma_init_q025          <- round(sigma_init_q025, 3)
  sigma_init_q975          <- round(sigma_init_q975, 3)
  thr_median               <- round(thr_median, 3)
  plots_declining_at_med   <- round(plots_declining_at_med, 3)
  draws_declining_med_plot <- round(draws_declining_med_plot, 3)
  frac_declining           <- round(frac_declining, 3)
}), row.names = FALSE)

cat("\nplots_declining_at_med   = share of plots whose pre-run declines at the posterior-median sigma_init\n")
cat("draws_declining_med_plot = share of posterior draws that decline for the median plot\n")
cat("frac_declining           = share of (draw x plot) pairs that decline\n\n")
