#!/usr/bin/env Rscript
# =============================================================================
# preview_marginal_style.R  --  preview/test the restyled marginal panels
# -----------------------------------------------------------------------------
# The cosmetic restyle of plot_one_marginal_honest() in
# Calibration_real_data/calibration_engine.R (light-grey filled prior backdrop +
# semi-transparent class-coloured posterior) only takes effect when a model is
# re-calibrated. This script previews the new look NOW, without re-running MCMC,
# by sourcing the engine's functions and feeding them the already-saved SP1
# posterior plus a reconstructed prior. SP1 is used because its transforms are
# simple (log / identity), so the prior pushforward is exact without the
# stick-break machinery.
#
# Run from the project root:  Rscript doublechecks/preview_marginal_style.R
# =============================================================================

OUT <- "doublechecks/figures"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# Source the engine for classify_param(), class_cols, plot_one_marginal_honest()
suppressMessages(source("Calibration_real_data/calibration_engine.R"))

post <- readRDS("Calibration_real_data_transient/runs/SP1_posterior_20260608_015026.rds")
pkg  <- readRDS("Data/model_inputs/SP1_inputs_20260608_015026.rds")
best_x <- pkg$best_x; sppm <- pkg$sigma_ppm

# SP1 transforms (from run_SP1_transient_calibration.R param_spec)
tf <- c(alpha = "log", beta1 = "log", beta2 = "id",
        gamma = "id", sigma_init = "log", sigma_input = "log")
set.seed(99)
n_prior <- 3000L
prior_phys <- sapply(names(best_x), function(nm) {
  z <- rnorm(n_prior, best_x[[nm]], sppm[[nm]])
  if (tf[[nm]] == "log") exp(z) else z
})

png(file.path(OUT, "SP1_marginals_RESTYLED_preview.png"),
    width = 6 * 150, height = 6 * 150, res = 150)
par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))
for (nm in colnames(post)) {
  plot_one_marginal_honest(post[, nm], prior_phys[, nm], nm, "SP1")
}
dev.off()
message("Preview written: ", file.path(OUT, "SP1_marginals_RESTYLED_preview.png"))
