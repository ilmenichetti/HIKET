source("manuscript/figures/run_ids.R")   # auto-selects current RUN_IDs
# S7 (appendix) -- posterior vs prior MARGINALS for the two most story-relevant
# auxiliary parameters: sigma_input (the input multiplier -- the Trickster) and
# sigma_init (the below-equilibrium 1917 start -- the protagonist). Both are clearly
# IDENTIFIED (posterior shifted/tightened off the prior) = the counterpoint to the
# non-identified transfer fractions (F8). Faithful prior via the engine setup; engine
# plot_one_marginal_honest styling.
#
# 2026-08-20: extended from Yasso20 alone to ALL THREE Yasso versions (6 panels,
# rows = model, columns = the two sigmas) and annotated with central tendency --
# MEAN as a solid rule, MEDIAN as a dashed rule, in the same colours as the
# distributions they belong to (class colour for posterior, grey for prior).
# Mean and median separate here because both sigmas are right-skewed on the
# natural scale, and the gap between the two rules IS that skew.

MODELS <- c("Yasso07", "Yasso15", "Yasso20")
PARAMS <- c("sigma_input", "sigma_init")

setup_model <- function(MODEL, RID) {
  sp  <- file.path("Calibration_real_data_transient", sprintf("run_%s_transient_calibration.R", MODEL))
  src <- readLines(sp, warn = FALSE); cut <- grep("^t_run <- system.time\\(\\{", src)[1]
  e   <- new.env(parent = globalenv())
  suppressWarnings(suppressMessages(
    source(textConnection(paste(src[seq_len(cut - 1L)], collapse = "\n")), local = e)))
  bx <- get("best_x", e); sp2 <- get("sigma_ppm", e); to <- get("to_original", e)
  set.seed(99); pr <- sapply(seq_along(bx), function(j) rnorm(3000L, bx[j], sp2[j]))
  colnames(pr) <- names(bx); prior <- t(apply(pr, 1, to))
  post <- readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_%s.rds", MODEL, RID))
  list(prior = prior, post = post,
       marg  = get("plot_one_marginal_honest", e),
       cols  = get("class_cols", e), cls = get("classify_param", e))
}

# Solid = mean, dashed = median, coloured to match the distribution it summarises.
add_central <- function(v, col, lwd = 2) {
  abline(v = mean(v),   col = col, lwd = lwd, lty = 1)
  abline(v = median(v), col = col, lwd = lwd, lty = 2)
}

M <- list()
for (m in MODELS) { cat("sourcing", m, "setup...\n"); M[[m]] <- setup_model(m, RID[[m]]) }

png("manuscript/figures/S7_param_marginals.png",
    width = 10, height = 12.6, units = "in", res = 200)
par(mfrow = c(3, 2), mar = c(4.6, 4.6, 3.2, 1), mgp = c(2.7, 0.8, 0),
    cex.axis = 1.05, cex.lab = 1.2)

summ <- data.frame()
for (m in MODELS) {
  for (nm in PARAMS) {
    pv <- M[[m]]$post [is.finite(M[[m]]$post [, nm]), nm]
    qv <- M[[m]]$prior[is.finite(M[[m]]$prior[, nm]), nm]
    M[[m]]$marg(pv, qv, nm, m)

    post_col  <- unname(M[[m]]$cols[M[[m]]$cls(nm)])
    add_central(qv, "grey55")     # prior  -- matches the grey backdrop
    add_central(pv, post_col)     # posterior -- class colour, drawn on top

    legend("topleft", legend = c("mean", "median"), lty = c(1, 2), lwd = 2,
           col = "grey25", bty = "n", cex = 0.7)

    summ <- rbind(summ, data.frame(
      model = m, param = nm,
      prior_mean = mean(qv),   prior_med = median(qv),
      post_mean  = mean(pv),   post_med  = median(pv),
      skew_gap   = mean(pv) - median(pv)))
  }
}
dev.off()

cat("\nwrote S7_param_marginals.png (", length(MODELS) * length(PARAMS), " panels)\n", sep = "")
print(format(summ, digits = 3), row.names = FALSE)
