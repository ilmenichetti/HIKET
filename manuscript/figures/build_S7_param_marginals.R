setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
# S7 (appendix) -- posterior vs prior MARGINALS for the two most story-relevant
# parameters (Yasso20): sigma_input (the input multiplier -- the Trickster) and
# sigma_init (the below-equilibrium 1917 start -- the protagonist). Both are clearly
# IDENTIFIED (posterior shifted/tightened off the prior) = the counterpoint to the
# non-identified transfer fractions (F8). Faithful prior via the engine setup; engine
# plot_one_marginal_honest styling. (More parameters can be added later.)
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
  list(prior = prior, post = post, marg = get("plot_one_marginal_honest", e))
}
cat("sourcing Yasso20 setup...\n")
M <- setup_model("Yasso20", "20260710_102431")

png("manuscript/figures/S7_param_marginals.png", width = 10, height = 4.6, units = "in", res = 200)
par(mfrow = c(1, 2), mar = c(4.6, 4.6, 3.2, 1), mgp = c(2.7, 0.8, 0), cex.axis = 1.05, cex.lab = 1.2)
for (nm in c("sigma_input", "sigma_init")) {
  pv <- M$post[is.finite(M$post[, nm]), nm]; qv <- M$prior[is.finite(M$prior[, nm]), nm]
  M$marg(pv, qv, nm, "Yasso20")
}
dev.off(); cat("wrote S7_param_marginals.png\n")
