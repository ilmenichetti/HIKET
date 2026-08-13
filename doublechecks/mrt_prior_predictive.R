# =============================================================================
# mrt_prior_predictive.R — is the short MRT LEARNED, or INHERITED from the prior?
#
# THE QUESTION. The posterior marginal peak is set by likelihood x prior x VOLUME.
# The third term has nothing to do with fit: intrinsic MRT is a many-to-one
# function of ~15 parameters, so if short MRT is simply reachable by more
# parameter combinations, the marginal peaks there even under a flat likelihood.
# In ~20 dimensions that effect is not small.
#
# THE TEST. Push the PRIOR draws through the same MRT function and compare:
#     published point   -- the single parameter set from the literature
#     PRIOR predictive  -- what MRT the prior alone implies, volume included
#     POSTERIOR         -- what we end up with
#
#   * prior predictive already at ~15-20 yr  => short MRT is INHERITED. The data
#     are not demanding fast turnover; our prior construction is. That would
#     undercut "the displacement is demanded by the data".
#   * prior predictive at ~33 yr, posterior pulled down to 15 => the DATA are
#     doing the work, and the displacement is real.
#
# Note the published POINT can sit far from the prior-predictive MEDIAN even
# though the prior is centred on the published parameters: MRT is a nonlinear
# many-to-one map, so centring the parameters does not centre the MRT.
#
# Prior recipe identical to run_diagnostics() and to S11/S12: seed 99,
# rnorm(best_x, sigma_ppm), pushed through to_original.
#
# Run from repo root:  Rscript doublechecks/mrt_prior_predictive.R
# =============================================================================

suppressMessages(library(BayesianTools))
source("doublechecks/intrinsic_mrt_lib.R")     # setup(), mrt_fun(), ref
options(width = 118)

RID <- c(Yasso07 = "20260812_080941", Yasso15 = "20260812_080940", Yasso20 = "20260812_080940")
PUB <- c(Yasso07 = 33.47, Yasso15 = 30.38, Yasso20 = 19.03)
N_PRIOR <- 3000L; N_POST <- 3000L
set.seed(2025)

OUT <- list()
for (M in names(RID)) {
  e <- setup(M); f <- mrt_fun(M, e)
  bx <- get("best_x", e); spm <- get("sigma_ppm", e)

  # --- prior predictive -------------------------------------------------------
  set.seed(99)
  pr <- sapply(seq_along(bx), function(j) rnorm(N_PRIOR, bx[j], spm[j]))
  colnames(pr) <- names(bx)
  mrt_pr <- vapply(seq_len(N_PRIOR), function(i)
    tryCatch(f(e$.to_original(pr[i, ])), error = function(z) NA_real_), numeric(1))
  mrt_pr <- mrt_pr[is.finite(mrt_pr) & mrt_pr > 0 & mrt_pr < 1e4]

  # --- posterior (chains are SAMPLING space; must go through to_original) -----
  s <- getSample(readRDS(sprintf("Calibration_real_data_transient/runs/%s_chains_%s.rds",
                                 M, RID[[M]])), parametersOnly = FALSE)
  i <- sample(nrow(s), min(N_POST, nrow(s)))
  mrt_po <- vapply(i, function(k)
    tryCatch(f(e$.to_original(s[k, names(bx)])), error = function(z) NA_real_), numeric(1))
  mrt_po <- mrt_po[is.finite(mrt_po)]

  OUT[[M]] <- list(prior = mrt_pr, post = mrt_po, pub = PUB[[M]])
  q <- function(v) quantile(v, c(0.05, 0.5, 0.95))
  cat(sprintf("\n=== %s ===\n", M))
  cat(sprintf("  published point     : %6.2f\n", PUB[[M]]))
  cat(sprintf("  PRIOR predictive    : %6.2f   90%% [%5.2f, %6.2f]   n=%d (%.0f%% of draws usable)\n",
              q(mrt_pr)[2], q(mrt_pr)[1], q(mrt_pr)[3], length(mrt_pr), 100*length(mrt_pr)/N_PRIOR))
  cat(sprintf("  POSTERIOR           : %6.2f   90%% [%5.2f, %6.2f]   n=%d\n",
              q(mrt_po)[2], q(mrt_po)[1], q(mrt_po)[3], length(mrt_po)))
  cat(sprintf("  -> prior median vs published : %+6.2f yr  (%s)\n",
              q(mrt_pr)[2] - PUB[[M]],
              if (abs(q(mrt_pr)[2] - PUB[[M]]) < 0.1 * PUB[[M]]) "prior is centred on published MRT"
              else "the MAP FROM PARAMETERS TO MRT ALREADY SHIFTS IT"))
  cat(sprintf("  -> posterior vs prior median : %+6.2f yr  (%s)\n",
              q(mrt_po)[2] - q(mrt_pr)[2],
              if (abs(q(mrt_po)[2] - q(mrt_pr)[2]) < 0.1 * q(mrt_pr)[2]) "DATA ADD ALMOST NOTHING -- inherited"
              else "the data move it"))
}

# --- figure ------------------------------------------------------------------
png("doublechecks/mrt_prior_predictive.png", width = 11.5, height = 4.2, units = "in", res = 190)
par(mfrow = c(1, 3), mar = c(4.2, 4.4, 3.2, 1.0), mgp = c(2.6, 0.7, 0), las = 1)
for (M in names(OUT)) {
  o <- OUT[[M]]
  xr <- range(c(quantile(o$prior, c(0.01, 0.99)), o$post, o$pub))
  dp <- density(o$prior, from = xr[1], to = xr[2]); dq <- density(o$post, from = xr[1], to = xr[2])
  plot(NA, xlim = xr, ylim = c(0, max(dp$y, dq$y) * 1.08), xlab = "intrinsic MRT (yr)",
       ylab = "density", main = M)
  polygon(c(dp$x, rev(dp$x)), c(dp$y, rep(0, length(dp$y))),
          col = adjustcolor("grey55", 0.45), border = "grey40")
  polygon(c(dq$x, rev(dq$x)), c(dq$y, rep(0, length(dq$y))),
          col = adjustcolor("#C1553B", 0.65), border = "#C1553B")
  abline(v = o$pub, lty = 2, lwd = 2, col = "grey15")
  legend("topright", bty = "n", cex = 0.78,
         fill = c(adjustcolor("grey55", 0.45), adjustcolor("#C1553B", 0.65)), border = NA,
         legend = c(sprintf("prior pred. (%.1f)", median(o$prior)),
                    sprintf("posterior (%.1f)", median(o$post))))
  text(o$pub, max(dp$y, dq$y), " published", srt = 90, adj = c(1, -0.4), cex = 0.72, col = "grey15")
}
dev.off()
cat("\nwrote doublechecks/mrt_prior_predictive.png\n")
saveRDS(OUT, "doublechecks/mrt_prior_predictive.rds")
