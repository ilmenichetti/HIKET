# =============================================================================
# equifinality_forecast.R   (2026-08-14)
#
# QUESTION. High-input/short-MRT and low-input/long-MRT fit the SOC data equally
# well -- the MRT x sigma_input ridge. Does that indifference PROPAGATE INTO THE
# FORECAST, or do the two ends of the ridge project the same future?
#
# THE PREDICTION BEING TESTED, in two parts:
#   (a) EQUILIBRIUM response is DEGENERATE. At steady state C = J/(k*xi), so a
#       warming that multiplies xi by f gives C' = C/f whatever the J/k split.
#       => the FRACTIONAL equilibrium change should be uncorrelated with
#          sigma_input, i.e. constant along the ridge.
#   (b) TRANSIENT response is NOT degenerate. MRT *is* the time constant: after
#       t years a step change is realised to 1 - exp(-t/MRT). Short-MRT draws
#       have nearly finished responding; long-MRT draws are still moving.
#       => the ABSOLUTE change realised by year t should vary along the ridge
#          even where the likelihood cannot tell the draws apart.
#
# METHOD. Per posterior draw, evaluate the model's own steady-state routine at a
# fixed reference condition (same construction as intrinsic_mrt.R) twice: at the
# dataset-mean climate, and at that climate + DT. Both MRTs are properties of the
# generator, uncontaminated by sigma_init or the fit. The draw's own likelihood
# comes from the chains, so "indistinguishable in fit" is measured, not assumed.
#
# ⚠ APPROXIMATIONS, stated because they matter for how far this can be pushed:
#   * the transient uses a SINGLE exponential with the bulk MRT as time constant.
#     A multi-pool system is a sum of exponentials, so this is a first-order
#     illustration of timing, not a substitute for a forward run.
#   * a step change in climate, not a trajectory.
#   * a common J_bar across models for the absolute scale.
#   It is deliberately a demonstration that the ridge maps onto the forecast --
#   NOT a projection to be quoted.
#
# Usage:  Rscript doublechecks/equifinality_forecast.R [N_DRAW] [DT]
# =============================================================================

suppressWarnings(suppressMessages({
  a <- commandArgs(trailingOnly = TRUE)
  N_DRAW <- if (length(a) >= 1) as.integer(a[[1]]) else 1500L
  DT     <- if (length(a) >= 2) as.numeric(a[[2]]) else 2.0     # warming, degrees C
  library(BayesianTools)
}))
set.seed(2025)
J_BAR <- 2.511      # common scale, tC/ha/yr (Boris tree-litter mean)
HORIZ <- 50         # years, for the transient comparison

RID <- c(Yasso07 = "20260812_080941", Yasso15 = "20260812_080940", Yasso20 = "20260812_080940")

setup <- function(M) {
  src <- readLines(sprintf("Calibration_real_data_transient/run_%s_transient_calibration.R", M),
                   warn = FALSE)
  cut <- grep("^t_run <- system.time", src)[1]
  e <- new.env(parent = globalenv())
  invisible(capture.output(suppressMessages(
    source(textConnection(paste(src[seq_len(cut-1)], collapse="\n")), local = e))))
  st <- grep("ll_fn <- make_likelihood", src)[1]; op <- 0L; en <- NA_integer_
  for (i in seq(st, length(src))) {
    ch <- strsplit(src[i], "")[[1]]; op <- op + sum(ch=="(") - sum(ch==")")
    if (op == 0L) { en <- i; break }
  }
  ml <- as.list(str2lang(sub("^\\s*ll_fn\\s*<-\\s*","",
                             paste(src[seq(st,en)], collapse="\n"))))[-1]
  e$.to_original <- eval(ml[["to_original"]], envir = e)
  e$.assemble    <- eval(ml[["assemble_params"]], envir = e)
  e
}

ref <- local({
  e <- setup("Yasso15")
  plots <- get("plots", e); cbp <- get("climate_by_plot", e); lms <- get("litter_means", e)
  clim <- data.frame(
    temp_mean      = mean(vapply(plots, function(p) mean(cbp[[p]]$temp_mean),      numeric(1))),
    temp_amplitude = mean(vapply(plots, function(p) mean(cbp[[p]]$temp_amplitude), numeric(1))),
    precip         = mean(vapply(plots, function(p) mean(cbp[[p]]$precip),         numeric(1))))
  gm <- function(f) rowMeans(vapply(plots, function(p) as.numeric(lms[[p]][[f]]), numeric(4)))
  nwl <- gm("nwl_mean"); fwl <- gm("fwl_mean"); cwl <- gm("cwl_mean")
  tot <- sum(nwl) + sum(fwl) + sum(cwl)
  list(clim = clim, nwl = nwl/tot, fwl = fwl/tot, cwl = cwl/tot)
})
cat(sprintf("Reference: T = %.2f C (warming +%.1f C) | precip = %.0f mm | horizon %d yr\n\n",
            ref$clim$temp_mean, DT, ref$clim$precip, HORIZ))

# MRT at an arbitrary climate -- same construction as intrinsic_mrt.R
mrt_fun <- function(M, e) {
  if (M == "Yasso07") {
    ss <- get("yasso07_steady_state", e); cxm <- get("compute_xi_mean_yasso07", e)
    function(p, clim) {
      mp <- e$.assemble(p)
      xi <- cxm(clim, mp[["beta1"]], mp[["beta2"]], mp[["gamma"]])
      sum(ss(mp, ref$nwl, ref$fwl, ref$cwl, xi))
    }
  } else {
    ss <- get("yasso15_steady_state", e); cxm <- get("compute_xi_mean_yasso15", e)
    .ypn <- sprintf("%s_PARAM_NAMES", toupper(M))
    YP <- if (exists(.ypn, envir = e, inherits = FALSE)) get(.ypn, envir = e) else NULL
    function(p, clim) {
      mp <- e$.assemble(p)
      xi <- cxm(clim_ss = clim, params = if (is.null(YP)) mp else mp[YP])
      sum(ss(mp, ref$nwl, ref$fwl, ref$cwl, xi, precip_mean = clim$precip))
    }
  }
}

clim_warm <- ref$clim; clim_warm$temp_mean <- clim_warm$temp_mean + DT

res <- list()
for (M in names(RID)) {
  e <- setup(M); f <- mrt_fun(M, e); p_def <- e$.to_original(get("best_x", e))

  # Per-draw parameters AND likelihood; start=2 drops the DEzs burn-in artefact.
  # ⚠ the CHAINS store parameters in the UNCONSTRAINED space (the posterior RDS
  # that intrinsic_mrt.R reads is already physical, but it carries no likelihood).
  # So every draw must go through to_original() -- skipping that silently yields
  # MRT = 0 and a sigma_input that is really a logit coordinate.
  ch <- readRDS(sprintf("Calibration_real_data_transient/runs/%s_chains_%s.rds", M, RID[[M]]))
  S  <- do.call(rbind, lapply(ch, function(s) getSample(s, parametersOnly = FALSE, start = 2)))
  xnames <- names(get("best_x", e))
  stopifnot(all(xnames %in% colnames(S)))
  idx <- round(seq(1, nrow(S), length.out = min(N_DRAW, nrow(S))))

  o <- t(vapply(idx, function(i) {
    p <- e$.to_original(S[i, xnames])
    m0 <- tryCatch(f(p, ref$clim),  error = function(z) NA_real_)
    m1 <- tryCatch(f(p, clim_warm), error = function(z) NA_real_)
    # unname(): a named scalar from S[] would turn the column into "ll.Llikelihood"
    c(mrt0 = m0, mrt1 = m1, si = unname(p[["sigma_input"]]),
      ll = unname(S[i, "Llikelihood"]))
  }, numeric(4)))
  d <- as.data.frame(o); d <- d[is.finite(d$mrt0) & is.finite(d$mrt1) & is.finite(d$ll), ]

  d$C0     <- d$si * J_BAR * d$mrt0                 # equilibrium stock, baseline
  d$C1     <- d$si * J_BAR * d$mrt1                 # equilibrium stock, warmed
  d$eqfrac <- d$C1 / d$C0 - 1                       # (a) fractional equilibrium change
  d$dC50   <- (d$C1 - d$C0) * (1 - exp(-HORIZ / d$mrt0))   # (b) realised by year HORIZ
  res[[M]] <- d

  # The posterior IS the set of draws the data accept, so its spread is the
  # answer. A "within 1 nat" band is too thin to be informative at this
  # subsample size, so also report a band that actually contains draws, and how
  # much of the forecast spread fit quality explains at all (R^2 of ll).
  band <- d[d$ll >= max(d$ll) - 5, ]
  r2   <- summary(lm(dC50 ~ ll, data = d))$r.squared
  cat(sprintf("=== %s  (n=%d draws; %d within 5 nats of the best) ===\n", M, nrow(d), nrow(band)))
  cat(sprintf("  MRT %.1f -> %.1f yr under +%.1f C   |  sigma_input %.2f-%.2f\n",
              median(d$mrt0), median(d$mrt1), DT, quantile(d$si,.025), quantile(d$si,.975)))
  cat(sprintf("  (a) EQUILIBRIUM fractional change : %+.1f%%  [%.1f, %.1f]   corr w/ log(sigma_input) = %+.3f\n",
              100*median(d$eqfrac), 100*quantile(d$eqfrac,.025), 100*quantile(d$eqfrac,.975),
              cor(log(d$si), d$eqfrac)))
  cat(sprintf("  (b) REALISED by yr %d (tC/ha)      : %+.2f  [%.2f, %.2f]   corr w/ log(sigma_input) = %+.3f\n",
              HORIZ, median(d$dC50), quantile(d$dC50,.025), quantile(d$dC50,.975),
              cor(log(d$si), d$dC50)))
  cat(sprintf("  posterior forecast SPREAD (95%%)   : %.2f tC/ha = %.0f%% of the median loss\n",
              diff(quantile(d$dC50, c(.025,.975))),
              100*abs(diff(quantile(d$dC50, c(.025,.975))))/abs(median(d$dC50))))
  if (nrow(band) > 5)
    cat(sprintf("  within 5 nats of best fit         : %.2f to %.2f tC/ha\n",
                min(band$dC50), max(band$dC50)))
  cat(sprintf("  R^2 of log-likelihood on forecast : %.3f   <- fraction of forecast spread the FIT explains\n\n",
              r2))
}

saveRDS(res, "doublechecks/equifinality_forecast.rds")

# --- figure ------------------------------------------------------------------
png("doublechecks/equifinality_forecast.png", width = 2000, height = 900, res = 190)
par(mfrow = c(1, 3), mar = c(4.4, 4.4, 2.6, 1.2), oma = c(0, 0, 2.6, 0))
cols <- c(Yasso07 = "#0072B2", Yasso15 = "#D55E00", Yasso20 = "#009E73")
for (M in names(res)) {
  d <- res[[M]]
  rng <- range(d$ll); keep <- d$ll > quantile(d$ll, 0.02)   # trim the far tail for readability
  d <- d[keep, ]
  plot(d$dC50, d$ll, pch = 16, cex = 0.35, col = adjustcolor(cols[M], 0.35),
       xlab = sprintf("SOC change realised by year %d  (tC/ha)", HORIZ),
       ylab = "log-likelihood", main = M)
  abline(h = max(d$ll) - 5, lty = 2, col = "grey30")
  band <- d[d$ll >= max(d$ll) - 5, ]
  if (nrow(band) > 2) {
    rug(band$dC50, col = cols[M], lwd = 1.2)
    mtext(sprintf("within 5 nats: %.2f to %.2f tC/ha", min(band$dC50), max(band$dC50)),
          side = 3, line = 0.2, cex = 0.62, col = "grey25")
  }
}
mtext(sprintf("Draws that fit equally well project differently (+%.1f C step, horizon %d yr) — DEMONSTRATION, not a projection",
              DT, HORIZ), side = 3, line = 0.7, outer = TRUE, cex = 0.72, col = "grey20")
dev.off()
cat("wrote doublechecks/equifinality_forecast.{rds,png}\n")
