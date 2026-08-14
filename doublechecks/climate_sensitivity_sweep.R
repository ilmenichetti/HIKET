# =============================================================================
# climate_sensitivity_sweep.R   (2026-08-14)
#
# WHY. The 2.3x spread in climate sensitivity across the Yasso family does not
# appear in ANY existing figure -- and equifinality_forecast.png actively hides
# it, because its three panels carry independent x-axes. This puts all models on
# one pair of shared axes, with the posterior propagated rather than a median, so
# the between-model divergence and the within-model uncertainty are visible
# together and can be compared by eye.
#
# It is also a prototype of the "climate reactivity" panel of the two-calibration
# experiment recorded in the storyline note (restrictive vs permissive
# sigma_input): the same sweep, run on two calibrations, is that comparison.
#
# WHAT IT SHOWS, per model, over a warming sweep:
#   (A) equilibrium fractional SOC change  -- where the system ends up
#   (B) change realised by year HORIZ      -- where it has GOT to by then;
#       this is the panel MRT governs, because MRT is the time constant
#
# ⚠ Same approximations as equifinality_forecast.R: single-exponential transient
# on the bulk MRT, step change rather than a trajectory, common J_bar for scale.
# A DEMONSTRATION of divergence, not a projection to quote.
#
# Usage:  Rscript doublechecks/climate_sensitivity_sweep.R [N_DRAW] [DTMAX]
# =============================================================================

suppressWarnings(suppressMessages({
  a <- commandArgs(trailingOnly = TRUE)
  N_DRAW <- if (length(a) >= 1) as.integer(a[[1]]) else 600L
  DTMAX  <- if (length(a) >= 2) as.numeric(a[[2]]) else 5.0
  library(BayesianTools)
}))
set.seed(2025)
J_BAR <- 2.511
HORIZ <- 50
DTS   <- seq(0, DTMAX, by = 0.5)

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

mrt_fun <- function(M, e) {
  if (M == "Yasso07") {
    ss <- get("yasso07_steady_state", e); cxm <- get("compute_xi_mean_yasso07", e)
    function(p, clim) {
      mp <- e$.assemble(p)
      sum(ss(mp, ref$nwl, ref$fwl, ref$cwl,
             cxm(clim, mp[["beta1"]], mp[["beta2"]], mp[["gamma"]])))
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

sweep <- list()
for (M in names(RID)) {
  e <- setup(M); f <- mrt_fun(M, e); p_def <- e$.to_original(get("best_x", e))
  ch <- readRDS(sprintf("Calibration_real_data_transient/runs/%s_chains_%s.rds", M, RID[[M]]))
  S  <- do.call(rbind, lapply(ch, function(s) getSample(s, parametersOnly = FALSE, start = 2)))
  xn <- names(get("best_x", e))
  idx <- round(seq(1, nrow(S), length.out = min(N_DRAW, nrow(S))))

  # MRT(draw, dT) and the draw's sigma_input
  P  <- lapply(idx, function(i) e$.to_original(S[i, xn]))
  si <- vapply(P, function(p) unname(p[["sigma_input"]]), numeric(1))
  Mm <- vapply(DTS, function(dt) {
    cl <- ref$clim; cl$temp_mean <- cl$temp_mean + dt
    vapply(P, function(p) tryCatch(f(p, cl), error = function(z) NA_real_), numeric(1))
  }, numeric(length(P)))                       # draws x length(DTS)

  ok <- apply(is.finite(Mm), 1, all)
  Mm <- Mm[ok, , drop = FALSE]; si <- si[ok]
  C  <- si * J_BAR * Mm                        # equilibrium stock at each dT
  eq <- C / C[, 1] - 1                         # (A) fractional equilibrium change
  tr <- (C - C[, 1]) * (1 - exp(-HORIZ / Mm[, 1]))   # (B) realised by year HORIZ
  sweep[[M]] <- list(eq = eq, tr = tr, mrt0 = Mm[, 1], n = nrow(Mm))
  cat(sprintf("%-8s n=%4d | MRT0 %.1f yr | at +2C: eq %+.1f%%  realised %+.2f tC/ha\n",
              M, nrow(Mm), median(Mm[, 1]),
              100*median(eq[, DTS == 2]), median(tr[, DTS == 2])))
}
saveRDS(list(DTS = DTS, sweep = sweep), "doublechecks/climate_sensitivity_sweep.rds")

cols <- c(Yasso07 = "#0072B2", Yasso15 = "#D55E00", Yasso20 = "#009E73")
png("doublechecks/climate_sensitivity_sweep.png", width = 1900, height = 900, res = 190)
par(mfrow = c(1, 2), mar = c(4.4, 4.6, 2.4, 1.0), oma = c(0, 0, 2.8, 0))

bandplot <- function(get, ylab, main, pct = FALSE) {
  ys <- lapply(sweep, function(s) apply(get(s), 2, quantile, c(.025,.5,.975), na.rm = TRUE))
  yl <- range(unlist(ys)) * if (pct) 100 else 1
  plot(NA, xlim = range(DTS), ylim = yl, xlab = "warming (°C)", ylab = ylab, main = main)
  abline(h = 0, col = "grey75")
  for (M in names(sweep)) {
    q <- ys[[M]] * if (pct) 100 else 1
    polygon(c(DTS, rev(DTS)), c(q[1, ], rev(q[3, ])), border = NA,
            col = adjustcolor(cols[M], 0.18))
    lines(DTS, q[2, ], col = cols[M], lwd = 2.4)
  }
  legend("bottomleft", names(sweep), col = cols[names(sweep)], lwd = 2.4, bty = "n", cex = 0.8)
}
bandplot(function(s) s$eq, "equilibrium SOC change (%)", "A. where it ends up", pct = TRUE)
bandplot(function(s) s$tr, sprintf("SOC change by year %d (tC/ha)", HORIZ),
         sprintf("B. where it has got to by year %d", HORIZ))
mtext(paste0("Climate sensitivity differs ~2x across the family, and the bands are posterior spread ",
             "— DEMONSTRATION, not a projection"),
      side = 3, line = 0.6, outer = TRUE, cex = 0.75, col = "grey20")
dev.off()
cat("wrote doublechecks/climate_sensitivity_sweep.{rds,png}\n")
