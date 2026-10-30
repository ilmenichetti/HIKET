# =============================================================================
# build_F_trends_by_region.R   (2026-10-30)
#
# Measured AND modelled SOC trends, South and North side by side, on the plots
# measured in all three campaigns. Decided with Samuli (2026-10-30): the regional
# split is shown for measured and modelled trends alike -- this reverses the
# 2026-10-26 "paper stays national" decision for the trend figures.
#
# Basis: identical to F3_mean_soc / S_soc_by_region -- balanced plot set, whole
# profile, unweighted within region (every plot in a region has the same design
# weight, so the plot mean is also the area mean), campaign markers at their mean
# true sampling year. "North" = design weight 3 (soc_homogenized_plot.csv).
# Model curves: production posterior predictive, per draw the mean over the region's
# plots, then posterior mean + 95% band. Rates: paired, each plot on its OWN
# campaign years (as eqinit_common.R), measured and modelled alike.
#
# Run from repo root:  Rscript manuscript/figures/build_F_trends_by_region.R
# Output: manuscript/figures/F_trends_by_region.png
#         manuscript/figures/F_trends_by_region_rates.csv
# =============================================================================
suppressMessages({
  source("manuscript/figures/run_ids.R"); source("manuscript/figures/obs_basis.R")
  source("manuscript/figures/model_palette.R")
})

# === Regions =================================================================
om  <- readRDS(sprintf("Data/model_inputs/Yasso20_inputs_%s.rds", RID[["Yasso20"]]))$obs_meta
BAL <- balanced_plots(om)
P   <- read.csv("Data/SOC_homogeneized/soc_homogenized_plot.csv")
WT  <- tapply(P$weight, P$plot_id, function(w) w[1])
REG <- list(South = BAL[WT[as.character(BAL)] == 1], North = BAL[WT[as.character(BAL)] == 3])
stopifnot(!anyNA(unlist(REG)), sum(lengths(REG)) == length(BAL))
message("F_trends_by_region ", basis_note(BAL), "; South ", length(REG$South),
        ", North ", length(REG$North))

# === Per-plot campaign years and measured stocks =============================
obs <- do.call(rbind, lapply(as.character(BAL), function(p) {
  z <- om[[p]]; data.frame(plot_id = p, year = 1984L + z$idx, soc = z$soc_obs) }))
obs$camp <- campaign_of(obs$year)
stopifnot(!anyDuplicated(obs[c("plot_id", "camp")]))
YRS <- reshape(obs[c("plot_id", "camp", "year")], idvar = "plot_id", timevar = "camp",
               direction = "wide")                       # year.1, year.2, year.3
SOC <- reshape(obs[c("plot_id", "camp", "soc")],  idvar = "plot_id", timevar = "camp",
               direction = "wide")
INTERVALS <- list("1985-2006" = c(1, 2), "2006-2024" = c(2, 3), "1985-2024" = c(1, 3))

# paired rate per plot; returns mean and 95% interval over plots (measured) or draws (model)
obs_rate <- function(pl) t(sapply(INTERVALS, function(ab) {
  i  <- match(as.character(pl), SOC$plot_id)
  d  <- (SOC[i, paste0("soc.", ab[2])] - SOC[i, paste0("soc.", ab[1])]) /
        (YRS[i, paste0("year.", ab[2])] - YRS[i, paste0("year.", ab[1])])
  se <- sd(d) / sqrt(length(d)); c(m = mean(d), lo = mean(d) - 1.96 * se, hi = mean(d) + 1.96 * se) }))

# === Model trajectories and rates, per region =================================
YR <- 1985:2024
model_region <- function(m) {
  b  <- readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_predictive_%s.rds", m, RID[[m]]))
  pp <- b$posterior_predictions[, c("plot_id", "year", "draw", "total_soc")]
  pp <- pp[pp$plot_id %in% as.character(BAL) & pp$year %in% YR, ]
  lapply(REG, function(pl) {
    x  <- pp[pp$plot_id %in% as.character(pl), ]
    H  <- tapply(x$total_soc, list(x$draw, x$year), mean)          # draw x year region mean
    traj <- data.frame(year = as.integer(colnames(H)), m = colMeans(H),
                       lo = apply(H, 2, quantile, 0.025), hi = apply(H, 2, quantile, 0.975))
    i  <- match(as.character(pl), YRS$plot_id)
    at <- function(k) {                                            # draw x plot at campaign k
      key <- paste(pl, YRS[i, paste0("year.", k)])
      s   <- x[paste(x$plot_id, x$year) %in% key, ]
      tapply(s$total_soc, list(s$draw, s$plot_id), mean)[, as.character(pl), drop = FALSE]
    }
    C <- lapply(1:3, at)
    rate <- t(sapply(INTERVALS, function(ab) {
      dt <- YRS[i, paste0("year.", ab[2])] - YRS[i, paste0("year.", ab[1])]
      r  <- rowMeans(sweep(C[[ab[2]]] - C[[ab[1]]], 2, dt, "/"))
      c(m = mean(r), lo = quantile(r, 0.025, names = FALSE), hi = quantile(r, 0.975, names = FALSE)) }))
    list(traj = traj, rate = rate)
  })
}
MOD <- setNames(lapply(MODEL_ORDER, model_region), MODEL_ORDER)

# === Rates table ==============================================================
tab <- do.call(rbind, lapply(names(REG), function(r) {
  o <- obs_rate(REG[[r]])
  rbind(data.frame(region = r, source = "measured", interval = rownames(o), o, row.names = NULL),
        do.call(rbind, lapply(MODEL_ORDER, function(m) {
          x <- MOD[[m]][[r]]$rate
          data.frame(region = r, source = m, interval = rownames(x), x, row.names = NULL) })))
}))
tab$n_plots <- lengths(REG)[tab$region]
write.csv(transform(tab, m = round(m, 3), lo = round(lo, 3), hi = round(hi, 3)),
          "manuscript/figures/F_trends_by_region_rates.csv", row.names = FALSE)
print(reshape(transform(tab[c("region", "source", "interval", "m")], m = round(m, 2)),
              idvar = c("region", "source"), timevar = "interval", direction = "wide"), row.names = FALSE)

# === Figure ===================================================================
# Print size: full text width (6.53 in), 9 pt. Shared y SPAN so slopes compare by eye.
CM  <- lapply(REG, function(pl) obs_campaigns(om, pl))
rng <- lapply(names(REG), function(r) range(CM[[r]]$lo, CM[[r]]$hi,
              sapply(MODEL_ORDER, function(m) range(MOD[[m]][[r]]$traj[c("lo", "hi")]))))
span <- max(sapply(rng, diff))
# Top row: trajectories. Bottom row: paired rates per interval, measured vs each model,
# on a y axis shared by both regions (rates, unlike levels, compare directly).
png("manuscript/figures/F_trends_by_region.png", width = 6.53, height = 5.6, units = "in",
    res = 300, pointsize = 9)
layout(matrix(1:4, 2, byrow = FALSE), heights = c(1.15, 1))
par(mar = c(3.2, 3.6, 1.6, 0.5), mgp = c(2.0, 0.5, 0), las = 1, tcl = -0.3)
par(cex = 1)
SRC  <- c("measured", MODEL_ORDER)
SCOL <- c(measured = "firebrick", MODEL_COL)
ryl  <- range(tab$lo, tab$hi)
for (r in names(REG)) {
  yl <- mean(rng[[match(r, names(REG))]]) + c(-0.5, 0.5) * span
  par(mar = c(3.2, 3.6, 1.6, 0.5))
  plot(NA, xlim = range(YR), ylim = yl, xlab = "Year",
       ylab = if (r == "South") "Mean SOC across plots (tC/ha)" else "")
  mtext(sprintf("%s (%d plots)", r, length(REG[[r]])), side = 3, line = 0.3, adj = 0, font = 2)
  for (m in MODEL_ORDER) { d <- MOD[[m]][[r]]$traj
    polygon(c(d$year, rev(d$year)), c(d$lo, rev(d$hi)), col = adjustcolor(MODEL_COL[m], 0.10), border = NA) }
  for (m in MODEL_ORDER) { d <- MOD[[m]][[r]]$traj; lines(d$year, d$m, col = MODEL_COL[m], lwd = 2) }
  cm <- CM[[r]]
  arrows(cm$year, cm$lo, cm$year, cm$hi, angle = 90, code = 3, length = 0.04, col = "black", lwd = 1.6)
  points(cm$year, cm$m, pch = 21, bg = "firebrick", col = "black", cex = 1.4, lwd = 1.1)
  if (r == "North") legend("bottomright", legend = MODEL_ORDER, col = MODEL_COL[MODEL_ORDER],
                           lwd = 2.2, bty = "n", cex = 0.85, ncol = 2)
  # --- rates ---
  par(mar = c(2.2, 3.6, 0.6, 0.5))
  plot(NA, xlim = c(0.5, 3.5), ylim = ryl, xaxt = "n", xlab = "",
       ylab = if (r == "South") "Rate of change (tC/ha/yr)" else "")
  axis(1, at = 1:3, labels = names(INTERVALS), tick = FALSE)
  abline(h = 0, col = "grey50"); abline(v = c(1.5, 2.5), col = "grey85")
  off <- seq(-0.33, 0.33, length.out = length(SRC))
  for (k in seq_along(INTERVALS)) for (j in seq_along(SRC)) {
    x <- tab[tab$region == r & tab$source == SRC[j] & tab$interval == names(INTERVALS)[k], ]
    segments(k + off[j], x$lo, k + off[j], x$hi, col = SCOL[SRC[j]], lwd = 1.6)
    points(k + off[j], x$m, pch = if (j == 1) 21 else 19, bg = "firebrick",
           col = if (j == 1) "black" else SCOL[SRC[j]], cex = if (j == 1) 1.3 else 1.0)
  }
}
dev.off(); cat("wrote F_trends_by_region.png\n")
