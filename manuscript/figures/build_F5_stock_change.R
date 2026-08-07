# =============================================================================
# build_F5_stock_change.R   (2026-08-07)
#
# THE INVENTORY FIGURE. Everything else in the set is about the SOC *stock*;
# a greenhouse-gas inventory reports the *change*. This is the only figure that
# shows it, and it is where the models fail hardest:
#
#   - all six UNDERestimate the observed sink
#   - Yasso15 and Yasso20 get the SIGN WRONG, reporting a source
#   - the flip is caused by sigma_init crossing the pre-run inversion threshold
#     (sigma_init > J_t0/J_full, median 0.818) -- see doublechecks/prerun_direction.R
#
# Paired plots only: a plot enters a period only if observed in BOTH campaigns,
# so changing plot composition cannot drive the comparison. Uncertainty on the
# observed change is the SE of the paired per-plot differences.
#
# Usage:  Rscript manuscript/figures/build_F5_stock_change.R
# =============================================================================

setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
source("manuscript/figures/run_ids.R")
source("manuscript/figures/model_palette.R")

campaign <- function(y) ifelse(y <= 1990, "c1985", ifelse(y <= 2015, "c2006", "c2024"))

# --- paired change per model, per period -------------------------------------
periods <- list("1985-2024" = c("c1985","c2024", 39),
                "2006-2024" = c("c2006","c2024", 18))

collect <- function(m) {
  f <- sort(list.files(file.path("Calibration_real_data_transient/diagnostics", m),
       pattern = sprintf("^%s_residuals_%s\\.csv$", m, RID[[m]]), full.names = TRUE),
       decreasing = TRUE)[1]
  if (is.na(f)) stop("no residuals CSV for ", m)
  d <- read.csv(f, stringsAsFactors = FALSE)
  d$plot_id <- as.character(d$plot_id); d$cp <- campaign(d$year)
  d <- d[!is.na(d$soc_obs_tCha), ]
  o <- reshape(d[, c("plot_id","cp","soc_obs_tCha")], idvar="plot_id", timevar="cp", direction="wide")
  p <- reshape(d[, c("plot_id","cp","soc_mean")],     idvar="plot_id", timevar="cp", direction="wide")
  merge(o, p, by = "plot_id")
}

res <- lapply(MODEL_ORDER, collect); names(res) <- MODEL_ORDER

stats <- do.call(rbind, lapply(names(periods), function(pn) {
  a <- periods[[pn]][1]; b <- periods[[pn]][2]; yrs <- as.numeric(periods[[pn]][3])
  do.call(rbind, lapply(MODEL_ORDER, function(m) {
    d  <- res[[m]]
    oa <- d[[paste0("soc_obs_tCha.",a)]]; ob <- d[[paste0("soc_obs_tCha.",b)]]
    pa <- d[[paste0("soc_mean.",a)]];     pb <- d[[paste0("soc_mean.",b)]]
    k  <- is.finite(oa)&is.finite(ob)&is.finite(pa)&is.finite(pb)
    od <- (ob[k]-oa[k])/yrs; pd <- (pb[k]-pa[k])/yrs
    data.frame(period = pn, model = m, n = sum(k),
               obs = mean(od), obs_se = sd(od)/sqrt(sum(k)),
               pred = mean(pd), pred_se = sd(pd)/sqrt(sum(k)),
               stringsAsFactors = FALSE)
  }))
}))

# --- plot ---------------------------------------------------------------------
png("manuscript/figures/F5_stock_change.png", width = 10, height = 5.6,
    units = "in", res = 200)
par(mfrow = c(1, 2), mar = c(4.4, 7.4, 2.6, 1.4), oma = c(0, 0, 2.6, 0),
    mgp = c(2.7, 0.7, 0), las = 1)

for (pn in names(periods)) {
  s  <- stats[stats$period == pn, ]
  s  <- s[match(MODEL_ORDER, s$model), ]
  ob <- s$obs[1]; ose <- s$obs_se[1]
  xr <- range(c(s$pred - 1.96*s$pred_se, s$pred + 1.96*s$pred_se,
                ob - 1.96*ose, ob + 1.96*ose, 0))
  xr <- xr + c(-0.06, 0.06) * diff(xr)
  yy <- rev(seq_along(MODEL_ORDER))

  plot(NA, xlim = xr, ylim = c(0.4, length(MODEL_ORDER) + 0.6),
       yaxt = "n", xlab = expression("Mean annual "*Delta*"SOC  (tC ha"^-1*" yr"^-1*")"),
       ylab = "", main = sprintf("%s   (n = %d paired plots)", pn, s$n[1]))
  # wrong-sign region: a model landing here reports a source, not a sink
  rect(xr[1] - 1, 0.2, 0, length(MODEL_ORDER) + 0.8,
       col = adjustcolor("#AA3333", 0.07), border = NA)
  # observed band, so points sit on top
  rect(ob - 1.96*ose, 0.2, ob + 1.96*ose, length(MODEL_ORDER) + 0.8,
       col = adjustcolor("#4477AA", 0.13), border = NA)
  abline(v = ob, col = "#4477AA", lwd = 2.2)
  abline(v = 0,  col = "grey35", lty = 2, lwd = 1.4)

  segments(s$pred - 1.96*s$pred_se, yy, s$pred + 1.96*s$pred_se, yy,
           col = MODEL_COL[s$model], lwd = 2.6)
  points(s$pred, yy, pch = 21, cex = 1.9, lwd = 1.6,
         bg = MODEL_COL[s$model], col = "grey20")
  axis(2, at = yy, labels = s$model, tick = FALSE, cex.axis = 1.05)

  # name the wrong-sign region rather than marking individual points
  if (any(s$pred < 0))
    mtext("reports a source", side = 3, line = -1.2, adj = 0.02,
          col = "#AA3333", cex = 0.78, font = 3)

  if (pn == names(periods)[1]) {
    legend("bottomright", bty = "n", cex = 0.82,
           legend = c("observed (paired)", "95% CI on observed", "zero"),
           col = c("#4477AA", adjustcolor("#4477AA", 0.35), "grey35"),
           lwd = c(2.2, 7, 1.4), lty = c(1, 1, 2))
  }
}
mtext("Modelled vs observed soil carbon stock change", outer = TRUE,
      line = 0.6, cex = 1.2, font = 2)
dev.off()

cat("\nWrote manuscript/figures/F5_stock_change.png\n\n")
print(within(stats, {obs <- round(obs,3); pred <- round(pred,3)
                     obs_se <- round(obs_se,3); pred_se <- round(pred_se,3)}),
      row.names = FALSE)
cat("\n✖ marks a model reporting a SOURCE where the observations show a SINK.\n")
