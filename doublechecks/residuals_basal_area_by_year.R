# =============================================================================
# residuals_basal_area_by_year.R — does the residual follow the stand of 1985
# or the stand of today?
#
# WHY. basal_area_85 is the top residual predictor in all six models, and in the
# broad RF (residuals_rf_everything.R) the 1985 stand size carries the signal
# while 1990/1995/2021-23 add nothing. Shown directly: one panel per inventory.
#
# WHAT. Plot-mean log residual (+ = model under-predicts) against tree-tally
# basal area of each inventory (ppa_kaikki_85/90/95/mustikka, one method for
# all four years, kasvillisuustk_kooste_13022026.xlsx). Points = mean of the
# Yasso07 and Yasso15 residuals (the reference target of the residual study);
# lines = loess of that mean and of each of the two models. Same axes in every panel.
#
# Run from repo root:  Rscript doublechecks/residuals_basal_area_by_year.R
# =============================================================================

suppressPackageStartupMessages({library(dplyr); library(readxl)})
source("manuscript/figures/run_ids.R")
cat("RUN_IDs:", paste(names(RID), RID, sep = "=", collapse = "  "), "\n\n")

TGT  <- c("Yasso07", "Yasso15")
COLS <- c(Yasso07 = "#eda100", Yasso15 = "#e87ba4")
YRS  <- c(`1985` = "ppa_kaikki_85", `1990` = "ppa_kaikki_90",
          `1995` = "ppa_kaikki_95", `2021-23` = "ppa_kaikki_mustikka")
OUT  <- "doublechecks/figures"; dir.create(OUT, showWarnings = FALSE)
INK  <- "#3d3d3a"; GRID <- "#e6e5df"; PT <- adjustcolor("#8a8983", 0.5)
XMAX <- 50

P <- bind_rows(lapply(TGT, function(m)
  as.data.frame(readRDS(sprintf(
    "Calibration_real_data_transient/runs/%s_posterior_predictive_%s.rds",
    m, RID[[m]]))$residuals_df) |>
    group_by(plot_id) |> summarise(r = mean(residual_log), .groups = "drop") |>
    mutate(model = m)))
w <- suppressWarnings(read_excel("Data/biomassatiedostot/kasvillisuustk_kooste_13022026.xlsx"))
ba <- function(ids, col) suppressWarnings(as.numeric(w[[col]][match(ids, as.numeric(w$koealatunnus_BIOSOIL))]))
ens <- P |> group_by(plot_id) |> summarise(r = mean(r), .groups = "drop")
for (y in names(YRS)) { ens[[y]] <- ba(ens$plot_id, YRS[[y]]); P[[y]] <- ba(P$plot_id, YRS[[y]]) }

cat("Spearman among inventories:\n")
print(round(cor(ens[names(YRS)], use = "pair", method = "spearman"), 2))

png(file.path(OUT, "residuals_basal_area_by_year.png"), 14, 4.4, units = "in", res = 150)
par(mfrow = c(1, 4), mar = c(4, 1, 2.4, 0.6), oma = c(0, 3.4, 2.2, 0),
    mgp = c(2.2, 0.6, 0), tcl = -0.3, col.axis = INK, col.lab = INK, fg = INK)
YL <- quantile(ens$r, c(.005, .995))
for (i in seq_along(YRS)) {
  y  <- names(YRS)[i]; x <- ens[[y]]; ok <- is.finite(x)
  rho <- cor(x[ok], ens$r[ok], method = "spearman")
  r2  <- summary(lm(ens$r[ok] ~ x[ok]))$r.squared
  plot(NA, xlim = c(0, XMAX), ylim = YL, bty = "l", yaxt = if (i == 1) "s" else "n",
       xlab = "Basal area (m2/ha)", ylab = "")
  abline(h = 0, col = GRID, lwd = 2)
  points(pmin(x, XMAX), ens$r, pch = 16, cex = 0.6, col = PT)
  for (m in TGT) {
    d <- P[P$model == m & is.finite(P[[y]]), ]
    d$x <- d[[y]]
    lo <- loess(r ~ x, d, span = 0.75)
    xs <- seq(quantile(d$x, .02), quantile(d$x, .98), length.out = 100)
    lines(xs, predict(lo, data.frame(x = xs)), col = COLS[[m]], lwd = 1.5)
  }
  lo <- loess(r ~ x, data.frame(x = x[ok], r = ens$r[ok]), span = 0.75)
  xs <- seq(quantile(x[ok], .02), quantile(x[ok], .98), length.out = 100)
  lines(xs, predict(lo, data.frame(x = xs)), col = INK, lwd = 2.5)
  mtext(y, 3, 0.9, adj = 0, font = 2, cex = 0.9)
  mtext(sprintf("ρ %.2f   R² %.2f   n %d", rho, r2, sum(ok)), 3, 0.1, adj = 0, cex = 0.7)
  if (i == 4) legend("topright", c("mean of Yasso07 + Yasso15", TGT), col = c(INK, COLS),
                     lwd = c(2.5, 1.5, 1.5), bty = "n", cex = 0.85)
  cat(sprintf("%-8s rho %.2f  R2 %.3f  n %d  (>%d m2/ha clipped: %d)\n",
              y, rho, r2, sum(ok), XMAX, sum(x > XMAX, na.rm = TRUE)))
}
mtext("Log residual (obs/pred), plot mean", 2, 2, outer = TRUE, cex = 0.75, col = INK)
mtext("Model residual against the basal area of each inventory; points = mean of Yasso07 and Yasso15 per plot, lines = loess",
      3, 0.6, outer = TRUE, cex = 0.8, col = INK)
dev.off()
cat("\nWrote", file.path(OUT, "residuals_basal_area_by_year.png"), "\n")
