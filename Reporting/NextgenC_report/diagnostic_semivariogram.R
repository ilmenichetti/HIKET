# === Diagnostic: empirical + fitted semivariograms (NOT part of Zenodo deposit) ===
# 2x3 panel, one per model: pooled empirical variogram of temporal-mean log-SOC
# with the fitted model overlaid. For inspection only -> diagnostic_semivariogram.png

suppressPackageStartupMessages({ library(gstat); library(sp) })

PROJ <- "/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling"
RUNS <- file.path(PROJ, "Calibration_real_data_transient", "runs")
OUT  <- file.path(PROJ, "Reporting", "NextgenC_report", "diagnostic_semivariogram.png")

site   <- read.csv(file.path(PROJ, "Data", "model_inputs", "site_raw.csv"))
coords <- site[, c("plot_id", "x_ETRS", "y_ETRS")]
coords <- coords[stats::complete.cases(coords), ]

bundles <- c(
  SP1     = "SP1_posterior_predictive_20260608_015026.rds",
  TP2     = "TP2_posterior_predictive_20260608_020212.rds",
  TP3     = "TP3_posterior_predictive_20260608_020212.rds",
  Yasso07 = "Yasso07_posterior_predictive_20260611_032825.rds",
  Yasso15 = "Yasso15_posterior_predictive_20260611_032825.rds",
  Yasso20 = "Yasso20_posterior_predictive_20260611_033140.rds")

png(OUT, width = 2400, height = 1500, res = 200)
par(mfrow = c(2, 3), mar = c(4.2, 4.4, 3, 1), mgp = c(2.4, 0.8, 0))

for (m in names(bundles)) {
  pp <- readRDS(file.path(RUNS, bundles[[m]]))
  s  <- merge(pp$posterior_summary[, c("plot_id", "year", "soc_mean")], coords, by = "plot_id")
  s  <- s[is.finite(s$soc_mean) & s$soc_mean > 0, ]
  s$zlog <- log(s$soc_mean)

  pooled <- aggregate(zlog ~ plot_id + x_ETRS + y_ETRS, data = s, FUN = mean)
  coordinates(pooled) <- ~ x_ETRS + y_ETRS
  v0    <- stats::var(pooled$zlog)
  vemp  <- variogram(zlog ~ 1, pooled, cutoff = 4e5)
  nug0  <- min(vemp$gamma); psill0 <- max(v0 - nug0, 0.05 * v0)
  vfit  <- suppressWarnings(fit.variogram(vemp,
            vgm(psill = psill0, model = "Exp", range = 8e4, nugget = nug0), fit.method = 7))
  if (isTRUE(attr(vfit, "singular")) || vfit$range[2] <= 0 || vfit$range[2] > 4e5)
    vfit <- vgm(psill = psill0, model = "Exp", range = 1e5, nugget = nug0)

  nug   <- vfit$psill[1]; sill <- sum(vfit$psill); rng <- vfit$range[2]
  nugfr <- nug / sill
  line  <- variogramLine(vfit, maxdist = max(vemp$dist))

  plot(vemp$dist/1000, vemp$gamma, pch = 19, col = "grey25",
       xlab = "separation distance (km)", ylab = "semivariance  γ(h)  [log-SOC]",
       main = m, ylim = c(0, max(vemp$gamma) * 1.15))
  lines(line$dist/1000, line$gamma, col = "firebrick", lwd = 2)
  abline(h = sill, lty = 2, col = "grey50")            # total sill (~field var)
  abline(h = nug,  lty = 3, col = "steelblue")          # nugget
  legend("bottomright", bty = "n", cex = 0.9,
         legend = c(sprintf("nugget = %.3f (%.0f%% of sill)", nug, 100*nugfr),
                    sprintf("sill = %.3f", sill),
                    sprintf("range = %.0f km", rng/1000)))
}
dev.off()
cat("wrote", OUT, "\n")
