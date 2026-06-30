# === NextGenC reporting: multi-panel SOC trajectory plot ===
# One panel per model; every plot's mean-SOC trajectory drawn as a
# transparent line (year on x, mean total SOC on y). Reads the *_SOC_mean.csv
# matrices produced by build_soc_matrices.R.

OUTDIR <- "/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling/Reporting/NextgenC_report"

models <- c("SP1", "TP2", "TP3", "Yasso07", "Yasso15", "Yasso20")

# load mean matrices (plot x year)
mats <- lapply(models, function(m) {
  d <- read.csv(file.path(OUTDIR, sprintf("%s_SOC_mean.csv", m)), check.names = FALSE)
  list(years = as.numeric(names(d)[-1]),
       mat   = as.matrix(d[, -1]))   # rows = plots, cols = years
})
names(mats) <- models

# common y-range across all models for comparability
ylim <- range(unlist(lapply(mats, function(x) range(x$mat))))

line_col <- rgb(0.10, 0.30, 0.65, alpha = 0.06)   # transparent blue

png(file.path(OUTDIR, "NextGenC_SOC_trajectories.png"),
    width = 2400, height = 1500, res = 200)
op <- par(mfrow = c(2, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))

for (m in models) {
  yrs <- mats[[m]]$years
  mat <- mats[[m]]$mat
  plot(NA, xlim = range(yrs), ylim = ylim,
       xlab = "Year", ylab = "Mean SOC (tC/ha)",
       main = sprintf("%s  (n = %d plots)", m, nrow(mat)))
  # one transparent line per plot
  for (i in seq_len(nrow(mat))) lines(yrs, mat[i, ], col = line_col, lwd = 1)
  # median trajectory across plots for reference
  lines(yrs, apply(mat, 2, median), col = "firebrick", lwd = 2)
}

mtext("Per-plot mean-SOC trajectories by model (red = across-plot median)",
      outer = TRUE, cex = 1.1, font = 2)
par(op)
dev.off()

cat("Wrote NextGenC_SOC_trajectories.png to", OUTDIR, "\n")
