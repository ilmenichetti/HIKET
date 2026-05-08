# =============================================================================
# diagnostics_inputs.R
#
# Sanity checks on litter inputs and SOC observations BEFORE MCMC.
# Reads input_raw_monthly.csv and site_raw.csv produced by Data_work.R.
#
# Run from the project root. Outputs go to ./diagnostics/inputs/.
#
# CHECKS:
#   1. MRT distribution: mean(SOC) / mean(annual total litter input) per plot.
#      Boreal expectation: 30-80 years for mineral SOC. Outliers worth checking.
#
#   2. Magnitude: observed SOC vs implied steady-state SOC = I * MRT_expected.
#      MRT_expected = 50 years as a sanity-check midpoint.
#      Ratio (obs / implied) should cluster around 1 if inputs are roughly right.
#
#   3. Oscillation analysis (3 time scales):
#        - dI/dt year-to-year:    I[t+1] - I[t]
#        - 2-year rolling slope:  linear slope across pairs of years
#        - 4-year rolling slope:  linear slope across 4-year windows
#      Histograms of each. CVs and max-absolute summaries per plot.
#      Big oscillations suggest data interpolation artefacts or NFI re-measurement
#      jumps that the smooth-input assumption of Yasso07 can't accommodate.
#
# OUTPUT FILES (./diagnostics/inputs/):
#   01_mrt_histogram.png
#   02_magnitude_check.png
#   03_oscillation_diff.png
#   04_oscillation_slope2.png
#   05_oscillation_slope4.png
#   06_per_plot_oscillation_summary.png
#   summary_stats.txt    (a short text report of headline numbers)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
})

DIR_OUT <- "./Calibration_real_data/diagnostics/inputs"
dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)

PX_PER_IN <- 300L

# ---- Load -------------------------------------------------------------------
input_raw <- read.csv("./Data/model_inputs/input_raw_monthly.csv")
site_raw  <- read.csv("./Data/model_inputs/site_raw.csv")

# Use calib-ready plots throughout.
calib_plots <- as.character(site_raw$plot_id[site_raw$calib_ready])
input_calib <- input_raw[as.character(input_raw$plot_id) %in% calib_plots, ]

cat("=============================================================\n")
cat("  Input diagnostics\n")
cat(sprintf("  %d calib-ready plots, %d input rows\n",
            length(calib_plots), nrow(input_calib)))
cat("=============================================================\n\n")


# ---- Helper: annual total litter (across NWL/FWL/CWL x AWEN) ----------------
litter_cols <- grep("^C_(nwl|fwl|cwl)_[AWEN]$", names(input_calib), value = TRUE)
stopifnot(length(litter_cols) == 12)

annual_lit <- input_calib %>%
  group_by(plot_id, year) %>%
  summarise(total_litter_yr = sum(across(all_of(litter_cols)), na.rm = TRUE),
            .groups = "drop")


# =============================================================================
# 1.  MRT distribution
# =============================================================================
# MRT = mean(SOC obs across years) / mean(annual total litter input across years)
# Expressed in years. With 1-2 SOC obs per plot this is a coarse estimate, but
# the distribution shape is informative.

soc_obs <- input_calib %>%
  filter(!is.na(soc_obs_tCha)) %>%
  group_by(plot_id) %>%
  summarise(mean_soc = mean(soc_obs_tCha, na.rm = TRUE),
            n_obs    = n(),
            .groups  = "drop")

mean_lit_per_plot <- annual_lit %>%
  group_by(plot_id) %>%
  summarise(mean_input = mean(total_litter_yr, na.rm = TRUE),
            .groups    = "drop")

mrt_df <- merge(soc_obs, mean_lit_per_plot, by = "plot_id")
mrt_df$MRT <- mrt_df$mean_soc / mrt_df$mean_input

# Drop plots with mean_input == 0 (would give Inf) and report them
n_zero_input <- sum(mrt_df$mean_input <= 0, na.rm = TRUE)
mrt_df_finite <- mrt_df[is.finite(mrt_df$MRT), ]
n_finite      <- nrow(mrt_df_finite)

cat(sprintf("[1] MRT: %d plots with finite MRT (%d dropped for zero input)\n",
            n_finite, n_zero_input))
cat(sprintf("    Median:  %.1f years\n", median(mrt_df_finite$MRT)))
cat(sprintf("    Mean:    %.1f years\n", mean(mrt_df_finite$MRT)))
cat(sprintf("    Range:   %.1f -- %.1f years\n",
            min(mrt_df_finite$MRT), max(mrt_df_finite$MRT)))
cat(sprintf("    >100y:   %d plots (%.1f%%)\n",
            sum(mrt_df_finite$MRT > 100),
            100 * sum(mrt_df_finite$MRT > 100) / n_finite))
cat(sprintf("    <20y:    %d plots (%.1f%%)\n",
            sum(mrt_df_finite$MRT < 20),
            100 * sum(mrt_df_finite$MRT < 20) / n_finite))
cat("\n")



png(file.path(DIR_OUT, "00_mean_inputs_histogram.png"),
    width = 9L * PX_PER_IN, height = 6L * PX_PER_IN, res = PX_PER_IN)
br <- pretty(c(0, quantile(mean_lit_per_plot$mean_input, 0.99, na.rm = TRUE)), n = 40)
hist(pmin(mean_lit_per_plot$mean_input, max(br)), breaks = br,
     col = "steelblue", border = "white",
     xlab = "Mean annual litter input",
     ylab = "Number of plots",
     main = sprintf("Plot-level mean litter input distribution  (n = %d)", n_finite))
dev.off()



png(file.path(DIR_OUT, "01_mrt_histogram.png"),
    width = 9L * PX_PER_IN, height = 6L * PX_PER_IN, res = PX_PER_IN)
par(mar = c(4, 4, 3, 1))
br <- pretty(c(0, quantile(mrt_df_finite$MRT, 0.99, na.rm = TRUE)), n = 40)
hist(pmin(mrt_df_finite$MRT, max(br)), breaks = br,
     col = "steelblue", border = "white",
     xlab = "MRT (years) = mean SOC / mean annual litter input",
     ylab = "Number of plots",
     main = sprintf("Plot-level MRT distribution  (n = %d)", n_finite))
abline(v = median(mrt_df_finite$MRT), col = "tomato", lwd = 2, lty = 2)
abline(v = c(30, 80), col = "darkgreen", lwd = 1, lty = 3)
legend("topright",
       legend = c(sprintf("Median: %.1f y", median(mrt_df_finite$MRT)),
                  "Boreal expected: 30-80 y"),
       col = c("tomato", "darkgreen"), lty = c(2, 3), lwd = c(2, 1),
       bty = "n")
mtext("Note: x-axis truncated at 99th percentile; right-tail outliers compressed",
      side = 1, line = 2.8, cex = 0.7, col = "grey40")
dev.off()


# =============================================================================
# 2.  Magnitude check
# =============================================================================
# Implied steady-state SOC under MRT_expected = 50 y:  SOC* ~ I * 50
# Compare to observed mean SOC. Ratio obs/implied:
#   ~1     -- inputs roughly consistent with observed SOC at typical boreal MRT
#   <0.5   -- inputs too high (or MRT shorter than expected)
#   >2     -- inputs too low (or MRT longer than expected)

MRT_EXPECTED <- 50L

mag <- mrt_df_finite
mag$soc_implied   <- mag$mean_input * MRT_EXPECTED
mag$ratio_obs_imp <- mag$mean_soc / mag$soc_implied

cat(sprintf("[2] Magnitude (assuming MRT = %d y):\n", MRT_EXPECTED))
cat(sprintf("    obs/implied ratio: median %.2f, mean %.2f\n",
            median(mag$ratio_obs_imp, na.rm = TRUE),
            mean(mag$ratio_obs_imp, na.rm = TRUE)))
cat(sprintf("    %d plots with ratio < 0.5  (inputs likely too high)\n",
            sum(mag$ratio_obs_imp < 0.5, na.rm = TRUE)))
cat(sprintf("    %d plots with ratio > 2    (inputs likely too low)\n",
            sum(mag$ratio_obs_imp > 2,   na.rm = TRUE)))
cat("\n")

png(file.path(DIR_OUT, "02_magnitude_check.png"),
    width = 12L * PX_PER_IN, height = 6L * PX_PER_IN, res = PX_PER_IN)
par(mfrow = c(1, 2), mar = c(4, 4.5, 3, 1))

# Left: scatter of observed vs implied SOC
ax_lim <- range(c(mag$mean_soc, mag$soc_implied), na.rm = TRUE)
plot(mag$soc_implied, mag$mean_soc,
     xlim = ax_lim, ylim = ax_lim,
     pch = 19, col = adjustcolor("steelblue", 0.5),
     xlab = sprintf("Implied SOC = mean_input * %d (tC/ha)", MRT_EXPECTED),
     ylab = "Observed mean SOC (tC/ha)",
     main = sprintf("Magnitude check: obs vs implied SOC at MRT = %dy", MRT_EXPECTED))
abline(0, 1, lty = 2, col = "grey40", lwd = 1.5)

# Right: histogram of obs/implied ratio
br <- pretty(log10(pmax(mag$ratio_obs_imp, 0.01)), n = 30)
hist(log10(mag$ratio_obs_imp), breaks = br,
     col = "steelblue", border = "white",
     xlab = "log10(obs / implied SOC)",
     ylab = "Number of plots",
     main = "Distribution of obs/implied ratio (log10)")
abline(v = 0, lty = 2, col = "tomato", lwd = 2)
abline(v = log10(c(0.5, 2)), lty = 3, col = "darkgreen")
legend("topright",
       legend = c("ratio = 1 (consistent)",
                  "0.5 - 2 (acceptable band)"),
       col = c("tomato", "darkgreen"), lty = c(2, 3),
       bty = "n", cex = 0.85)
dev.off()


# =============================================================================
# 3.  Oscillation analysis (three time scales)
# =============================================================================
# Three measures of input variability per plot:
#   diff:    I[t+1] - I[t]                       year-to-year
#   slope2:  lm(I ~ t) over a 2-year rolling window
#   slope4:  lm(I ~ t) over a 4-year rolling window
#
# Each measure has units of tC/ha/yr per year. The slope quantities are
# numerically equivalent to the diff for window=2; for window=4 they
# average over more time and reveal slower trends.

# Helper: rolling slope using least-squares formula (no list-apply needed)
rolling_slope <- function(y, w) {
  n <- length(y)
  if (n < w) return(rep(NA_real_, n))
  x_centered <- seq(0, w - 1) - (w - 1) / 2
  denom      <- sum(x_centered^2)
  out        <- rep(NA_real_, n)
  for (i in seq_len(n - w + 1)) {
    yi    <- y[i:(i + w - 1)]
    yc    <- yi - mean(yi)
    out[i + w - 1] <- sum(x_centered * yc) / denom   # slope per yr
  }
  out
}

osc <- annual_lit %>%
  arrange(plot_id, year) %>%
  group_by(plot_id) %>%
  mutate(
    diff_yr  = c(NA_real_, diff(total_litter_yr)),
    slope_2  = rolling_slope(total_litter_yr, 2L),
    slope_4  = rolling_slope(total_litter_yr, 4L)
  ) %>%
  ungroup()

osc_summary <- osc %>%
  group_by(plot_id) %>%
  summarise(
    mean_input  = mean(total_litter_yr, na.rm = TRUE),
    sd_input    = sd(total_litter_yr,   na.rm = TRUE),
    cv_input    = sd_input / mean_input,
    max_abs_diff   = max(abs(diff_yr),  na.rm = TRUE),
    max_abs_slope2 = max(abs(slope_2),  na.rm = TRUE),
    max_abs_slope4 = max(abs(slope_4),  na.rm = TRUE),
    .groups     = "drop"
  )

cat("[3] Oscillation summary across plots:\n")
cat(sprintf("    Per-plot CV of annual inputs: median %.3f, max %.3f\n",
            median(osc_summary$cv_input, na.rm = TRUE),
            max(osc_summary$cv_input,    na.rm = TRUE)))
cat(sprintf("    Year-to-year |diff|:    median %.3f, max %.3f tC/ha\n",
            median(abs(osc$diff_yr),  na.rm = TRUE),
            max(abs(osc$diff_yr),     na.rm = TRUE)))
cat(sprintf("    2-year |slope|:         median %.3f, max %.3f tC/ha/yr\n",
            median(abs(osc$slope_2),  na.rm = TRUE),
            max(abs(osc$slope_2),     na.rm = TRUE)))
cat(sprintf("    4-year |slope|:         median %.3f, max %.3f tC/ha/yr\n",
            median(abs(osc$slope_4),  na.rm = TRUE),
            max(abs(osc$slope_4),     na.rm = TRUE)))
cat("\n")

# Helper: histogram with median + percentile annotations
plot_osc_hist <- function(values, title, xlab, file) {
  v <- values[is.finite(values)]
  if (length(v) == 0) return(invisible(NULL))
  q <- quantile(abs(v), 0.99, na.rm = TRUE)
  v_clip <- pmax(pmin(v, q), -q)   # clip extreme tails
  png(file, width = 9L * PX_PER_IN, height = 6L * PX_PER_IN, res = PX_PER_IN)
  par(mar = c(4, 4, 3, 1))
  hist(v_clip, breaks = 60, col = "steelblue", border = "white",
       xlab = xlab, ylab = "Frequency",
       main = sprintf("%s  (n = %d)", title, length(v)))
  abline(v = 0, lty = 2, col = "grey40")
  med <- median(v, na.rm = TRUE)
  abline(v = med, lty = 1, col = "tomato", lwd = 1.5)
  legend("topright",
         legend = c(sprintf("Median: %+.3f", med),
                    sprintf("99%%ile |v|: %.3f", q),
                    "x-axis clipped at 99%ile"),
         col = c("tomato", NA, NA), lty = c(1, NA, NA),
         lwd = c(1.5, NA, NA), bty = "n", cex = 0.8)
  dev.off()
}

plot_osc_hist(osc$diff_yr,
              "Year-to-year change in annual litter input",
              "I[t+1] - I[t]  (tC/ha/yr)",
              file.path(DIR_OUT, "03_oscillation_diff.png"))

plot_osc_hist(osc$slope_2,
              "2-year rolling slope of annual litter input",
              "Slope (tC/ha/yr per year), window = 2",
              file.path(DIR_OUT, "04_oscillation_slope2.png"))

plot_osc_hist(osc$slope_4,
              "4-year rolling slope of annual litter input",
              "Slope (tC/ha/yr per year), window = 4",
              file.path(DIR_OUT, "05_oscillation_slope4.png"))


# Per-plot oscillation summary: 3-panel histogram of max-abs-* per plot
png(file.path(DIR_OUT, "06_per_plot_oscillation_summary.png"),
    width = 14L * PX_PER_IN, height = 5L * PX_PER_IN, res = PX_PER_IN)
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))

plot_per_plot <- function(x, title, xlab) {
  v <- x[is.finite(x)]
  q <- quantile(v, 0.99)
  hist(pmin(v, q), breaks = 40, col = "steelblue", border = "white",
       xlab = xlab, ylab = "Number of plots",
       main = title)
  abline(v = median(v), col = "tomato", lwd = 2, lty = 2)
  legend("topright",
         legend = sprintf("Median: %.3f", median(v)),
         col = "tomato", lty = 2, lwd = 2, bty = "n", cex = 0.85)
}

plot_per_plot(osc_summary$max_abs_diff,
              "Max |yr-to-yr diff| per plot",
              "tC/ha/yr")
plot_per_plot(osc_summary$max_abs_slope2,
              "Max |2-yr slope| per plot",
              "tC/ha/yr per year")
plot_per_plot(osc_summary$max_abs_slope4,
              "Max |4-yr slope| per plot",
              "tC/ha/yr per year")

mtext("Per-plot maximum-absolute oscillation across the time series",
      side = 3, outer = TRUE, line = -1.5, cex = 1.0, font = 2)
dev.off()


# =============================================================================
# Text summary
# =============================================================================

sink(file.path(DIR_OUT, "summary_stats.txt"))
cat("=============================================================\n")
cat(sprintf("  Input diagnostics summary -- %s\n", format(Sys.time())))
cat(sprintf("  %d calib-ready plots\n", length(calib_plots)))
cat("=============================================================\n\n")

cat("[1] MRT distribution\n")
cat(sprintf("    Plots with finite MRT:  %d  (%d dropped for zero input)\n",
            n_finite, n_zero_input))
cat(sprintf("    Median MRT:             %.1f years\n", median(mrt_df_finite$MRT)))
cat(sprintf("    Mean MRT:               %.1f years\n", mean(mrt_df_finite$MRT)))
cat(sprintf("    Range:                  %.1f -- %.1f years\n",
            min(mrt_df_finite$MRT), max(mrt_df_finite$MRT)))
cat(sprintf("    Boreal expected:        ~30-80 years\n"))
cat(sprintf("    Plots with MRT > 100y:  %d (%.1f%%)\n",
            sum(mrt_df_finite$MRT > 100),
            100 * sum(mrt_df_finite$MRT > 100) / n_finite))
cat(sprintf("    Plots with MRT < 20y:   %d (%.1f%%)\n",
            sum(mrt_df_finite$MRT < 20),
            100 * sum(mrt_df_finite$MRT < 20) / n_finite))
cat("\n")

cat(sprintf("[2] Magnitude (assuming MRT_expected = %dy)\n", MRT_EXPECTED))
cat(sprintf("    obs/implied ratio: median %.2f, mean %.2f\n",
            median(mag$ratio_obs_imp, na.rm = TRUE),
            mean(mag$ratio_obs_imp,    na.rm = TRUE)))
cat(sprintf("    Plots with ratio < 0.5: %d  (inputs likely too high)\n",
            sum(mag$ratio_obs_imp < 0.5, na.rm = TRUE)))
cat(sprintf("    Plots with ratio > 2:   %d  (inputs likely too low)\n",
            sum(mag$ratio_obs_imp > 2,   na.rm = TRUE)))
cat("\n")

cat("[3] Oscillation (across all plot-years)\n")
cat(sprintf("    Year-to-year |diff|:  median %.3f, max %.3f tC/ha\n",
            median(abs(osc$diff_yr),  na.rm = TRUE),
            max(abs(osc$diff_yr),     na.rm = TRUE)))
cat(sprintf("    2-year |slope|:       median %.3f, max %.3f tC/ha/yr\n",
            median(abs(osc$slope_2),  na.rm = TRUE),
            max(abs(osc$slope_2),     na.rm = TRUE)))
cat(sprintf("    4-year |slope|:       median %.3f, max %.3f tC/ha/yr\n",
            median(abs(osc$slope_4),  na.rm = TRUE),
            max(abs(osc$slope_4),     na.rm = TRUE)))
cat("\n")

cat("    Per-plot CV of annual inputs:\n")
cat(sprintf("        Median:  %.3f\n", median(osc_summary$cv_input)))
cat(sprintf("        Max:     %.3f\n", max(osc_summary$cv_input)))
cat(sprintf("    Per-plot max |yr-to-yr diff|:\n"))
cat(sprintf("        Median:  %.3f tC/ha\n", median(osc_summary$max_abs_diff)))
cat(sprintf("        Max:     %.3f tC/ha\n", max(osc_summary$max_abs_diff)))

cat("\n")
cat(sprintf("    Output PNGs written to %s\n", DIR_OUT))
sink()


cat("=============================================================\n")
cat(sprintf("  Done. PNGs and summary in %s\n", DIR_OUT))
cat("=============================================================\n")






# =============================================================================
# 7.  Mean litter input trajectories by region (North vs South)
# =============================================================================

library(ggplot2)

# North/South split via kasvyo_syke: zones 1-2 = South, 3-5 = North.
# Fall back to latitude median split if kasvyo_syke is missing.
region_lookup <- site_raw[, c("plot_id",
                              if ("kasvyo_syke" %in% names(site_raw)) "kasvyo_syke" else NULL,
                              if ("lat_WGS84"   %in% names(site_raw)) "lat_WGS84"   else NULL)]
region_lookup$plot_id <- as.character(region_lookup$plot_id)

if ("kasvyo_syke" %in% names(region_lookup)) {
  region_lookup$region <- ifelse(region_lookup$kasvyo_syke %in% c(1, 2),
                                 "South Finland", "North Finland")
} else {
  lat_median <- median(region_lookup$lat_WGS84, na.rm = TRUE)
  region_lookup$region <- ifelse(region_lookup$lat_WGS84 >= lat_median,
                                 "North Finland", "South Finland")
}

annual_lit_region <- annual_lit %>%
  mutate(plot_id = as.character(plot_id)) %>%
  left_join(region_lookup[, c("plot_id", "region")], by = "plot_id") %>%
  filter(!is.na(region))

# Per-year quantile bands across plots within each region
traj_region <- annual_lit_region %>%
  group_by(region, year) %>%
  summarise(
    mean  = mean(total_litter_yr, na.rm = TRUE),
    q25   = quantile(total_litter_yr, 0.25, na.rm = TRUE),
    q75   = quantile(total_litter_yr, 0.75, na.rm = TRUE),
    q025  = quantile(total_litter_yr, 0.025, na.rm = TRUE),
    q975  = quantile(total_litter_yr, 0.975, na.rm = TRUE),
    n     = sum(!is.na(total_litter_yr)),
    .groups = "drop"
  )

region_colours <- c("South Finland" = "#2b4f9e",
                    "North Finland" = "#9e2b2b")

p_region <- ggplot(traj_region, aes(x = year, colour = region, fill = region)) +
  geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.18, colour = NA) +
  geom_ribbon(aes(ymin = q25,  ymax = q75),  alpha = 0.35, colour = NA) +
  geom_line(aes(y = mean), linewidth = 1.0) +
  scale_colour_manual(values = region_colours, name = "Region") +
  scale_fill_manual(  values = region_colours, name = "Region") +
  labs(
    title    = "Annual litter inputs by region: cross-plot uncertainty",
    subtitle = "Solid line = mean across plots  |  Dark band = 50% interval  |  Light band = 95% interval",
    x        = "Year",
    y        = "Total annual litter input (tC/ha/yr)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    plot.subtitle    = element_text(size = 10, colour = "grey40")
  )

ggsave(file.path(DIR_OUT, "07_litter_trajectories_by_region.png"),
       p_region,
       width = 10, height = 6, dpi = PX_PER_IN)

cat(sprintf("\nRegion trajectory plot: %s\n",
            file.path(DIR_OUT, "07_litter_trajectories_by_region.png")))



# =============================================================================
# 8.  Average ACF across all plots
# =============================================================================
# Computes the ACF of annual litter input for each plot separately, then
# averages the autocorrelation coefficients across plots at each lag.
# The piecewise-linear interpolation between NFI anchor years (~5y spacing)
# should produce strong positive autocorrelation decaying over ~5 lags.
# A spike in the mean ACF at lag 5 would confirm the NFI cycle is the
# dominant source of temporal structure in the inputs.

max_lag <- 10L

plot_ids_acf <- unique(annual_lit$plot_id)

acf_list <- lapply(plot_ids_acf, function(pid) {
  y <- annual_lit$total_litter_yr[annual_lit$plot_id == pid]
  y <- y[order(annual_lit$year[annual_lit$plot_id == pid])]
  if (length(y) < max_lag + 2L || sd(y, na.rm = TRUE) < 1e-10)
    return(NULL)   # flat series: ACF undefined, skip
  a <- acf(y, lag.max = max_lag, plot = FALSE, na.action = na.pass)
  a$acf[, 1, 1]   # vector of length max_lag + 1 (lag 0 ... max_lag)
})

acf_mat   <- do.call(rbind, acf_list[!sapply(acf_list, is.null)])
n_acf     <- nrow(acf_mat)
mean_acf  <- colMeans(acf_mat, na.rm = TRUE)
se_acf    <- apply(acf_mat, 2, sd, na.rm = TRUE) / sqrt(n_acf)
lags      <- 0:max_lag

# Approximate 95% pointwise CI on the mean ACF (not the individual-plot CI)
ci_lo <- mean_acf - 1.96 * se_acf
ci_hi <- mean_acf + 1.96 * se_acf

# Significance threshold for a single white-noise series (for reference)
wn_threshold <- qnorm(0.975) / sqrt(mean(table(annual_lit$plot_id)))

cat(sprintf("[8] Average ACF: %d plots included (%d skipped — flat or too short)\n",
            n_acf, length(plot_ids_acf) - n_acf))
cat(sprintf("    Lag-1 mean ACF:  %.3f\n", mean_acf[2]))
cat(sprintf("    Lag-5 mean ACF:  %.3f\n", mean_acf[6]))
cat(sprintf("    Lag-10 mean ACF: %.3f\n", mean_acf[11]))
cat("\n")

png(file.path(DIR_OUT, "08_mean_acf.png"),
    width = 9L * PX_PER_IN, height = 6L * PX_PER_IN, res = PX_PER_IN)
par(mar = c(4, 4.5, 3, 1))

plot(lags, mean_acf,
     type = "n",
     xlim = c(0, max_lag), ylim = range(c(ci_lo, ci_hi, -0.15, 1.05)),
     xlab = "Lag (years)",
     ylab = "Mean autocorrelation across plots",
     main = sprintf("Mean ACF of annual litter inputs (n = %d plots)", n_acf),
     xaxt = "n")
axis(1, at = 0:max_lag)

# Shaded CI band on the mean
polygon(c(lags, rev(lags)), c(ci_hi, rev(ci_lo)),
        col = adjustcolor("steelblue", alpha.f = 0.2), border = NA)

# White-noise reference lines
abline(h =  wn_threshold, lty = 3, col = "grey50")
abline(h = -wn_threshold, lty = 3, col = "grey50")
abline(h = 0,             lty = 1, col = "grey80")

# Stems + points
segments(lags, 0, lags, mean_acf, col = "steelblue", lwd = 2)
points(lags, mean_acf, pch = 19, col = "steelblue", cex = 0.9)

# Mark lag 5 explicitly
points(5, mean_acf[6], pch = 19, col = "tomato", cex = 1.3)
text(5, mean_acf[6], labels = sprintf("lag 5\n%.3f", mean_acf[6]),
     pos = 4, cex = 0.8, col = "tomato")

legend("topright",
       legend = c("Mean ACF", "±1.96 SE on mean", "White-noise threshold"),
       col    = c("steelblue", adjustcolor("steelblue", 0.3), "grey50"),
       lty    = c(1, NA, 3), lwd = c(2, NA, 1),
       pch    = c(19, 15, NA), pt.cex = c(0.9, 2, NA),
       bty    = "n", cex = 0.85)
dev.off()


# =============================================================================
# 9.  Average periodogram across all plots
# =============================================================================
# Computes the raw (unnormalised) periodogram for each plot via spec.pgram(),
# then averages power across plots at matching frequencies.
# Because plots have different lengths, we interpolate each periodogram onto a
# common frequency grid before averaging.
# Expected: a dominant peak at frequency 1/5 yr^-1 (period = 5 years) from
# the NFI measurement cycle. Secondary peaks at harmonics (1/10, 1/2.5) are
# also possible.

# Common frequency grid: 0.05 to 0.5 in steps of 0.01 (Nyquist at 0.5 yr^-1)
freq_grid <- seq(0.05, 0.50, by = 0.01)

pgram_list <- lapply(plot_ids_acf, function(pid) {
  y <- annual_lit$total_litter_yr[annual_lit$plot_id == pid]
  y <- y[order(annual_lit$year[annual_lit$plot_id == pid])]
  y <- y[!is.na(y)]
  if (length(y) < 10L || sd(y) < 1e-10) return(NULL)
  
  # Detrend and taper before FFT to reduce spectral leakage
  sp <- spec.pgram(y, detrend = TRUE, taper = 0.1, plot = FALSE)
  
  # Interpolate onto common frequency grid (log-power for stability)
  log_power <- approx(sp$freq, log(sp$spec), xout = freq_grid,
                      rule = 2)$y
  log_power
})

pgram_mat   <- do.call(rbind, pgram_list[!sapply(pgram_list, is.null)])
n_pgram     <- nrow(pgram_mat)
mean_lpower <- colMeans(pgram_mat, na.rm = TRUE)
se_lpower   <- apply(pgram_mat, 2, sd, na.rm = TRUE) / sqrt(n_pgram)

cat(sprintf("[9] Average periodogram: %d plots included\n", n_pgram))
peak_idx    <- which.max(mean_lpower)
peak_freq   <- freq_grid[peak_idx]
cat(sprintf("    Dominant frequency: %.3f yr^-1  (period = %.1f years)\n",
            peak_freq, 1 / peak_freq))
cat("\n")

png(file.path(DIR_OUT, "09_mean_periodogram.png"),
    width = 10L * PX_PER_IN, height = 6L * PX_PER_IN, res = PX_PER_IN)
par(mar = c(4.5, 4.5, 3, 1))

y_lo <- mean_lpower - 1.96 * se_lpower
y_hi <- mean_lpower + 1.96 * se_lpower

plot(freq_grid, mean_lpower,
     type = "n",
     xlim = c(0.05, 0.50),
     ylim = range(c(y_lo, y_hi)),
     xlab = "Frequency (yr⁻¹)",
     ylab = "Mean log spectral power",
     main = sprintf("Mean periodogram of annual litter inputs (n = %d plots)", n_pgram))

# Period labels on top axis
period_ticks <- c(2, 2.5, 3, 4, 5, 7, 10, 20)
freq_ticks   <- 1 / period_ticks
axis(3, at = freq_ticks, labels = paste0(period_ticks, "y"), cex.axis = 0.8)
mtext("Period", side = 3, line = 2.2, cex = 0.85)

# CI band
polygon(c(freq_grid, rev(freq_grid)), c(y_hi, rev(y_lo)),
        col = adjustcolor("steelblue", alpha.f = 0.2), border = NA)

lines(freq_grid, mean_lpower, col = "steelblue", lwd = 1.8)

# Mark NFI-cycle frequency (1/5) and harmonics
nfi_freqs   <- c(1/5, 1/10, 2/5)
nfi_labels  <- c("NFI cycle\n(5y)", "2nd harmonic\n(10y)", "1st harmonic\n(2.5y)")
for (j in seq_along(nfi_freqs)) {
  abline(v = nfi_freqs[j], lty = 2, col = "tomato", lwd = 1.2)
  text(nfi_freqs[j], max(y_hi) * 0.97 - (j - 1) * diff(range(y_hi)) * 0.06,
       labels = nfi_labels[j], col = "tomato", cex = 0.72, adj = c(-0.1, 1))
}

legend("topright",
       legend = c("Mean log power", "±1.96 SE on mean", "NFI cycle & harmonics"),
       col    = c("steelblue", adjustcolor("steelblue", 0.3), "tomato"),
       lty    = c(1, NA, 2), lwd = c(1.8, NA, 1.2),
       pch    = c(NA, 15, NA), pt.cex = c(NA, 2, NA),
       bty    = "n", cex = 0.85)
dev.off()

cat("=============================================================\n")
cat(sprintf("  Sections 8-9 done. PNGs written to %s\n", DIR_OUT))
cat("=============================================================\n")


# =============================================================================
# 10.  Climate space by site fertility type (kasvup_tyyppi)
# =============================================================================
# Scatter plot of plot-level mean annual temperature vs mean annual
# precipitation for all calibration-ready plots, coloured and shaped by
# kasvup_tyyppi — the full 8-class Finnish NFI site fertility classification
# recorded at each plot in the 1985 kuvio data.
#
# kasvup_tyyppi coding (from Data_work.R):
#   1 = Lehdot (OMaT)           -- most productive, herb-rich
#   2 = Lehtomainen kangas (OMT)
#   3 = Tuorekangas (MT)        -- mesic heath, most common
#   4 = Kuivahko kangas (VT)    -- sub-xeric heath
#   5 = Kuiva kangas (CT)       -- xeric heath
#   6 = Karukkokangas (ClT)     -- very xeric, lichen-dominated
#   7 = Kalliot ja hietikot     -- rocky / sandy soils
#   8 = Lakimetsät ja tunturit  -- fell / subalpine
#
# This is finer than KA (which collapses to 4 mineral classes) and more
# meaningful for characterising climate-productivity co-variation, because
# it preserves the distinction between the xeric end-members (5-8) that KA
# lumps together.
#
# PURPOSE: check whether the 8 fertility classes occupy distinct climate
# niches across the calibration set, and whether the set covers Finland's
# full temperature-precipitation gradient without strong class-specific gaps.
#
# Data: mean_temp and mean_precip_annual from site_raw (gridded climate
# normals pre-computed in Data_work.R; same variables used in the residual
# random-forest). kasvup_tyyppi from the 1985 NFI kuvio merge in site_raw.
#
# WHY site_raw RATHER THAN input_raw:
#   The site_raw climate variables are long-run normals covering the full
#   model run period — they characterise "which climate space does this plot
#   occupy". Computing a temporal mean from the monthly temp_air / precip
#   columns in input_raw would give the same numbers (same gridded source).
# =============================================================================

# Identify which of the needed columns are actually present. The guard below
# lets the script run and emit a warning if site_raw is an older version
# missing one of the climate columns.
climate_site_cols <- c("plot_id", "mean_temp", "mean_precip_annual",
                       "kasvup_tyyppi")
available_cols    <- intersect(climate_site_cols, names(site_raw))
missing_cols      <- setdiff(climate_site_cols, names(site_raw))

if (length(missing_cols) > 0)
  warning(sprintf("[10] Column(s) not found in site_raw, output may be incomplete: %s",
                  paste(missing_cols, collapse = ", ")))

if (!all(c("mean_temp", "mean_precip_annual") %in% names(site_raw))) {
  message("[10] Skipping climate-space plot: mean_temp or mean_precip_annual absent from site_raw")
} else {
  
  clim_site <- site_raw[as.character(site_raw$plot_id) %in% calib_plots,
                        available_cols, drop = FALSE]
  
  # Human-readable labels for the 8 NFI site fertility codes.
  # Codes are integers 1-8; any value outside this range (or NA) gets "Unknown".
  ktp_labels <- c(
    "1" = "1 Lehdot (OMaT)",
    "2" = "2 Lehtomainen (OMT)",
    "3" = "3 Tuorekangas (MT)",
    "4" = "4 Kuivahko (VT)",
    "5" = "5 Kuiva (CT)",
    "6" = "6 Karukkokangas (ClT)",
    "7" = "7 Kalliot/hietikot",
    "8" = "8 Lakimetsät/tunturit"
  )
  
  if ("kasvup_tyyppi" %in% names(clim_site)) {
    n_ktp_na <- sum(is.na(clim_site$kasvup_tyyppi))
    if (n_ktp_na > 0)
      message(sprintf("[10] %d plots with NA kasvup_tyyppi shown as 'Unknown'",
                      n_ktp_na))
    # Map codes to labels; unmapped codes fall through to NA -> "Unknown"
    raw_codes   <- as.character(clim_site$kasvup_tyyppi)
    label_vec   <- ktp_labels[raw_codes]
    label_vec[is.na(label_vec)] <- "Unknown"
    clim_site$ktp_factor <- factor(label_vec,
                                   levels = c(ktp_labels, "Unknown"))
    clim_site$ktp_factor <- droplevels(clim_site$ktp_factor)
  } else {
    clim_site$ktp_factor <- factor(rep("Unknown", nrow(clim_site)))
  }
  
  ktp_levels <- levels(clim_site$ktp_factor)
  n_ktp      <- length(ktp_levels)
  
  # Divergent palette: MetBrewer::Austria via paletteer (8 colours, one per
  # named fertility class in order). Unknown gets grey regardless.
  # paletteer returns a character vector of hex codes; naming them to match
  # the ordered factor levels makes the mapping explicit and robust to absent
  # levels in the calibration subset.
  named_levels <- c("1 Lehdot (OMaT)", "2 Lehtomainen (OMT)",
                    "3 Tuorekangas (MT)", "4 Kuivahko (VT)",
                    "5 Kuiva (CT)", "6 Karukkokangas (ClT)",
                    "7 Kalliot/hietikot", "8 Lakimetsät/tunturit")
  austria_8   <- as.character(paletteer::paletteer_d("colorblindr::OkabeIto_black",
                                                     n = 8L))
  all_colours <- setNames(c(austria_8, "#999999"),
                          c(named_levels, "Unknown"))
  
  # pch: solid/open variants of circle/triangle/square/diamond — distinguishable
  # in greyscale and for red-green colourblind readers.
  all_pch <- c(
    "1 Lehdot (OMaT)"        = 19L,
    "2 Lehtomainen (OMT)"    = 17L,
    "3 Tuorekangas (MT)"     = 15L,
    "4 Kuivahko (VT)"        = 18L,
    "5 Kuiva (CT)"           = 1L,
    "6 Karukkokangas (ClT)"  = 2L,
    "7 Kalliot/hietikot"     = 0L,
    "8 Lakimetsät/tunturit"  = 5L,
    "Unknown"                = 4L
  )
  
  ktp_colours <- all_colours[ktp_levels]
  ktp_pch     <- all_pch[ktp_levels]
  
  col_vec <- ktp_colours[as.integer(clim_site$ktp_factor)]
  pch_vec <- ktp_pch[as.integer(clim_site$ktp_factor)]
  
  n_climate_plot <- sum(
    !is.na(clim_site$mean_temp) & !is.na(clim_site$mean_precip_annual))
  
  cat(sprintf("[10] Climate space by kasvup_tyyppi: %d plots plotted, %d site-type levels\n",
              n_climate_plot, n_ktp))
  
  png(file.path(DIR_OUT, "10_climate_space_site_type.png"),
      width = 10L * PX_PER_IN, height = 8L * PX_PER_IN, res = PX_PER_IN)
  par(mar = c(4.5, 4.5, 3.5, 1))
  
  plot(clim_site$mean_temp, clim_site$mean_precip_annual,
       pch = pch_vec, col = col_vec,
       cex = 0.9, lwd = 1.2,
       xlab = "Mean annual temperature (\u00b0C)",
       ylab = "Mean annual precipitation (mm yr\u207b\u00b9)",
       main = sprintf(
         "Climate space of calibration plots by site fertility type  (n = %d)",
         n_climate_plot))
  
  legend("topleft",
         legend = ktp_levels,
         col    = ktp_colours,
         pch    = ktp_pch,
         pt.cex = 1.1, pt.lwd = 1.2,
         bty    = "n", cex  = 0.82,
         title  = "kasvup_tyyppi")
  
  dev.off()
  
  cat(sprintf("    Written: %s\n\n",
              file.path(DIR_OUT, "10_climate_space_site_type.png")))
  
}  # end climate-space guard



# =============================================================================
# 11.  SOC trajectory across measurement years (1985 / 2006 / 2024)
# =============================================================================
# Boxplot of 1m-extrapolated soc_obs_tCha by observation year, restricted to
# calib-ready plots. Two panels:
#   Left  — all calib-ready plots with an SOC observation in each year
#   Right — plots present in ALL THREE years (tightest comparison)
#
# Uses soc_obs_tCha (placed in June of each observation year in input_raw_monthly)
# so values are fully consistent with what enters the calibration likelihood.
#
# The 1985→2006 jump is expected to be large and likely reflects a measurement
# protocol shift between survey waves rather than a real sink signal. Document
# in M&M if the jump is confirmed here.
# =============================================================================

soc_traj <- input_calib[!is.na(input_calib$soc_obs_tCha) &
                          input_calib$month == 6L,
                        c("plot_id", "year", "soc_obs_tCha")]

obs_years   <- sort(unique(soc_traj$year))
plots_by_yr <- lapply(obs_years, function(y)
  unique(soc_traj$plot_id[soc_traj$year == y]))
names(plots_by_yr) <- as.character(obs_years)

common_plots_traj <- Reduce(intersect, plots_by_yr)

cat(sprintf("[11] SOC trajectory: %d obs years detected: %s\n",
            length(obs_years), paste(obs_years, collapse = ", ")))
cat(sprintf("     Plots per year: %s\n",
            paste(sapply(plots_by_yr, length), collapse = " / ")))
cat(sprintf("     Plots in ALL years: %d\n", length(common_plots_traj)))

yr_cols <- c("grey70", "steelblue", "tomato")
names(yr_cols) <- as.character(obs_years)

png(file.path(DIR_OUT, "11_soc_trajectory_by_year.png"),
    width = 14L * PX_PER_IN, height = 6L * PX_PER_IN, res = PX_PER_IN)
par(mfrow = c(1, 2), mar = c(4, 4.5, 3.5, 1))

# Left panel — all calib-ready plots per year
soc_traj$year_f <- factor(soc_traj$year, levels = obs_years)
boxplot(soc_obs_tCha ~ year_f, data = soc_traj,
        col    = yr_cols[as.character(obs_years)],
        border = "grey30", outline = TRUE,
        xlab   = "Observation year",
        ylab   = "SOC 1 m extrapolated (Mg C/ha)",
        main   = sprintf("All calib-ready plots per year\n(n = %s)",
                         paste(sapply(plots_by_yr, length), collapse = " / ")))

# Right panel — common plots only
soc_common <- soc_traj[soc_traj$plot_id %in% common_plots_traj, ]
boxplot(soc_obs_tCha ~ year_f, data = soc_common,
        col    = yr_cols[as.character(obs_years)],
        border = "grey30", outline = TRUE,
        xlab   = "Observation year",
        ylab   = "SOC 1 m extrapolated (Mg C/ha)",
        main   = sprintf("Common plots across all years\n(n = %d each)",
                         length(common_plots_traj)))

dev.off()

cat(sprintf("    Written: %s\n\n",
            file.path(DIR_OUT, "11_soc_trajectory_by_year.png")))