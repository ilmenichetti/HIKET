# =============================================================================
# residuals_events_over_time.R — does management BEFORE a soil sampling show in
# that sampling's residual?
#
# WHY. The per-plot GIS summaries collapse 1997-2026 into one number and set it
# against the plot-MEAN residual, so a harvest in 2010 is compared with the
# 1985/2006 misfit too, and declarations after 2024 are counted. Here each SOIL
# OBSERVATION is one point, with its own residual and only the events dated
# strictly BEFORE its sampling year (the sampling date inside the year is unknown,
# so same-year events are left out).
#
# WHAT. Residual = mean of Yasso07 and Yasso15 log residuals for that plot-year.
# Events: Metsäkeskus forest use declarations (harvest INTENTIONS, filing date;
# Data/GIS_points/events/harvest_declaration_events.csv) and KEMERA completed works
# (kemera_events.csv), both from GIS/. Records start 1997 (declarations dense
# from 2004), so only observations from 2000 on are shown.
# Panels: years since last regeneration felling / thinning (none = own group),
# declarations in the 10 years before sampling, damage-driven cutting before
# sampling, years since last young-stand tending, fertilised before sampling.
#
# Run from repo root:  Rscript doublechecks/residuals_events_over_time.R
# Writes doublechecks/figures/residuals_events_over_time.png
# =============================================================================

suppressPackageStartupMessages(library(dplyr))
source("manuscript/figures/run_ids.R")
OUT <- "doublechecks/figures"
TGT <- c("Yasso07", "Yasso15")
INK <- "#3d3d3a"; GRID <- "#e6e5df"
CAMP_COL <- c(`2006` = "#2a78d6", `2024` = "#eb6834")

# --- one row per soil observation ---------------------------------------------
R <- bind_rows(lapply(TGT, function(m)
  as.data.frame(readRDS(sprintf(
    "Calibration_real_data_transient/runs/%s_posterior_predictive_%s.rds",
    m, RID[[m]]))$residuals_df)[, c("plot_id", "year", "residual_log")])) |>
  group_by(plot_id, year) |> summarise(r = mean(residual_log), n_mod = n(), .groups = "drop") |>
  filter(n_mod == length(TGT), year >= 2000)
R$campaign <- ifelse(R$year < 2015, "2006", "2024")
cat("Observations from 2000 on:", nrow(R), " by sampling year:\n"); print(table(R$year))

D <- read.csv("Data/GIS_points/events/harvest_declaration_events.csv")
K <- read.csv("Data/GIS_points/events/kemera_events.csv")

before <- function(ev, pid, yr, f) {   # apply f to the event years before yr
  e <- ev$year[ev$plot_id == pid & ev$year < yr]; f(e)
}
last_gap <- function(type) mapply(function(p, y)
  before(D[D$type == type, ], p, y, function(e) if (length(e)) y - max(e) else NA), R$plot_id, R$year)
R$yrs_since_regen  <- last_gap("regeneration")
R$yrs_since_thin   <- last_gap("thinning")
R$n_decl_10yr      <- mapply(function(p, y) before(D, p, y, function(e) sum(e >= y - 10)), R$plot_id, R$year)
R$damage_before    <- mapply(function(p, y) before(D[D$type == "damage", ], p, y, length) > 0, R$plot_id, R$year)
R$yrs_since_tend   <- mapply(function(p, y) before(K[K$group == "tending", ], p, y,
                               function(e) if (length(e)) y - max(e) else NA), R$plot_id, R$year)
R$fert_before      <- mapply(function(p, y) before(K[K$group == "fertilise", ], p, y, length) > 0, R$plot_id, R$year)

# --- drawing -------------------------------------------------------------------
YL <- quantile(R$r, c(.005, .995))
panel_gap <- function(x, lab) {                 # continuous gap + a "none" column
  none <- is.na(x); xmax <- max(x, na.rm = TRUE); xn <- xmax + 4
  plot(NA, xlim = c(0, xn + 1.5), ylim = YL, bty = "l", xaxt = "n",
       xlab = lab, ylab = "Log residual (obs/pred), this sampling")
  axis(1, pretty(c(0, xmax))); axis(1, xn, "none", tick = FALSE)
  abline(h = 0, col = GRID, lwd = 2)
  xx <- ifelse(none, xn + runif(length(x), -0.8, 0.8), x + runif(length(x), -0.3, 0.3))
  points(xx, R$r, pch = 16, cex = 0.6, col = adjustcolor(CAMP_COL[R$campaign], 0.55))
  if (sum(!none) >= 15) {
    lo <- loess(R$r[!none] ~ x[!none], span = 0.9)
    xs <- seq(min(x, na.rm = TRUE), xmax, length.out = 60)
    lines(xs, predict(lo, xs), col = INK, lwd = 2.5)
    rho <- cor(x[!none], R$r[!none], method = "spearman")
  } else rho <- NA
  segments(xn - 1, mean(R$r[none]), xn + 1, mean(R$r[none]), lwd = 3, col = INK)
  mtext(sprintf("with event n %d, Spearman ρ %.2f | none n %d, mean %.2f vs event mean %.2f",
                sum(!none), rho, sum(none), mean(R$r[none]), mean(R$r[!none])), 3, 0.3, adj = 0, cex = 0.62)
}
panel_group <- function(g, lab, levs) {          # discrete: jittered columns + means
  g <- factor(g, levels = levs); gi <- as.integer(g)
  plot(NA, xlim = c(0.5, nlevels(g) + 0.5), ylim = YL, bty = "l", xaxt = "n",
       xlab = lab, ylab = "Log residual (obs/pred), this sampling")
  axis(1, seq_along(levs), levs); abline(h = 0, col = GRID, lwd = 2)
  points(gi + runif(length(gi), -0.25, 0.25), R$r, pch = 16, cex = 0.6,
         col = adjustcolor(CAMP_COL[R$campaign], 0.55))
  m <- tapply(R$r, g, mean); n <- tapply(R$r, g, length)
  segments(seq_along(levs) - 0.3, m, seq_along(levs) + 0.3, m, lwd = 3, col = INK)
  p <- if (all(n >= 3, na.rm = TRUE)) tryCatch(summary(lm(R$r ~ g))$coefficients[2, 4], error = function(e) NA) else NA
  mtext(sprintf("%s | linear p %.2g", paste(sprintf("%s n %d mean %.2f", levs, n, m), collapse = "; "), p),
        3, 0.3, adj = 0, cex = 0.62)
}

set.seed(2025)
png(file.path(OUT, "residuals_events_over_time.png"), 14, 8.8, units = "in", res = 150)
par(mfrow = c(2, 3), mar = c(4.2, 4, 2.6, 1), oma = c(0, 0, 2.4, 0), mgp = c(2.3, 0.6, 0),
    tcl = -0.3, col.axis = INK, col.lab = INK, fg = INK)
panel_gap(R$yrs_since_regen, "Years since last regeneration-felling declaration")
panel_gap(R$yrs_since_thin, "Years since last thinning declaration")
panel_group(pmin(R$n_decl_10yr, 3), "Declarations in the 10 years before sampling (3 = 3+)", 0:3)
panel_group(ifelse(R$damage_before, "yes", "no"), "Damage-driven cutting declared before sampling", c("no", "yes"))
panel_gap(R$yrs_since_tend, "Years since last young-stand tending (KEMERA, completed)")
panel_group(ifelse(R$fert_before, "yes", "no"), "Remedial fertilisation before sampling (KEMERA)", c("no", "yes"))
legend("topright", c("2006 sampling", "2024 sampling"), pch = 16, col = CAMP_COL, bty = "n", cex = 0.9)
mtext(sprintf("Each point = one soil observation (n %d); residual = mean of Yasso07 and Yasso15; only events dated BEFORE the sampling year",
              nrow(R)), 3, 0.7, outer = TRUE, cex = 0.85, col = INK)
dev.off()

cat("\nObservations with an event before sampling:\n")
print(sapply(R[, c("yrs_since_regen", "yrs_since_thin", "yrs_since_tend")], function(x) sum(!is.na(x))))
print(table(damage = R$damage_before, campaign = R$campaign)); print(table(fert = R$fert_before, campaign = R$campaign))
cat("Wrote", file.path(OUT, "residuals_events_over_time.png"), "\n")
