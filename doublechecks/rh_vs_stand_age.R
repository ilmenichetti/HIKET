# =============================================================================
# rh_vs_stand_age.R — modelled heterotrophic respiration against stand age,
# laid out like Samuli's flux-site figure (Rh vs Age, colour = air temperature,
# size = fertility), for comparison
#
# WHAT. All calibrated plots, 1986-2024. Per model, plot and year:
#   Rh = J * median(sigma_input) - (C_t - C_{t-1})   [g C m-2 yr-1]
#   (carbon balance on the posterior-mean trajectory; Yasso15/20 include leaching)
# Stand age = tally stand age in 1985 + years since 1985, reset to 0 at a dated
#   regeneration cut after 1985 (NFI 1990/1995, MUSTIKKA 2021-23, declarations;
#   as in residuals_time_since_cut.R).
# Colour = plot mean annual air temperature (scale -2...6 °C, our range); size = fertility, high = NFI site
#   type 1-3 (fertile, mesic), low = 4-5. All plots are mineral soils.
# Panels: a) mean of Yasso07 and Yasso15 (the residual-study reference), one
#   point per plot every 5 years.
#   b) binned medians, six models.
#
# Run from repo root, after residuals_rf_everything.R:
#   Rscript doublechecks/rh_vs_stand_age.R
# Writes doublechecks/figures/rh_vs_stand_age.png
# =============================================================================

suppressPackageStartupMessages(library(dplyr))
source("manuscript/figures/run_ids.R")
OUT <- "doublechecks/figures"
num <- function(x) suppressWarnings(as.numeric(as.character(x)))

# === 1. Plot attributes, stand age =============================================
B <- readRDS(file.path(OUT, "residuals_rf_everything_Yasso07_Yasso15.rds"))
X <- B$X; X$plot_id <- as.numeric(B$plot_id)
REGEN <- c("7", "8", "P"); ev <- list()
add <- function(pid, yr) { ok <- !is.na(yr); if (any(ok)) ev[[length(ev) + 1]] <<- data.frame(plot_id = pid[ok], year = round(yr[ok])) }
g <- as.character(X$k90.hakk_laatu) %in% REGEN; a <- num(X$k90.hakk_aika[g]); a[is.na(a)] <- 3; add(X$plot_id[g], 1990 - a)
g <- as.character(X$k95.tehd_hakkuut) %in% REGEN; a <- num(X$k95.tehd_hakk_aika[g]); a[is.na(a)] <- 2; add(X$plot_id[g], 1995 - a)
for (k in 1:2) { g <- as.character(X[[paste0("m23.hakkuu", k)]]) %in% REGEN
  add(X$plot_id[g], 2022 - num(X[[paste0("m23.hakkuu_aika", k)]][g])) }
D <- read.csv("Data/GIS_points/events/harvest_declaration_events.csv"); D <- D[D$type == "regeneration", ]
add(as.numeric(D$plot_id), D$year)
EV <- do.call(rbind, ev); EV <- EV[EV$year > 1985, ]

lm_ <- readRDS(sprintf("Data/model_inputs/Yasso15_inputs_%s.rds", RID[["Yasso15"]]))
R1 <- as.data.frame(readRDS(sprintf("Calibration_real_data_transient/runs/Yasso15_posterior_predictive_%s.rds",
                                    RID[["Yasso15"]]))$residuals_df)
PA <- R1 |> group_by(plot_id) |> summarise(T = first(mean_temp), .groups = "drop") |> mutate(plot_id = as.numeric(plot_id))
PA$age85 <- num(X$tal.stand_age_85_tally)[match(PA$plot_id, X$plot_id)]
ft <- num(X$cov.kasvup_tyyppi)[match(PA$plot_id, X$plot_id)]
PA$fert <- ifelse(ft <= 3, "high", ifelse(ft <= 5, "low", NA))

J0 <- bind_rows(lapply(lm_$inputs_by_plot, function(z)
  data.frame(plot_id = z$plot_id, year = z$year, J = rowSums(z[, grep("^(nwl|fwl|cwl)_", names(z))]))))
age_of <- function(pid, yr, a85) {           # vectorised over one plot's years
  e <- EV$year[EV$plot_id == pid[1]]
  sapply(yr, function(y) { c0 <- e[e <= y]; if (length(c0)) y - max(c0) else a85[1] + (y - 1985) })
}

# === 2. Rh per model, plot and year ============================================
RH <- bind_rows(lapply(FIG_MODELS, function(m) {
  si <- median(readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_%s.rds", m, RID[[m]]))[, "sigma_input"])
  P <- readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_predictive_%s.rds",
                       m, RID[[m]]))$posterior_summary[, c("plot_id", "year", "soc_mean")]
  P$plot_id <- as.numeric(P$plot_id)
  P |> inner_join(J0, by = c("plot_id", "year")) |> arrange(plot_id, year) |> group_by(plot_id) |>
    mutate(Rh = 100 * (J * si - (soc_mean - lag(soc_mean)))) |> ungroup() |>
    filter(!is.na(Rh)) |> mutate(model = m)
})) |> inner_join(PA, by = "plot_id") |> filter(!is.na(age85))
AG <- RH |> distinct(plot_id, year, age85) |> group_by(plot_id) |>
  mutate(age = age_of(plot_id, year, age85)) |> ungroup()
RH <- RH |> inner_join(AG[, c("plot_id", "year", "age")], by = c("plot_id", "year"))
cat("Plots:", n_distinct(RH$plot_id), " plot-years per model:", nrow(RH) / length(FIG_MODELS), "\n")

REF <- RH |> filter(model %in% c("Yasso07", "Yasso15")) |> group_by(plot_id, year) |>
  summarise(Rh = mean(Rh), age = first(age), T = first(T), fert = first(fert), .groups = "drop")
BRK <- c(0, 5, 10, 15, 20, 30, 40, 50, 60, 80, 100, 130, 300)
BIN <- RH |> mutate(b = cut(age, BRK, right = FALSE)) |> group_by(model, b) |>
  summarise(age = median(age), Rh = median(Rh), q1 = quantile(Rh, .25), q3 = quantile(Rh, .75), n = n(), .groups = "drop") |>
  filter(n >= 30)
cat("\nBinned median Rh (g C m-2 yr-1) by stand age:\n")
print(tidyr::pivot_wider(BIN[, c("model", "b", "Rh")], names_from = model, values_from = Rh) |>
        mutate(across(-b, round)) |> as.data.frame())
y20 <- REF |> filter(age <= 20)
cat(sprintf("\nYasso07+15, first 20 yr: linear slope %.1f g C m-2 yr-1 per yr (n %d plot-years)\n",
            coef(lm(Rh ~ age, y20))[2], nrow(y20)))
for (m in FIG_MODELS) { z <- RH |> filter(model == m, age <= 20)
  cat(sprintf("  %-8s slope %.1f\n", m, coef(lm(Rh ~ age, z))[2])) }

# === 3. Figure =================================================================
INK <- "#3d3d3a"; GRID <- "#e6e5df"; MUTED <- "#8a8983"
tpal <- colorRampPalette(c("#3a5fcd", "#b9c9e8", "#f1d6c9", "#d9603b", "#b2182b"))(101)
TR <- c(-2, 6)                                            # our plots: -1.9 ... 5.7 °C (Samuli's sites run to 10)
tcol <- function(t) tpal[pmin(101, pmax(1, round(100 * (t - TR[1]) / diff(TR)) + 1))]
MCOL <- c(SP1 = "#8a8983", TP2 = "#1baf7a", TP3 = "#0b6a48",
          Yasso07 = "#86b6ef", Yasso15 = "#2a78d6", Yasso20 = "#0d366b")
pt <- REF |> filter(year %% 5 == 0, !is.na(fert))
set.seed(2025); pt <- pt[sample(nrow(pt)), ]
YL <- c(0, 1300)

png(file.path(OUT, "rh_vs_stand_age.png"), 13.5, 5.6, units = "in", res = 150)
layout(matrix(1:2, 1), widths = c(1.15, 1))
par(mar = c(4.2, 4.6, 3, 1), oma = c(1.6, 0, 0, 0), mgp = c(2.5, 0.6, 0), tcl = -0.3,
    col.axis = INK, col.lab = INK, fg = INK)
plot(NA, xlim = c(0, 210), ylim = YL, bty = "l", xlab = "Stand age (yr)",
     ylab = expression(R[h]~"(g C"~m^-2~a^-1*")"))
points(pt$age, pt$Rh, pch = 21, bg = adjustcolor(tcol(pt$T), 0.85), col = adjustcolor(INK, 0.5),
       cex = ifelse(pt$fert == "high", 1.15, 0.7), lwd = 0.4)
mtext("a) Mean of Yasso07 and Yasso15, one point per plot every 5 years (1990-2020)", 3, 0.6, adj = 0, cex = 0.8, font = 2)
legend("topright", c("-2", "0", "2", "4", "6"), pt.bg = tcol(c(-2, 0, 2, 4, 6)), pch = 21,
       title = "Tair (°C)", bty = "n", cex = 0.8, col = INK)
legend(160, 860, c("High (1-3)", "Low (4-5)"), pch = 21, pt.cex = c(1.15, 0.7), pt.bg = "grey80", col = INK,
       title = "Fertility", bty = "n", cex = 0.8)
plot(NA, xlim = c(0, 210), ylim = YL, bty = "l", xlab = "Stand age (yr)",
     ylab = expression(R[h]~"(g C"~m^-2~a^-1*")"))
for (m in FIG_MODELS) { s <- BIN[BIN$model == m, ]
  polygon(c(s$age, rev(s$age)), c(s$q1, rev(s$q3)), col = adjustcolor(MCOL[[m]], 0.08), border = NA)
  lines(s$age, s$Rh, col = MCOL[[m]], lwd = 2.2) }
legend("topright", FIG_MODELS, col = MCOL, lwd = 2.2, bty = "n", cex = 0.8, ncol = 2)
mtext("b) Six models: binned median and interquartile range", 3, 0.6, adj = 0, cex = 0.8, font = 2)
mtext(paste("Modelled Rh = litter input x posterior median sigma_input - annual change in posterior-mean SOC (Yasso15/20 include leaching).",
            "Stand age = tally age 1985 + years, reset at dated regeneration cuts after 1985. Mineral soils only."),
      1, 0.4, outer = TRUE, cex = 0.62, col = MUTED)
dev.off()
cat("Wrote", file.path(OUT, "rh_vs_stand_age.png"), "\n")
