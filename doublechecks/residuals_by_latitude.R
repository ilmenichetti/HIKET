# =============================================================================
# residuals_by_latitude.R — does ONE calibration fit North and South alike?
#
# WHY. The calibration is unweighted: each plot is one observation, so the North
# (12% of plots, ~30% of forest area) pulls less than its area. That is correct
# for the likelihood IF the parameters are common to both regions. If they are
# not, the posterior settles on southern behaviour and the North is fitted
# systematically worse. Cheap check: residuals by region and latitude band.
#
# WHAT. On the posterior-mean log residual (log obs - log pred), per model:
#   LEVEL  plot-mean residual, North vs South, and slope on latitude
#   TREND  per-plot change in residual first -> last campaign per decade
#          (a trend misfit that differs by region = rate fitted from the South)
# Plot means first: one value per plot, so plot-years do not pseudo-replicate.
#
# Run from repo root:  Rscript doublechecks/residuals_by_latitude.R
# =============================================================================

suppressPackageStartupMessages(library(dplyr))
source("manuscript/figures/run_ids.R")
cat("RUN_IDs:", paste(names(RID), RID, sep = "=", collapse = "  "), "\n\n")

rd <- function(m) {
  f <- sprintf("Calibration_real_data_transient/runs/%s_posterior_predictive_%s.rds", m, RID[[m]])
  as.data.frame(readRDS(f)$residuals_df) |>
    transmute(model = m, plot_id, year, r = residual_log, lat = lat_WGS84,
              north = region == 2, holdout = is_holdout)
}
R <- bind_rows(lapply(FIG_MODELS, rd))

# --- LEVEL: one mean residual per plot ---------------------------------------
L <- R |> group_by(model, plot_id, north, lat) |> summarise(r = mean(r), .groups = "drop")
lev <- L |> group_by(model) |> group_modify(function(d, k) {
  t <- t.test(d$r[d$north], d$r[!d$north])
  s <- summary(lm(r ~ lat, d))$coefficients
  tibble(n_S = sum(!d$north), n_N = sum(d$north),
         res_S = mean(d$r[!d$north]), res_N = mean(d$r[d$north]),
         N_minus_S = diff(rev(t$estimate)), p = t$p.value,
         slope_per_deg = s["lat", 1], p_slope = s["lat", 4])
}) |> ungroup()

# --- TREND: residual change first -> last observation, per decade ------------
Tr <- R |> group_by(model, plot_id, north) |> filter(n() >= 2) |>
  arrange(year, .by_group = TRUE) |>
  summarise(dr = (last(r) - first(r)) / (last(year) - first(year)) * 10, .groups = "drop")
tre <- Tr |> group_by(model) |> group_modify(function(d, k) {
  t <- t.test(d$dr[d$north], d$dr[!d$north])
  tibble(n_S = sum(!d$north), n_N = sum(d$north),
         dres_S = mean(d$dr[!d$north]), dres_N = mean(d$dr[d$north]),
         N_minus_S = diff(rev(t$estimate)), p = t$p.value)
}) |> ungroup()

# --- 8 equal-count latitude bands (as tau_R), plot-mean level residual -------
L$band <- cut(L$lat, quantile(unique(L[c("plot_id","lat")])$lat, 0:8/8),
              include.lowest = TRUE, labels = 1:8)
bands <- L |> group_by(model, band) |> summarise(r = mean(r), .groups = "drop") |>
  tidyr::pivot_wider(names_from = band, values_from = r, names_prefix = "b")

options(width = 140, pillar.sigfig = 3)
cat("LEVEL residual (log scale; + = model under-predicts), plot means\n"); print(lev, n = 6)
cat("\nTREND residual change per decade (+ = model under-predicts the gain)\n"); print(tre, n = 6)
cat("\nLEVEL residual by latitude band (1 = south ... 8 = north)\n"); print(bands, n = 6)
# --- IS THE NORTHERN OFFSET BASAL AREA? (added 2026-09-30) -------------------
# The North has sparser stands, and sparse stands are underpredicted everywhere
# (manuscript sec:plotscale). If the offset is that gradient, basal area absorbs it
# and the SOUTHERN slope alone predicts it; temperature should not absorb it.
BA <- bind_rows(lapply(FIG_MODELS, function(m) {
  f <- sprintf("Calibration_real_data_transient/runs/%s_posterior_predictive_%s.rds", m, RID[[m]])
  as.data.frame(readRDS(f)$residuals_df) |> group_by(plot_id) |>
    summarise(model = m, r = mean(residual_log), north = region[1] == 2,
              ba = basal_area_85[1], T = mean_temp[1], .groups = "drop")
})) |> filter(is.finite(ba))
ba_tab <- BA |> group_by(model) |> group_modify(function(d, k) {
  S <- d[!d$north, ]; N <- d[d$north, ]
  bS <- coef(lm(r ~ ba, S))[["ba"]]
  cN <- function(f) coef(summary(lm(f, d)))["northTRUE", c(1, 4)]
  a <- cN(r ~ north); b <- cN(r ~ north + ba); t <- cN(r ~ north + T)
  Sm <- S[S$ba <= quantile(N$ba, 0.9), ]
  tibble(ba_N = mean(N$ba), ba_S = mean(S$ba), slope_S = bS,
         NmS = a[1], p = a[2], pred_from_ba = bS * (mean(N$ba) - mean(S$ba)),
         NmS_given_ba = b[1], p_ba = b[2], NmS_given_T = t[1], p_T = t[2],
         res_S_sparse = mean(Sm$r), res_N = mean(N$r))
}) |> ungroup()
cat("\nNORTH-SOUTH offset vs basal area (plot means; + = under-prediction)\n"); print(ba_tab, n = 6, width = 200)

saveRDS(list(run_ids = RID, level = lev, trend = tre, bands = bands, basal_area = ba_tab),
        "doublechecks/residuals_by_latitude.rds")
