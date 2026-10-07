# =============================================================================
# harvest_dynamics.R — what the models do to soil C and heterotrophic
# respiration around a regeneration cut (Samuli's question, 2026-10-07)
#
# WHY. "For sites that are harvested, do you get credible dynamics of soil C and
# Rh during early stand development?" The residual analysis looks only at the
# sampling dates; this looks at the modelled annual trajectories.
#
# WHAT. Plots with a DATED regeneration cut in 1988-2015 (NFI 1990/1995 records,
# MUSTIKKA 2021-23, Metsäkeskus declarations; same sources and codes as
# residuals_time_since_cut.R), aligned on the cut year (first such cut per plot).
# Per model, per plot and year (posterior mean trajectory, 100 draws):
#   J   litter input = Tupek input (all AWEN x size classes) x posterior median
#       sigma_input  [tC/ha/yr]
#   C   modelled total SOC (posterior mean)  [tC/ha]
#   Rh  = J - (C_t - C_{t-1})  [tC/ha/yr]; carbon balance, so for Yasso15/20 it
#       also contains the small leaching term
# Averaged over plots by years since the cut, -4 ... +25.
# Observed: soil observations of the same plots, residual by years since cut.
#
# Run from repo root:  Rscript doublechecks/harvest_dynamics.R
# Writes doublechecks/figures/harvest_dynamics.png
# =============================================================================

suppressPackageStartupMessages(library(dplyr))
source("manuscript/figures/run_ids.R")
OUT <- "doublechecks/figures"

# === 1. Dated regeneration cuts (as in residuals_time_since_cut.R) ============
B <- readRDS(file.path(OUT, "residuals_rf_everything_Yasso07_Yasso15.rds"))
X <- B$X; X$plot_id <- as.numeric(B$plot_id)
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
REGEN <- c("7", "8", "P"); ev <- list()
add <- function(pid, yr, src) { ok <- !is.na(yr)
  if (any(ok)) ev[[length(ev) + 1]] <<- data.frame(plot_id = pid[ok], year = round(yr[ok]), src = src) }
g <- as.character(X$k90.hakk_laatu) %in% REGEN; a <- num(X$k90.hakk_aika[g]); a[is.na(a)] <- 3
add(X$plot_id[g], 1990 - a, "NFI 1990")
g <- as.character(X$k95.tehd_hakkuut) %in% REGEN; a <- num(X$k95.tehd_hakk_aika[g]); a[is.na(a)] <- 2
add(X$plot_id[g], 1995 - a, "NFI 1995")
for (k in 1:2) { g <- as.character(X[[paste0("m23.hakkuu", k)]]) %in% REGEN
  add(X$plot_id[g], 2022 - num(X[[paste0("m23.hakkuu_aika", k)]][g]), "MUSTIKKA 2021-23") }
D <- read.csv("Data/GIS_points/events/harvest_declaration_events.csv")
D <- D[D$type == "regeneration", ]; add(as.numeric(D$plot_id), D$year, "declaration")
EV <- do.call(rbind, ev)
CUT <- EV |> filter(year >= 1988, year <= 2015) |> group_by(plot_id) |>
  summarise(cut = min(year), src = src[which.min(year)], .groups = "drop") |>
  filter(cut >= 1989, cut <= 2004)                     # so -3 ... +20 lies inside 1985-2024
cat("Plots with a dated regeneration cut 1989-2004:", nrow(CUT), "\n"); print(table(CUT$src))

# === 2. Input, SOC and Rh per model ============================================
lm_ <- readRDS(sprintf("Data/model_inputs/Yasso15_inputs_%s.rds", RID[["Yasso15"]]))
J0 <- bind_rows(lapply(lm_$inputs_by_plot, function(z)
  data.frame(plot_id = z$plot_id, year = z$year, J = rowSums(z[, grep("^(nwl|fwl|cwl)_", names(z))]))))
REL <- -3:20                      # balanced window: every plot covers all of it
traj <- bind_rows(lapply(FIG_MODELS, function(m) {
  s <- readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_%s.rds", m, RID[[m]]))
  si <- median(s[, "sigma_input"])
  P <- readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_predictive_%s.rds",
                       m, RID[[m]]))$posterior_summary[, c("plot_id", "year", "soc_mean")]
  P$plot_id <- as.numeric(P$plot_id)
  P |> inner_join(J0, by = c("plot_id", "year")) |> inner_join(CUT, by = "plot_id") |>
    arrange(plot_id, year) |> group_by(plot_id) |>
    mutate(J = J * si, Rh = J - (soc_mean - lag(soc_mean)), rel = year - cut,
           C_pre = mean(soc_mean[rel %in% -3:-1]), J_pre = mean(J[rel %in% -3:-1]),
           Rh_pre = mean(Rh[rel %in% -2:-1])) |> ungroup() |>
    filter(rel %in% REL, is.finite(C_pre)) |> group_by(plot_id) |>
    filter(n() == length(REL)) |> ungroup() |> mutate(model = m)
}))
S <- traj |> group_by(model, rel) |>
  summarise(n = n_distinct(plot_id), J = mean(J), J_rel = mean(J / J_pre), C_rel = mean(soc_mean / C_pre),
            dC = mean(soc_mean - C_pre), C = mean(soc_mean), Rh = mean(Rh, na.rm = TRUE),
            J_pre = mean(J_pre), Rh_pre = mean(Rh_pre), .groups = "drop")
cat("\nPlots in the average at selected years since the cut:\n")
print(S |> filter(model == "Yasso15", rel %in% c(-3, 0, 10, 20)) |> select(rel, n))
cat("\nPre-cut levels (mean of years -3..-1; Rh -2..-1):\n")
print(S |> filter(rel == 0) |> transmute(model, J_pre = round(J_pre, 2), Rh_pre = round(Rh_pre, 2),
                                         C_pre = round(C / C_rel, 1)) |> as.data.frame())
cat("\nBy model (years since cut 0, 1, 5, 10, 20):\n")
print(S |> filter(rel %in% c(0, 1, 5, 10, 20)) |>
        transmute(model, rel, J_rel = round(J_rel, 2), C_rel = round(C_rel, 3), dC = round(dC, 2),
                  Rh = round(Rh, 2)) |> as.data.frame())
mn <- S |> group_by(model) |> summarise(C_min_rel = round(min(C_rel), 3), yr_min = rel[which.min(C_rel)],
                                        C_25 = round(C_rel[rel == max(rel)], 3), .groups = "drop")
cat("\nModelled SOC minimum after the cut:\n"); print(as.data.frame(mn))

# Observed: residual of these plots' soil observations by years since the cut
R <- bind_rows(lapply(c("Yasso07", "Yasso15"), function(m)
  as.data.frame(readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_predictive_%s.rds",
                                m, RID[[m]]))$residuals_df)[, c("plot_id", "year", "residual_log")])) |>
  group_by(plot_id, year) |> summarise(r = mean(residual_log), .groups = "drop") |>
  mutate(plot_id = as.numeric(plot_id)) |> inner_join(CUT, by = "plot_id") |> mutate(rel = year - cut)
cat("\nObserved residual (Yasso07+15) on these plots, by years since cut:\n")
print(R |> mutate(b = cut(rel, c(-40, -1, 5, 10, 20, 40))) |> group_by(b) |>
        summarise(n = n(), se = round(sd(r) / sqrt(n()), 3), r = round(mean(r), 3)))

# === 3. Figure =================================================================
INK <- "#3d3d3a"; GRID <- "#e6e5df"; MUTED <- "#8a8983"
MCOL <- c(SP1 = "#8a8983", TP2 = "#1baf7a", TP3 = "#0b6a48",
          Yasso07 = "#86b6ef", Yasso15 = "#2a78d6", Yasso20 = "#0d366b")
panel <- function(v, ylab, h = NULL, title) {
  yl <- range(S[[v]]); plot(NA, xlim = range(REL), ylim = yl, bty = "l", xlab = "Years since regeneration cut",
                            ylab = ylab)
  abline(v = 0, col = GRID, lwd = 2); if (!is.null(h)) abline(h = h, col = GRID, lwd = 2)
  for (m in FIG_MODELS) { s <- S[S$model == m, ]; lines(s$rel, s[[v]], col = MCOL[[m]], lwd = 2.2) }
  mtext(title, 3, 0.4, adj = 0, cex = 0.8, font = 2)
}
png(file.path(OUT, "harvest_dynamics.png"), 13, 4.6, units = "in", res = 150)
par(mfrow = c(1, 3), mar = c(4.2, 4.4, 2.6, 1), oma = c(1.6, 0, 1.6, 0), mgp = c(2.4, 0.6, 0),
    tcl = -0.3, col.axis = INK, col.lab = INK, fg = INK)
panel("J", "Litter input (tC/ha/yr)", title = "Litter input (Tupek x sigma_input)")
panel("C_rel", "SOC relative to the 3 pre-cut years", h = 1, title = "Modelled soil carbon")
legend("bottomleft", FIG_MODELS, col = MCOL, lwd = 2.2, bty = "n", cex = 0.85, ncol = 2)
panel("Rh", "Rh = input - stock change (tC/ha/yr)", title = "Derived heterotrophic respiration")
mtext(sprintf("Mean over the same %d calibrated plots in every year (dated regeneration cut 1989-2004); posterior-mean trajectories",
              S$n[1]), 3, 0.2, outer = TRUE, cex = 0.8)
mtext("Sources: NFI 1990/1995 stand records, MUSTIKKA 2021-23, Metsäkeskus declarations (filing year). Rh includes leaching for Yasso15/20.",
      1, 0.4, outer = TRUE, cex = 0.65, col = MUTED)
dev.off()
cat("Wrote", file.path(OUT, "harvest_dynamics.png"), "\n")
