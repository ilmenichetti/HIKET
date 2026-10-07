# =============================================================================
# residuals_time_since_cut.R — the residual against years since the last
# regeneration cut, at each soil sampling
#
# WHY. Two readings of the local residual compete. (a) It follows the 1985 STAND
# (residuals_stand_and_soil.R: basal area 1985 alone R² 0.42). (b) It follows the
# DISTURBANCE HISTORY: plots open in 1985 were cut shortly before, and are the
# most under-predicted; mature 1985 stands were cut later and are over-predicted
# (residuals_events_over_time.R: −0.21 in the 25 years after a regeneration
# felling). If (b), a transient input after cutting (pioneer vegetation) that no
# tree-based litter product records is a candidate for the missing flux. One
# variable per reading, same observations, and see which survives the other.
#
# WHAT. One row per soil observation; residual = mean of Yasso07 and Yasso15 log
# residuals for that plot-year. T = years from the last regeneration cut to the
# sampling year, using only events BEFORE the sampling year:
#   DATED events (regeneration cut = codes 7, 8; 2021-23 also P, small clearcut)
#     NFI 1985 stand record: hakk_laatu + hakk_aika class (this summer 0,
#       previous season 1, 2-5 yr 3.5, 6-10 yr 8) before 1985
#     NFI 1990 / 1995: hakk_laatu / tehd_hakkuut + years before 1990 / 1995
#     MUSTIKKA 2021-23: hakkuu1/2 + hakkuu_aika1/2 years before 2022
#     Metsäkeskus regeneration-felling declarations (filing year, from 1997)
#   INFERRED, when no dated event precedes the sampling: stand age 1985 from the
#     tally sample trees + years since 1985, i.e. time since stand establishment.
# Models (lme4, plot random intercept; ML for AIC):
#   ba    r ~ ns(basal area 1985, 3)
#   T     r ~ ns(log1p(T), 3)
#   both  r ~ ns(basal area 1985, 3) + ns(log1p(T), 3)
#   each also + fertility + texture + stoniness, as in residuals_stand_and_soil.R
# Figure (kept SEPARATE from the stand-and-soil figure):
#   top    field-layer cover 2023 (herbs + grasses + raspberry), 2024 samplings
#   middle raw residual vs T; the two post-1985 samplings of a plot joined
#   bottom residual after the 1985 basal-area curve vs T
#
# Run from repo root, after residuals_rf_everything.R:
#   Rscript doublechecks/residuals_time_since_cut.R
# Writes doublechecks/figures/residuals_time_since_cut.png
# =============================================================================

suppressPackageStartupMessages({library(dplyr); library(splines); library(lme4)})
source("manuscript/figures/run_ids.R")
OUT <- "doublechecks/figures"
TGT <- c("Yasso07", "Yasso15")
B <- readRDS(file.path(OUT, "residuals_rf_everything_Yasso07_Yasso15.rds"))
X <- B$X; X$plot_id <- as.numeric(B$plot_id)
num <- function(x) suppressWarnings(as.numeric(as.character(x)))

# === 1. Observations ===========================================================
R <- bind_rows(lapply(TGT, function(m)
  as.data.frame(readRDS(sprintf(
    "Calibration_real_data_transient/runs/%s_posterior_predictive_%s.rds",
    m, RID[[m]]))$residuals_df)[, c("plot_id", "year", "residual_log")])) |>
  group_by(plot_id, year) |> summarise(r = mean(residual_log), n_mod = n(), .groups = "drop") |>
  filter(n_mod == length(TGT)) |> mutate(plot_id = as.numeric(plot_id))
R$campaign <- ifelse(R$year < 2000, "1985", ifelse(R$year < 2015, "2006", "2024"))

# === 2. Dated regeneration cuts ================================================
REGEN <- c("7", "8", "P")
ev <- list()
add <- function(pid, yr, src) {
  ok <- !is.na(yr); if (any(ok)) ev[[length(ev) + 1]] <<- data.frame(plot_id = pid[ok], year = yr[ok], src = src)
}
cls85 <- c(`1` = 0, `2` = 1, `3` = 3.5, `4` = 8)
g <- as.character(X$k85.hakk_laatu) %in% REGEN
add(X$plot_id[g], 1985 - cls85[as.character(X$k85.hakk_aika[g])], "NFI 1985")
g <- as.character(X$k90.hakk_laatu) %in% REGEN
a <- num(X$k90.hakk_aika[g]); a[is.na(a)] <- 3                       # timing missing: mid-period
add(X$plot_id[g], 1990 - a, "NFI 1990")
g <- as.character(X$k95.tehd_hakkuut) %in% REGEN
a <- num(X$k95.tehd_hakk_aika[g]); a[is.na(a)] <- 2
add(X$plot_id[g], 1995 - a, "NFI 1995")
for (k in 1:2) {
  g <- as.character(X[[paste0("m23.hakkuu", k)]]) %in% REGEN
  add(X$plot_id[g], 2022 - num(X[[paste0("m23.hakkuu_aika", k)]][g]), "MUSTIKKA 2021-23")
}
D <- read.csv("Data/GIS_points/events/harvest_declaration_events.csv")
D <- D[D$type == "regeneration", ]
add(as.numeric(D$plot_id), D$year, "declaration")
EV <- do.call(rbind, ev)
cat("Dated regeneration cuts by source:\n"); print(table(EV$src))

# === 3. Years since the last regeneration at each sampling =====================
age85 <- setNames(num(X$tal.stand_age_85_tally), X$plot_id)
R$T <- NA_real_; R$T_src <- NA_character_
for (i in seq_len(nrow(R))) {
  e <- EV$year[EV$plot_id == R$plot_id[i] & EV$year < R$year[i]]
  if (length(e)) { R$T[i] <- R$year[i] - max(e); R$T_src[i] <- "dated" }
  else if (!is.na(a0 <- age85[as.character(R$plot_id[i])])) {
    R$T[i] <- a0 + (R$year[i] - 1985); R$T_src[i] <- "inferred" }
}
cat("\nT source by campaign:\n"); print(table(R$campaign, R$T_src, useNA = "a"))

# === 4. Covariates and models ==================================================
ft <- num(X$cov.kasvup_tyyppi); ml <- num(X$k85.maap_laatu)
P <- data.frame(plot_id = X$plot_id, ba = X$wb.ppa_kaikki_85,
  fert  = factor(c("fertile", "fertile", "mesic", "subxeric", "xeric")[ft]),
  tex   = factor(ifelse(ml %in% c(3, 5, 6), "coarse", ifelse(ml == 4, "fine moraine",
                 ifelse(ml %in% 7:9, "fine", NA)))),
  stone = cut(X$cov.CoarseFragments, c(-Inf, 20, 50, Inf)),
  field = X$und.herbs + X$und.graminoids + X$und.raspberry)
d <- merge(R, P, by = "plot_id")
d$lT <- log1p(d$T)
dm <- d[complete.cases(d[, c("r", "ba", "lT", "fert", "tex", "stone")]), ]
cat(sprintf("\nModel sample: %d observations on %d plots\n", nrow(dm), length(unique(dm$plot_id))))

fitm <- function(f) lmer(f, dm, REML = FALSE)
r2m <- function(m) { v <- var(predict(m, re.form = NA)); v / (v + sum(as.data.frame(VarCorr(m))$vcov)) }
soil <- "+ fert + tex + stone"
M <- list(
  ba        = fitm(r ~ ns(ba, 3) + (1 | plot_id)),
  T         = fitm(r ~ ns(lT, 3) + (1 | plot_id)),
  both      = fitm(r ~ ns(ba, 3) + ns(lT, 3) + (1 | plot_id)),
  ba_soil   = fitm(as.formula(paste("r ~ ns(ba, 3)", soil, "+ (1 | plot_id)"))),
  T_soil    = fitm(as.formula(paste("r ~ ns(lT, 3)", soil, "+ (1 | plot_id)"))),
  both_soil = fitm(as.formula(paste("r ~ ns(ba, 3) + ns(lT, 3)", soil, "+ (1 | plot_id)"))))
tab <- data.frame(model = names(M), AIC = round(sapply(M, AIC), 1),
                  marginal_R2 = round(sapply(M, r2m), 3))
tab$dAIC <- round(tab$AIC - min(tab$AIC), 1); rownames(tab) <- NULL
cat("\nPer-observation mixed models (plot random intercept):\n"); print(tab)
lr <- function(a, b) anova(M[[a]], M[[b]])[2, "Pr(>Chisq)"]
cat(sprintf("\nT beyond basal area: p %.2g (no soil), %.2g (with soil)\n", lr("ba", "both"), lr("ba_soil", "both_soil")))
cat(sprintf("Basal area beyond T: p %.2g (no soil), %.2g (with soil)\n", lr("T", "both"), lr("T_soil", "both_soil")))
cat(sprintf("cor(basal area 1985, log T) over observations: %.2f\n", cor(dm$ba, dm$lT)))
cat("Same, dated only:\n")
dd <- dm[dm$T_src == "dated", ]
if (nrow(dd) > 30) {
  m1 <- lmer(r ~ ns(ba, 3) + (1 | plot_id), dd, REML = FALSE)
  m2 <- lmer(r ~ ns(ba, 3) + ns(lT, 3) + (1 | plot_id), dd, REML = FALSE)
  cat(sprintf("  %d obs: T beyond basal area p %.2g\n", nrow(dd), anova(m1, m2)[2, "Pr(>Chisq)"]))
}

# Residual after the basal-area curve (fixed part of the 'ba' model)
d$r_ba <- NA
ok <- !is.na(d$ba); d$r_ba[ok] <- d$r[ok] - predict(M$ba, newdata = d[ok, ], re.form = NA)

# === 5. Figure =================================================================
INK <- "#3d3d3a"; GRID <- "#e6e5df"; MUTED <- "#8a8983"
CAMP_COL <- c(`1985` = "#7a5cc4", `2006` = "#2a78d6", `2024` = "#eb6834")
BRK <- c(0, 5, 10, 20, 30, 45, 60, 90, 400)
tx  <- function(t) log1p(t)                                  # x axis: log(1 + T)
TCK <- c(0, 2, 5, 10, 20, 40, 80, 150)
binned <- function(x, y) {
  b <- cut(x, BRK, right = FALSE)
  s <- data.frame(b, x, y)[!is.na(x) & !is.na(y), ] |> group_by(b) |>
    summarise(xm = median(x), m = mean(y), se = sd(y) / sqrt(n()), n = n(), .groups = "drop")
  s[s$n >= 8, ]
}
draw_bins <- function(s) {
  arrows(tx(s$xm), s$m - 1.96 * s$se, tx(s$xm), s$m + 1.96 * s$se, angle = 90, code = 3,
         length = 0.04, col = INK, lwd = 1.6)
  lines(tx(s$xm), s$m, col = INK, lwd = 2.2); points(tx(s$xm), s$m, pch = 21, bg = "white", col = INK, cex = 1.3, lwd = 1.6)
  text(tx(s$xm), par("usr")[4], s$n, pos = 1, cex = 0.7, col = MUTED)
}
xax <- function() { axis(1, tx(TCK), TCK); abline(v = tx(TCK), col = GRID, lwd = 0.6) }
pts <- function(y, lab) {
  yl <- quantile(y, c(.005, .995), na.rm = TRUE)
  plot(NA, xlim = tx(c(0, 200)), ylim = yl, bty = "l", xaxt = "n", xlab = "", ylab = lab)
  xax(); abline(h = 0, col = GRID, lwd = 2)
}

set.seed(2025)
png(file.path(OUT, "residuals_time_since_cut.png"), 9.5, 11, units = "in", res = 150)
layout(matrix(1:3), heights = c(0.75, 1, 1))
par(mar = c(3.2, 4.4, 2.8, 1), oma = c(3.2, 0, 2.2, 0), mgp = c(2.4, 0.6, 0), tcl = -0.3,
    col.axis = INK, col.lab = INK, fg = INK)

# top: field layer cover 2023 on the 2024 samplings
c24 <- d[d$campaign == "2024" & !is.na(d$T) & !is.na(d$field), ]
pts(c24$field, "Field-layer cover 2023 (%)")
points(tx(c24$T) + rnorm(nrow(c24), 0, 0.02), c24$field, pch = ifelse(c24$T_src == "dated", 16, 1),
       cex = 0.6, col = adjustcolor(CAMP_COL["2024"], 0.6))
draw_bins(binned(c24$T, c24$field))
mtext("Herbs + grasses + raspberry, 2024 samplings only; binned means ± 95% CI, n per bin on top", 3, 0.4, adj = 0, cex = 0.7)

# middle: raw residual, post-1985 pairs joined
pts(d$r, "Log residual (obs / pred)")
pp <- d[d$campaign %in% c("2006", "2024") & !is.na(d$T), ] |> group_by(plot_id) |> filter(n() == 2) |> arrange(year)
for (p in split(pp, pp$plot_id)) lines(tx(p$T), p$r, col = adjustcolor(MUTED, 0.18), lwd = 0.7)
dt <- d[!is.na(d$T), ]
points(tx(dt$T), dt$r, pch = ifelse(dt$T_src == "dated", 16, 1), cex = 0.55,
       col = adjustcolor(CAMP_COL[dt$campaign], 0.6))
draw_bins(binned(dt$T, dt$r))
mtext("Raw residual; grey lines join the 2006 and 2024 samplings of one plot", 3, 0.4, adj = 0, cex = 0.7)
legend("bottomleft", c(paste(names(CAMP_COL), "sampling"), "dated cut", "inferred from stand age"),
       pch = c(16, 16, 16, 16, 1), col = c(CAMP_COL, INK, INK), bty = "n", cex = 0.8, ncol = 2)

# bottom: residual after the basal-area curve
db <- d[!is.na(d$T) & !is.na(d$r_ba), ]
pts(db$r_ba, "Residual after 1985 basal area (log)")
points(tx(db$T), db$r_ba, pch = ifelse(db$T_src == "dated", 16, 1), cex = 0.55,
       col = adjustcolor(CAMP_COL[db$campaign], 0.6))
draw_bins(binned(db$T, db$r_ba))
mtext(sprintf("After removing the 1985 basal-area curve (mixed model, plot intercept). T beyond basal area: p %.2g; with soil classes: p %.2g",
              lr("ba", "both"), lr("ba_soil", "both_soil")), 3, 0.4, adj = 0, cex = 0.7)
mtext("Years since the last regeneration cut, at sampling (log scale)", 1, 2.3, cex = 0.85)

mtext("Residual against time since the last regeneration cut (mean of Yasso07 and Yasso15; + = models predict too little SOC)",
      3, 0.6, outer = TRUE, cex = 0.85)
mtext(paste("Dated cuts: NFI stand records 1985/1990/1995 (codes 7-8), MUSTIKKA 2021-23 (7, 8, P), Metsäkeskus regeneration declarations (filing year).",
            "\nInferred: no dated cut before the sampling, so stand age in 1985 (tally sample trees) + years since 1985 = time since stand establishment."),
      1, 1.6, outer = TRUE, cex = 0.62, col = MUTED, adj = 0, at = 0.02)
dev.off()
cat("Wrote", file.path(OUT, "residuals_time_since_cut.png"), "\n")

# === 6. Summary figure: the time course, overall and by 1985 stand ============
# Solid: all observations (raw residual). Dashed: the same residual split by
# basal area 1985 (< 10, 10-20, > 20 m2/ha, about thirds of the observations),
# coarser bins so each point has >= 18 observations. All lines are the SAME
# quantity (mean log obs/pred), so they share the axis; the classes differ in
# level. A class line that still rises toward 20-30 yr = time since cut matters
# at a given 1985 stand. Right axis in % SOC (exp(r) - 1).
BRK_CL <- c(0, 10, 20, 30, 50, 90, 400)
binned_cl <- function(x, y) { BRK_old <- BRK; BRK <<- BRK_CL; on.exit(BRK <<- BRK_old); binned(x, y) }
dcl <- d[!is.na(d$T) & !is.na(d$ba), ]
dcl$cl <- cut(dcl$ba, c(-Inf, 10, 20, Inf), labels = c("< 10", "10-20", "> 20"))
S_all <- binned(dt$T, dt$r)
S_cl  <- lapply(split(dcl, dcl$cl), function(z) binned_cl(z$T, z$r))
COL_CL <- c(`< 10` = "#1baf7a", `10-20` = "#2a78d6", `> 20` = "#eb6834")
LTY_CL <- c(`< 10` = 2, `10-20` = 4, `> 20` = 5)
png(file.path(OUT, "residuals_time_since_cut_summary.png"), 9, 6, units = "in", res = 150)
par(mar = c(5.2, 4.6, 4.2, 4.6), mgp = c(2.5, 0.6, 0), tcl = -0.3,
    col.axis = INK, col.lab = INK, fg = INK)
all_lo <- unlist(lapply(c(list(S_all), S_cl), function(s) s$m - 1.96 * s$se))
all_hi <- unlist(lapply(c(list(S_all), S_cl), function(s) s$m + 1.96 * s$se))
yl <- range(all_lo, all_hi) + c(-0.03, 0.1)
plot(NA, xlim = tx(c(1.5, 160)), ylim = yl, bty = "u", xaxt = "n",
     xlab = "Years since the last regeneration cut, at sampling (log scale)",
     ylab = "Mean log residual (obs / pred)")
rect(tx(20), yl[1], tx(30), yl[2], col = adjustcolor(INK, 0.05), border = NA)
xax(); abline(h = 0, col = GRID, lwd = 2)
pc <- c(-40, -30, -20, -10, 0, 10, 20, 30, 50, 100, 150)
pc <- pc[log(1 + pc / 100) > yl[1] & log(1 + pc / 100) < yl[2]]
axis(4, log(1 + pc / 100), sprintf("%+d%%", pc), las = 1)
mtext("Observed SOC relative to predicted", 4, 3.2, col = INK)
dodge <- c(`< 10` = -0.03, `10-20` = 0, `> 20` = 0.03)
for (k in names(S_cl)) {
  s <- S_cl[[k]]; x <- tx(s$xm) + dodge[[k]]
  arrows(x, s$m - 1.96 * s$se, x, s$m + 1.96 * s$se, angle = 90, code = 3, length = 0.03,
         col = adjustcolor(COL_CL[[k]], 0.6), lwd = 1.2)
  lines(x, s$m, col = COL_CL[[k]], lwd = 2, lty = LTY_CL[[k]])
  points(x, s$m, pch = 21, bg = COL_CL[[k]], col = "white", cex = 1.1)
}
arrows(tx(S_all$xm), S_all$m - 1.96 * S_all$se, tx(S_all$xm), S_all$m + 1.96 * S_all$se,
       angle = 90, code = 3, length = 0.04, col = INK, lwd = 1.6)
lines(tx(S_all$xm), S_all$m, col = INK, lwd = 2.8)
points(tx(S_all$xm), S_all$m, pch = 21, bg = "white", col = INK, cex = 1.3, lwd = 1.8)
legend("topright", c("all observations", paste("basal area 1985", names(S_cl), "m2/ha")),
       col = c(INK, COL_CL[names(S_cl)]), lty = c(1, LTY_CL[names(S_cl)]), lwd = c(2.8, 2, 2, 2),
       pch = 21, pt.bg = c("white", COL_CL[names(S_cl)]), bty = "n", cex = 0.82, bg = "white")
mtext("Residual over time since the last regeneration cut, overall and by 1985 stand", 3, 2.4, adj = 0, font = 2, cex = 1.05)
mtext(sprintf("Binned means ± 95%% CI, %d soil observations, mean of Yasso07 and Yasso15. Shaded: 20-30 yr.",
              nrow(dt)), 3, 1.25, adj = 0, cex = 0.72, col = MUTED)
mtext(sprintf("Time since cut beyond basal area: p %.1g; with soil classes p %.1g; dated cuts only p %.1g.",
              lr("ba", "both"), lr("ba_soil", "both_soil"), anova(m1, m2)[2, "Pr(>Chisq)"]),
      3, 0.35, adj = 0, cex = 0.72, col = MUTED)
mtext(paste("Time since cut = dated cut where recorded, else stand age in 1985 + years since 1985.",
            "Class lines: bins 0-10, 10-20, 20-30, 30-50, 50-90, > 90 yr."),
      1, 3.3, adj = 0, cex = 0.68, col = MUTED, at = par("usr")[1])
dev.off()
cat("\nClass means by time bin:\n")
for (k in names(S_cl)) { cat(k, ": "); cat(sprintf("%.0f:%+.2f(n%d) ", S_cl[[k]]$xm, S_cl[[k]]$m, S_cl[[k]]$n), "\n") }
cat("Wrote", file.path(OUT, "residuals_time_since_cut_summary.png"), "\n")
