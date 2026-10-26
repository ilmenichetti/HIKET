# =============================================================================
# ghgi_steady_state_vs_obs.R   (2026-10-26)
#
# Aleksi's question (email 2026-10-26): if Yasso runs with the inventory input and
# DEFAULT parameters from a steady state in the 1970s -- as the GHG inventory does --
# does it fall outside the uncertainty of the measured change?
#
# No new run needed: the inventory's own Yasso07 run is in the LUKE package
# (Data/GHG_datasets/GHGI_mineral_soils_data/, received 2026-10-12, gitignored):
#   mineral_stocks_and_changes.csv  stock per region (1 South, 2 North) x year x pool
#                                   (soc, litter, deadwood); steady state at the first
#                                   year (South 1972, North 1975; mineral_spinup.csv)
#   mineral_area.csv                mineral forest area per region (weights)
#
# COMPARABILITY. Inventory stocks have no depth basis, observations run to 1 m (or
# refusal) -> LEVELS are not comparable, only CHANGES. Everything is shown as change
# since the first campaign, and as rates over the campaign intervals.
#   weighting  both regions weighted by inventory mineral area (1990: 55/45) -- the
#              average hectare, i.e. what the inventory reports. Observations and our
#              models are averaged per region first, then weighted the same way.
#   plot set   balanced (all three campaigns), true observation years (obs_basis.R).
#   timing     inventory read at each region's mean campaign years.
#   regions    shown SEPARATELY: the national agreement is a CANCELLATION -- the South
#              gains faster than the inventory run, the North (38 plots) loses carbon
#              while the inventory run gains. A national-only figure would hide that.
# Our calibrated models are left out on purpose (they are in F17); the question here is
# the inventory's own run.
#   dead wood  observations exclude coarse dead wood; the inventory's 'deadwood' pool
#              is drawn both in and out (dashed = without).
# ASSUMPTION: inventory region code == site_attributes$region (as in
# litter_reconciliation_ghgi.R). Confirm with LUKE.
# Uncertainty: observed = 95% CI over plots (per-plot paired changes); models = 90%
# over posterior draws; inventory = none given.
#
# Run from repo root:  Rscript doublechecks/ghgi_steady_state_vs_obs.R
# Output: doublechecks/figures/ghgi_steady_state_vs_obs.{png,csv}
# =============================================================================
suppressMessages({
  source("manuscript/figures/run_ids_eqinit.R"); source("manuscript/figures/obs_basis.R")
})
G   <- "Data/GHG_datasets/GHGI_mineral_soils_data/"
IVS <- list("1985-2006" = c(1L, 2L), "2006-2024" = c(2L, 3L), "1985-2024" = c(1L, 3L))

# === Inventory: total and without dead wood, per region x year ===
d  <- read.csv(paste0(G, "mineral_stocks_and_changes.csv"))
ar <- read.csv(paste0(G, "mineral_area.csv"))
W  <- with(ar[ar$year == 1990, ], setNames(area_ha / sum(area_ha), region))
inv_tot <- with(d, tapply(stock_tC_ha, list(region, year), sum))
inv_nod <- with(d[d$crt_pool != "deadwood", ], tapply(stock_tC_ha, list(region, year), sum))
inv_at  <- function(M, r, y) {                       # linear in year (mean campaign years are fractional)
  yy <- as.numeric(colnames(M)); approx(yy, M[as.character(r), ], y)$y }

# === Observations and plots, by region ===
sa  <- read.csv("Data/model_inputs/site_attributes.csv")
om  <- readRDS(arm_files("Yasso07", "prod")$inp)$obs_meta
bal <- balanced_plots(om)
reg <- setNames(sa$region[match(bal, sa$plot_id)], bal)
stopifnot(!anyNA(reg))
cy <- do.call(rbind, lapply(as.character(bal), function(p) {   # per plot: campaign years + obs
  yr <- 1984L + om[[p]]$idx; cp <- campaign_of(yr); s <- om[[p]]$soc_obs
  data.frame(plot_id = p, region = reg[[p]],
             y1 = yr[cp == 1][1], y2 = yr[cp == 2][1], y3 = yr[cp == 3][1],
             s1 = s[cp == 1][1], s2 = s[cp == 2][1], s3 = s[cp == 3][1])
}))
cat(sprintf("balanced plots: %d (South %d, North %d)\n", nrow(cy), sum(cy$region == 1), sum(cy$region == 2)))
Ymean <- sapply(1:2, function(r) colMeans(cy[cy$region == r, c("y1", "y2", "y3")]))  # 3 x 2

# --- per-plot paired quantities -> South, North, and area-weighted national ---
summ <- function(x) {
  m  <- tapply(x, cy$region, mean); se <- tapply(x, cy$region, function(v) sd(v) / sqrt(length(v)))
  mn <- sum(W[names(m)] * m); sen <- sqrt(sum((W[names(se)] * se)^2))
  data.frame(area = c("South", "North", "National"), med = c(m, mn),
             lo = c(m, mn) - 1.96 * c(se, sen), hi = c(m, mn) + 1.96 * c(se, sen))
}
GROUPS <- c("South", "North", "National")
obs_rate <- do.call(rbind, lapply(names(IVS), function(iv) { ab <- IVS[[iv]]
  x <- (cy[[paste0("s", ab[2])]] - cy[[paste0("s", ab[1])]]) /
       (cy[[paste0("y", ab[2])]] - cy[[paste0("y", ab[1])]])
  cbind(interval = iv, summ(x)) }))
obs_chg <- lapply(1:2, function(r) rbind(c(Ymean[1, r], 0, 0, 0),
  t(sapply(2:3, function(k) { x <- (cy[[paste0("s", k)]] - cy$s1)[cy$region == r]
    se <- sd(x) / sqrt(length(x)); c(Ymean[k, r], mean(x), mean(x) - 1.96 * se, mean(x) + 1.96 * se) }))))

# --- inventory rates on the same timing ---
inv_rate <- function(M) do.call(rbind, lapply(names(IVS), function(iv) { ab <- IVS[[iv]]
  r <- sapply(1:2, function(r) (inv_at(M, r, Ymean[ab[2], r]) - inv_at(M, r, Ymean[ab[1], r])) /
                               (Ymean[ab[2], r] - Ymean[ab[1], r]))
  data.frame(interval = iv, area = GROUPS, med = c(r, sum(W * r))) }))
it <- inv_rate(inv_tot); ino <- inv_rate(inv_nod)

tab <- merge(merge(obs_rate, setNames(it, c("interval", "area", "inv_total"))),
             setNames(ino, c("interval", "area", "inv_no_deadwood")))
tab$inv_inside_CI <- tab$inv_total >= tab$lo & tab$inv_total <= tab$hi
tab <- tab[order(match(tab$interval, names(IVS)), match(tab$area, GROUPS)), ]
write.csv(tab, "doublechecks/figures/ghgi_steady_state_vs_obs.csv", row.names = FALSE)
num <- sapply(tab, is.numeric); tab[num] <- round(tab[num], 3); print(tab, row.names = FALSE)

# === Figure (print size, 9 pt) ===
png("doublechecks/figures/ghgi_steady_state_vs_obs.png", width = 6.5, height = 5.6, units = "in",
    res = 300, pointsize = 9)
layout(rbind(c(1, 1, 1, 2, 2, 2), c(3, 3, 4, 4, 5, 5)), heights = c(1.15, 1), widths = c(1.25, 1, 1, 1, 1, 1))
par(cex = 1); par(mgp = c(2.0, 0.5, 0), las = 1, tcl = -0.3)
OCOL <- "#4477AA"; ICOL <- "grey15"
yl <- range(-8, sapply(obs_chg, function(o) o[, 3:4]),
            sapply(1:2, function(r) inv_tot[r, ] - inv_at(inv_tot, r, Ymean[1, r])), na.rm = TRUE)

# (a, b) change since the first campaign, per region
for (r in 1:2) {
  yy <- as.numeric(colnames(inv_tot)); ok <- !is.na(inv_tot[r, ])
  y0 <- yy[ok][1]
  par(mar = c(2.6, if (r == 1) 4.2 else 2.4, 1.8, 0.6))
  plot(NA, xlim = c(1970, 2026), ylim = yl, xlab = "",
       ylab = if (r == 1) expression("Change since first campaign  (tC ha"^-1*")") else "")
  mtext(sprintf("(%s) %s, %d plots", letters[r], c("South", "North")[r], sum(cy$region == r)),
        3, 0.4, adj = 0, font = 2)
  abline(h = 0, col = "grey60", lty = 2); abline(v = y0, col = "grey75", lty = 3)
  text(y0, yl[2], sprintf("inventory\nsteady state\n%d", y0), adj = c(-0.08, 1), cex = 0.8, col = "grey40")
  lines(yy[ok], inv_tot[r, ok] - inv_at(inv_tot, r, Ymean[1, r]), col = ICOL, lwd = 2)
  lines(yy[ok], inv_nod[r, ok] - inv_at(inv_nod, r, Ymean[1, r]), col = ICOL, lwd = 1.3, lty = 2)
  o <- obs_chg[[r]]
  segments(o[-1, 1], o[-1, 3], o[-1, 1], o[-1, 4], col = OCOL, lwd = 2)
  points(o[, 1], o[, 2], pch = 21, bg = OCOL, col = "grey15", cex = 1.3)
  if (r == 1) legend("bottomright", bty = "n", cex = 0.85,
       legend = c("inventory run, total", "inventory run, without dead wood", "measured (95% CI)"),
       col = c(ICOL, ICOL, OCOL), lty = c(1, 2, NA), lwd = c(2, 1.3, NA), pch = c(NA, NA, 21),
       pt.bg = OCOL, pt.cex = 1.2)
}

# (c-e) rates per interval: measured interval vs inventory, South / North / national
ypos <- setNames(3:1, GROUPS)
xr <- range(tab$lo, tab$hi, tab$inv_total, tab$inv_no_deadwood, 0); xr <- xr + c(-0.04, 0.04) * diff(xr)
for (k in seq_along(IVS)) {
  s <- tab[tab$interval == names(IVS)[k], ]
  par(mar = c(3.3, if (k == 1) 5.2 else 0.6, 1.8, 0.5))
  plot(NA, xlim = xr, ylim = c(0.5, 3.75), yaxt = "n", ylab = "",
       xlab = expression("Change  (tC ha"^-1*" yr"^-1*")"))
  mtext(sprintf("(%s) %s", letters[k + 2], sub("-", " → ", names(IVS)[k])), 3, 0.4, adj = 0, font = 2)
  abline(v = 0, col = "grey35", lty = 2)
  y <- ypos[s$area]
  rect(s$lo, y - 0.22, s$hi, y + 0.22, col = adjustcolor(OCOL, 0.25), border = OCOL)
  segments(s$med, y - 0.22, s$med, y + 0.22, col = OCOL, lwd = 2)
  points(s$inv_total, y, pch = 15, cex = 1.4, col = ICOL)
  points(s$inv_no_deadwood, y, pch = 0, cex = 1.3, col = ICOL)
  if (k == 1) axis(2, at = ypos, labels = c("South", "North", "National\n(area-wt.)"), tick = FALSE)
  if (k == 2) legend("topright", bty = "n", cex = 0.8, bg = "white", legend = c("inventory, total", "without dead wood"),
                     pch = c(15, 0), col = ICOL, pt.cex = 1.2)
}
dev.off()
cat("\nWrote doublechecks/figures/ghgi_steady_state_vs_obs.png\n")
