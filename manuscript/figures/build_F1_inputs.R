setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
# F1 (redesign) -- INPUTS-FIRST. The litter input is the protagonist of the opening;
# climate is demoted to a supporting context strip. Three views of the inputs:
#   (a) temporal composition  -- what rises, and the rise-then-plateau shape
#   (b) spatial gradient       -- south->north productivity decline + plot spread
#   (c) climate context        -- interannual only, no accumulation trend (support)
# Data: Yasso20 bundle (only the Yasso bundles carry the nwl/fwl/cwl size-class split).

pkg <- readRDS("Data/model_inputs/Yasso20_inputs_20260710_102431.rds")
ib  <- do.call(rbind, pkg$inputs_by_plot); ib <- ib[ib$year >= 1986, ]
cb  <- do.call(rbind, pkg$climate_by_plot); cb <- cb[cb$year >= 1986, ]
pi  <- pkg$plot_info
gf  <- function(f) sapply(pi, function(p) if (is.null(p[[f]])) NA else p[[f]])
lat <- setNames(as.numeric(gf("lat_WGS84")), gf("plot_id"))

# --- size-class time series (cross-plot mean) ---
nw <- rowSums(ib[, grep("^nwl_", names(ib))]); fw <- rowSums(ib[, grep("^fwl_", names(ib))])
cw <- rowSums(ib[, grep("^cwl_", names(ib))])
ym <- function(v) tapply(v, ib$year, mean, na.rm = TRUE)
NW <- ym(nw); FW <- ym(fw); CW <- ym(cw); yr <- as.integer(names(NW))
peak_y <- yr[which.max(NW + FW + CW)]

# --- spatial: per-plot mean litter + latitude band, band time series (IQR) ---
tot_py <- rowSums(ib[, c(grep("^nwl_|^fwl_|^cwl_", names(ib)))])
band   <- cut(lat[as.character(ib$plot_id)],
              breaks = quantile(lat, c(0, 1/3, 2/3, 1), na.rm = TRUE),
              labels = c("South", "Central", "North"), include.lowest = TRUE)
bcol <- c(South = "#c62828", Central = "#ef8f00", North = "#1565c0")
qtile <- function(b, p) tapply(tot_py[band == b], ib$year[band == b], quantile, probs = p, na.rm = TRUE)

# --- climate context (monthly -> annual per plot-year, then cross-plot mean) ---
py <- interaction(cb$plot_id, cb$year, drop = TRUE)
aT <- tapply(cb$temp_air, py, mean, na.rm = TRUE)     # annual mean temperature
aP <- tapply(cb$precip,   py, sum,  na.rm = TRUE)      # annual precip total
yof <- as.integer(sub(".*\\.", "", names(aT)))
Tm <- tapply(aT, yof, mean, na.rm = TRUE)
Pm <- tapply(aP, yof, mean, na.rm = TRUE)
cyr <- as.integer(names(Tm))

png("manuscript/figures/F1_inputs.png", width = 8.6, height = 9.6, units = "in", res = 200)
layout(matrix(c(1,1, 2,2, 3,4), nrow = 3, byrow = TRUE), heights = c(1.15, 1.0, 0.62))
par(mar = c(3.4, 4.6, 2.4, 1.2), mgp = c(2.5, 0.7, 0), las = 1, cex.axis = 1.05, cex.lab = 1.2, cex.main = 1.15)

## (a) temporal composition -- stacked size classes -------------------------
ccol <- c(nonwoody = "#2e7d32", finewoody = "#8d6e63", coarsewoody = "#4e342e")
b1 <- NW; b2 <- NW + FW; b3 <- NW + FW + CW
plot(NA, xlim = range(yr), ylim = c(0, max(b3) * 1.08), xlab = "", ylab = expression("Litter input (tC ha"^-1*" yr"^-1*")"),
     main = "(a)  Litter input: composition and the rise-then-plateau")
polyb <- function(lo, hi, col) polygon(c(yr, rev(yr)), c(lo, rev(hi)), col = col, border = NA)
polyb(rep(0, length(yr)), b1, ccol["nonwoody"])
polyb(b1, b2, ccol["finewoody"])
polyb(b2, b3, ccol["coarsewoody"])
lines(yr, b3, lwd = 2.2, col = "grey15")
abline(v = peak_y, lty = 3, col = "grey40")
text(peak_y, max(b3) * 1.05, sprintf("peak %d", peak_y), pos = 2, cex = 0.8, col = "grey30")
text(yr[2], b1[2] * 0.5, "non-woody\n(foliage + fine root)", pos = 4, cex = 0.78, col = "white", font = 2)
text(2018, (b1[yr==2018]+b2[yr==2018])/2, "fine woody", pos = 3, cex = 0.72, col = "white", font = 2)
legend("bottomright", inset = c(0.01, 0.04), bty = "n", cex = 0.8, fill = ccol, border = NA,
       legend = c("non-woody litter", "fine woody", "coarse woody"))
mtext(sprintf("total  %.1f -> %.1f (peak %d) -> %.1f tC/ha/yr;  non-woody drives ~55%% of the rise",
              b3[1], max(b3), peak_y, b3[length(b3)]), side = 3, line = -0.1, cex = 0.72, col = "grey35", adj = 0.02)

## (b) spatial gradient -- litter by latitude band (median + IQR) -----------
plot(NA, xlim = range(yr), ylim = c(1.2, 3.8), xlab = "", ylab = expression("Litter input (tC ha"^-1*" yr"^-1*")"),
     main = "(b)  A strong south–north productivity gradient (447 plots, latitude terciles)")
for (b in c("South","Central","North")) {          # all bands first
  lo <- qtile(b, .25); hi <- qtile(b, .75); by <- as.integer(names(lo))
  polygon(c(by, rev(by)), c(lo, rev(hi)), col = adjustcolor(bcol[b], 0.13), border = NA)
}
for (b in c("North","Central","South")) {          # then medians, on top
  md <- qtile(b, .5); by <- as.integer(names(md)); lines(by, md, col = bcol[b], lwd = 2.6)
}
mn <- tapply(tot_py, band, mean, na.rm = TRUE)
legend("topright", bty = "n", cex = 0.82, lwd = 2.6, col = bcol, bg = "white", box.col = NA,
       legend = sprintf("%s  (mean %.2f)", names(bcol), mn[names(bcol)]))
mtext("shaded = interquartile spread across plots within band (the input heterogeneity that sigma_input engages)",
      side = 3, line = -0.1, cex = 0.72, col = "grey35", adj = 0.02)

## (c)(d) climate context -- demoted --------------------------------------
par(mar = c(3.4, 4.4, 2.0, 1.0))
plot(cyr, Tm, type = "l", col = "firebrick", lwd = 1.6, xlab = "Year",
     ylab = expression("T ("*degree*"C)"), main = "(c)  Temperature — context")
plot(cyr, Pm, type = "l", col = "steelblue", lwd = 1.6, xlab = "Year",
     ylab = "Precip (mm)", main = "(d)  Precipitation — context")
dev.off()
cat("Wrote manuscript/figures/F1_inputs.png\n")
