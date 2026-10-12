# =============================================================================
# build_S18_litter_vs_inventory.R   (2026-10-12)
#
# SUPPLEMENTARY FIGURE S18: the litter input of the plot product used for calibration
# (Tupek et al.) against the national GHG inventory's estimate for mineral forest soils,
# on the scope the two share (tree turnover + harvest residues), South and North.
# Cited in the draft's Methods, subsection "Litter input".
#
#   inventory  Data/GHG_datasets/GHGI_mineral_soils_data/mineral_litter_input_by_type.csv
#              (received 2026-10-12; gitignored): tC/ha/yr per region and year, by source.
#              Scope used: living_trees + harvest_residues (natural mortality and ground
#              vegetation excluded, because the plot product carries neither).
#   product    Data/model_inputs/input_raw_monthly.csv, summed to the year; calibration-
#              ready plots, unweighted mean within each region.
# ASSUMPTION: inventory region codes = site_raw$region (1 South, 2 North). Confirm with LUKE.
# The change-by-interval diagnostic (with natural mortality) is
# doublechecks/litter_reconciliation_ghgi.R.
#
# Run from repo root:  Rscript manuscript/figures/build_S18_litter_vs_inventory.R
# =============================================================================
L  <- read.csv("Data/GHG_datasets/GHGI_mineral_soils_data/mineral_litter_input_by_type.csv")
ir <- read.csv("Data/model_inputs/input_raw_monthly.csv")
st <- read.csv("Data/model_inputs/site_raw.csv")
REG <- c("1" = "Southern Finland", "2" = "Northern Finland")
Y0 <- 1985L; Y1 <- 2024L

G <- aggregate(C_total ~ region + year, L[L$source %in% c("living_trees", "harvest_residues"), ], sum)
names(G)[3] <- "J"
# years where the living-tree term exists (absent 1970-71 in the South, which made a false dip)
lt <- unique(L[L$source == "living_trees", c("region", "year")])
G  <- G[paste(G$region, G$year) %in% paste(lt$region, lt$year), ]
cal <- st[st$calib_ready, c("plot_id", "region")]
ir <- ir[ir$plot_id %in% cal$plot_id, ]
lc <- grep("^C_(nwl|fwl|cwl)_", names(ir), value = TRUE)
py <- aggregate(list(J = rowSums(ir[, lc])), list(plot_id = ir$plot_id, year = ir$year), sum)
py <- merge(py, cal, by = "plot_id")
T  <- aggregate(J ~ region + year, py, mean)
n_plots <- table(cal$region)

C_G <- "#2a78d6"; C_T <- "#eb6834"; INK <- "#52514e"
png("manuscript/figures/S18_litter_vs_inventory.png", width = 11, height = 4.8, units = "in", res = 200)
par(mfrow = c(1, 2), mar = c(4, 4.6, 2.6, 1), oma = c(2, 0, 0, 0), mgp = c(2.5, 0.6, 0), las = 1,
    col.axis = INK, fg = "grey55")
yl <- c(0, max(c(G$J, T$J)) * 1.05)              # from zero (Lorenzo, 2026-10-12)
for (r in names(REG)) {
  plot(NA, xlim = c(1970, 2025), ylim = yl, xlab = "Year",
       ylab = expression("Litter input  (tC ha"^-1*" yr"^-1*")"), yaxs = "i",
       main = sprintf("(%s) %s", if (r == "1") "a" else "b", REG[[r]]))
  abline(h = pretty(yl), col = "grey92"); rect(Y0, -1, Y1, yl[2] + 1, col = "#00000008", border = NA)
  x <- G[G$region == as.integer(r), ]; lines(x$year, x$J, col = C_G, lwd = 2.4)
  x <- T[T$region == as.integer(r), ]; lines(x$year, x$J, col = C_T, lwd = 2.4)
  legend("bottomright", bty = "n", cex = 0.82, text.col = "grey15", y.intersp = 1.4,
         legend = c("national inventory, per hectare of mineral forest land",
                    sprintf("plot product used here, mean of %d study plots", n_plots[[r]])),
         col = c(C_G, C_T), lwd = 2.4)
}
mtext(paste0("Both series: litter from the turnover of living trees and from harvest residues. ",
             "The inventory's natural-mortality and ground-vegetation terms are excluded. Shaded: 1985-2024."),
      side = 1, outer = TRUE, line = 0.6, cex = 0.66, col = INK, adj = 0.01)
dev.off()
cat("saved manuscript/figures/S18_litter_vs_inventory.png\n")
