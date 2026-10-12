# =============================================================================
# litter_reconciliation_ghgi.R -- national inventory litter vs the Tupek plot product
# (2026-10-12)
# -----------------------------------------------------------------------------
# WHY. The trilemma (storyline v4, Intro 6-7) rests on the SHAPE of the litter series:
# the plot product (Tupek et al., Zenodo 10.5281/zenodo.19736499) rises ~45% to 2006 and
# is flat after, while the national GHG inventory's litter grows over the whole period.
# This script puts the two on a comparable footing and shows where they differ.
#
# SOURCES.
#   GHGI  Data/GHG_datasets/GHGI_mineral_soils_data/ (received 2026-10-12, gitignored):
#         mineral_litter_input_by_type.csv -- tC/ha/yr per region (1 South, 2 North) and
#         year, by source (living_trees, harvest_residues, natural_mortality,
#         ground_vegetation) and litter type (nwl, fwl, cwl). Per hectare of mineral
#         forest land in the region (mineral_area.csv gives the areas).
#   Tupek Data/model_inputs/input_raw_monthly.csv -- per plot, monthly, summed to the year;
#         scope = tree turnover + harvest residues (no natural mortality, no understorey).
#         Calibration-ready plots, unweighted mean within each region.
# ASSUMPTION: the GHGI region code matches site_raw$region (1 South, 2 North; our region 2
#   plots lie above 66 N). Confirm with LUKE.
#
# SCOPES compared:
#   GHGI total            all four sources
#   GHGI trees            minus ground vegetation
#   GHGI trees, no mort.  living trees + harvest residues = the Tupek scope
#
# Run from repo root:  Rscript doublechecks/litter_reconciliation_ghgi.R
# Output: doublechecks/figures/litter_reconciliation_ghgi.{png,csv}
# =============================================================================
D <- "Data/GHG_datasets/GHGI_mineral_soils_data"
L <- read.csv(file.path(D, "mineral_litter_input_by_type.csv"))
A <- read.csv(file.path(D, "mineral_area.csv"))
REG <- c("1" = "South", "2" = "North")
Y0 <- 1985L; Y1 <- 2024L; BASE <- 1986:1990

# --- GHGI by scope, region, year ---------------------------------------------
g_scope <- function(src) {
  x <- aggregate(C_total ~ region + year, L[L$source %in% src, ], sum)
  setNames(x, c("region", "year", "J"))
}
G <- list(
  total   = g_scope(c("living_trees", "harvest_residues", "natural_mortality", "ground_vegetation")),
  trees   = g_scope(c("living_trees", "harvest_residues", "natural_mortality")),
  nomort  = g_scope(c("living_trees", "harvest_residues")))
# years where all four sources exist (ground vegetation and living trees start later)
ok_years <- with(aggregate(source ~ region + year, L, function(s) length(unique(s))), cbind(region, year)[source == 4, ])
okk <- paste(ok_years[, 1], ok_years[, 2])
G <- lapply(G, function(x) x[paste(x$region, x$year) %in% okk, ])
by_src <- aggregate(C_total ~ region + year + source, L, sum)

# --- Tupek: plot means by region ---------------------------------------------
ir <- read.csv("Data/model_inputs/input_raw_monthly.csv")
st <- read.csv("Data/model_inputs/site_raw.csv")
cal <- st[st$calib_ready, c("plot_id", "region")]
ir <- ir[ir$plot_id %in% cal$plot_id, ]
lc <- grep("^C_(nwl|fwl|cwl)_", names(ir), value = TRUE)
py <- aggregate(ir[, lc], list(plot_id = ir$plot_id, year = ir$year), sum)
py$J <- rowSums(py[, lc])
for (k in c("nwl", "fwl", "cwl")) py[[k]] <- rowSums(py[, grep(paste0("^C_", k, "_"), lc, value = TRUE)])
py <- merge(py, cal, by = "plot_id")
T <- aggregate(cbind(J, nwl, fwl, cwl) ~ region + year, py, mean)
n_plots <- table(cal$region)

# --- numbers: levels and changes ---------------------------------------------
idx <- function(x, yrs) { b <- mean(x$J[x$year %in% BASE]); x$J / b }
rate <- function(x, a, b) { s <- x[x$year >= a & x$year <= b, ]; unname(coef(lm(J ~ year, s))[2]) }
summ <- do.call(rbind, lapply(names(REG), function(r) {
  rows <- list(GHGI_total = G$total, GHGI_trees = G$trees, GHGI_trees_no_mortality = G$nomort, Tupek = T)
  do.call(rbind, lapply(names(rows), function(n) {
    x <- rows[[n]]; x <- x[x$region == as.integer(r), ]
    b <- mean(x$J[x$year %in% BASE])
    data.frame(region = REG[[r]], series = n,
               J_1986_90 = b, J_2004_08 = mean(x$J[x$year %in% 2004:2008]),
               J_2020_24 = mean(x$J[x$year %in% 2020:2024]),
               change_to_2006_pct = 100 * (mean(x$J[x$year %in% 2004:2008]) / b - 1),
               change_to_2024_pct = 100 * (mean(x$J[x$year %in% 2020:2024]) / b - 1),
               slope_1986_2006 = rate(x, 1986, 2006), slope_2006_2024 = rate(x, 2006, 2024))
  }))
}))
num <- sapply(summ, is.numeric); summ[num] <- lapply(summ[num], round, 3)
print(summ, row.names = FALSE)
dir.create("doublechecks/figures", showWarnings = FALSE)
write.csv(summ, "doublechecks/figures/litter_reconciliation_ghgi.csv", row.names = FALSE)

# --- figure --------------------------------------------------------------------
C_G <- "#2a78d6"; C_T <- "#eb6834"                     # inventory blue, plot product orange
SRC <- c(living_trees = "#1baf7a", harvest_residues = "#eda100",
         natural_mortality = "#6250d6", ground_vegetation = "#008300")
SRC_LAB <- c(living_trees = "living trees (turnover)", harvest_residues = "harvest residues",
             natural_mortality = "natural mortality", ground_vegetation = "ground vegetation")
LT <- c(total = 1, trees = 2, nomort = 3)
INK <- "#52514e"
png("doublechecks/figures/litter_reconciliation_ghgi.png", width = 11, height = 5, units = "in", res = 170)
par(mfrow = c(1, 2), mar = c(4, 4.6, 3, 1), oma = c(2.4, 0, 3.2, 0), mgp = c(2.5, 0.6, 0), las = 1,
    col.axis = INK, fg = "grey55")
# (a, b) of the old 2x2 = the levels on the same scope -> manuscript S18
# (manuscript/figures/build_S18_litter_vs_inventory.R). This figure keeps the change bars.
# (c, d) change since 1986-90, side by side (no stacking): the product, the inventory on the
# same scope, and the inventory with natural mortality added (the source the product lacks)
chg <- function(x, yrs) mean(x$J[x$year %in% yrs]) - mean(x$J[x$year %in% BASE])
C_G2 <- "#8fb6e6"
for (r in names(REG)) {
  ri <- as.integer(r)
  ser <- list(T[T$region == ri, ], G$nomort[G$nomort$region == ri, ], G$trees[G$trees$region == ri, ])
  M <- sapply(ser, function(x) c(chg(x, 2004:2008), chg(x, 2020:2024)))
  colnames(M) <- c("plot product", "inventory, same scope", "inventory + mortality")
  bp <- barplot(t(M), beside = TRUE, col = c(C_T, C_G, C_G2), border = "white", ylim = c(0, 1.15),
                names.arg = c("1986-90 to 2004-08", "1986-90 to 2020-24"),
                ylab = expression("Change in litter input  (tC ha"^-1*" yr"^-1*")"),
                main = sprintf("(%s) %s: how much the litter input changed", if (r == "1") "a" else "b", REG[[r]]))
  abline(h = 0, col = "grey40")
  text(bp, t(M) + 0.03, sprintf("%+.2f", t(M)), cex = 0.8, col = "grey15")
  if (r == "1") legend("topright", bty = "n", cex = 0.85, text.col = "grey15", border = "white",
                       legend = c("plot product (Tupek)", "inventory: tree turnover + harvest residues (same scope)",
                                  "inventory: the same + natural mortality"),
                       fill = c(C_T, C_G, C_G2))
}
mtext("Litter input change: national GHG inventory against the plot product (levels: manuscript S18)", outer = TRUE,
      line = 1.3, cex = 1.1, font = 2, col = "grey10")
mtext(paste0("Inventory: per hectare of mineral forest land in each region (mineral_litter_input_by_type.csv); ",
             "change in 5-yr means; same scope = tree turnover + harvest residues (ground vegetation: constant, omitted)."),
      side = 1, outer = TRUE, line = 0.2, cex = 0.62, col = INK, adj = 0.01)
mtext(paste0("Plot product: unweighted mean of the calibration plots in the region. ",
             "Region codes assumed to match (1 South, 2 North)."),
      side = 1, outer = TRUE, line = 1.2, cex = 0.62, col = INK, adj = 0.01)
dev.off()
cat("saved doublechecks/figures/litter_reconciliation_ghgi.{png,csv}\n")
