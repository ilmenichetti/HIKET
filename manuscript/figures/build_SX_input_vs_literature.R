# =============================================================================
# SX (APPENDIX / DIAGNOSTIC) -- litter input: our calibration vs the inventory
#                              products vs the independent NPP literature
#
# ⚠ NOT FOR PRODUCTION. Built 2026-08-14 to support the decision on the
# `flux_pair` input-flux window. It is a working diagnostic: several of the
# comparisons on it are NOT like-for-like and are drawn precisely to show that.
# Do NOT cite it, put it in the manuscript, or quote numbers off it without
# reading the caveats in the M&M working document, section "Prior specification:
# the litter-input flux window".
#
# WHAT IT SHOWS. Three families of number that all look like "litter input to
# Finnish forest soil" but are not the same quantity:
#   (1) INVENTORY PRODUCTS (NFI biomass x turnover, the LUKE lineage)
#       - Boris/Tupek J_bar: TREE litter only (incl. harvest residues and
#         natural mortality; understorey EXCLUDED)  -- Zenodo 10.5281/zenodo.19736499
#       - Lehtonen & Heikkinen 2015: TOTAL litter (tree + understorey)
#   (2) INDEPENDENT NPP ESTIMATES (an upper bound: litter <= NPP)
#       - Gower et al. 2001, Nordic Class I stands, and all Class I evergreens
#       - Zheng et al. 2004, gridded, Finland+Sweden  [DRY MATTER -> C x 0.5]
#   (3) OUR CALIBRATED effective flux, sigma_input x J_bar, per model
#
# THE UNIT TRAP THIS FIGURE EXISTS TO RECORD. Gower reports g C m-2 yr-1
# explicitly; Zheng reports DRY MATTER. Zheng cites Gower's world-boreal TNPP as
# 109-1827 (mean 892) where Gower's own table says 218-912 (mean 424) g C -- a
# ratio of exactly 2.00-2.10, which is how the dry-matter reading is established.
# Mixing the two inflates the ceiling by 2x and was the original error.
# =============================================================================

set.seed(2025)
REPO <- "/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling"
setwd(REPO)
OUT <- "manuscript/figures/SX_input_vs_literature.png"

J_BAR <- 2.511   # Boris tree-litter cross-plot mean, corrected-target run 20260813

# --- our calibrated effective flux, per model -------------------------------
source("manuscript/figures/run_ids.R")   # auto-selects current RUN_IDs
rid <- as.list(RID)
ours <- do.call(rbind, lapply(names(rid), function(m) {
  f <- sprintf("Calibration_real_data_transient/runs/%s_posterior_%s.rds", m, rid[[m]])
  if (!file.exists(f)) return(NULL)
  x <- readRDS(f); v <- x[, grep("sigma_input", colnames(x))[1]] * J_BAR
  data.frame(model = m, lo = quantile(v, .025), mid = median(v), hi = quantile(v, .975))
}))

# --- literature -------------------------------------------------------------
# lo/hi are RANGES where reported (not CIs) -- see notes column.
lit <- data.frame(
  label = c("Plot product, trees only",
            "Lehtonen & Heikkinen 2015",
            "Zheng et al. 2004, NPP",
            "Gower et al. 2001, Nordic NPP",
            "Gower et al. 2001, evergreen NPP"),
  mid = c(2.511, 2.70, 2.81, 3.21, 3.87),
  lo  = c(NA,    2.43, 1.26, 2.15, 2.14),
  hi  = c(NA,    2.97, 7.13, 4.62, 9.12),
  kind = c("inv", "inv", "npp", "npp", "npp"),
  stringsAsFactors = FALSE)

# --- candidate window ceilings ----------------------------------------------
# Only the upper bound of the prior window on the flux is drawn (Methods, input prior);
# the Nordic 4.62 and Zheng-mean lines belonged to the abandoned hard window (2026-08-31).
ceil <- data.frame(label = "upper bound of the prior", x = 8.70, col = "grey35", stringsAsFactors = FALSE)

# Print size (0.95 x 6.53 in text width, 9 pt), so nothing prints below ~7 pt.
png(OUT, width = 6.2, height = 4.96, units = "in", res = 300, pointsize = 9)
op <- par(mar = c(3.4, 12.5, 1.0, 0.6), xaxs = "i", mgp = c(2.0, 0.5, 0), tcl = -0.3)
XMAX <- 9.6
n <- nrow(lit) + nrow(ours) + 1
plot(NA, xlim = c(0, XMAX), ylim = c(0.3, n + 0.7), yaxt = "n", bty = "n",
     xlab = expression("Litter input or NPP  (tC ha"^-1*" yr"^-1*")"), ylab = "")

# ceilings first, behind everything
for (i in seq_len(nrow(ceil))) {
  abline(v = ceil$x[i], col = ceil$col[i], lwd = 1.5, lty = 3)
  text(ceil$x[i], n + 0.6, ceil$label[i], srt = 90, adj = c(1, -0.4), cex = 0.85, col = ceil$col[i], xpd = NA)
}

ypos <- n; labs <- c(); cols <- c()
draw <- function(y, lo, mid, hi, col, pch = 19) {
  if (!is.na(lo)) segments(lo, y, hi, y, col = col, lwd = 3, lend = 1)
  points(mid, y, pch = pch, col = col, cex = 1.2)
}
# literature
for (i in seq_len(nrow(lit))) {
  col <- if (lit$kind[i] == "inv") "#009E73" else "#56B4E9"
  draw(ypos, lit$lo[i], lit$mid[i], lit$hi[i], col)
  labs <- c(labs, lit$label[i]); cols <- c(cols, col); ypos <- ypos - 1
}
# separator
abline(h = ypos + 0.5, col = "grey80", lty = 1); ypos <- ypos - 0.6
# ours
for (i in seq_len(nrow(ours))) {
  draw(ypos, ours$lo[i], ours$mid[i], ours$hi[i], "#CC79A7", pch = 18)
  labs <- c(labs, ours$model[i])
  cols <- c(cols, "#CC79A7"); ypos <- ypos - 1
}
axis(2, at = seq(n, by = -1, length.out = nrow(lit)), labels = labs[seq_len(nrow(lit))],
     las = 1, cex.axis = 0.9, tick = FALSE, line = -0.5)
axis(2, at = seq(n - nrow(lit) - 0.6, by = -1, length.out = nrow(ours)),
     labels = labs[-seq_len(nrow(lit))], las = 1, cex.axis = 0.9, tick = FALSE, line = -0.5)

legend("bottomright", bty = "n", cex = 0.85, inset = c(0.06, 0),
       legend = c("inventory-based litter", "independent NPP", "calibrated effective input"),
       col = c("#009E73", "#56B4E9", "#CC79A7"), pch = c(19, 19, 18), lwd = 3, lty = 1)
par(op); dev.off()
cat("wrote", OUT, "\n")
print(ours, row.names = FALSE)
