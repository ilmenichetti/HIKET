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
  label = c("Boris / Tupek  J̄  (tree only, no understorey)",
            "Lehtonen & Heikkinen 2015  (tree + understorey)",
            "Zheng 2004  gridded TNPP  (DM→C)",
            "Gower 2001  Nordic  TNPP",
            "Gower 2001  Class I evergreen  TNPP"),
  mid = c(2.511, 2.70, 2.81, 3.21, 3.87),
  lo  = c(NA,    2.43, 1.26, 2.15, 2.14),
  hi  = c(NA,    2.97, 7.13, 4.62, 9.12),
  kind = c("inv", "inv", "npp", "npp", "npp"),
  stringsAsFactors = FALSE)

# --- candidate window ceilings ----------------------------------------------
ceil <- data.frame(
  label = c("CURRENT window (Gower global max)", "Gower Nordic max", "Zheng gridded mean"),
  x     = c(8.70, 4.62, 2.81),
  col   = c("grey35", "#0072B2", "#D55E00"), stringsAsFactors = FALSE)

png(OUT, width = 2100, height = 1500, res = 200)
op <- par(mar = c(5.5, 16, 5.2, 1.5), xaxs = "i")
XMAX <- 9.6
n <- nrow(lit) + nrow(ours) + 1
plot(NA, xlim = c(0, XMAX), ylim = c(0.3, n + 0.7), yaxt = "n", bty = "n",
     xlab = expression("litter input to soil / NPP  (tC ha"^-1*" yr"^-1*")"), ylab = "")

# ceilings first, behind everything
for (i in seq_len(nrow(ceil))) {
  abline(v = ceil$x[i], col = ceil$col[i], lwd = 2, lty = c(3, 2, 4)[i])
  text(ceil$x[i], n + 0.55, ceil$label[i], srt = 90, adj = c(1, -0.35),
       cex = 0.62, col = ceil$col[i], xpd = NA)
}

ypos <- n; labs <- c(); cols <- c()
draw <- function(y, lo, mid, hi, col, pch = 19) {
  if (!is.na(lo)) segments(lo, y, hi, y, col = col, lwd = 3, lend = 1)
  points(mid, y, pch = pch, col = col, cex = 1.4)
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
  labs <- c(labs, paste0("OURS: ", ours$model[i], "  (σ_input × J̄)"))
  cols <- c(cols, "#CC79A7"); ypos <- ypos - 1
}
axis(2, at = seq(n, by = -1, length.out = nrow(lit)), labels = labs[seq_len(nrow(lit))],
     las = 1, cex.axis = 0.72, tick = FALSE, line = -0.5)
axis(2, at = seq(n - nrow(lit) - 0.6, by = -1, length.out = nrow(ours)),
     labels = labs[-seq_len(nrow(lit))], las = 1, cex.axis = 0.72, tick = FALSE, line = -0.5)

title_x <- grconvertX(0, "ndc", "user")   # left device edge, so long titles are not clipped
mtext("Litter input to Finnish forest soil: inventory products vs independent NPP vs our calibration",
      side = 3, line = 3.6, at = title_x, adj = 0, cex = 0.92, font = 2, xpd = NA)
mtext("Bars are reported ranges (Gower, Zheng) or 95% intervals (L&H, ours).",
      side = 3, line = 2.4, at = title_x, adj = 0, cex = 0.6, col = "grey30", xpd = NA)
mtext("litter = NPP − biomass increment − stemwood removed − exudates.",
      side = 3, line = 0.8, at = title_x, adj = 0, cex = 0.6, col = "grey30", xpd = NA)
legend("bottomright", bty = "n", cex = 0.66,
       legend = c("inventory product (NFI biomass × turnover)", "independent NPP (upper bound)",
                  "our posterior effective flux"),
       col = c("#009E73", "#56B4E9", "#CC79A7"), pch = c(19, 19, 18), lwd = 3, lty = 1)
par(op); dev.off()
cat("wrote", OUT, "\n")
print(ours, row.names = FALSE)
