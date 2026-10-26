# =============================================================================
# build_S15S17_attribution.R   (2026-10-08)
#
# ATTRIBUTION FIGURES (methodology: manuscript/appendices/appendix_attribution.tex;
# numbers: doublechecks/attribution_decomposition.R -> doublechecks/attribution.rds).
#
#   S15_attribution_by_model.png   six panels, one per model: the annual sink, 1985-2024,
#                                  split into historical carbon, litter inputs and climate
#                                  (stacked bars; positive up, negative down; the line is
#                                  the sink). Common y axis so the models compare.
#   S16_attribution_ensemble.png   the six-model average: the annual sink as above.
#   S17_attribution_cumulative.png the six-model average, cumulative: carbon gained since
#                                  1985 by each component (tC/ha); the line is the stock
#                                  change. HELD BACK for now (see the S17 section).
#   *_smooth10.png                 the same with a 10-yr centred running mean of the annual
#                                  values, pre-run and 1985 onwards smoothed SEPARATELY
#                                  with a window that shrinks at the edges (linear, so the
#                                  parts still add up exactly). S17 needs no smoothing.
#
# Components (all exact, see the appendix):
#   litter rise before 1985 (initialisation; "historical C")
#                  1917-1984: the pre-run sink (soil catching up with the 1917 -> 1985 rise
#                  of litter; climate held at the 1985-2004 mean). From 1985: the sink owed
#                  to the deficit carried into 1985 = production start minus equilibrium
#                  start under the same forcing. THE INITIALISATION STORY.
#   litter change since 1985   litter departing from its 1985-89 level.
#   climate        each year's climate departing from the 1985-2004 mean (one band, all
#                  of a model's modifiers together).
# Split: SEQUENTIAL, three runs (history = 1-2, climate = 2-3, litter = 3); see the appendix.
# Per model: POSTERIOR MEAN over the draws (means add up). Average: unweighted over models.
#
# Usage:  Rscript manuscript/figures/build_S15S17_attribution.R
#         HIKET_ATTR_RECAL=1 -> the RECALIBRATED decomposition (every scenario its own
#         calibration; doublechecks/attribution_recalibrated.rds); outputs carry "_recal".
# =============================================================================

setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
source("manuscript/figures/model_palette.R")
RECAL <- identical(Sys.getenv("HIKET_ATTR_RECAL"), "1")
Z <- readRDS(if (RECAL) "doublechecks/attribution_recalibrated.rds" else "doublechecks/attribution.rds")
if (!Z$version %in% c("level_v3_sequential", "recalibrated")) stop("unknown attribution version -- rerun the decomposition")
if (RECAL != identical(Z$version, "recalibrated")) stop("HIKET_ATTR_RECAL does not match the file read")
MS <- intersect(MODEL_ORDER, names(Z$models))
COMP     <- c("history", "inputs", "climate")
# Labels (Lorenzo, 2026-10-08): "historical C" IS an input effect -- the litter rise before
# 1985 and its delayed legacy -- so both input parts are green, split by era.
COMP_LAB <- c(history = "Inherited state (litter rise before 1985)", inputs = "Litter change since 1985",
              climate = "Climate")
COMP_COL <- c(history = "#2F6B3B", inputs = "#9CCB86", climate = "#C44E52")   # climate: red (Lorenzo)
HATCH <- c(history = FALSE, inputs = FALSE, climate = FALSE)
# JOINT VARIANT (2026-10-12, Lorenzo): HIKET_ATTR_JOINT=1 draws the two joint effects of
# climate as their own HATCHED bands, in the colour of the part they are assigned to in the
# sequential split (inheritance x climate -> history, litter x climate -> climate). The main
# effects are then measured at the reference corner: history = run 1' - 3 (reference
# climate), climate = run 4 - 5 (reference litter). Same total; outputs carry "_joint".
JOINT <- identical(Sys.getenv("HIKET_ATTR_JOINT"), "1")
if (JOINT) {
  COMP     <- c("history", "hist_x_clim", "inputs", "climate", "inputs_x_clim")
  COMP_LAB <- c(COMP_LAB[c("history")], hist_x_clim = "Litter rise before 1985 \u00d7 climate (joint)",
                COMP_LAB[c("inputs", "climate")], inputs_x_clim = "Litter change since 1985 \u00d7 climate (joint)")
  COMP_COL <- c(COMP_COL[c("history")], hist_x_clim = COMP_COL[["history"]], COMP_COL[c("inputs", "climate")],
                inputs_x_clim = COMP_COL[["climate"]])
  HATCH    <- setNames(COMP %in% c("hist_x_clim", "inputs_x_clim"), COMP)
}
SFX_J <- paste0(if (JOINT) "_joint" else "", if (RECAL) "_recal" else "")
Y0 <- 1985; Y_PROJ <- 2024

# Projection CUT (Lorenzo, 2026-10-08): the figures stop at 2024. The decomposition
# still computes 2025-2084 (attribution.rds); only the display is restricted.
# Display window CUT to the simulation window 1985-2024 (Lorenzo, 2026-10-09): every part is
# anchored at 1985, so the pre-run adds nothing the reading needs; its total is reported in
# the table only. The decomposition still holds 1917-2084 (attribution.rds).
Y_START <- Y0; Y_END <- Y_PROJ
mean_arr_full <- function(m) apply(Z$models[[m]]$arr, c(1, 2), mean)
mean_arr <- function(m) { a <- mean_arr_full(m); a <- a[a[, "year"] >= Y_START & a[, "year"] <= Y_END, , drop = FALSE]
  if (JOINT) {
    if (anyNA(a[, c("inter_hist_clim", "interaction")]))
      stop("five-term figure: the joint effects need arms S4 (_eqinit_litheld) and S1' (_climheld)")
    a <- cbind(a, hist_x_clim = a[, "inter_hist_clim"], inputs_x_clim = a[, "interaction"])
    a[, "history"] <- a[, "history"] - a[, "inter_hist_clim"]   # run 1' - 3 (reference climate)
    a[, "climate"] <- a[, "climate"] - a[, "interaction"]       # run 4 - 5 (reference litter)
    a <- a[, c("year", "F", COMP)]
  }
  a }   # year x var
yr <- mean_arr(MS[1])[, "year"]
P_all <- lapply(MS, function(m) mean_arr(m)[, COMP, drop = FALSE]); names(P_all) <- MS
F_all <- lapply(MS, function(m) mean_arr(m)[, "F"]);               names(F_all) <- MS
Pm <- Reduce(`+`, P_all) / length(MS); Fm <- Reduce(`+`, F_all) / length(MS)
stopifnot(max(abs(rowSums(Pm) - Fm)) < 1e-9)

# Running mean applied SEPARATELY to the pre-run and to 1985 onwards: a window straddling
# 1985 would leak post-1985 inputs/climate into pre-run years, where both are zero by
# construction.
sm <- function(P, W) {
  if (W == 1L) return(P)
  P <- as.matrix(P); out <- P
  for (k in list(yr < Y0, yr >= Y0)) if (any(k) && !anyNA(P[k, ]))
    out[k, ] <- apply(P[k, , drop = FALSE], 2, rmean, W = W)
  out
}
# Centred running mean that SHRINKS at the edges (no gaps): same weights for every
# component, so the smoothed parts still add up exactly to the smoothed sink.
rmean <- function(v, W) {
  n <- length(v); lo <- (W - 1L) %/% 2L; hi <- W - 1L - lo
  vapply(seq_len(n), function(i) mean(v[max(1L, i - lo):min(n, i + hi)]), numeric(1))
}
stack_bars <- function(x, P) {
  for (part in list(pmax(P, 0), pmin(P, 0))) {
    base <- rep(0, nrow(part))
    for (k in COMP) { top <- base + part[, k]
      if (HATCH[[k]]) {
        rect(x - 0.45, base, x + 0.45, top, col = adjustcolor(COMP_COL[[k]], 0.30), border = NA)
        rect(x - 0.45, base, x + 0.45, top, density = 60, angle = 45, col = COMP_COL[[k]], border = NA, lwd = 0.6)
      } else rect(x - 0.45, base, x + 0.45, top, col = COMP_COL[[k]], border = NA)
      base <- top }
  }
  lines(x, rowSums(P), col = "grey10", lwd = 1.5)
}
stack_area <- function(x, P) {
  for (part in list(pmax(P, 0), pmin(P, 0))) {
    base <- rep(0, nrow(part))
    for (k in COMP) { top <- base + part[, k]
      polygon(c(x, rev(x)), c(base, rev(top)), col = adjustcolor(COMP_COL[[k]], if (HATCH[[k]]) 0.30 else 0.85), border = NA)
      if (HATCH[[k]]) polygon(c(x, rev(x)), c(base, rev(top)), density = 25, angle = 45, col = COMP_COL[[k]], border = NA)
      base <- top }
  }
  lines(x, rowSums(P), col = "grey10", lwd = 2.2)
}
periods <- function(yl) abline(h = 0, col = "grey40")
period_labels <- function(yl, cex = 0.72) invisible(NULL)        # one period only (1985-2024)
legend_top <- function(ncol = 4) {                                # ncol = 2 on single-panel figures
  op <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE); on.exit(par(op))
  plot.new()
  legend("top", inset = 0.005, ncol = ncol, bty = "n", cex = 0.9, x.intersp = 0.5,
         legend = c(COMP_LAB[COMP], "total"), border = NA,
         fill = c(COMP_COL[COMP], NA), density = c(ifelse(HATCH[COMP], 40, -1), -1), angle = 45,
         col = c(rep(NA, length(COMP)), "grey10"), lty = c(rep(NA, length(COMP)), 1),
         lwd = c(rep(NA, length(COMP)), 1.8))
}
note <- if (RECAL) sprintf("RECALIBRATED: every scenario its own calibration; balanced plot set (n = %d), unweighted national mean; %d draws per arm (mean); runs %s",
                Z$models[[MS[1]]]$n_plots, Z$n_draw, paste(unique(na.omit(substr(unlist(Z$rid[MS]), 1, 8))), collapse = "/")) else
  sprintf("balanced plot set (n = %d), unweighted national mean; %d posterior draws per model (mean); production runs %s",
                Z$models[[MS[1]]]$n_plots, Z$n_draw, paste(unique(substr(Z$rid[MS], 1, 8)), collapse = "/"))
ylab_sink <- expression("Sink  (tC ha"^-1*" yr"^-1*")")

draw <- function(W) {
  sfx <- if (W == 1L) "" else sprintf("_smooth%d", W)
  wtxt <- if (W == 1L) "annual values" else sprintf("%d-yr centred running mean", W)
  # --- S15: six panels -------------------------------------------------------------
  PS <- lapply(P_all, sm, W = W); ok <- stats::complete.cases(PS[[1]])
  yl <- range(unlist(lapply(PS, function(P) c(rowSums(pmax(P[ok, ], 0)), rowSums(pmin(P[ok, ], 0))))))
  yl <- yl + c(0, 0.08) * diff(yl)
  # Print size (0.95 x 6.53 in, 9 pt); title and basis note live in the caption (fig:S_attr_models).
  png(sprintf("manuscript/figures/S15_attribution_by_model%s%s.png", sfx, SFX_J), width = 6.2, height = 4.8, units = "in", res = 300, pointsize = 9)
  par(mfrow = c(2, 3), mar = c(3.0, 3.6, 1.6, 0.4), oma = c(0, 0, 2.4, 0), mgp = c(2.0, 0.5, 0), las = 1, tcl = -0.3)
  for (m in MODEL_ORDER) {
    plot(NA, xlim = range(yr), ylim = yl, xlab = "Year", ylab = if (match(m, MODEL_ORDER) %in% c(1, 4)) ylab_sink else "")
    mtext(sprintf("(%s) %s", letters[match(m, MODEL_ORDER)], m), 3, 0.3, adj = 0, font = 2)
    if (!m %in% MS) { text(mean(range(yr)), 0, "not computed", col = "grey40"); next }
    periods(yl); stack_bars(yr[ok], PS[[m]][ok, , drop = FALSE])
    if (m == MS[1]) period_labels(yl, 0.62)
  }
  message("S15 ", wtxt, "; ", note)
  legend_top(); dev.off()
  # --- S16: average, annual sink ----------------------------------------------------
  Ps <- sm(Pm, W); okm <- stats::complete.cases(Ps)
  # Print size (0.85 x 6.53 in, 9 pt); title and basis note live in the caption (fig:S16).
  png(sprintf("manuscript/figures/S16_attribution_ensemble%s%s.png", sfx, SFX_J), width = 5.55, height = 3.6, units = "in", res = 300, pointsize = 9)
  par(mar = c(3.2, 3.8, 0.6, 0.6), oma = c(0, 0, 3.0, 0), mgp = c(2.2, 0.5, 0), las = 1, tcl = -0.3)
  yl1 <- range(c(rowSums(pmax(Ps[okm, ], 0)), rowSums(pmin(Ps[okm, ], 0)))); yl1 <- yl1 + c(0, 0.1) * diff(yl1)
  plot(NA, xlim = range(yr), ylim = yl1, xlab = "Year", ylab = ylab_sink, main = "")
  periods(yl1); stack_bars(yr[okm], Ps[okm, , drop = FALSE]); period_labels(yl1)
  legend_top(2); dev.off()
  cat(sprintf("Wrote S15_attribution_by_model%s.png and S16_attribution_ensemble%s.png\n", sfx, sfx))
}
for (W in c(1L, 10L)) draw(W)

# --- S17: average, carbon gained since 1917 (cumulative; HELD BACK for now) ----------------
# Held back (Lorenzo, 2026-10-08): the climate band keeps growing through the projection,
# which readers will take for future warming. The projection has NO climate change (it
# repeats 2005-2024): that band is the delayed response to warming already observed,
# fading as the soils settle. Labelled as such on the figure.
Cm <- apply(Pm, 2, cumsum)
png(sprintf("manuscript/figures/S17_attribution_cumulative%s.png", SFX_J), width = 9, height = 5.6, units = "in", res = 200)
par(mar = c(4, 4.6, 2.4, 1), oma = c(1.6, 0, 5, 0), mgp = c(2.5, 0.6, 0), las = 1)
yl2 <- range(c(rowSums(pmax(Cm, 0)), rowSums(pmin(Cm, 0)))); yl2 <- yl2 + c(0, 0.1) * diff(yl2)
plot(NA, xlim = range(yr), ylim = yl2, xlab = "Year", main = "Carbon gained since 1985",
     ylab = expression("Stock change since 1985  (tC ha"^-1*")"))
periods(yl2); stack_area(yr, Cm); period_labels(yl2)
if (FALSE) text(Y_PROJ + 2, yl2[2] - 0.065 * diff(yl2), "climate: response to warming already\nobserved, no further warming assumed",
     adj = c(0, 1), cex = 0.66, col = "#8E2F33", font = 3)
mtext(sprintf("%d-model average (%s)", length(MS), paste(MS, collapse = ", ")), outer = TRUE, line = 3.4, cex = 1.05, font = 2)
mtext(note, side = 1, outer = TRUE, line = 0.4, cex = 0.55, col = "grey35", adj = 0.01)
legend_top(2); dev.off()
cat("Wrote S17_attribution_cumulative.png\n")

# --- numbers: carbon gained by component over each period (tC/ha) ------------------------
# 1985-2024 by component, plus the pre-run gain 1917-1984 (not shown in the figures)
pre <- sapply(MS, function(m) { a <- mean_arr_full(m); sum(a[a[, "year"] < Y0, "F"]) })
tab <- do.call(rbind, lapply(c(MS, "average"), function(m) {
  P <- if (m == "average") Pm else P_all[[m]]
  data.frame(model = m, period = "1985-2024", t(round(colSums(P), 2)), total = round(sum(P), 2),
             prerun_gain_1917_1984 = round(if (m == "average") mean(pre) else pre[[m]], 2))
}))
cat("\nCarbon gained by component (tC/ha):\n"); print(tab, row.names = FALSE)
write.csv(tab, sprintf("manuscript/figures/S17_attribution_by_period%s.csv", SFX_J), row.names = FALSE)
