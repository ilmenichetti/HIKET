# =============================================================================
# build_S_saturation.R   (2026-10-26)
#
# S_saturation -- SATURATION, POSTPONED BUT NOT REMOVED (R3, side arc "inputs and
# saturation"; Conclusion 1).
# No new runs: the predictive stage already projects every draw twice from 2025
# (Calibration_real_data_transient/forward_scenarios.R):
#   total_soc      litter held at 2024, the 2005-2024 climate recycled, 60 years
#   total_soc_up   the same with the litter x HIKET_INPUT_SCEN_MULT (1.2) from 2025
#   C_eq           equilibrium stock at the 2024 litter and the recycled climate
# The models are linear in their input, so the equilibrium of the +20% arm is
# MULT x C_eq exactly.
#   (a-f) national mean stock 1985-2084 per model: history, held litter (solid),
#         +20% (dashed), with the two equilibria (dotted) and the campaign means;
#         text = distance below balance in 2024 and 2084, and the share of the
#         +20% gain realised by 2084
#   (g)   the sink over each WHOLE 20-yr climate cycle of the projection (2025-44,
#         2045-64, 2065-84). The recycled climate alternates by decade, so shorter
#         windows carry its phase (doublechecks/yasso20_late_source_check.R); a
#         trailing 20-yr window was tried first and REJECTED -- across the step in
#         input at 2025 it draws a false peak at 2044.
#   "headroom" = C_eq / C - 1, the definition of build_forward_scenarios.R and F19.
# Basis: balanced plot set, unweighted national mean; 90% over posterior draws.
#
# Run from repo root:  Rscript manuscript/figures/build_S_saturation.R
# =============================================================================
suppressMessages({
  source("manuscript/figures/run_ids.R"); source("manuscript/figures/obs_basis.R")
  source("manuscript/figures/model_palette.R")
})
cat("RUN_IDs:", paste(names(RID), RID, collapse = " | "), "\n")
MULT <- suppressWarnings(as.numeric(Sys.getenv("HIKET_INPUT_SCEN_MULT", "1.2")))
W    <- 20L                                              # sink window = one climate cycle
q90  <- function(M) apply(M, 2, function(x) c(med = median(x), lo = unname(quantile(x, .05)),
                                              hi = unname(quantile(x, .95))))

# === Extract national means per draw ===
Z <- list()
for (m in FIG_MODELS) {
  pk <- readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_predictive_%s.rds", m, RID[[m]]))
  om <- readRDS(sprintf("Data/model_inputs/%s_inputs_%s.rds", m, RID[[m]]))$obs_meta
  bal <- balanced_plots(om)
  pp <- pk$projection_predictions; pp <- pp[as.integer(pp$plot_id) %in% bal, ]
  ph <- pk$posterior_predictions;  ph <- ph[as.integer(ph$plot_id) %in% bal, ]
  if (!all(c("total_soc_up", "C_eq") %in% names(pp)))
    stop(m, ": projection_predictions lacks the side-analysis columns -- rerun the predictive stage")
  H <- tapply(ph$total_soc,    list(ph$draw, ph$year), mean)
  S <- tapply(pp$total_soc,    list(pp$draw, pp$year), mean)
  U <- tapply(pp$total_soc_up, list(pp$draw, pp$year), mean)
  first <- !duplicated(pp[c("draw", "plot_id")])
  CE <- tapply(pp$C_eq[first], pp$draw[first], mean)
  dr <- intersect(rownames(H), rownames(S)); dr <- intersect(dr, names(CE))
  held <- cbind(H[dr, ], S[dr, ]); up <- cbind(H[dr, ], U[dr, ])
  yr <- as.integer(colnames(held)); CE <- CE[dr]
  ends <- c(2044L, 2064L, 2084L)                         # whole climate cycles
  rate <- function(M) { R <- (M[, match(ends, yr)] - M[, match(ends - W, yr)]) / W
    colnames(R) <- ends; R }
  ann <- function(M) { py <- 2025:2084                     # annual sink, median over draws
    apply(M[, match(py, yr)] - M[, match(py - 1L, yr)], 2, median) }
  Z[[m]] <- list(yr = yr, held = held, up = up, CE = CE, Rh = rate(held), Ru = rate(up),
                 Ah = ann(held), Au = ann(up),
                 obs = obs_campaigns(om, bal), n = length(bal))
}

# === Numbers (printed, and used for the panel text) ===
num <- do.call(rbind, lapply(FIG_MODELS, function(m) {
  z <- Z[[m]]; at <- function(M, y) M[, match(y, z$yr)]
  data.frame(model = m,
    head_2024    = median(z$CE / at(z$held, 2024) - 1),
    head_2084    = median(z$CE / at(z$held, 2084) - 1),
    sink_25_44   = median(z$Rh[, "2044"]),
    sink_65_84   = median(z$Rh[, "2084"]),
    up_25_44     = median(z$Ru[, "2044"]),
    up_65_84     = median(z$Ru[, "2084"]),
    up_realised  = median((at(z$up, 2084) - at(z$held, 2084)) / ((MULT - 1) * z$CE)))
}))
print(num, digits = 3, row.names = FALSE)
cat("annual sink range (medians):", round(range(unlist(lapply(Z, function(z) c(z$Ah, z$Au)))), 3), "\n")

# === Figure ===
# Drawn at the size it is printed (0.98 x the 6.53 in text width), so the type is the
# type: pointsize 9, nothing below cex 0.8 (~7 pt). Title and basis note live in the
# caption, not in the figure.
png("manuscript/figures/S_saturation.png", width = 6.5, height = 8.4, units = "in", res = 300, pointsize = 9)
layout(rbind(c(1,1,2,2,3,3), c(4,4,5,5,6,6), c(7,7,7,8,8,8)), heights = c(1, 1, 1.15))
par(cex = 1)   # layout/mfrow shrink text to 0.66 (or 0.83 in 2x2); keep the 9 pt
par(oma = c(0, 0, 0, 0), mgp = c(2.0, 0.5, 0), las = 1, tcl = -0.3)
band <- function(x, Q, col) polygon(c(x, rev(x)), c(Q["lo", ], rev(Q["hi", ])),
                                    col = adjustcolor(col, 0.22), border = NA)

for (k in seq_along(FIG_MODELS)) {
  m <- FIG_MODELS[k]; z <- Z[[m]]; col <- MODEL_COL[[m]]
  Qh <- q90(z$held); Qu <- q90(z$up); ce <- median(z$CE)
  pj <- z$yr >= 2024
  par(mar = c(2.4, if (k %in% c(1, 4)) 3.4 else 2.4, 1.8, 0.4))
  plot(NA, xlim = c(1985, 2084), ylim = c(0, max(Qh, Qu[, pj], ce * MULT, z$obs$hi) * 1.12),
       xlab = "", ylab = if (k %in% c(1, 4)) expression("Stock  (tC ha"^-1*")") else "",
       main = sprintf("(%s) %s", letters[k], m), font.main = 1, cex.main = 1.1, yaxs = "i",
       xaxt = "n")
  axis(1, at = c(1990, 2030, 2070))
  rect(2024.5, -1e3, 2100, 1e3, col = "grey94", border = NA)
  abline(h = c(ce, ce * MULT), col = "grey30", lty = 3, lwd = 1.1)
  band(z$yr[pj], Qu[, pj], "grey40"); band(z$yr, Qh, col)
  lines(z$yr[pj], Qu["med", pj], col = "grey15", lwd = 1.4, lty = 2)
  lines(z$yr, Qh["med", ], col = col, lwd = 2)
  segments(z$obs$year, z$obs$lo, z$obs$year, z$obs$hi, col = CAMPAIGN_COL[z$obs$camp], lwd = 1.4)
  points(z$obs$year, z$obs$m, pch = 21, bg = CAMPAIGN_COL[z$obs$camp], col = "grey15", cex = 1.1)
  nk <- num[num$model == m, ]
  text(1986, c(ce * MULT, ce), c("balance, +20%", "balance"), adj = c(0, -0.3), cex = 0.8, col = "grey25")
  text(2083, 0, adj = c(1, -0.25), cex = 0.8, col = "grey15",
       sprintf("to gain: %.0f%% → %.0f%%\n+20%% realised: %.0f%%",
               round(100 * nk$head_2024) + 0, round(100 * nk$head_2084) + 0, 100 * nk$up_realised))
}

# (g)/(h): the sink as grouped bars -- one group per 20-year repetition of the climate
# (whole repetitions, because the repeated climate alternates by decade), one bar per model,
# 90% interval as whiskers. Tried before: dodged points joined by lines (read as a time
# series) and steps on a time axis (cluttered); bars say "one value per block" directly.
S90 <- lapply(Z, function(z) list(h = q90(z$Rh), u = q90(z$Ru)))
yl  <- range(-0.1, 0.55, unlist(lapply(S90, function(s) c(s$h[c("lo", "hi"), ], s$u[c("lo", "hi"), ]))))
BENCH <- c("SP1", "TP2", "TP3")
BAR_FILL   <- ifelse(FIG_MODELS %in% BENCH, adjustcolor(MODEL_COL[FIG_MODELS], 0.5), MODEL_COL[FIG_MODELS])
BAR_BORDER <- ifelse(FIG_MODELS %in% BENCH, MODEL_COL[FIG_MODELS], NA)
sink_panel <- function(a, tag, first) {
  M  <- sapply(FIG_MODELS, function(m) S90[[m]][[a]]["med", ])          # 3 blocks x 6 models
  LO <- sapply(FIG_MODELS, function(m) S90[[m]][[a]]["lo", ])
  HI <- sapply(FIG_MODELS, function(m) S90[[m]][[a]]["hi", ])
  # the simple benchmarks (SP1/TP2/TP3) at half opacity with a solid border, so they read as
  # benchmarks, not as alternative estimates of the Finnish sink
  bp <- barplot(t(M), beside = TRUE, col = BAR_FILL, border = BAR_BORDER, lwd = 1, ylim = yl,
                names.arg = c("2025-44", "2045-64", "2065-84"), space = c(0, 1.2),
                ylab = if (first) expression("Sink  (tC ha"^-1*" yr"^-1*")") else "")
  abline(h = 0, col = "grey35")
  abline(h = c(0.26, 0.12), col = "#AA3333", lty = 3)
  arrows(bp, t(LO), bp, t(HI), angle = 90, code = 3, length = 0.015, col = "grey25", lwd = 0.9)
  mtext(tag, 3, 0.4, adj = 0)
  bp
}
par(mar = c(2.6, 3.8, 2.0, 0.4))
bp <- sink_panel("h", "(g) Sink, litter held at 2024", TRUE)
text(max(bp) + 0.6, c(0.26, 0.12), c("observed 1985-2024", "observed 2006-2024"),
     adj = c(1, -0.4), cex = 0.8, col = "#AA3333")
par(mar = c(2.6, 2.4, 2.0, 0.6))
sink_panel("u", "(h) Sink, litter +20% from 2025", FALSE)
legend("topright", bty = "n", cex = 0.85, ncol = 2, legend = FIG_MODELS,
       fill = BAR_FILL, border = ifelse(is.na(BAR_BORDER), "white", BAR_BORDER))
dev.off()
cat("\nWrote manuscript/figures/S_saturation.png\n")
