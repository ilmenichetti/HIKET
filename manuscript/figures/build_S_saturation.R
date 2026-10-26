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
  Z[[m]] <- list(yr = yr, held = held, up = up, CE = CE, Rh = rate(held), Ru = rate(up),
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

# === Figure ===
png("manuscript/figures/S_saturation.png", width = 14, height = 7.2, units = "in", res = 200)
layout(rbind(c(1, 2, 3, 7), c(4, 5, 6, 7)), widths = c(1, 1, 1, 1.45))
par(oma = c(1.6, 0, 2.6, 0), mgp = c(2.3, 0.55, 0), las = 1, cex.axis = 0.85)
band <- function(x, Q, col) polygon(c(x, rev(x)), c(Q["lo", ], rev(Q["hi", ])),
                                    col = adjustcolor(col, 0.22), border = NA)

for (k in seq_along(FIG_MODELS)) {
  m <- FIG_MODELS[k]; z <- Z[[m]]; col <- MODEL_COL[[m]]
  Qh <- q90(z$held); Qu <- q90(z$up); ce <- median(z$CE)
  pj <- z$yr >= 2024
  par(mar = c(3.4, 3.9, 2.2, 0.6))
  plot(NA, xlim = c(1985, 2084), ylim = range(0, Qh, Qu[, pj], ce * MULT, z$obs$lo, z$obs$hi) * c(0.97, 1.08),
       xlab = "", ylab = if (k %in% c(1, 4)) expression("Stock  (tC ha"^-1*")") else "",
       main = sprintf("(%s) %s", letters[k], m), font.main = 1, yaxs = "i")
  rect(2024.5, -1e3, 2100, 1e3, col = "grey95", border = NA)
  abline(h = c(ce, ce * MULT), col = "grey30", lty = 3, lwd = 1.2)
  band(z$yr[pj], Qu[, pj], "grey40"); band(z$yr, Qh, col)
  lines(z$yr[pj], Qu["med", pj], col = "grey15", lwd = 1.8, lty = 2)
  lines(z$yr, Qh["med", ], col = col, lwd = 2.4)
  segments(z$obs$year, z$obs$lo, z$obs$year, z$obs$hi, col = CAMPAIGN_COL[z$obs$camp], lwd = 1.6)
  points(z$obs$year, z$obs$m, pch = 21, bg = CAMPAIGN_COL[z$obs$camp], col = "grey15", cex = 1.2)
  text(1987, ce, "balance, held", adj = c(0, -0.35), cex = 0.68, col = "grey25")
  text(1987, ce * MULT, sprintf("balance, +%d%%", round(100 * (MULT - 1))), adj = c(0, -0.35),
       cex = 0.68, col = "grey25")
  nk <- num[num$model == m, ]
  legend("bottomright", bty = "n", cex = 0.7, text.col = "grey20",
         legend = c(sprintf("to gain before balance: %.0f%% (2024) -> %.0f%% (2084)",
                            round(100 * nk$head_2024) + 0, round(100 * nk$head_2084) + 0),
                    sprintf("+%d%% gain realised by 2084: %.0f%%", round(100 * (MULT - 1)),
                            100 * nk$up_realised)))
}

par(mar = c(3.4, 4.4, 2.2, 1))
OFF <- setNames(seq(-2.5, 2.5, length.out = 6), FIG_MODELS)  # dodge, years
cx  <- c(2034.5, 2054.5, 2074.5)                              # cycle mid-years
S90 <- lapply(Z, function(z) list(h = q90(z$Rh), u = q90(z$Ru)))
yl  <- range(0, 0.27, unlist(lapply(S90, function(s) c(s$h[c("lo", "hi"), ], s$u[c("lo", "hi"), ]))))
plot(NA, xlim = c(2024, 2085), ylim = yl + c(-0.01, 0.02), xaxt = "n", xlab = "",
     ylab = expression("Sink over the climate cycle  (tC ha"^-1*" yr"^-1*")"),
     main = "(g) The sink fades with either input", font.main = 1)
axis(1, at = cx, labels = c("2025-44", "2045-64", "2065-84"), tick = FALSE)
abline(h = 0, col = "grey35", lty = 2)
abline(h = c(0.26, 0.12), col = "#AA3333", lty = 3)
text(2024, c(0.26, 0.12), c("observed 1985-2024", "observed 2006-2024"), adj = c(0, -0.4),
     cex = 0.7, col = "#AA3333")
for (m in FIG_MODELS) for (a in c("h", "u")) {
  Q <- S90[[m]][[a]]; x <- cx + OFF[[m]] + if (a == "u") 0.4 else -0.4
  lines(x, Q["med", ], col = MODEL_COL[[m]], lwd = 1.6, lty = if (a == "u") 2 else 1)
  segments(x, Q["lo", ], x, Q["hi", ], col = adjustcolor(MODEL_COL[[m]], 0.6), lwd = 1.4)
  points(x, Q["med", ], pch = 21, cex = 1.15, col = "grey20",
         bg = if (a == "h") MODEL_COL[[m]] else "white")
}
legend("topright", bty = "n", cex = 0.78, ncol = 2,
       legend = c(FIG_MODELS, "litter held at 2024", sprintf("litter +%d%% from 2025", round(100 * (MULT - 1)))),
       col = c(MODEL_COL[FIG_MODELS], "grey20", "grey20"), lty = c(rep(1, 6), 1, 2),
       pch = c(rep(NA, 6), 21, 21), pt.bg = c(rep(NA, 6), "grey50", "white"),
       lwd = c(rep(2.2, 6), 1.6, 1.6))

mtext("Saturation, postponed but not removed: litter held at 2024 against +20% from 2025 (2005-2024 climate recycled)",
      outer = TRUE, line = 0.7, cex = 1.0, font = 2)
mtext(sprintf("%s; bands: 90%% over posterior draws of the national mean; points: campaign means (95%% CI) at their mean sampling year; dotted: equilibrium stock of each arm",
              basis_note(seq_len(Z[[1]]$n))),
      side = 1, outer = TRUE, line = 0.4, cex = 0.62, col = "grey35", adj = 0.01)
dev.off()
cat("\nWrote manuscript/figures/S_saturation.png\n")
