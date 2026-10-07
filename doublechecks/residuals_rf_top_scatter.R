# =============================================================================
# residuals_rf_top_scatter.R — scatterplots of the top predictors of the
# extended residual RF (residuals_rf_everything.R)
#
# WHY. Permutation importance ranks one signal many times (1985 basal area,
# volume and five biomass components are the same stand). Plotting the top six
# by rank would show six views of one thing. Instead walk down the ranking and
# skip any predictor correlated |Spearman| > RHO_MAX with one already chosen,
# so each panel is a distinct signal.
#
# WHAT. Target = mean plot residual of the RF's target models (default
# Yasso07 + Yasso15). Numeric predictors: points + loess of the target, thin
# loess per target model. Coded predictors: class means (n >= 10) over points.
#
# Run from repo root, after residuals_rf_everything.R:
#   Rscript doublechecks/residuals_rf_top_scatter.R
#   (HIKET_RF_MODELS as in the RF script selects which .rds to read)
# =============================================================================

OUT <- "doublechecks/figures"
TGT <- strsplit(Sys.getenv("HIKET_RF_MODELS", "Yasso07,Yasso15"), ",")[[1]]
TAG <- if (identical(TGT, "all")) "six" else paste(TGT, collapse = "_")
f   <- file.path(OUT, sprintf("residuals_rf_everything_%s.rds", TAG))
if (!file.exists(f)) stop("Run residuals_rf_everything.R first: ", f, call. = FALSE)
B   <- readRDS(f)
source("doublechecks/residuals_rf_common.R")
cat("Target:", paste(B$target_models, collapse = " + "),
    " RUN_IDs:", paste(B$run_ids, collapse = " "), sprintf(" OOB R² %.3f\n", B$r2))

N_PANEL <- 6; RHO_MAX <- 0.7
INK <- "#3d3d3a"; GRID <- "#e6e5df"; PT <- adjustcolor("#8a8983", 0.5)
MCOL <- c(Yasso07 = "#eda100", Yasso15 = "#e87ba4", SP1 = "#2a78d6",
          TP2 = "#eb6834", TP3 = "#1baf7a", Yasso20 = "#008300")



# --- choose distinct predictors ----------------------------------------------
X <- B$X; y <- X$resid
num_of <- function(x) if (is.numeric(x)) x else suppressWarnings(as.numeric(as.character(x)))
cand <- sub("_missing$", "", names(B$perm))
cand <- unique(ifelse(cand %in% names(PREFER), PREFER[cand], cand))
cand <- cand[cand %in% names(X) & cand != "noise.random"]
pick <- character(0)
for (v in cand) {
  xv <- num_of(X[[v]])
  dup <- any(vapply(pick, function(p) {
    ok <- is.finite(xv) & is.finite(num_of(X[[p]]))
    sum(ok) > 30 && abs(cor(xv[ok], num_of(X[[p]])[ok], method = "spearman")) > RHO_MAX
  }, logical(1)))
  if (!dup) pick <- c(pick, v)
  if (length(pick) == N_PANEL) break
}
cat("Panels (permutation rank):\n")
for (v in pick) cat(sprintf("  %-28s family rank %3d  %s\n", v, match(v, cand), lab(v)))

# --- draw ----------------------------------------------------------------------
loess_line <- function(x, r, col, lwd) {
  ok <- is.finite(x) & is.finite(r)
  lo <- loess(r ~ x, data.frame(x = x[ok], r = r[ok]), span = 0.75)
  xs <- seq(quantile(x[ok], .02), quantile(x[ok], .98), length.out = 100)
  lines(xs, predict(lo, data.frame(x = xs)), col = col, lwd = lwd)
}
YL <- quantile(y, c(.005, .995), na.rm = TRUE)
set.seed(2025)
png(file.path(OUT, sprintf("residuals_rf_top_scatter_%s.png", TAG)), 13, 8.6, units = "in", res = 150)
par(mfrow = c(2, 3), mar = c(7.4, 4, 2.4, 1), oma = c(0, 0, 2.2, 0), mgp = c(2.2, 0.6, 0),
    tcl = -0.3, col.axis = INK, col.lab = INK, fg = INK)
stats <- NULL
for (v in pick) {
  xr <- X[[v]]; xn <- num_of(xr)
  discrete <- !is.numeric(xr) || length(unique(na.omit(xn))) <= 10
  ok <- is.finite(y) & (if (discrete) !is.na(xr) & as.character(xr) != "NA" else is.finite(xn))
  if (discrete) {
    g  <- factor(as.character(xr[ok])); gi <- as.integer(g)
    r2 <- summary(lm(y[ok] ~ g))$r.squared
    plot(NA, xlim = c(0.5, nlevels(g) + 0.5), ylim = YL, xaxt = "n", bty = "l",
         xlab = "", ylab = "Log residual (obs/pred), plot mean")
    mtext(lab(v), 1, if (v %in% names(CODE_LABELS)) 6.3 else 2.2, cex = 0.7, col = INK)
    coded <- v %in% names(CODE_LABELS)
    tl <- if (coded) CODE_LABELS[[v]][levels(g)] else levels(g)
    tl[is.na(tl)] <- levels(g)[is.na(tl)]
    axis(1, seq_len(nlevels(g)), tl, las = if (coded) 2 else 1,
         cex.axis = if (coded) 0.75 else 1)
    abline(h = 0, col = GRID, lwd = 2)
    points(gi + runif(sum(ok), -.25, .25), y[ok], pch = 16, cex = 0.6, col = PT)
    m <- tapply(y[ok], g, mean); n <- tapply(y[ok], g, length); k <- n >= 10
    lines(which(k), m[k], type = "o", pch = 16, lwd = 2.5, col = INK)
    sub_txt <- sprintf("R² %.2f   n %d   (coded; means for classes n ≥ 10)", r2, sum(ok))
    rho <- NA
  } else {
    rho <- cor(xn[ok], y[ok], method = "spearman")
    r2  <- summary(lm(y[ok] ~ xn[ok]))$r.squared
    plot(NA, xlim = quantile(xn[ok], c(0, 1)), ylim = YL, bty = "l",
         xlab = lab(v), ylab = "Log residual (obs/pred), plot mean")
    abline(h = 0, col = GRID, lwd = 2)
    points(xn[ok], y[ok], pch = 16, cex = 0.6, col = PT)
    for (m in B$target_models) {
      rm <- B$per_model[[m]]
      loess_line(xn, rm, adjustcolor(MCOL[[m]], 0.9), 1.5)
    }
    loess_line(xn, y, INK, 2.5)
    sub_txt <- sprintf("Spearman ρ %.2f   R² %.2f   n %d", rho, r2, sum(ok))
  }
  mtext(sub_txt, 3, 0.3, adj = 0, cex = 0.72)
  stats <- rbind(stats, data.frame(variable = v, label = lab(v), rho = round(rho, 3),
                                   R2 = round(r2, 3), n = sum(ok)))
}
legend("topright", c(sprintf("mean of %s", paste(B$target_models, collapse = " + ")),
                     B$target_models), col = c(INK, MCOL[B$target_models]),
       lwd = c(2.5, rep(1.5, length(B$target_models))), bty = "n", cex = 0.85)
mtext(sprintf("Top distinct predictors of the residual (extended RF, OOB R² %.2f; |ρ| > %.1f with a higher-ranked predictor = skipped)",
              B$r2, RHO_MAX), 3, 0.6, outer = TRUE, cex = 0.8, col = INK)
dev.off()
print(stats, row.names = FALSE)
cat("\nWrote", file.path(OUT, sprintf("residuals_rf_top_scatter_%s.png", TAG)), "\n")
