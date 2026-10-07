# =============================================================================
# build_F11_alt_rf_heatmap.R — ALTERNATIVE F11: residual-predictor importance
# from the EXTENDED predictor set (residuals_rf_everything.R)
#
# WHY. F11 (manuscript/figures/build_F11_rf_heatmap.R) uses the pipeline RF:
# plot-years, calibration set, 57 covariates, the 1985 stand only. This version
# uses ~340 predictors from all four inventories, one mean residual per plot,
# all plots. ⚠ F11 itself is ON HOLD until the coauthors agree what it should
# show; this file only writes to doublechecks/figures/.
#
# WHAT. Columns: each model's own plot-mean residual, plus the mean of Yasso07
# and Yasso15 (the reference target of the residual study, outlined).
# Importance: permutation, 5 seeds. Near-identical predictors are one FAMILY
# (PREFER in residuals_rf_common.R; "_missing" flags join their variable);
# family importance = max over members, since correlated members share the
# permutation importance between them. Rows = top N_ROWS distinct families of
# the reference column, skipping any |Spearman| > RHO_MAX with a higher row.
# Relative importance within column (max over ALL families = 1), as in F11.
#
# Run from repo root, after residuals_rf_everything.R:
#   Rscript doublechecks/build_F11_alt_rf_heatmap.R
# Importances are cached in doublechecks/figures/F11_alt_importance.rds
# (keyed on the RUN_IDs); delete it to recompute.
# =============================================================================

suppressPackageStartupMessages({library(dplyr); library(ranger)})
source("manuscript/figures/run_ids.R")
source("doublechecks/residuals_rf_common.R")
OUT <- "doublechecks/figures"
SEEDS <- 2025 + 0:4; N_ROWS <- 15; RHO_MAX <- 0.7
REF <- c("Yasso07", "Yasso15"); REF_COL <- "Yasso07+15"

B <- readRDS(file.path(OUT, "residuals_rf_everything_Yasso07_Yasso15.rds"))
X <- B$X; X$resid <- NULL

# --- per-plot residual of each column -----------------------------------------
res <- sapply(FIG_MODELS, function(m) {
  d <- as.data.frame(readRDS(sprintf(
    "Calibration_real_data_transient/runs/%s_posterior_predictive_%s.rds",
    m, RID[[m]]))$residuals_df)
  tapply(d$residual_log, d$plot_id, mean)[as.character(B$plot_id)]
})
res <- cbind(res, rowMeans(res[, REF]))
colnames(res)[ncol(res)] <- REF_COL
COLS <- colnames(res)

# --- importance per column (cached) -------------------------------------------
cache <- file.path(OUT, "F11_alt_importance.rds")
key   <- paste(RID[FIG_MODELS], collapse = "|")
C <- if (file.exists(cache)) readRDS(cache) else NULL
if (is.null(C) || !identical(C$key, key) || !setequal(names(C$imp), COLS)) {
  C <- list(key = key, imp = list(), r2 = c())
  for (k in COLS) {
    D <- prep(cbind(resid = res[, k], X)); D <- D[is.finite(D$resid), ]
    fits <- lapply(SEEDS, function(s)
      ranger(resid ~ ., D, num.trees = 1000, importance = "permutation", seed = s,
             respect.unordered.factors = "order"))
    C$imp[[k]] <- rowMeans(sapply(fits, `[[`, "variable.importance"))
    C$r2[k]    <- mean(sapply(fits, `[[`, "r.squared"))
    cat(sprintf("%-11s OOB R² %.3f\n", k, C$r2[k]))
  }
  saveRDS(C, cache)
} else cat("Using cached importances (", cache, ")\n")

# --- collapse to families -----------------------------------------------------
fam_of <- function(v) { b <- sub("_missing$", "", v); ifelse(b %in% names(PREFER), PREFER[b], b) }
FAMS <- sort(unique(unlist(lapply(C$imp, function(x) fam_of(names(x))))))
fam <- sapply(COLS, function(k) {
  i <- C$imp[[k]]; tapply(i, fam_of(names(i)), max)[FAMS]
})
rownames(fam) <- FAMS
fam[is.na(fam)] <- 0
fam <- fam[rownames(fam) != "noise.random", ]
rel <- sweep(pmax(fam, 0), 2, apply(fam, 2, max), "/")

# --- distinct rows from the reference column ----------------------------------
num_of <- function(x) if (is.numeric(x)) x else suppressWarnings(as.numeric(as.character(x)))
ord <- rownames(fam)[order(-fam[, REF_COL])]
rows <- character(0)
for (v in ord) {
  if (!v %in% names(X)) next
  xv <- num_of(X[[v]])
  dup <- any(vapply(rows, function(p) {
    xp <- num_of(X[[p]]); ok <- is.finite(xv) & is.finite(xp)
    sum(ok) > 30 && abs(cor(xv[ok], xp[ok], method = "spearman")) > RHO_MAX
  }, logical(1)))
  if (!dup) rows <- c(rows, v)
  if (length(rows) == N_ROWS) break
}
mat <- rel[rev(rows), , drop = FALSE]        # most important at the top
write.csv(data.frame(variable = rows, label = sapply(rows, lab), round(rel[rows, ], 3),
                     check.names = FALSE),
          file.path(OUT, "F11_alt_rf_heatmap.csv"), row.names = FALSE)

# --- render (F11 style) -------------------------------------------------------
np <- nrow(mat); nm <- ncol(mat)
pal <- colorRampPalette(c("#f7fbff", "#deebf7", "#9ecae1", "#4292c6", "#08306b"))(100)
collab <- sprintf("%s\nR2=%.2f", COLS, C$r2[COLS])
png(file.path(OUT, "F11_alt_rf_heatmap.png"), width = 10.5, height = max(6, np * 0.34 + 2.4),
    units = "in", res = 200)
layout(matrix(c(1, 2), nrow = 1), widths = c(1, 0.15))
par(mar = c(4.4, 17, 4.4, 0.6), mgp = c(2.4, 0.6, 0))
image(1:nm, 1:np, t(mat), col = pal, zlim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
axis(2, at = 1:np, labels = sapply(rownames(mat), lab), las = 1, tick = FALSE, cex.axis = 0.84)
text(x = 1:nm, y = np + 0.62, labels = collab, xpd = NA, cex = 0.8, adj = c(0.5, 0),
     font = ifelse(COLS == REF_COL, 2, 1))
abline(h = seq(0.5, np + 0.5, 1), v = seq(0.5, nm + 0.5, 1), col = "white", lwd = 1.2)
rect(nm - 0.5, 0.5, nm + 0.5, np + 0.5, border = "#3d3d3a", lwd = 2.2, xpd = NA)
mtext("RF residual-predictor importance, extended predictors (relative, within-column max = 1)",
      side = 1, line = 1.4, cex = 0.9, font = 2)
mtext(sprintf("One residual per plot (n %d). Near-identical predictors merged (1985 basal area = basal area, volume, biomass). Outlined: mean of Yasso07 and Yasso15.",
              nrow(X)), side = 1, line = 2.7, cex = 0.62)
par(mar = c(4.4, 0.5, 3.6, 3.2))
leg <- seq(0, 1, length.out = 100)
image(1, leg, matrix(leg, nrow = 1), col = pal, axes = FALSE, xlab = "", ylab = "")
axis(4, at = seq(0, 1, 0.25), las = 1, cex.axis = 0.8); mtext("relative importance", side = 4, line = 2.0, cex = 0.8)
dev.off()
cat("Rows:\n"); print(data.frame(label = sapply(rows, lab), round(rel[rows, ], 2)), right = FALSE)
cat("\nWrote", file.path(OUT, "F11_alt_rf_heatmap.png"), "\n")
