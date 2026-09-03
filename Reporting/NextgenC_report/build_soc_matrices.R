# === NextGenC reporting: per-model SOC matrices ===
# For each of the 6 models, build two plot x year matrices from the
# posterior-predictive bundle:
#   - mean  : soc_mean (mean total SOC across MCMC posterior-predictive draws)
#   - sd    : soc_sd   (SD across draws = "deviation")
# Rows = simulated points (plots), columns = years.
# Output: 12 CSVs + one multi-sheet ODS workbook (12 sheets).

suppressPackageStartupMessages({
  library(tidyr)
  library(readODS)
})

# --- paths ---
PROJ   <- "/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling"
RUNS   <- file.path(PROJ, "Calibration_real_data_transient", "runs")
OUTDIR <- file.path(PROJ, "Reporting", "NextgenC_report")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# --- model -> predictive bundle (AUTO-SELECTED, never hardcoded) --------------
# Sources manuscript/figures/run_ids.R, the single source of truth for which
# calibration everything is built from. It picks the newest production posterior,
# FAILS LOUDLY if the companion predictive/input files are missing, and echoes the
# RUN_IDs it chose. Rewired 2026-09-03: these four NextGenC scripts had carried
# hardcoded 20260710_* bundles since July and were five production runs stale --
# they did not error, they rebuilt happily from old posteriors.
# Pin an older run for comparison with HIKET_FIG_RID (see run_ids.R).
local({ owd <- setwd(PROJ); on.exit(setwd(owd)); source("manuscript/figures/run_ids.R") })
bundles <- setNames(sprintf(
  "%s_posterior_predictive_%s.rds", FIG_MODELS, RID[FIG_MODELS]),
  FIG_MODELS)
# --- pivot posterior_summary -> plot x year matrix for one value column ---
make_matrix <- function(summary_df, value_col) {
  w <- pivot_wider(
    summary_df[, c("plot_id", "year", value_col)],
    names_from = "year", values_from = dplyr::all_of(value_col)
  )
  w <- w[order(suppressWarnings(as.numeric(w$plot_id)), w$plot_id), ]
  # order year columns numerically
  yr_cols <- setdiff(names(w), "plot_id")
  yr_cols <- yr_cols[order(as.numeric(yr_cols))]
  w[, c("plot_id", yr_cols)]
}

sheets <- list()
for (m in names(bundles)) {
  pp <- readRDS(file.path(RUNS, bundles[[m]]))
  s  <- pp$posterior_summary

  mean_mat <- make_matrix(s, "soc_mean")
  sd_mat   <- make_matrix(s, "soc_sd")

  write.csv(mean_mat, file.path(OUTDIR, sprintf("%s_SOC_mean.csv", m)), row.names = FALSE)
  write.csv(sd_mat,   file.path(OUTDIR, sprintf("%s_SOC_sd.csv",   m)), row.names = FALSE)

  sheets[[paste0(m, "_mean")]] <- as.data.frame(mean_mat)
  sheets[[paste0(m, "_sd")]]   <- as.data.frame(sd_mat)

  cat(sprintf("%-8s  plots=%d  years=%d (%s..%s)\n",
              m, nrow(mean_mat), ncol(mean_mat) - 1L,
              names(mean_mat)[2], names(mean_mat)[ncol(mean_mat)]))
}

# --- single ODS workbook, one sheet per matrix (12 sheets) ---
ods_path <- file.path(OUTDIR, "NextGenC_SOC_matrices.ods")
write_ods(sheets, path = ods_path)

cat("\nWrote 12 CSVs + ", basename(ods_path), " to:\n", OUTDIR, "\n", sep = "")
