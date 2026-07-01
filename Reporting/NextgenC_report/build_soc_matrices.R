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

# --- model -> predictive bundle (latest production RUN_IDs) ---
# TP3 is the EXACT-integrator re-calibration (RUN 20260630_090644, synced from
# Puhti 2026-07-01). Its mean-SOC trajectory now shows the genuine (physical)
# interannual oscillation of the fast active pool, not the former Euler ringing.
bundles <- c(
  SP1     = "SP1_posterior_predictive_20260608_015026.rds",
  TP2     = "TP2_posterior_predictive_20260608_020212.rds",
  TP3     = "TP3_posterior_predictive_20260630_090644.rds",
  Yasso07 = "Yasso07_posterior_predictive_20260611_032825.rds",
  Yasso15 = "Yasso15_posterior_predictive_20260611_032825.rds",
  Yasso20 = "Yasso20_posterior_predictive_20260611_033140.rds"
)

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
