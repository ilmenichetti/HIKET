# =============================================================================
# quarantine_ablation_runs.R   (2026-08-04)
#
# WHY THIS EXISTS.
# run_*_transient_predictive.R selects its posterior by listing
# Calibration_real_data_transient/runs/ and taking the LAST filename after a
# decreasing sort -- i.e. the newest RUN_ID wins, with no notion of what kind of
# run produced it:
#
#     rns    <- sort(rns, decreasing = TRUE)
#     RUN_ID <- sub("^TP2_posterior_(.+)\\.rds$", "\\1", rns[1])
#
# The ablation suite writes its posteriors into that same directory. Those runs
# use 3 short chains and deliberately non-production settings (C5 off, C3 off,
# ...), so if one of them ever sorts newest, the whole downstream stage would
# silently analyse an ablation as though it were the production calibration.
#
# This script moves ablation posteriors (and their paired input bundles) out of
# the production directories into doublechecks/ablation_runs/, so runs/ holds
# ONLY production calibrations. Nothing is deleted; the ablation index CSVs are
# rewritten to point at the new location, and summarise_ablation.R looks here
# first.
#
# Run AFTER an ablation suite completes and BEFORE any downstream/predictive
# stage or any rsync of runs/ to or from the cluster.
#
# Usage:  Rscript doublechecks/quarantine_ablation_runs.R
# =============================================================================

DIR_RUNS  <- "Calibration_real_data_transient/runs"
DIR_INP   <- "Data/model_inputs"
DIR_QUAR  <- "doublechecks/ablation_runs"
DIR_IDX   <- "doublechecks/ablation_logs"

dir.create(DIR_QUAR, showWarnings = FALSE, recursive = TRUE)

idx_files <- list.files(DIR_IDX, pattern = "_ablation_index\\.csv$", full.names = TRUE)
if (!length(idx_files)) {
  message("No ablation index files in ", DIR_IDX, " -- nothing to quarantine.")
  quit(status = 0L)
}

moved <- 0L; missing <- 0L
for (f in idx_files) {
  idx   <- read.csv(f, stringsAsFactors = FALSE)
  MODEL <- sub("_ablation_index\\.csv$", "", basename(f))
  for (i in seq_len(nrow(idx))) {
    rid <- idx$run_id[i]
    if (is.na(rid) || !nzchar(rid)) next
    # posterior + its paired input bundle (predictive keys the bundle by RUN_ID)
    src <- c(file.path(DIR_RUNS, sprintf("%s_posterior_%s.rds", MODEL, rid)),
             file.path(DIR_RUNS, sprintf("%s_chains_%s.rds",    MODEL, rid)),
             file.path(DIR_INP,  sprintf("%s_inputs_%s.rds",    MODEL, rid)))
    for (s in src) {
      if (!file.exists(s)) { missing <- missing + 1L; next }
      dst <- file.path(DIR_QUAR, basename(s))
      if (file.rename(s, dst)) {
        moved <- moved + 1L
        message(sprintf("  moved %-46s -> %s", basename(s), DIR_QUAR))
      }
    }
  }
  idx$quarantined_to <- DIR_QUAR
  write.csv(idx, f, row.names = FALSE)
}

cat(sprintf("\nQuarantined %d files into %s (%d expected files absent).\n",
            moved, DIR_QUAR, missing))

# --- report what remains in runs/, so production runs stay visible -----------
cat("\nProduction posteriors still in runs/ (newest first, per model):\n")
for (M in c("SP1","TP2","TP3","Yasso07","Yasso15","Yasso20")) {
  r <- list.files(DIR_RUNS, pattern = sprintf("^%s_posterior_[0-9]{8}_[0-9]{6}\\.rds$", M))
  r <- sort(r, decreasing = TRUE)
  cat(sprintf("  %-9s %s\n", M,
              if (length(r)) paste(head(sub(".*_posterior_|\\.rds", "", r), 3), collapse = "  ")
              else "(none)"))
}
cat("\nThe first RUN_ID on each line is what run_*_predictive.R will auto-select.\n")
