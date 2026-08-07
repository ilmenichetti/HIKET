# =============================================================================
# run_ids.R   (2026-08-07)
#
# Single source of truth for which calibration the manuscript figures are built
# from. Every build_*.R sources this and interpolates RID[["<MODEL>"]] instead of
# hardcoding a RUN_ID.
#
# WHY THIS EXISTS. The builders used to hardcode the 20260710 RUN_IDs. After a
# re-calibration they did not error -- the old posteriors are still on disk, so
# they rebuilt happily from stale data and produced figures that looked current.
# The same failure was found the same day in two doublechecks scripts and four
# NextGenC scripts. Auto-selecting removes the whole class.
#
# Ablation posteriors are moved out of runs/ by
# doublechecks/quarantine_ablation_runs.R, so "newest" here is a production run.
# Run that first if an ablation suite has just finished.
# =============================================================================

FIG_MODELS <- c("SP1", "TP2", "TP3", "Yasso07", "Yasso15", "Yasso20")

.fig_runs   <- "Calibration_real_data_transient/runs"
.fig_inputs <- "Data/model_inputs"

RID <- vapply(FIG_MODELS, function(m) {
  fs <- list.files(.fig_runs,
                   pattern = sprintf("^%s_posterior_[0-9]{8}_[0-9]{6}\\.rds$", m))
  if (!length(fs))
    stop("run_ids.R: no posterior found for ", m, " in ", .fig_runs, call. = FALSE)
  sub(sprintf("^%s_posterior_(.+)\\.rds$", m), "\\1", sort(fs, decreasing = TRUE)[1])
}, character(1))

# Fail LOUDLY if the companion files a figure needs are absent, rather than
# letting a builder fall back to whatever else is lying around.
local({
  missing <- character(0)
  for (m in FIG_MODELS) {
    want <- c(file.path(.fig_runs,   sprintf("%s_posterior_predictive_%s.rds", m, RID[[m]])),
              file.path(.fig_inputs, sprintf("%s_inputs_%s.rds",              m, RID[[m]])))
    missing <- c(missing, want[!file.exists(want)])
  }
  if (length(missing))
    stop("run_ids.R: RUN_ID selected but companion files missing:\n  ",
         paste(missing, collapse = "\n  "),
         "\nRun stages 2-4 (run_hiket_pipeline.R --skip-calibration) first.",
         call. = FALSE)
})

message("Figures build from: ",
        paste(sprintf("%s=%s", FIG_MODELS, RID[FIG_MODELS]), collapse = "  "))
