# =============================================================================
# run_hiket_pipeline.R
#
# Top-level orchestrator for the HIKET Bayesian SOC model intercomparison.
# Runs the full pipeline for all three models in sequence:
#   Calibration -> Predictive -> Residual analysis -> Multi-model comparison
#
# USAGE:
#   Rscript run_hiket_pipeline.R [--skip-calibration]
#
#   --skip-calibration  : skip the (slow) MCMC stage and use existing
#                         posterior files. Useful for re-running downstream
#                         steps after a completed calibration.
#
# DESIGN:
#   Each stage runs in its own clean R process via Rscript. This gives
#   environment isolation (no object shadowing between models) and means any
#   stage can be re-run independently if it fails, without re-running the
#   expensive calibration. RUN_IDs are auto-detected from the filesystem
#   after each calibration completes.
#
# =============================================================================
# PUHTI WORKFLOW  (command by command)
# =============================================================================
#
# --- ONE-TIME SETUP (only needed after first clone or after git pull) --------
#
#   [on Puhti, from /scratch/project_2019134/HIKET/]
#
#   # Recompile Fortran .so files (Mac binaries won't run on Linux)
#   cd Model_functions_real_data/Decomposition_functions/Yasso/
#   rm -f yasso07.so yasso07.o yasso07_mod.mod yasso15.so yasso15.o yasso15_mod.mod
#   apptainer_wrapper exec R CMD SHLIB yasso07.f90 yasso15.f90
#
#   cd ../RothC/
#   rm -f rothc_step_f.so rothc_step_f.o
#   apptainer_wrapper exec R CMD SHLIB rothc_step_f.f90
#
#   cd /scratch/project_2019134/HIKET/
#
# --- STAGE 1: CALIBRATION (one SLURM job per model) -------------------------
#
#   [on Puhti, from /scratch/project_2019134/HIKET/]
#
#   # Submit all three calibration jobs simultaneously
#   sbatch Calibration_real_data/hiket_sp1.sh
#   sbatch Calibration_real_data/hiket_tp2.sh
#   sbatch Calibration_real_data/hiket_yasso07.sh
#   sbatch Calibration_real_data/hiket_yasso15.sh
#   sbatch Calibration_real_data/hiket_yasso20.sh
#
#   # Monitor job queue
#   squeue -u menichet
#
#   # Monitor chain progress (replace RUNID with the timestamp from the .out log)
#   tail -f Calibration_real_data/progress_logs/Yasso07_chain_01_<RUNID>.log
#
#   # Check output and error logs (replace JOBID with the SLURM job number)
#   cat Calibration_real_data/progress_logs/yasso07_<JOBID>.out
#   cat Calibration_real_data/progress_logs/yasso07_<JOBID>.err
#
# --- STAGE 2-4: PREDICTIVE + RESIDUALS + COMPARISON -------------------------
#
#   [on Puhti, from /scratch/project_2019134/HIKET/, after all calibrations done]
#
#   # Run downstream pipeline (fast enough for a login node, ~10-30 min total)
#   module load r-env
#   apptainer_wrapper exec Rscript --no-save \
#     Calibration_real_data/run_hiket_pipeline.R --skip-calibration
#
# --- SYNC RESULTS BACK TO MAC -----------------------------------------------
#
#   [on Mac, from your local HIKET folder]
#
#   # Pull posterior RDS files
#   rsync -av menichet@puhti.csc.fi:/scratch/project_2019134/HIKET/Calibration_real_data/runs/ \
#     ./Calibration_real_data/runs/
#
#   # Pull diagnostic PNGs and CSVs
#   rsync -av menichet@puhti.csc.fi:/scratch/project_2019134/HIKET/Calibration_real_data/diagnostics/ \
#     ./Calibration_real_data/diagnostics/
#
# --- SYNC CODE UPDATES TO PUHTI ---------------------------------------------
#
#   [on Mac]   git add ... && git commit -m "..." && git push
#   [on Puhti] cd /scratch/project_2019134/HIKET/ && git pull
#              # Then recompile .so if any .f90 files changed (see above)
#
# =============================================================================
#
# OUTPUTS:
#   All outputs are written to the same directories as if the scripts were
#   run individually. This script adds no new output files of its own.
# =============================================================================

# =============================================================================
# 0.  Configuration
# =============================================================================

SCRIPT_DIR <- "./Calibration_real_data"
RUNS_DIR   <- "./Calibration_real_data/runs"

MODELS <- c("SP1", "TP2", "Yasso07", "Yasso15", "Yasso20")

# Parse --skip-calibration flag
args <- commandArgs(trailingOnly = TRUE)
SKIP_CALIBRATION <- "--skip-calibration" %in% args

# Helper: run a script as a subprocess and stop the pipeline if it fails.
# stdout/stderr are inherited so all progress messages appear in the console.
run_stage <- function(script, args_str = "", label = NULL) {
  cmd   <- sprintf("Rscript %s/%s %s", SCRIPT_DIR, script, args_str)
  label <- if (!is.null(label)) label else script
  message(sprintf("\n>>> %s", label))
  message(sprintf("    %s", cmd))
  t_start <- proc.time()["elapsed"]
  status  <- system(cmd)
  elapsed <- proc.time()["elapsed"] - t_start
  if (status != 0)
    stop(sprintf("FAILED: %s  (exit status %d)", label, status))
  message(sprintf("    OK  (%.1f min)", elapsed / 60))
  invisible(status)
}

# Helper: auto-detect the most recently modified posterior RDS for a model.
detect_run_id <- function(model) {
  pat <- sprintf("^%s_posterior_[0-9]{8}_[0-9]{6}\\.rds$", model)
  fns <- list.files(RUNS_DIR, pattern = pat, full.names = TRUE)
  if (length(fns) == 0)
    stop(sprintf("No posterior found for %s in %s", model, RUNS_DIR))
  fn <- fns[order(file.info(fns)$mtime, decreasing = TRUE)[1]]
  sub(sprintf("^.*%s_posterior_(.+)\\.rds$", model), "\\1", fn)
}

# =============================================================================
# 1.  Calibration (skipped with --skip-calibration)
# =============================================================================

if (!SKIP_CALIBRATION) {
  message("\n=============================================================")
  message("  STAGE 1: Calibration (MCMC)")
  message("=============================================================")
  message("NOTE: Running all three models sequentially. On Mac this will")
  message("  take ~1.5-3 hours. On Roihu, use --skip-calibration and")
  message("  submit calibration scripts as separate SLURM jobs instead.")
  
  for (m in MODELS) {
    run_stage(sprintf("run_%s_calibration.R", m),
              label = sprintf("%s calibration", m))
  }
} else {
  message("\n[Calibration skipped -- using existing posteriors]\n")
}

# =============================================================================
# 2.  Detect RUN_IDs
# =============================================================================

message("\n=============================================================")
message("  Detecting RUN_IDs from filesystem...")
message("=============================================================")

run_ids <- setNames(lapply(MODELS, detect_run_id), MODELS)

for (m in MODELS)
  message(sprintf("  %-8s  %s", m, run_ids[[m]]))

# =============================================================================
# 3.  Posterior predictive
# =============================================================================

message("\n=============================================================")
message("  STAGE 2: Posterior predictive")
message("=============================================================")

for (m in MODELS) {
  run_stage(sprintf("run_%s_predictive.R", m),
            args_str = run_ids[[m]],
            label    = sprintf("%s predictive", m))
}

# =============================================================================
# 4.  Per-model residual analysis
# =============================================================================

message("\n=============================================================")
message("  STAGE 3: Per-model residual analysis")
message("=============================================================")

for (m in MODELS) {
  run_stage("run_residual_analysis.R",
            args_str = sprintf("%s %s", m, run_ids[[m]]),
            label    = sprintf("%s residual analysis", m))
}

# =============================================================================
# 5.  Multi-model comparison
# =============================================================================

message("\n=============================================================")
message("  STAGE 4: Multi-model comparison")
message("=============================================================")

comparison_args <- paste(
  sapply(MODELS, function(m) sprintf("%s:%s", m, run_ids[[m]])),
  collapse = " "
)

run_stage("run_multimodel_comparison.R",
          args_str = comparison_args,
          label    = "Multi-model comparison")

# =============================================================================
# 6.  Summary
# =============================================================================

message("\n=============================================================")
message("  HIKET Pipeline Complete")
message("=============================================================")
for (m in MODELS)
  message(sprintf("  %-8s  %s", m, run_ids[[m]]))
message(sprintf("\n  Outputs: %s/diagnostics/", SCRIPT_DIR))
message(sprintf("  Multi-model: %s/diagnostics/multimodel/", SCRIPT_DIR))
message("\nNext steps:")
message("  1. Inspect multimodel_metrics_*.txt for headline results.")
message("  2. Inspect per-model residuals for structure to motivate extensions.")
message("  3. Scale to full ~400 plots on Roihu with N_PLOTS_TEST = NA.")