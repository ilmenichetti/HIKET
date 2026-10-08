# =============================================================================
# run_ids_eqinit.R   (2026-10-08)
#
# Run selection for the equilibrium-init counterfactual figures (F16-F19, T_eqinit;
# design in NEXT_RUN_equilibrium_init.md). Two arms:
#   RID     production (transient start)  -- from run_ids.R, unchanged
#   RID_EQ  equilibrium start             -- newest <MODEL>_eqinit_posterior_*.rds
#
# The eqinit outputs carry MODEL_NAME = "<MODEL>_eqinit", so run_ids.R can never
# select them; this file is the ONLY place that does.
#
# Partial landing is allowed: EQ_MODELS lists the models whose equilibrium arm has
# BOTH the posterior and the predictive bundle; figures draw those and say which
# are missing. Overrides:
#   HIKET_FIG_RID_EQ="SP1=<id>,..."   pin equilibrium RUN_IDs
#   HIKET_FIG_EQ_RUNS / HIKET_FIG_EQ_INPUTS   other directories (used for testing)
# =============================================================================

source("manuscript/figures/run_ids.R")

.eq_runs   <- Sys.getenv("HIKET_FIG_EQ_RUNS",   .fig_runs)
.eq_inputs <- Sys.getenv("HIKET_FIG_EQ_INPUTS", .fig_inputs)

RID_EQ <- vapply(FIG_MODELS, function(m) {
  fs <- list.files(.eq_runs, pattern = sprintf("^%s_eqinit_posterior_[0-9]{8}_[0-9]{6}\\.rds$", m))
  if (!length(fs)) return(NA_character_)
  sub(sprintf("^%s_eqinit_posterior_(.+)\\.rds$", m), "\\1", sort(fs, decreasing = TRUE)[1])
}, character(1))

local({
  ov <- Sys.getenv("HIKET_FIG_RID_EQ", "")
  if (nzchar(ov)) {
    kv <- strsplit(strsplit(ov, ",")[[1]], "=")
    for (e in kv) RID_EQ[[trimws(e[1])]] <<- trimws(e[2])
    message("run_ids_eqinit.R: RID_EQ OVERRIDDEN by HIKET_FIG_RID_EQ")
  }
})

# File paths for one arm of one model.
arm_files <- function(m, arm = c("prod", "eq")) {
  arm <- match.arg(arm)
  if (arm == "prod")
    return(list(post  = file.path(.fig_runs,   sprintf("%s_posterior_%s.rds", m, RID[[m]])),
                pred  = file.path(.fig_runs,   sprintf("%s_posterior_predictive_%s.rds", m, RID[[m]])),
                chain = file.path(.fig_runs,   sprintf("%s_chains_%s.rds", m, RID[[m]])),
                inp   = file.path(.fig_inputs, sprintf("%s_inputs_%s.rds", m, RID[[m]])),
                rid   = RID[[m]]))
  tag <- paste0(m, "_eqinit")
  list(post  = file.path(.eq_runs,   sprintf("%s_posterior_%s.rds", tag, RID_EQ[[m]])),
       pred  = file.path(.eq_runs,   sprintf("%s_posterior_predictive_%s.rds", tag, RID_EQ[[m]])),
       chain = file.path(.eq_runs,   sprintf("%s_chains_%s.rds", tag, RID_EQ[[m]])),
       inp   = file.path(.eq_inputs, sprintf("%s_inputs_%s.rds", tag, RID_EQ[[m]])),
       rid   = RID_EQ[[m]])
}

EQ_MODELS <- FIG_MODELS[vapply(FIG_MODELS, function(m) {
  if (is.na(RID_EQ[[m]])) return(FALSE)
  f <- arm_files(m, "eq")
  all(file.exists(c(f$post, f$pred, f$inp)))
}, logical(1))]

EQ_MISSING <- setdiff(FIG_MODELS, EQ_MODELS)
if (!length(EQ_MODELS))
  stop("run_ids_eqinit.R: no equilibrium arm with posterior + predictive bundle in ", .eq_runs,
       "\nRun run_<MODEL>_transient_predictive.R with HIKET_EQUILIBRIUM_INIT=1 first.", call. = FALSE)
message("Equilibrium arm: ", paste(sprintf("%s=%s", EQ_MODELS, RID_EQ[EQ_MODELS]), collapse = "  "),
        if (length(EQ_MISSING)) paste0("\n  MISSING (not drawn): ", paste(EQ_MISSING, collapse = ", ")) else "")
