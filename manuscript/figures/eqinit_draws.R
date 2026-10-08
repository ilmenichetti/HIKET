# =============================================================================
# eqinit_draws.R   (2026-10-08)
#
# Per-draw quantities for F18 and T_eqinit that need each model's own transforms,
# i.e. a re-sourced calibration script (as F14/S13 do), once per arm:
#   ll_data   the DATA log-likelihood of every retained draw. The sampler's stored
#             Llikelihood includes the log-Jacobian, which DIFFERS between the arms
#             (production carries the extra sigma_init coordinate), so it is
#             subtracted before any comparison. Comparable across arms: same data,
#             same correlated error model, same sigma.
#   mtt       intrinsic mean transit time (unit input, dataset-mean reference),
#             computed exactly as F14 (Yasso) and S13 (SP1/TP2/TP3), on a thinned
#             subsample of N_MTT draws. Independent of sigma_input and sigma_init.
#   si        sigma_input of the same draws.
# Chains are extracted per sampler with start = 2 (the first retained iteration of
# each internal DEzs chain sits 100-200 ll below the bulk -- CLAUDE.md, F14).
#
# Cached in manuscript/figures/eqinit_draws.rds, stamped with both arms' RUN_IDs.
# Re-sourcing a calibration script writes an input bundle, a sanity PNG and, for the
# equilibrium arm, diagnostics/<M>_eqinit/; everything created here is removed.
#
# Usage:  Rscript manuscript/figures/eqinit_draws.R [N_MTT]     (~minutes per model)
# Testing only: HIKET_FIG_EQ_SETUP_PROD=1 sources the eq arm WITHOUT the switch
# (for the symlinked prod-as-eq test, whose chains have production's dimension).
# =============================================================================

setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
a <- commandArgs(trailingOnly = TRUE)
N_MTT <- if (length(a)) as.integer(a[[1]]) else 1500L
set.seed(2025)

Sys.setenv(HIKET_EQUILIBRIUM_INIT = "0")            # ref (inside the lib) is built on production
snap <- function() c(list.files("Data/model_inputs", full.names = TRUE),
                     list.files("Calibration_real_data_transient/diagnostics", recursive = TRUE,
                                full.names = TRUE, include.dirs = TRUE))
before <- snap()
source("doublechecks/intrinsic_mrt_lib.R")          # setup(), ref, mrt_fun() (Yasso)
# AFTER the lib: it defines its own Yasso-only RID, which run_ids.R must overwrite.
source("manuscript/figures/run_ids_eqinit.R")
suppressMessages(library(BayesianTools))

# Simple models: as build_S13_mrt_ridge_benchmark.R (unit input, sigma_input -> 1,
# PURE steady state at the reference climate -- never the transient-init binding).
mrt_fun_simple <- function(M, e) {
  ss  <- get(sprintf("%s_steady_state", tolower(M)), envir = e)
  cxm <- get(sprintf("compute_xi_mean_%s_engine", tolower(M)), envir = e)
  function(p) {
    mp <- e$.assemble(p); mp["sigma_input"] <- 1
    sum(ss(mp, list(J_total_mean = 1), cxm(clim_ss = ref$clim, model_params = mp)))
  }
}

arm_draws <- function(M, arm) {
  Sys.setenv(HIKET_EQUILIBRIUM_INIT =
               if (arm == "eq" && Sys.getenv("HIKET_FIG_EQ_SETUP_PROD") != "1") "1" else "0")
  e <- setup(M)
  f <- if (M %in% c("SP1", "TP2", "TP3")) mrt_fun_simple(M, e) else mrt_fun(M, e)
  ch <- readRDS(arm_files(M, arm)$chain)
  s  <- do.call(rbind, lapply(ch, function(z) getSample(z, parametersOnly = FALSE, start = 2)))
  free <- names(get("best_x", e))
  if (!all(free %in% colnames(s)))
    stop(M, " ", arm, ": chain columns do not match the script's free parameters -- wrong arm?")
  X  <- s[, free, drop = FALSE]
  lj <- vapply(seq_len(nrow(X)), function(k) { x <- X[k, ]; e$log_jacobian(x, e$.to_original(x)) }, numeric(1))
  si <- vapply(seq_len(nrow(X)), function(k) e$.to_original(X[k, ])[["sigma_input"]], numeric(1))
  ll <- s[, "Llikelihood"] - lj
  k  <- sort(sample(nrow(X), min(N_MTT, nrow(X))))
  mtt <- vapply(k, function(i) tryCatch(f(e$.to_original(X[i, ])), error = function(z) NA_real_), numeric(1))
  cat(sprintf("  %-8s %-4s %6d draws | data ll max %.2f median %.2f | MTT median %.1f | sigma_input %.3f\n",
              M, arm, nrow(X), max(ll), median(ll), median(mtt, na.rm = TRUE), median(si)))
  list(ll = ll, si = si, thin = data.frame(mtt = mtt, si = si[k], ll = ll[k]),
       J_bar = get("J_bar", e), n_free = length(free))
}

stamp <- list(prod = RID[EQ_MODELS], eq = RID_EQ[EQ_MODELS])
out <- list(stamp = stamp, models = EQ_MODELS, d = list())
for (M in EQ_MODELS) out$d[[M]] <- list(prod = arm_draws(M, "prod"), eq = arm_draws(M, "eq"))
Sys.setenv(HIKET_EQUILIBRIUM_INIT = "0")

new <- setdiff(snap(), before)
unlink(new[!dir.exists(new)]); unlink(new[dir.exists(new)], recursive = TRUE)
cat(sprintf("removed %d setup side-effect files\n", length(new)))
saveRDS(out, "manuscript/figures/eqinit_draws.rds")
cat("saved manuscript/figures/eqinit_draws.rds\n")
