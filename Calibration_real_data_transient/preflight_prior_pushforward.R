# =============================================================================
# Prior-pushforward pre-flight  (local, no Puhti / no fairshare)
# =============================================================================
# Purpose: before spending a Puhti queue cycle, confirm on the laptop that the
# homogenized priors keep the FULL forward pipeline (steady state + transient
# pre-run + forward + multiplicative likelihood) finite across the prior mass,
# at FULL N. This is the faithful, conservative proxy for the MCMC `inf_rate`:
# the prior is wider than the proposals the sampler will use, so a low blowup
# rate here implies a healthy run.
#
# Faithfulness: we do NOT re-implement the likelihood. We source each model's
# real calibration script up to (but not including) the `run_mcmc_chains(...)`
# launch, which leaves the genuine `ll_fn`, `prior`, `best_x`, `sigma_ppm`,
# `FREE_NAMES` in scope (and runs the script's own pre-MCMC forward sanity).
#
# We separate two non-finite sources:
#   - prior CONSTRAINT rejections (recycling/simplex -> -Inf): legitimate, not
#     a pathology; the sampler simply never visits them.
#   - forward-model BLOWUPS (constraints pass, but ll_fn is non-finite): the
#     pathology we are hunting (the old beta2=0.05 xi detonation).
#
# Usage:
#   Rscript Calibration_real_data_transient/preflight_prior_pushforward.R \
#     [MODEL] [K] [SCALE]
#   MODEL : SP1 TP2 TP3 Yasso07 Yasso15 Yasso20  (default: Yasso07)
#   K     : prior draws (default 300)
#   SCALE : new | old | both   (default both -> before/after contrast)
# =============================================================================

suppressWarnings(suppressMessages({
  args  <- commandArgs(trailingOnly = TRUE)
  MODEL <- if (length(args) >= 1) args[[1]] else "Yasso07"
  K     <- if (length(args) >= 2) as.integer(args[[2]]) else 300L
  SCALE <- if (length(args) >= 3) args[[3]] else "both"
}))

stopifnot(MODEL %in% c("SP1","TP2","TP3","Yasso07","Yasso15","Yasso20"))
set.seed(2025)

script_path <- file.path("Calibration_real_data_transient",
                         sprintf("run_%s_transient_calibration.R", MODEL))
if (!file.exists(script_path)) stop("Cannot find ", script_path,
                                    " (run from project root).")

# --- source the real calibration script up to the MCMC launch -----------------
src   <- readLines(script_path, warn = FALSE)
cutix <- grep("^t_run <- system.time\\(\\{", src)[1]   # line that launches MCMC
if (is.na(cutix)) stop("Could not locate MCMC launch in ", script_path)
message(sprintf("[%s] sourcing setup (lines 1..%d of %d) ...",
                MODEL, cutix - 1L, length(src)))

e <- new.env(parent = globalenv())
ok <- tryCatch({
  source(textConnection(paste(src[seq_len(cutix - 1L)], collapse = "\n")),
         local = e)
  TRUE
}, error = function(err) { message("SETUP ERROR: ", conditionMessage(err)); FALSE })
if (!isTRUE(ok)) quit(status = 1L)

for (nm in c("ll_fn","prior","best_x","sigma_ppm","FREE_NAMES"))
  if (!exists(nm, envir = e)) stop("Setup did not define `", nm, "`.")

ll_fn     <- get("ll_fn",     e)
prior     <- get("prior",     e)
best_x    <- get("best_x",    e)
sigma_new <- get("sigma_ppm", e)
FREE      <- get("FREE_NAMES", e)
N_FREE    <- length(FREE)
n_plots   <- if (exists("plots", e)) length(get("plots", e)) else NA_integer_

# --- reconstruct the pre-fix ("old") sigma for the before/after contrast ------
old_overrides <- c(beta1 = 0.20, beta2 = 0.05, gamma = 0.30,
                   delta1 = 0.15, delta2 = 0.10, r = 0.015)
sigma_old <- sigma_new
for (nm in names(old_overrides))
  if (nm %in% names(sigma_old)) sigma_old[nm] <- old_overrides[[nm]]
frac_ix <- grepl("^p_", names(sigma_old))      # transfer fractions / p_H / p_S
sigma_old[frac_ix] <- 1.0                       # pre-fix fraction width

# --- one timing call so we can size things honestly --------------------------
t1 <- system.time(ll0 <- ll_fn(best_x))[["elapsed"]]
message(sprintf("[%s] N_plots=%s, N_free=%d, 1 likelihood eval = %.3fs, ll@defaults=%.2f",
                MODEL, n_plots, N_FREE, t1, ll0))

pushforward <- function(sigma_vec, label) {
  prior_reject <- 0L; blowup <- 0L; finite_ll <- numeric(0)
  for (k in seq_len(K)) {
    x  <- rnorm(N_FREE, mean = best_x, sd = sigma_vec)   # full prior width
    pd <- tryCatch(prior$density(x), error = function(z) -Inf)
    if (!is.finite(pd)) { prior_reject <- prior_reject + 1L; next }
    ll <- tryCatch(ll_fn(x), error = function(z) NA_real_)
    if (is.finite(ll)) finite_ll <- c(finite_ll, ll) else blowup <- blowup + 1L
  }
  accepted <- K - prior_reject
  data.frame(
    scale            = label,
    draws            = K,
    prior_rejects    = prior_reject,
    constraint_ok    = accepted,
    fwd_blowups      = blowup,
    blowup_rate_pct  = if (accepted > 0) round(100 * blowup / accepted, 1) else NA,
    finite_ll_median = if (length(finite_ll)) round(median(finite_ll), 1) else NA,
    stringsAsFactors = FALSE)
}

res <- list()
if (SCALE %in% c("old","both"))  res[["old"]] <- pushforward(sigma_old, "OLD (pre-fix)")
if (SCALE %in% c("new","both"))  res[["new"]] <- pushforward(sigma_new, "NEW (homogenized)")
out <- do.call(rbind, res)

cat("\n=====================================================================\n")
cat(sprintf("PRIOR-PUSHFORWARD PRE-FLIGHT  ::  %s  (K=%d draws, full N=%s)\n",
            MODEL, K, n_plots))
cat("blowup_rate_pct = forward-model non-finite among CONSTRAINT-PASSING draws\n")
cat("(this is the MCMC inf_rate proxy; lower is better)\n")
cat("---------------------------------------------------------------------\n")
print(out, row.names = FALSE)
cat("=====================================================================\n")
