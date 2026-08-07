# =============================================================================
# preflight_naive_priors.R   (2026-08-07)
#
# QUESTION: can the flux_pair bounded transform be replaced by plain lognormal
# priors on the two auxiliary parameters, without the forward model blowing up?
#
# WHY IT MIGHT BE POSSIBLE. flux_pair was introduced against a real pathology:
# under the old Tier-3 prior (log SD 0.50) TP2/TP3 drove sigma_input to 13-20x,
# i.e. effective litter fluxes of 34-50 tC/ha/yr, above boreal NPP. But that SD
# 0.50 prior puts ~2.4% of its mass above the NPP ceiling, whereas SD 0.30 puts
# ~0%. So the containment may be achievable with a tighter ORDINARY prior, and
# the bounded transform -- with its logit-coordinate priors that cannot be
# stated in physical units -- may no longer be earning its complexity.
#
# WHAT THIS DOES. Same machinery as preflight_prior_pushforward.R: sources the
# real calibration script up to the MCMC launch so the genuine ll_fn is used.
# But it first PATCHES the sourced text to swap flux_pair for two "log" entries
# and to set the naive prior widths. Then it pushes prior draws through and
# reports the forward blow-up rate, exactly as the production pre-flight does.
#
# This isolates ONE change (the transform). Anchoring, centres and everything
# else are left as they are in production, so a difference in blow-up rate is
# attributable to removing the bound and nothing else.
#
# Usage:
#   Rscript doublechecks/preflight_naive_priors.R [MODEL] [K] [SD_INPUT] [SD_INIT]
#   defaults: TP2 300 0.30 0.20
# =============================================================================

suppressWarnings(suppressMessages({
  args     <- commandArgs(trailingOnly = TRUE)
  MODEL    <- if (length(args) >= 1) args[[1]] else "TP2"
  K        <- if (length(args) >= 2) as.integer(args[[2]]) else 300L
  SD_INPUT <- if (length(args) >= 3) as.numeric(args[[3]]) else 0.30
  SD_INIT  <- if (length(args) >= 4) as.numeric(args[[4]]) else 0.20
}))
stopifnot(MODEL %in% c("SP1","TP2","TP3","Yasso07","Yasso15","Yasso20"))
set.seed(2025)

script_path <- file.path("Calibration_real_data_transient",
                         sprintf("run_%s_transient_calibration.R", MODEL))
src   <- readLines(script_path, warn = FALSE)
cutix <- grep("^t_run <- system.time\\(\\{", src)[1]
src   <- src[seq_len(cutix - 1L)]

# ---------------------------------------------------------------------------
# PATCH 1: flux_pair -> two independent log-transformed parameters.
# The spec spans two lines; collapse and rewrite.
# ---------------------------------------------------------------------------
fp_line <- grep('type = "flux_pair"', src)
stopifnot(length(fp_line) == 1L)
src[fp_line] <- '  list(names = c("sigma_input", "sigma_init"), type = "log")'
# the continuation line carries "window = ..., J_bar = ...)" -- blank it, but
# keep the closing paren of param_spec intact.
cont <- fp_line + 1L
src[cont] <- sub('^\\s*window\\s*=.*?\\)\\s*$', '', src[cont])

# ---------------------------------------------------------------------------
# PATCH 2: the J_bar injection block assumes a flux_pair group exists.
# Neutralise the injection but KEEP the transform rebuild and best_x.
# ---------------------------------------------------------------------------
fpidx <- grep('^fp_idx <- which', src)
if (length(fpidx)) src[fpidx] <- 'fp_idx <- integer(0)'
inj <- grep('^param_spec\\[\\[fp_idx\\]\\]\\$J_bar', src)
if (length(inj)) src[inj] <- '# (J_bar injection disabled: naive priors, no flux_pair)'
msg <- grep('^message\\(sprintf\\("flux_pair J_bar', src)
if (length(msg)) { src[msg] <- 'message("NAIVE PRIORS: flux_pair replaced by log-normal")'
                   src[msg + 1L] <- '' }

# ---------------------------------------------------------------------------
# PATCH 3: the naive prior widths, applied after sigma_ppm is built.
# ---------------------------------------------------------------------------
src <- c(src, sprintf(
  'sigma_ppm["sigma_input"] <- %f; sigma_ppm["sigma_init"] <- %f', SD_INPUT, SD_INIT))

message(sprintf("[%s] sourcing PATCHED setup (naive priors: sd_input=%.2f, sd_init=%.2f) ...",
                MODEL, SD_INPUT, SD_INIT))
e  <- new.env(parent = globalenv())
ok <- tryCatch({ source(textConnection(paste(src, collapse = "\n")), local = e); TRUE },
               error = function(err) { message("SETUP ERROR: ", conditionMessage(err)); FALSE })
if (!isTRUE(ok)) quit(status = 1L)

ll_fn  <- get("ll_fn", e);  prior <- get("prior", e)
best_x <- get("best_x", e); sigma_vec <- get("sigma_ppm", e)
FREE   <- get("FREE_NAMES", e); N_FREE <- length(FREE)
n_plots <- if (exists("plots", e)) length(get("plots", e)) else NA_integer_
to_original <- get("to_original", e)

t1 <- system.time(ll0 <- ll_fn(best_x))[["elapsed"]]
message(sprintf("[%s] N_plots=%s N_free=%d  1 eval=%.3fs  ll@defaults=%.2f",
                MODEL, n_plots, N_FREE, t1, ll0))

# ---------------------------------------------------------------------------
# Pushforward, additionally recording the implied effective flux so we can see
# how much prior mass lands above the NPP ceiling WITHOUT a bound to stop it.
# ---------------------------------------------------------------------------
J_bar <- if (exists("J_bar", e)) get("J_bar", e) else NA_real_
prior_reject <- 0L; blowup <- 0L; fin <- numeric(0); si <- numeric(0)
for (k in seq_len(K)) {
  x  <- rnorm(N_FREE, mean = best_x, sd = sigma_vec)
  pd <- tryCatch(prior$density(x), error = function(z) -Inf)
  if (!is.finite(pd)) { prior_reject <- prior_reject + 1L; next }
  p  <- tryCatch(to_original(x), error = function(z) NULL)
  if (!is.null(p) && "sigma_input" %in% names(p)) si <- c(si, unname(p["sigma_input"]))
  ll <- tryCatch(ll_fn(x), error = function(z) NA_real_)
  if (is.finite(ll)) fin <- c(fin, ll) else blowup <- blowup + 1L
}
acc <- K - prior_reject

cat("\n=====================================================================\n")
cat(sprintf("NAIVE-PRIOR PRE-FLIGHT :: %s  (K=%d, full N=%s)\n", MODEL, K, n_plots))
cat(sprintf("flux_pair REMOVED; sigma_input ~ logN(., %.2f), sigma_init ~ logN(., %.2f)\n",
            SD_INPUT, SD_INIT))
cat("---------------------------------------------------------------------\n")
cat(sprintf("  prior constraint rejects : %d / %d\n", prior_reject, K))
cat(sprintf("  forward blow-ups         : %d / %d  (%.1f%% of constraint-passing)\n",
            blowup, acc, if (acc > 0) 100*blowup/acc else NA))
cat(sprintf("  median finite ll         : %.1f\n", if (length(fin)) median(fin) else NA))
if (length(si) && is.finite(J_bar)) {
  cat(sprintf("  sigma_input prior draws  : median %.2f, 97.5%% %.2f\n",
              median(si), quantile(si, .975)))
  cat(sprintf("  implied effective flux   : median %.2f, 97.5%% %.2f tC/ha/yr\n",
              median(si)*J_bar, quantile(si, .975)*J_bar))
  cat(sprintf("  prior mass above NPP 8.7 : %.1f%%\n", 100*mean(si*J_bar > 8.7)))
}
cat("=====================================================================\n")
cat("Compare the blow-up rate with preflight_prior_pushforward.R (same model,\n")
cat("same K) which runs the production flux_pair configuration.\n")
