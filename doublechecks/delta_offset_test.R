# =============================================================================
# delta_offset_test.R   (2026-08-05)
#
# QUESTION: can the concern behind C5 be expressed FAITHFULLY instead of by proxy?
#
# C5 inflates the per-observation SD of the 1985 campaign. But the concern it
# encodes -- that the VMI8 methodology differs from 2006/2024 and could be biased
# -- is a SYSTEMATIC, SHARED offset, not extra independent noise. Those behave
# differently: independent noise averages away over 441 plots (the campaign mean's
# SE scales as sigma/sqrt(n)), while a shared offset does not shrink with n at all.
# So C5 inflates precisely the component that averages out.
#
# The faithful representation is a campaign-level offset:
#     obs_1985  ~  Normal( (1 + delta) * SOC_hat ,  (1 + delta) * SOC_hat * sigma_obs )
# with delta < 0 meaning "VMI8 reads LOW relative to the model".
#
# The revision plan (removed; NEXT_SESSION.md §5) ruled this out as "non-identifiable vs sigma_init". That may be too
# strong: sigma_init sets C_init and propagates through ALL THREE campaigns (humus
# MRT ~165 yr, so most of it persists to 2024), whereas delta shifts ONLY the 1985
# comparison. Different signatures => 2006/2024 should partially separate them.
#
# This script tests that directly. It adds delta as a free parameter with a weak
# prior and asks:
#   1. Is delta IDENTIFIED?  (posterior SD << prior SD => the data inform it)
#   2. Where does it sit?    (delta < 0 would support "VMI8 reads low")
#   3. How confounded is it with sigma_init? (posterior correlation)
#   4. What does admitting delta do to sigma_init itself?
#
# It touches NO production files: the forward model is rebuilt here from the
# sourced calibration setup and VALIDATED against the engine's own ll_fn at
# delta = 0 before anything is sampled.
#
# Usage:  Rscript doublechecks/delta_offset_test.R [MODEL] [TEST_ITER] [N_CHAINS] [TEST_CORES]
# =============================================================================

suppressWarnings(suppressMessages({
  library(BayesianTools)
  args    <- commandArgs(trailingOnly = TRUE)
  MODEL   <- if (length(args) >= 1) args[[1]] else "TP2"
  TEST_ITER   <- if (length(args) >= 2) as.integer(args[[2]]) else 4000L
  TEST_CHAINS <- if (length(args) >= 3) as.integer(args[[3]]) else 3L
  TEST_CORES  <- if (length(args) >= 4) as.integer(args[[4]]) else 4L
}))
set.seed(2025)

# NB: the run script sources calib_config.R WITHOUT local=, so its N_ITER / N_CHAINS land in
# the global environment and would silently overwrite same-named variables here. Hence the
# TEST_* prefixes -- a plain N_ITER was clobbered to 50000 on the first attempt.

DELTA_SD <- 0.25   # weak prior: +-25% systematic offset at 1 sd. Deliberately WIDE --
                   # a tight prior would return itself and prove nothing about identifiability.

script <- sprintf("Calibration_real_data_transient/run_%s_transient_calibration.R", MODEL)
src    <- readLines(script, warn = FALSE)
cutix  <- grep("^t_run <- system.time\\(\\{", src)[1]
message(sprintf("[%s] sourcing setup ...", MODEL))
e <- new.env(parent = globalenv())
invisible(capture.output(suppressMessages(
  source(textConnection(paste(src[seq_len(cutix - 1L)], collapse = "\n")), local = e))))

# --- engine bindings (parse the real make_likelihood call; args may be expressions)
start <- grep("ll_fn <- make_likelihood\\(", src)[1]
open <- 0L; end <- NA_integer_
for (i in seq(start, length(src))) {
  ch <- strsplit(src[i], "")[[1]]
  open <- open + sum(ch == "(") - sum(ch == ")")
  if (open == 0L) { end <- i; break }
}
ml <- as.list(str2lang(sub("^\\s*ll_fn\\s*<-\\s*", "",
                           paste(src[seq(start, end)], collapse = "\n"))))[-1]
argof <- function(nm, d = NULL) if (is.null(ml[[nm]])) d else eval(ml[[nm]], envir = e)

to_original     <- argof("to_original")
log_jacobian    <- argof("log_jacobian")
assemble_params <- argof("assemble_params")
compute_xi      <- argof("compute_xi")
compute_xi_mean <- argof("compute_xi_mean")
steady_state    <- argof("steady_state")
run_model       <- argof("run_model")
STEADY_N        <- argof("steady_state_n", NULL)
sigma_obs_fixed <- argof("sigma_obs_fixed")

plots           <- get("plots",           e)
climate_by_plot <- get("climate_by_plot", e)
inputs_by_plot  <- get("inputs_by_plot",  e)
litter_means    <- get("litter_means",    e)
obs_meta        <- get("obs_meta",        e)
prior           <- get("prior",           e)
best_x          <- get("best_x",          e)
FREE_NAMES      <- get("FREE_NAMES",      e)
ll_fn           <- get("ll_fn",           e)
SOC_obs_all     <- get("SOC_obs_all",     e)
N_FREE          <- length(FREE_NAMES)

# which observations belong to the 1985 campaign, per plot
is85 <- lapply(plots, function(pid)
  SOC_obs_all$year[as.character(SOC_obs_all$plot_id) == pid] == 1985L)
names(is85) <- plots

# =============================================================================
# Likelihood with a campaign-level offset.
#   use_c5 = TRUE, delta = 0  -> must reproduce the engine's ll_fn EXACTLY (guard)
#   use_c5 = FALSE, delta free -> the model under test
# =============================================================================
ll_delta <- function(x, delta = 0, use_c5 = FALSE) {
  p_free <- to_original(x)
  mp     <- assemble_params(p_free)
  lj     <- log_jacobian(x, p_free)

  ll <- unlist(parallel::mclapply(plots, function(pid) {
    clim <- climate_by_plot[[pid]]; inputs <- inputs_by_plot[[pid]]
    lm   <- litter_means[[pid]];    meta   <- obs_meta[[pid]]
    if (any(is.na(meta$idx))) return(-Inf)
    xa <- tryCatch(compute_xi(clim, mp), error = function(z) NULL)
    if (is.null(xa)) return(-Inf)
    n_ss <- if (is.null(STEADY_N)) nrow(clim) else min(STEADY_N, nrow(clim))
    xs <- tryCatch(compute_xi_mean(clim[seq_len(n_ss), , drop = FALSE], mp),
                   error = function(z) NULL)
    if (is.null(xs)) return(-Inf)
    ci <- tryCatch(steady_state(mp, lm, xs), error = function(z) NULL)
    if (is.null(ci) || any(!is.finite(ci)) || any(ci < 0)) return(-Inf)
    ro <- tryCatch(run_model(inputs, mp, ci, xa), error = function(z) NULL)
    if (is.null(ro)) return(-Inf)
    SOC_hat <- ro$total_soc[meta$idx]
    if (any(!is.finite(SOC_hat)) || any(SOC_hat <= 0)) return(-Inf)

    # campaign offset applied to the 1985 predictions only
    mu <- SOC_hat
    f  <- is85[[pid]]
    if (any(f)) mu[f] <- mu[f] * (1 + delta)
    if (any(mu <= 0)) return(-Inf)

    sd_vec <- mu * sigma_obs_fixed
    if (use_c5 && !is.null(meta$sigma_infl)) sd_vec <- sd_vec * meta$sigma_infl
    sum(dnorm(meta$soc_obs, mean = mu, sd = sd_vec, log = TRUE))
  }, mc.cores = TEST_CORES))

  if (any(!is.finite(ll))) return(-Inf)
  sum(ll) + lj
}

# --- GUARD: at delta = 0 with C5 on, this must BE the engine's likelihood -------
ll_engine <- ll_fn(best_x)
ll_mine   <- ll_delta(best_x, delta = 0, use_c5 = TRUE)
cat(sprintf("\nFaithfulness check: engine %.4f | reconstructed %.4f | |diff| %.2e\n",
            ll_engine, ll_mine, abs(ll_engine - ll_mine)))
if (abs(ll_engine - ll_mine) > 1e-6) {
  cat("MISMATCH -- not the engine's likelihood. Aborting.\n"); quit(status = 1L)
}
cat("Reconstructed likelihood == engine likelihood. Proceeding.\n")

# =============================================================================
# Sample over (free params, delta). C5 is OFF: delta REPLACES it.
# =============================================================================
target <- function(z) ll_delta(z[seq_len(N_FREE)], delta = z[N_FREE + 1L], use_c5 = FALSE)

dens <- function(z) {
  pd <- tryCatch(prior$density(z[seq_len(N_FREE)]), error = function(q) -Inf)
  if (!is.finite(pd)) return(-Inf)
  pd + dnorm(z[N_FREE + 1L], 0, DELTA_SD, log = TRUE)
}
samp <- function(n = 1) {
  s <- prior$sampler(n)
  if (is.null(dim(s))) s <- matrix(s, nrow = 1)
  cbind(s, rnorm(nrow(s), 0, DELTA_SD))
}
low <- c(rep(-Inf, N_FREE), -0.9)
upp <- c(rep( Inf, N_FREE),  3.0)

setup <- createBayesianSetup(
  likelihood = target,
  prior      = createPrior(density = dens, sampler = samp, lower = low, upper = upp),
  names      = c(FREE_NAMES, "delta"))

cat(sprintf("\nSampling: %d chains x %d iter, %d free + delta, %d cores, prior delta ~ N(0, %.2f)\n",
            TEST_CHAINS, TEST_ITER, N_FREE, TEST_CORES, DELTA_SD))
t0 <- Sys.time()
chains <- lapply(seq_len(TEST_CHAINS), function(k) {
  message(sprintf("  chain %d/%d ...", k, TEST_CHAINS))
  out <- runMCMC(setup, sampler = "DEzs",
                 settings = list(iterations = TEST_ITER, message = FALSE,
                                 startValue = 3, consoleUpdates = 1e9))
  getSample(out, start = floor(TEST_ITER / 5), coda = FALSE)
})
cat(sprintf("done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

post <- do.call(rbind, chains)
colnames(post) <- c(FREE_NAMES, "delta")

# --- transform the model parameters back to physical scale for reporting -------
phys <- t(apply(post[, seq_len(N_FREE), drop = FALSE], 1, to_original))
colnames(phys) <- FREE_NAMES
d <- post[, "delta"]

# --- convergence: simple between/within on delta -------------------------------
dl <- lapply(chains, function(m) m[, N_FREE + 1L])
W  <- mean(vapply(dl, var, numeric(1)))
B  <- var(vapply(dl, mean, numeric(1))) * length(dl[[1]])
n  <- length(dl[[1]])
rhat_d <- sqrt(((n - 1)/n * W + B/n) / W)

cat("\n=====================================================================\n")
cat(sprintf("DELTA-OFFSET TEST  ::  %s   (C5 OFF; delta replaces it)\n", MODEL))
cat("=====================================================================\n")
cat(sprintf("\ndelta   prior  N(0, %.3f)\n", DELTA_SD))
cat(sprintf("        post   median %+0.4f  [2.5%% %+0.4f, 97.5%% %+0.4f]  sd %.4f  R-hat %.3f\n",
            median(d), quantile(d, .025), quantile(d, .975), sd(d), rhat_d))
cat(sprintf("\n1. IDENTIFIED?   posterior sd / prior sd = %.3f   -> %s\n",
            sd(d)/DELTA_SD,
            if (sd(d)/DELTA_SD < 0.7) "YES, the data inform delta" else
              "NO, delta is essentially prior-driven"))
cat(sprintf("2. SIGN:         P(delta < 0) = %.3f   -> %s\n", mean(d < 0),
            if (mean(d < 0) > 0.95) "VMI8 reads LOW (supports the C5 premise)" else
            if (mean(d > 0) > 0.95) "VMI8 reads HIGH" else
              "no clear direction -- premise NOT supported"))
cat(sprintf("3. CONFOUNDING:  cor(delta, sigma_init) = %+0.3f\n",
            cor(d, phys[, "sigma_init"])))
cat(sprintf("4. sigma_init:   median %.4f  [%.4f, %.4f]  (C5=2.0 ref: 0.285; C5 off ref: see A1)\n",
            median(phys[, "sigma_init"]),
            quantile(phys[, "sigma_init"], .025), quantile(phys[, "sigma_init"], .975)))
cat(sprintf("   sigma_input:  median %.4f  [%.4f, %.4f]\n",
            median(phys[, "sigma_input"]),
            quantile(phys[, "sigma_input"], .025), quantile(phys[, "sigma_input"], .975)))
cat("\nInterpretation: an identified delta near 0 means the 1985 campaign shows no\n")
cat("systematic offset once the model is free to find one -- i.e. C5 had nothing to fix.\n")
cat("=====================================================================\n")

saveRDS(list(post = post, phys = phys, delta = d, DELTA_SD = DELTA_SD),
        file.path("doublechecks", sprintf("delta_offset_%s.rds", MODEL)))
