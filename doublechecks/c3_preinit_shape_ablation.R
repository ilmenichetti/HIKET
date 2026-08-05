# =============================================================================
# c3_preinit_shape_ablation.R   (C3 attribution, 2026-08-04)
#
# Question: how much of the modelled 1985 stock -- and of the "+26 tC/ha above
# observed 1985" gap -- is caused by the LINEAR pre-run litter ramp, as opposed
# to anything the calibration does?
#
# Why a forward-only test: C3 changes the 1917->1985 spin-up FORCING, so its
# mechanical effect is visible with the parameters HELD FIXED. No MCMC needed.
# That makes it the one factor we can attribute cleanly before spending a
# six-model Roihu run in which C1/C2/C3/C4b/C5 + the new SOC target all move at
# once.
#
# Method: source the real calibration script up to the MCMC launch (same trick
# as preflight_prior_pushforward.R, so we use the GENUINE wrappers, priors and
# data), then run the engine's own forward path twice per plot at the SAME
# parameters -- once with the growing-stock shape, once with it stripped so the
# wrappers fall back to the linear ramp.
#
# Faithfulness guard: we re-derive the total log-likelihood from our own forward
# pass and check it against the engine's ll_fn(best_x). If those agree, our
# forward IS the engine's forward. The script refuses to report if they diverge.
#
# Usage:
#   Rscript doublechecks/c3_preinit_shape_ablation.R [MODEL]
#   MODEL : SP1 TP2 TP3 Yasso07 Yasso15 Yasso20   (default: TP2)
# =============================================================================

suppressWarnings(suppressMessages({
  args  <- commandArgs(trailingOnly = TRUE)
  MODEL <- if (length(args) >= 1) args[[1]] else "TP2"
}))
stopifnot(MODEL %in% c("SP1","TP2","TP3","Yasso07","Yasso15","Yasso20"))
set.seed(2025)

script_path <- file.path("Calibration_real_data_transient",
                         sprintf("run_%s_transient_calibration.R", MODEL))
if (!file.exists(script_path)) stop("Cannot find ", script_path, " (run from project root).")

# --- source the real calibration script up to the MCMC launch ----------------
src   <- readLines(script_path, warn = FALSE)
cutix <- grep("^t_run <- system.time\\(\\{", src)[1]
if (is.na(cutix)) stop("Could not locate MCMC launch in ", script_path)
message(sprintf("[%s] sourcing setup (lines 1..%d of %d) ...", MODEL, cutix - 1L, length(src)))

e <- new.env(parent = globalenv())
source(textConnection(paste(src[seq_len(cutix - 1L)], collapse = "\n")), local = e)

# --- recover the engine bindings actually passed to make_likelihood ----------
# The argument names are model-specific (compute_xi_tp2_engine,
# steady_state_tp2_engine, ...) AND some arguments are EXPRESSIONS, not bare
# symbols: Yasso20 passes `steady_state_n = STEADY_STATE_YEARS * 12L` because its
# climate_by_plot holds MONTHLY rows. So we parse the call properly and evaluate
# each argument in the sourced environment -- a symbol-only regex silently got
# Yasso20's steady-state window wrong by a factor of 12 (caught by the guard).
start <- grep("ll_fn <- make_likelihood\\(", src)[1]
open <- 0L; end <- NA_integer_
for (i in seq(start, length(src))) {
  chars <- strsplit(src[i], "")[[1]]
  open <- open + sum(chars == "(") - sum(chars == ")")
  if (open == 0L) { end <- i; break }
}
if (is.na(end)) stop("Could not find the end of the make_likelihood call.")
call_txt <- paste(src[seq(start, end)], collapse = "\n")
call_txt <- sub("^\\s*ll_fn\\s*<-\\s*", "", call_txt)
ml_args  <- as.list(str2lang(call_txt))[-1]
argof <- function(nm, default = NULL) {
  if (is.null(ml_args[[nm]])) return(default)
  eval(ml_args[[nm]], envir = e)
}
to_original     <- argof("to_original")
assemble_params <- argof("assemble_params")
compute_xi      <- argof("compute_xi")
compute_xi_mean <- argof("compute_xi_mean")
steady_state    <- argof("steady_state")
run_model       <- argof("run_model")
transient_init  <- isTRUE(argof("transient_init", FALSE))

plots           <- get("plots",           e)
climate_by_plot <- get("climate_by_plot", e)
inputs_by_plot  <- get("inputs_by_plot",  e)
litter_means    <- get("litter_means",    e)
obs_meta        <- get("obs_meta",        e)
sigma_obs_fixed <- get("sigma_obs_fixed", e)
best_x          <- get("best_x",          e)
ll_fn           <- get("ll_fn",           e)
STEADY_N        <- argof("steady_state_n", NULL)   # NB Yasso20: monthly rows (x12)

# Observation years, recovered per plot (obs_meta carries sigma_infl, not year).
SOC_obs_all <- get("SOC_obs_all", e)
obs_years   <- lapply(plots, function(pid)
  SOC_obs_all$year[as.character(SOC_obs_all$plot_id) == pid])
names(obs_years) <- plots

# =============================================================================
# Forward pass -- mirrors calibration_engine_transient.R make_likelihood()
# =============================================================================
forward_one <- function(pid, model_params, sigma_init, use_shape) {
  clim   <- climate_by_plot[[pid]]
  inputs <- inputs_by_plot[[pid]]
  lm     <- litter_means[[pid]]
  meta   <- obs_meta[[pid]]
  if (any(is.na(meta$idx))) return(NULL)

  # THE TOGGLE: strip the growing-stock shape -> wrappers fall back to
  # frac <- (i-1)/(n_pre-1), i.e. the old linear 1917->1985 ramp.
  if (!use_shape) lm$preinit_shape <- NULL

  xi_array <- tryCatch(compute_xi(clim, model_params), error = function(z) NULL)
  if (is.null(xi_array)) return(NULL)
  n_ss <- if (is.null(STEADY_N)) nrow(clim) else min(STEADY_N, nrow(clim))
  xi_ss <- tryCatch(compute_xi_mean(clim[seq_len(n_ss), , drop = FALSE], model_params),
                    error = function(z) NULL)
  if (is.null(xi_ss)) return(NULL)

  C_init <- tryCatch(steady_state(model_params, lm, xi_ss), error = function(z) NULL)
  if (is.null(C_init) || any(!is.finite(C_init)) || any(C_init < 0)) return(NULL)

  run_out <- tryCatch(run_model(inputs, model_params, C_init, xi_array),
                      error = function(z) NULL)
  if (is.null(run_out)) return(NULL)

  SOC_hat <- run_out$total_soc[meta$idx]
  if (any(!is.finite(SOC_hat)) || any(SOC_hat <= 0)) return(NULL)

  # Engine's observation-error model under transient_init = TRUE (all six models):
  # sigma_init is propagated PHYSICALLY through the 68-yr pre-run, so it is NOT
  # added in quadrature here -- doing so would double-count it.
  sd_vec <- SOC_hat * sigma_obs_fixed
  # C5: down-weight the suspect 1985 (VMI8) campaign.
  if (!is.null(meta$sigma_infl)) sd_vec <- sd_vec * meta$sigma_infl

  list(C_init_tot = sum(C_init),
       SOC_hat    = SOC_hat,
       soc_obs    = meta$soc_obs,
       years      = obs_years[[pid]],
       ll         = sum(dnorm(meta$soc_obs, SOC_hat, sd_vec, log = TRUE)))
}

run_all <- function(use_shape, sigma_init_set = NA_real_) {
  p_free <- to_original(best_x)
  if (is.finite(sigma_init_set)) p_free["sigma_init"] <- sigma_init_set
  mp     <- assemble_params(p_free)
  si     <- p_free["sigma_init"]
  res    <- lapply(plots, forward_one, model_params = mp, sigma_init = si,
                   use_shape = use_shape)
  names(res) <- plots
  res[!vapply(res, is.null, logical(1))]
}

# -----------------------------------------------------------------------------
# WHY sigma_init must be swept, not left at its prior centre.
#
#   J_1917 = lm$J_full_mean * sigma_init * sigma_input     (whole-record mean)
#   J_1985 = lm$J_t0_mean   * sigma_input                  (1985 value)
#
# Litter RISES after 1985 (population mean ~2.1 -> ~3.0 by the mid-2000s), so
# J_full_mean > J_t0_mean. At the prior centre sigma_init = 1 this makes
# J_1917 > J_1985: the pre-run litter DECLINES from 1917 to 1985, the reverse of
# the depleted-forest history C3 is meant to encode. sigma_init is precisely the
# parameter that fixes this -- the old posteriors put it at 0.19-0.72, which
# pulls J_1917 well below J_1985 and restores the recovering-forest shape.
#
# So C3's leverage is a FUNCTION of sigma_init, and evaluating only at the prior
# centre would understate it (and even flip its sign). We sweep the posterior
# range instead.
# -----------------------------------------------------------------------------
SIGMA_INIT_GRID <- c(0.20, 0.30, 0.50, 0.70, 1.00)

message(sprintf("[%s] forward pass: growing-stock shape (C3 ON) ...", MODEL))
on  <- run_all(TRUE)
message(sprintf("[%s] forward pass: linear ramp (C3 OFF) ...", MODEL))
off <- run_all(FALSE)

# --- faithfulness guard ------------------------------------------------------
ll_engine <- ll_fn(best_x)
ll_mine   <- sum(vapply(on, `[[`, numeric(1), "ll")) +
  get("log_jacobian", e)(best_x, to_original(best_x))
d <- abs(ll_engine - ll_mine)
cat(sprintf("\nFaithfulness check: engine ll = %.4f | reconstructed = %.4f | |diff| = %.2e\n",
            ll_engine, ll_mine, d))
if (d > 1e-6) {
  cat("MISMATCH -- reconstructed forward is NOT the engine's forward. Not reporting.\n")
  quit(status = 1L)
}
cat("Reconstructed forward == engine forward. Proceeding.\n")

# =============================================================================
# Report
# =============================================================================
gather <- function(res, yr) {
  unlist(lapply(res, function(r) r$SOC_hat[r$years == yr]))
}
gather_obs <- function(res, yr) {
  unlist(lapply(res, function(r) r$soc_obs[r$years == yr]))
}
common <- intersect(names(on), names(off))
on <- on[common]; off <- off[common]

cat(sprintf("\n=====================================================================\n"))
cat(sprintf("C3 PRE-RUN SHAPE ABLATION  ::  %s   (%d plots, parameters FIXED at prior centre)\n",
            MODEL, length(common)))
cat(sprintf("=====================================================================\n"))

ci_on  <- vapply(on,  `[[`, numeric(1), "C_init_tot")
ci_off <- vapply(off, `[[`, numeric(1), "C_init_tot")
cat(sprintf("\nInitial state at 1985 (C_init, tC/ha, cross-plot median):\n"))
cat(sprintf("  growing-stock shape : %7.2f\n", median(ci_on)))
cat(sprintf("  linear ramp         : %7.2f\n", median(ci_off)))
cat(sprintf("  C3 effect           : %+7.2f  (%+.1f%%)\n",
            median(ci_on) - median(ci_off),
            100 * (median(ci_on) / median(ci_off) - 1)))

cat(sprintf("\nModelled vs observed SOC by campaign (tC/ha, cross-plot median):\n"))
cat(sprintf("  %-6s %10s %10s %10s %10s\n", "year", "observed", "C3 ON", "C3 OFF", "C3 effect"))
for (yr in c(1985L, 2006L, 2024L)) {
  o <- gather_obs(on, yr); a <- gather(on, yr); b <- gather(off, yr)
  if (!length(a)) next
  cat(sprintf("  %-6d %10.2f %10.2f %10.2f %+10.2f\n",
              yr, median(o), median(a), median(b), median(a) - median(b)))
}

cat(sprintf("\nGap to observed (modelled - observed, tC/ha, median of per-obs gaps):\n"))
for (yr in c(1985L, 2006L, 2024L)) {
  o <- gather_obs(on, yr); a <- gather(on, yr); b <- gather(off, yr)
  if (!length(a)) next
  cat(sprintf("  %-6d  C3 ON %+8.2f   C3 OFF %+8.2f   shrinkage %+8.2f\n",
              yr, median(a - o), median(b - o), median(a - o) - median(b - o)))
}

ll_on  <- sum(vapply(on,  `[[`, numeric(1), "ll"))
ll_off <- sum(vapply(off, `[[`, numeric(1), "ll"))
cat(sprintf("\nLog-likelihood at the SAME parameters (higher = better):\n"))
cat(sprintf("  C3 ON  %12.2f\n  C3 OFF %12.2f\n  delta  %+12.2f  (%s)\n",
            ll_on, ll_off, ll_on - ll_off,
            if (ll_on > ll_off) "growing-stock shape fits better" else
              "linear ramp fits better"))

# =============================================================================
# sigma_init sweep -- where C3 actually has leverage
# =============================================================================
cat(sprintf("\n=====================================================================\n"))
cat(sprintf("C3 EFFECT vs sigma_init   (%s)\n", MODEL))
cat(sprintf("sigma_init sets J_1917/J_full; below ~0.8 the pre-run rises (physical),\n"))
cat(sprintf("at 1.0 it falls (unphysical). Old posteriors sat at 0.19-0.72.\n"))
cat(sprintf("=====================================================================\n"))
cat(sprintf("  %-10s %11s %11s %11s %13s %11s\n",
            "sigma_init", "C_init ON", "C_init OFF", "C3 effect", "1985 gap ON", "gap OFF"))
for (si in SIGMA_INIT_GRID) {
  a <- run_all(TRUE,  si); b <- run_all(FALSE, si)
  cm <- intersect(names(a), names(b)); a <- a[cm]; b <- b[cm]
  cia <- median(vapply(a, `[[`, numeric(1), "C_init_tot"))
  cib <- median(vapply(b, `[[`, numeric(1), "C_init_tot"))
  oa  <- gather_obs(a, 1985L); ga <- gather(a, 1985L); gb <- gather(b, 1985L)
  cat(sprintf("  %-10.2f %11.2f %11.2f %+11.2f %+13.2f %+11.2f\n",
              si, cia, cib, cia - cib,
              median(ga - oa), median(gb - oa)))
}
cat(sprintf("=====================================================================\n"))
