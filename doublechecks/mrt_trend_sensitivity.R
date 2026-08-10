# =============================================================================
# mrt_trend_sensitivity.R   (2026-08-10)
#
# Question: how much of the 2006->2024 trend gap is model SPEED, as opposed to
# the shape of the litter driver?
#
# NOT a test of "should MRT equal the ICBM anchor". Bulk MRT in a forest is an
# EMERGENT quantity and the arable Ultuna anchor is not a target -- C1 leaves the
# slow rate free precisely so its displacement is a transferability diagnostic.
# So this sweeps MRT and reports the TREND RESPONSE, leaving where MRT should sit
# to physical judgement.
#
# Method, per model, at the posterior MEDIAN parameter vector:
#   - scale the effective decomposition rates by 1/f (implemented as a uniform
#     scaling of xi, which multiplies the rates in every model) -> MRT scales ~f
#   - RE-SOLVE sigma_input so the modelled SOC level still matches observed.
#     SOC is exactly linear in sigma_input (it scales the inputs, and the pre-run
#     and forward system are both linear), so one rescale is exact.
#   - report the ACHIEVED bulk MRT measured from the output (C/J), not the
#     nominal f -- so the x-axis stays honest even where a model does not apply
#     xi to every pool (e.g. a climate-independent humus rate).
#
# Holding the stock fit fixed is what isolates the question: any trend change is
# then attributable to turnover speed alone, not to a better/worse level fit.
#
# Usage:  Rscript doublechecks/mrt_trend_sensitivity.R [MODEL] [f1,f2,...]
# =============================================================================

suppressWarnings(suppressMessages({
  args  <- commandArgs(trailingOnly = TRUE)
  MODEL <- if (length(args) >= 1) args[[1]] else "SP1"
  GRID  <- if (length(args) >= 2) as.numeric(strsplit(args[[2]], ",")[[1]])
           else c(0.7, 1.0, 1.4, 1.8, 2.2, 2.8, 3.5)
  library(BayesianTools)
}))
stopifnot(MODEL %in% c("SP1","TP2","TP3","Yasso07","Yasso15","Yasso20"))
set.seed(2025)

script_path <- file.path("Calibration_real_data_transient",
                         sprintf("run_%s_transient_calibration.R", MODEL))
src   <- readLines(script_path, warn = FALSE)
cutix <- grep("^t_run <- system.time\\(\\{", src)[1]
if (is.na(cutix)) stop("Could not locate MCMC launch in ", script_path)
message(sprintf("[%s] sourcing setup ...", MODEL))
e <- new.env(parent = globalenv())
source(textConnection(paste(src[seq_len(cutix - 1L)], collapse = "\n")), local = e)

# --- engine bindings (args may be expressions, e.g. Yasso20's monthly x12) ----
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
assemble_params <- argof("assemble_params")
compute_xi      <- argof("compute_xi")
compute_xi_mean <- argof("compute_xi_mean")
steady_state    <- argof("steady_state")
run_model       <- argof("run_model")
STEADY_N        <- argof("steady_state_n", NULL)

plots           <- get("plots",           e)
climate_by_plot <- get("climate_by_plot", e)
inputs_by_plot  <- get("inputs_by_plot",  e)
litter_means    <- get("litter_means",    e)
obs_meta        <- get("obs_meta",        e)
sigma_obs_fixed <- get("sigma_obs_fixed", e)
best_x          <- get("best_x",          e)
ll_fn           <- get("ll_fn",           e)
log_jacobian    <- get("log_jacobian",    e)

SOC_obs_all <- get("SOC_obs_all", e)
obs_years   <- lapply(plots, function(pid)
  SOC_obs_all$year[as.character(SOC_obs_all$plot_id) == pid])
names(obs_years) <- plots

# --- forward pass: mirrors calibration_engine_transient.R make_likelihood() ---
# XI_SCALE divides xi -> divides every rate xi multiplies -> multiplies MRT.
# Same default rule as the engine: log-normal unless HIKET_LOGNORMAL_LIK=0.
LOGNORM <- !identical(Sys.getenv("HIKET_LOGNORMAL_LIK"), "0")
# Honour the total-sigma override so the faithfulness guard stays meaningful.
.sig_tot <- suppressWarnings(as.numeric(Sys.getenv("HIKET_SIGMA_TOTAL", NA)))
if (is.finite(.sig_tot)) { sigma_obs_fixed <- .sig_tot
  message(sprintf("[error model] total sigma overridden: %.3f", .sig_tot)) }
message("[error model] ", if (LOGNORM) "LOG-NORMAL (default)" else "multiplicative normal")

# xi comes in two shapes across the six models: a numeric vector (SP1/TP2/TP3,
# Yasso07) or a named LIST of per-pool-group modifiers (Yasso15/20:
# xi_awe, xi_n, xi_h). Scaling every element scales every rate xi multiplies,
# which for Yasso15/20 includes the humus pool -- so the MRT scaling stays
# uniform across pools rather than sparing the slow one.
scale_xi <- function(x, f) if (is.list(x)) lapply(x, function(z) z / f) else x / f

# Raw litter flux per plot, also two shapes: a scalar J_total_mean (simple
# models) or AWEN component vectors nwl/fwl/cwl_mean (Yasso).
raw_J <- function(lm) {
  if (!is.null(lm$J_total_mean)) return(unname(lm$J_total_mean))
  comp <- c("nwl_mean", "fwl_mean", "cwl_mean")
  sum(unlist(lm[intersect(comp, names(lm))]))
}

forward_one <- function(pid, model_params, xi_scale, sigma_init) {
  clim <- climate_by_plot[[pid]]; inputs <- inputs_by_plot[[pid]]
  lm <- litter_means[[pid]];      meta   <- obs_meta[[pid]]
  if (any(is.na(meta$idx))) return(NULL)

  xi_array <- tryCatch(compute_xi(clim, model_params), error = function(z) NULL)
  if (is.null(xi_array)) return(NULL)
  n_ss  <- if (is.null(STEADY_N)) nrow(clim) else min(STEADY_N, nrow(clim))
  xi_ss <- tryCatch(compute_xi_mean(clim[seq_len(n_ss), , drop = FALSE], model_params),
                    error = function(z) NULL)
  if (is.null(xi_ss)) return(NULL)

  xi_array <- scale_xi(xi_array, xi_scale)
  xi_ss    <- scale_xi(xi_ss,    xi_scale)

  C_init <- tryCatch(steady_state(model_params, lm, xi_ss), error = function(z) NULL)
  if (is.null(C_init) || any(!is.finite(C_init)) || any(C_init < 0)) return(NULL)
  run_out <- tryCatch(run_model(inputs, model_params, C_init, xi_array),
                      error = function(z) NULL)
  if (is.null(run_out)) return(NULL)

  SOC_hat <- run_out$total_soc[meta$idx]
  if (any(!is.finite(SOC_hat)) || any(SOC_hat <= 0)) return(NULL)

  # Mirror the engine's ERROR-MODEL SWITCH exactly (calibration_engine_transient.R).
  # Log-normal is the DEFAULT since 2026-08-07 and uses sigma_obs_fixed * infl
  # DIRECTLY -- it neither scales the sd by SOC_hat nor adds the is_first
  # sigma_init quadrature term that the multiplicative-normal branch uses.
  infl <- if (!is.null(meta$sigma_infl)) meta$sigma_infl else 1
  ll <- if (LOGNORM) {
    sum(dnorm(log(meta$soc_obs), mean = log(SOC_hat),
              sd = sigma_obs_fixed * infl, log = TRUE))
  } else {
    sd_vec <- if (!is.null(meta$is_first))
      ifelse(meta$is_first, SOC_hat * sqrt(sigma_obs_fixed^2 + sigma_init^2),
                            SOC_hat * sigma_obs_fixed)
    else SOC_hat * sigma_obs_fixed
    sum(dnorm(meta$soc_obs, mean = SOC_hat, sd = sd_vec * infl, log = TRUE))
  }

  list(SOC_hat = SOC_hat, soc_obs = meta$soc_obs, years = obs_years[[pid]],
       J_raw   = raw_J(lm), ll = ll)
}

run_all <- function(p_free, xi_scale) {
  mp  <- assemble_params(p_free)
  res <- lapply(plots, forward_one, model_params = mp, xi_scale = xi_scale,
                sigma_init = unname(p_free["sigma_init"]))
  names(res) <- plots
  res[!vapply(res, is.null, logical(1))]
}

# --- faithfulness guard: reconstruct the engine's own ll at best_x ------------
p_best <- to_original(best_x)
base   <- run_all(p_best, 1)
ll_eng <- ll_fn(best_x)
ll_mine <- sum(vapply(base, `[[`, numeric(1), "ll")) + log_jacobian(best_x, p_best)
d <- abs(ll_eng - ll_mine)
cat(sprintf("\nFaithfulness: engine ll = %.4f | reconstructed = %.4f | |diff| = %.2e\n",
            ll_eng, ll_mine, d))
if (d > 1e-6) { cat("MISMATCH -- not reporting.\n"); quit(status = 1L) }
cat("Reconstructed forward == engine forward. Proceeding.\n\n")

# --- posterior MEDIAN parameter vector ---------------------------------------
rid <- sub(sprintf("^%s_posterior_(.+)\\.rds$", MODEL), "\\1",
           sort(list.files("Calibration_real_data_transient/runs",
                pattern = sprintf("^%s_posterior_[0-9]{8}_[0-9]{6}\\.rds$", MODEL)),
                decreasing = TRUE)[1])
post   <- readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_%s.rds", MODEL, rid))
smp    <- getSample(post)
p_med  <- p_best
for (nm in intersect(names(p_med), colnames(smp))) p_med[nm] <- median(smp[, nm])
message(sprintf("[%s] posterior median from run %s", MODEL, rid))

# --- helpers ------------------------------------------------------------------
paired <- function(res, y0, y1) {
  ch <- vapply(res, function(r) {
    a <- r$SOC_hat[r$years == y0]; b <- r$SOC_hat[r$years == y1]
    if (length(a) != 1 || length(b) != 1) return(NA_real_)
    (b - a) / (y1 - y0)
  }, numeric(1))
  mean(ch, na.rm = TRUE)
}
paired_obs <- function(res, y0, y1) {
  ch <- vapply(res, function(r) {
    a <- r$soc_obs[r$years == y0]; b <- r$soc_obs[r$years == y1]
    if (length(a) != 1 || length(b) != 1) return(NA_real_)
    (b - a) / (y1 - y0)
  }, numeric(1))
  mean(ch, na.rm = TRUE)
}
soc_med <- function(res) median(unlist(lapply(res, `[[`, "SOC_hat")))
obs_med <- function(res) median(unlist(lapply(res, `[[`, "soc_obs")))
bulk_mrt <- function(res, si) {
  median(vapply(res, function(r) median(r$SOC_hat) / (si * r$J_raw), numeric(1)))
}

target <- obs_med(base)
cat(sprintf("Target SOC level (observed median) = %.2f tC/ha\n\n", target))

# --- sweep --------------------------------------------------------------------
out <- do.call(rbind, lapply(GRID, function(f) {
  p <- p_med
  r1 <- run_all(p, f)                                  # pass 1
  if (!length(r1)) return(NULL)
  # SOC is exactly linear in sigma_input -> one analytic rescale hits the target
  p["sigma_input"] <- p["sigma_input"] * target / soc_med(r1)
  r2 <- run_all(p, f)
  if (!length(r2)) return(NULL)
  data.frame(f = f, sigma_input = unname(p["sigma_input"]),
             MRT = bulk_mrt(r2, unname(p["sigma_input"])),
             SOC_med = soc_med(r2),
             trend_0624 = paired(r2, 2006, 2024),
             trend_8524 = paired(r2, 1985, 2024),
             obs_0624   = paired_obs(r2, 2006, 2024),
             obs_8524   = paired_obs(r2, 1985, 2024),
             n = length(r2))
}))

cat(sprintf("=== %s : MRT sweep at fixed SOC level ===\n", MODEL))
print(out, row.names = FALSE, digits = 3)
cat("\ntrend_* = mean paired stock change (tC/ha/yr); obs_* = same plots, observed.\n")
cat("MRT is MEASURED from the output (C/J), not the nominal scaling f.\n")

dir.create("doublechecks/mrt_sweep", showWarnings = FALSE)
write.csv(out, sprintf("doublechecks/mrt_sweep/%s_mrt_sweep.csv", MODEL), row.names = FALSE)
cat(sprintf("\nwrote doublechecks/mrt_sweep/%s_mrt_sweep.csv\n", MODEL))
