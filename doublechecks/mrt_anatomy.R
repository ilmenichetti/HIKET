# =============================================================================
# mrt_anatomy.R   (2026-08-10)
#
# Two questions, one harness:
#
# (A) REACHABILITY. The MRT sweep scaled xi artificially. The calibration can
#     only move xi through the climate parameters (beta1, beta2, gamma), which
#     carry informative Tier-1 priors. How far would they have to move -- in
#     PRIOR SIGMA -- to halve xi (i.e. double MRT)? If that is many sigma, the
#     slow-MRT solution is unreachable and the reported MRT is a design
#     consequence, not a calibration outcome.
#
# (B) ANATOMY. Yasso's decomposition rates are FIXED. So what actually sets its
#     bulk MRT? Decompose:
#       MRT_bulk  =  MRT_intrinsic (fixed a-vector + fitted fractions)  /  xi
#     and report where the carbon actually sits at steady state. If almost no
#     carbon reaches the slow humus pool, bulk MRT is governed by the fast pools
#     and is a structural constant, not something the data chose.
#
# Usage:  Rscript doublechecks/mrt_anatomy.R [MODEL]
# =============================================================================

suppressWarnings(suppressMessages({
  args  <- commandArgs(trailingOnly = TRUE)
  MODEL <- if (length(args) >= 1) args[[1]] else "Yasso15"
  library(BayesianTools)
}))
set.seed(2025)

script_path <- file.path("Calibration_real_data_transient",
                         sprintf("run_%s_transient_calibration.R", MODEL))
src   <- readLines(script_path, warn = FALSE)
cutix <- grep("^t_run <- system.time\\(\\{", src)[1]
message(sprintf("[%s] sourcing setup ...", MODEL))
e <- new.env(parent = globalenv())
source(textConnection(paste(src[seq_len(cutix - 1L)], collapse = "\n")), local = e)

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
STEADY_N        <- argof("steady_state_n", NULL)

plots <- get("plots", e); climate_by_plot <- get("climate_by_plot", e)
litter_means <- get("litter_means", e); best_x <- get("best_x", e)

# posterior median parameter vector
rid <- sub(sprintf("^%s_posterior_(.+)\\.rds$", MODEL), "\\1",
           sort(list.files("Calibration_real_data_transient/runs",
                pattern = sprintf("^%s_posterior_[0-9]{8}_[0-9]{6}\\.rds$", MODEL)),
                decreasing = TRUE)[1])
smp   <- getSample(readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_%s.rds",
                                   MODEL, rid)))
p_med <- to_original(best_x)
for (nm in intersect(names(p_med), colnames(smp))) p_med[nm] <- median(smp[, nm])
mp <- assemble_params(p_med)

cat(sprintf("\n=====================  %s  (run %s)  =====================\n", MODEL, rid))

# --- the fixed kinetic constants ---------------------------------------------
FIXED <- if (exists("FIXED_RATE_NAMES", e)) get("FIXED_RATE_NAMES", e) else character(0)
cat("\n--- fixed kinetic constants (NOT calibrated) ---\n")
if (length(FIXED)) {
  for (nm in FIXED) if (!is.na(mp[nm])) cat(sprintf("  %-10s %12.6g", nm, mp[nm]),
                                            sprintf("   (1/rate = %.1f yr)\n", 1/mp[nm]))
} else cat("  (none declared)\n")

# --- xi at the posterior median ----------------------------------------------
xis <- unlist(lapply(plots, function(pid) {
  x <- tryCatch(compute_xi(climate_by_plot[[pid]], mp), error = function(z) NULL)
  if (is.null(x)) return(NULL)
  if (is.list(x)) vapply(x, function(z) median(z), numeric(1)) else median(x)
}))
cat("\n--- climate modifier xi at posterior median ---\n")
if (!is.null(names(xis)) && length(unique(names(xis))) > 1) {
  for (g in unique(names(xis)))
    cat(sprintf("  %-8s median %.4f   [%.3f, %.3f]\n", g,
                median(xis[names(xis)==g]),
                quantile(xis[names(xis)==g], .05), quantile(xis[names(xis)==g], .95)))
} else {
  cat(sprintf("  xi median %.4f   [%.3f, %.3f]\n", median(xis),
              quantile(xis,.05), quantile(xis,.95)))
}
xi_bar <- median(xis)
cat(sprintf("  => intrinsic MRT (xi=1) is  1/xi = %.2f x  the realised bulk MRT\n", 1/xi_bar))

# --- where does the carbon sit at steady state? -------------------------------
pool <- lapply(plots, function(pid) {
  clim <- climate_by_plot[[pid]]
  n_ss <- if (is.null(STEADY_N)) nrow(clim) else min(STEADY_N, nrow(clim))
  xs <- tryCatch(compute_xi_mean(clim[seq_len(n_ss), , drop=FALSE], mp), error=function(z) NULL)
  if (is.null(xs)) return(NULL)
  tryCatch(steady_state(mp, litter_means[[pid]], xs), error = function(z) NULL)
})
pool <- pool[!vapply(pool, is.null, logical(1))]
if (length(pool)) {
  M <- do.call(rbind, lapply(pool, function(v) v / sum(v)))
  cat("\n--- steady-state pool shares (median across plots) ---\n")
  sh <- apply(M, 2, median)
  nm <- if (!is.null(colnames(M))) colnames(M) else paste0("pool", seq_along(sh))
  for (i in seq_along(sh)) cat(sprintf("  %-8s %6.2f %%\n", nm[i], 100*sh[i]))
}

# --- (A) how far must the climate params move to halve xi? --------------------
cat("\n--- reachability: displacement needed to HALVE xi (= double MRT) ---\n")
psp <- if (exists("PRIOR_SPECS", e)) get("PRIOR_SPECS", e) else
       if (exists("param_spec", e)) get("param_spec", e) else NULL
clim_pars <- intersect(c("beta1","beta2","gamma"), names(p_med))

med_xi_at <- function(p) {
  m <- assemble_params(p)
  v <- unlist(lapply(plots[seq_len(min(120, length(plots)))], function(pid) {
    x <- tryCatch(compute_xi(climate_by_plot[[pid]], m), error = function(z) NULL)
    if (is.null(x)) return(NULL)
    if (is.list(x)) median(unlist(x)) else median(x)
  }))
  if (!length(v)) return(NA_real_) else median(v)
}
base_xi <- med_xi_at(p_med)
cat(sprintf("  baseline median xi = %.4f ; target = %.4f\n", base_xi, base_xi/2))

for (nm in clim_pars) {
  f <- function(val) { p <- p_med; p[nm] <- val; med_xi_at(p) - base_xi/2 }
  lo <- p_med[[nm]] - 50*abs(p_med[[nm]]) - 5
  hi <- p_med[[nm]] + 50*abs(p_med[[nm]]) + 5
  sol <- tryCatch(uniroot(f, c(lo, hi))$root, error = function(z) NA_real_)
  # prior sd on this parameter, if recoverable
  sd_nm <- NA_real_
  if (!is.null(psp) && !is.null(psp[[nm]]))
    sd_nm <- tryCatch(psp[[nm]]$sd, error = function(z) NA_real_)
  cat(sprintf("  %-7s posterior %+.5g -> needs %+.5g   (delta %+.5g%s)\n",
              nm, p_med[[nm]], sol, sol - p_med[[nm]],
              if (is.finite(sd_nm)) sprintf(" = %.1f prior SD", abs(sol-p_med[[nm]])/sd_nm) else ""))
}
cat("\nNOTE: single-parameter moves; the calibration could combine them, so these\n")
cat("are upper bounds on the displacement needed from any ONE parameter.\n")
