# =============================================================================
# scenario_switch_tests.R   (2026-10-13)
#
# Local tests for the scenario switches (Calibration_real_data_transient/
# scenario_switches.R) BEFORE any Roihu launch. Each test evaluates the REAL
# likelihood (the calibration script sourced up to the MCMC launch) at a fixed
# point in sampling space.
#
#   baseline  ll at best_x with every switch off, written to a file BEFORE the
#             engine is edited
#   noop      the same after the edit: must be IDENTICAL (the switches are inert)
#   held      the climate-held and litter-held hooks against an INDEPENDENT
#             construction: the same likelihood built with the switch OFF but
#             with compute_xi / run_model patched to apply the scenario definition
#             of doublechecks/attribution_decomposition.R (pre-run climate modifier
#             repeated every year; litter at its 1985-89 level). Must be IDENTICAL.
#
# Usage:  Rscript doublechecks/scenario_switch_tests.R baseline [MODELS...]
#         Rscript doublechecks/scenario_switch_tests.R noop     [MODELS...]
#         Rscript doublechecks/scenario_switch_tests.R held     [MODELS...]
# Setup side effects (input bundle, sanity PNG) are removed at the end.
# =============================================================================

setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
a     <- commandArgs(trailingOnly = TRUE)
MODE  <- if (length(a)) a[[1]] else "noop"
MODELS <- if (length(a) >= 2) a[-1] else c("SP1", "Yasso15")
BASE  <- "doublechecks/scenario_switch_baseline.rds"
N_TEST_PLOTS <- 40L                                  # enough to touch every code path

snap <- function() c(list.files("Data/model_inputs", full.names = TRUE),
                     list.files("Calibration_real_data_transient/diagnostics", recursive = TRUE,
                                full.names = TRUE, include.dirs = TRUE))
before <- snap()

setup_model <- function(M) {
  src <- readLines(sprintf("Calibration_real_data_transient/run_%s_transient_calibration.R", M), warn = FALSE)
  cut <- grep("^t_run <- system.time", src)[1]
  e <- new.env(parent = globalenv())
  suppressMessages(capture.output(
    source(textConnection(paste(src[seq_len(cut - 1)], collapse = "\n")), local = e)))
  st <- grep("ll_fn <- make_likelihood", src)[1]; op <- 0L; en <- NA_integer_
  for (i in seq(st, length(src))) {
    ch <- strsplit(src[i], "")[[1]]; op <- op + sum(ch == "(") - sum(ch == ")")
    if (op == 0L) { en <- i; break }
  }
  e$.ml <- as.list(str2lang(sub("^\\s*ll_fn\\s*<-\\s*", "", paste(src[seq(st, en)], collapse = "\n"))))[-1]
  e
}
# Build a likelihood on a plot subset, optionally with patched engine functions.
build_ll <- function(e, patch = list()) {
  args <- lapply(e$.ml, eval, envir = e)
  args$plots <- head(args$plots, N_TEST_PLOTS)
  for (k in names(patch)) args[[k]] <- patch[[k]]
  do.call(make_likelihood, args)
}
x0 <- function(e) get("best_x", e)

res <- list()
for (M in MODELS) {
  e  <- setup_model(M)
  x  <- x0(e)
  if (MODE %in% c("baseline", "noop")) {
    ll <- build_ll(e)(x)
    res[[M]] <- ll
    cat(sprintf("%-8s ll(best_x) = %.10f\n", M, ll))
  }
  if (MODE == "held") {
    stopifnot(exists("hiket_xi_held"), exists("hiket_litter_held"))
    cbp <- args_cbp <- eval(e$.ml$climate_by_plot, envir = e)
    ibp <- eval(e$.ml$inputs_by_plot, envir = e)
    lms <- eval(e$.ml$litter_means, envir = e)
    nss <- eval(e$.ml$steady_state_n, envir = e)
    cxi <- eval(e$.ml$compute_xi, envir = e); cxm <- eval(e$.ml$compute_xi_mean, envir = e)
    rmd <- eval(e$.ml$run_model, envir = e)
    pid_of <- function(obj, lst) { for (p in names(lst)) if (identical(lst[[p]], obj)) return(p); stop("plot not found") }
    # independent constructions: the scenario definitions of attribution_decomposition.R
    xi_rep <- function(z, n) if (is.list(z)) lapply(z, function(v) rep(unname(v)[1], n)) else rep(unname(z)[1], n)
    cxi_held <- function(clim, mp) {
      pid <- pid_of(clim, cbp); n <- length(if (is.list(xa <- cxi(clim, mp))) xa[[1]] else xa)
      xi_rep(cxm(hiket_ss_clim(pid, clim, min(nss, nrow(clim))), mp), n)
    }
    LIT <- c(paste0("nwl_", c("A","W","E","N")), paste0("fwl_", c("A","W","E","N")), paste0("cwl_", c("A","W","E","N")))
    rmd_held <- function(inputs, mp, C, xi) {
      lm <- lms[[pid_of(inputs, ibp)]]; r <- inputs
      if ("J_total" %in% names(r)) r$J_total <- lm$J_t0_mean
      else for (s in c("nwl", "fwl", "cwl")) { v <- lm[[paste0(s, "_t0_mean")]]
        for (k in c("A","W","E","N")) r[[paste0(s, "_", k)]] <- unname(v[[paste0(s, "_", k)]]) }
      rmd(r, mp, C, xi)
    }
    for (sw in c("climate", "litter")) {
      flag <- if (sw == "climate") ".hiket_climate_held" else ".hiket_litter_held"
      assign(flag, TRUE, envir = globalenv());  llA <- build_ll(e)(x)
      assign(flag, FALSE, envir = globalenv())
      llB <- build_ll(e, if (sw == "climate") list(compute_xi = cxi_held) else list(run_model = rmd_held))(x)
      ll0 <- build_ll(e)(x)
      ok  <- isTRUE(all.equal(llA, llB, tolerance = 1e-12)) && abs(llA - ll0) > 1e-6
      cat(sprintf("%-8s %-7s held: hook %.10f | independent %.10f | off %.10f  %s\n",
                  M, sw, llA, llB, ll0, if (ok) "PASS" else "FAIL"))
      res[[paste(M, sw)]] <- ok
    }
  }
}

if (MODE == "baseline") { saveRDS(res, BASE); cat("baseline saved ->", BASE, "\n") }
if (MODE == "noop") {
  b <- readRDS(BASE)
  for (M in names(res)) cat(sprintf("%-8s no-op: %s (baseline %.10f, now %.10f)\n", M,
      if (identical(b[[M]], res[[M]])) "PASS (identical)" else "FAIL", b[[M]], res[[M]]))
}
if (MODE == "held" && !all(unlist(res))) stop("held-switch test FAILED")

new <- setdiff(snap(), before)
unlink(new[!dir.exists(new)]); unlink(new[dir.exists(new)], recursive = TRUE)
cat(sprintf("removed %d setup side-effect files\n", length(new)))
