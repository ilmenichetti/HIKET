# =============================================================================
# init_state_plausibility.R   (2026-08-19)
#
# WHY THIS EXISTS -- and what it REPLACES.
#
# It supersedes the *bound* argument in sigma_init_vs_growing_stock.R. That
# script judged sigma_init against the NFI growing-stock ratio,
# (V_1917/V_1985)^eps = (1400/1775)^0.43 = 0.90. The arithmetic is fine but the
# TARGET was wrong: sigma_init does TWO jobs at once and the NFI constrains only
# one of them.
#
#   job 1  the 1917 LITTER FLUX      J_1917 = J_t0_mean * sigma_init * sigma_input
#          -> the growing-stock record speaks to this directly.
#   job 2  the 1917 SOIL STOCK       the pre-run STARTS at equilibrium of J_1917
#          -> the NFI is silent on whether 1917 soils had caught up with their
#             flux. After a century of slash-and-burn, litter raking and heavy
#             cutting they plausibly had not, and the one-parameter transient
#             init cannot express "high flux, disequilibrated soil".
#
# Using an EQUILIBRIUM-INIT benchmark to judge our own initial state also
# re-imports the very convention the manuscript argues against (Lehtonen 2016,
# Palosuo 2008, Peltoniemi 2004). Decision 2026-08-19 (Lorenzo): RETIRE the flux
# bound for now, run the correlated likelihood without it, and judge the result
# on BIOLOGICAL PLAUSIBILITY of the derived state instead. It may come back.
#
# WHAT THIS SCRIPT REPORTS, per model, at the posterior median:
#   * the 1917 and 1985 litter fluxes it implies (tC/ha/yr)
#   * the 1917 SOC STOCK it implies (tC/ha)          <- the plausibility target
#   * the modelled 1985 stock, and the pre-run rate  (C_1985 - C_1917)/68
#
# HOW C_1917 IS OBTAINED without reimplementing six initialisers: the pre-run
# starts at steady state for J_1917 and ramps to J_1985 along lm$preinit_shape.
# Setting that shape to ZERO holds the flux at J_1917 for all 68 years, so the
# routine returns its own starting state. Identical trick works for both families.
#
# THE TWO CRITERIA (Lorenzo, 2026-08-19):
#   (1) 1917 STOCK. 20-30 tC/ha is too little for a forest. The anchor is to be
#       sourced from the boreal-forest SOC literature as a MINIMUM -- a floor,
#       not a central estimate, so a model that fails it fails conservatively.
#       ⚠ Screen on a COMPARABLE basis: our target is the whole profile to a
#       per-plot z_cap (organic + mineral + modelled deep tail), NOT "0-30 cm"
#       and not mineral-only. STOCK_FLOOR below is a PLACEHOLDER pending that.
#   (2) PRE-RUN RATE, which is independent of the equilibrium confound because
#       it reads the trajectory rather than the initial condition. Context from
#       Korhonen 2024: growing stock rose +5.5 Mm3/yr over 1917-1985 and
#       +21.1 Mm3/yr over 1985-2024, so a pre-run rate ABOVE the observed modern
#       +0.259 tC/ha/yr has soil carbon accumulating faster while its driver grew
#       ~4x slower. ⚠ Not decisive on its own: a depleted soil relaxes fastest
#       early whatever the driver does, so this bounds how large the implied
#       depletion is, it does not by itself refute it.
#
# ⚠ TWO EXTRACTION TRAPS, both already documented and both hit while writing this:
#   * the posterior .rds is in PHYSICAL space, the *_chains_*.rds in SAMPLING space.
#     Running to_original on the former double-transforms (SP1 beta1 0.099 -> 1.104,
#     xi -> 1.2e7). This script reads the CHAINS and transforms, as build_F14 does.
#   * getSample() on the saved list returns each internal DEzs chain's FIRST retained
#     iteration -- a start-up state, not a posterior sample. Extract per sampler with
#     start = 2, which also avoids the 1-in-3 thinning getSample applies to a list.
#   Quantities are computed PER DRAW (a coherent parameter vector) and summarised
#   afterwards -- never from a vector of per-parameter medians, which the stick-
#   breaking fractions and the coupled sigma_input/sigma_init would not respect.
#
# Usage:  Rscript doublechecks/init_state_plausibility.R [N_DRAW]
# =============================================================================

suppressWarnings(suppressMessages({
  a <- commandArgs(trailingOnly = TRUE)
  N_DRAW <- if (length(a) >= 1) as.integer(a[[1]]) else 25L
  library(BayesianTools)
}))
set.seed(2025)

MODELS      <- strsplit(Sys.getenv("HIKET_MODELS", "SP1,TP2,TP3,Yasso07,Yasso15,Yasso20"), ",")[[1]]
N_PRE       <- 68L          # 1917 -> 1985
OBS_RATE    <- 0.259        # observed 1985->2024, balanced set (obs_basis.R)
STOCK_FLOOR <- 40           # PLACEHOLDER, tC/ha whole profile -- to be sourced
DIR_RUNS    <- "Calibration_real_data_transient/runs"
NCORE       <- max(1L, parallel::detectCores() - 1L)

latest_run_id <- function(model) {
  fs <- list.files(DIR_RUNS,
                   pattern = sprintf("^%s_posterior_[0-9]{8}_[0-9]{6}\\.rds$", model))
  if (!length(fs)) return(NA_character_)
  sub(sprintf("^%s_posterior_(.+)\\.rds$", model), "\\1", sort(fs, decreasing = TRUE)[1])
}

# Source a calibration script up to the MCMC launch, then pull the four pieces
# the likelihood was built from. Same harness as equifinality_forecast.R.
setup <- function(M) {
  src <- readLines(sprintf("Calibration_real_data_transient/run_%s_transient_calibration.R", M),
                   warn = FALSE)
  cut <- grep("^t_run <- system.time", src)[1]
  e <- new.env(parent = globalenv())
  invisible(capture.output(suppressMessages(
    source(textConnection(paste(src[seq_len(cut - 1)], collapse = "\n")), local = e))))
  st <- grep("ll_fn <- make_likelihood", src)[1]; op <- 0L; en <- NA_integer_
  for (i in seq(st, length(src))) {
    ch <- strsplit(src[i], "")[[1]]; op <- op + sum(ch == "(") - sum(ch == ")")
    if (op == 0L) { en <- i; break }
  }
  ml <- as.list(str2lang(sub("^\\s*ll_fn\\s*<-\\s*", "",
                             paste(src[seq(st, en)], collapse = "\n"))))[-1]
  for (k in c("to_original", "assemble_params", "compute_xi_mean", "steady_state"))
    assign(paste0(".", k), eval(ml[[k]], envir = e), envir = e)
  e
}

# Total litter flux carried by a plot's bundle, before the two sigmas.
j_t0 <- function(x) {
  if (!is.null(x$J_t0_mean)) return(sum(x$J_t0_mean))
  sum(x$nwl_t0_mean, x$fwl_t0_mean, x$cwl_t0_mean)
}

# One draw -> national means of (C_1917, C_1985). Plots that fail are dropped,
# which is the same treatment the likelihood gives them.
state_for_draw <- function(e, p_free, plots, cbp, lms) {
  mp <- e$.assemble_params(p_free)
  out <- do.call(rbind, parallel::mclapply(plots, function(pid) {
    clim <- cbp[[pid]]; lm <- lms[[pid]]
    xi <- tryCatch(e$.compute_xi_mean(clim, mp), error = function(z) NULL)
    if (is.null(xi)) return(NULL)
    lm0 <- lm; lm0$preinit_shape <- rep(0, N_PRE)   # hold J at J_1917 -> start state
    c17 <- tryCatch(sum(e$.steady_state(mp, lm0, xi)), error = function(z) NA_real_)
    c85 <- tryCatch(sum(e$.steady_state(mp, lm,  xi)), error = function(z) NA_real_)
    if (!is.finite(c17) || !is.finite(c85)) return(NULL)
    c(C_1917 = c17, C_1985 = c85, J_t0 = j_t0(lm))
  }, mc.cores = NCORE))
  if (is.null(out) || !nrow(out)) return(NULL)
  colMeans(out)
}

cat("\n=== Derived 1917 state: is it biologically plausible? ===\n")
cat(sprintf("stock floor (PLACEHOLDER) %.0f tC/ha whole profile | observed modern rate %+.3f tC/ha/yr\n",
            STOCK_FLOOR, OBS_RATE))
cat("NFI context: growing stock +5.5 Mm3/yr over 1917-1985 vs +21.1 Mm3/yr over 1985-2024\n\n")
cat(sprintf("%-8s %-16s %7s %8s %8s %8s %8s %8s   %s\n", "model", "run_id",
            "s_init", "J_1917", "J_1985", "C_1917", "C_1985", "rate", "verdict"))

rows <- list()
for (M in MODELS) {
  rid <- latest_run_id(M)
  if (is.na(rid)) { cat(sprintf("%-8s  (no posterior)\n", M)); next }
  e <- tryCatch(setup(M), error = function(z) NULL)
  if (is.null(e)) { cat(sprintf("%-8s  (setup failed)\n", M)); next }

  ch <- file.path(DIR_RUNS, sprintf("%s_chains_%s.rds", M, rid))
  if (!file.exists(ch)) { cat(sprintf("%-8s  (no chains file)\n", M)); next }
  s_ <- do.call(rbind, lapply(readRDS(ch),
                function(z) getSample(z, parametersOnly = FALSE, start = 2)))
  free  <- names(get("best_x", e))          # order matters for to_original
  plots <- get("plots", e); cbp <- get("climate_by_plot", e); lms <- get("litter_means", e)

  idx <- sample(nrow(s_), min(N_DRAW, nrow(s_)))
  st <- lapply(idx, function(k) {
    pf <- tryCatch(e$.to_original(s_[k, free]), error = function(z) NULL)
    if (is.null(pf)) return(NULL)
    v <- state_for_draw(e, pf, plots, cbp, lms)
    if (is.null(v)) return(NULL)
    c(v, sigma_init = unname(pf["sigma_init"]), sigma_input = unname(pf["sigma_input"]))
  })
  st <- st[!vapply(st, is.null, logical(1))]
  if (!length(st)) { cat(sprintf("%-8s  (no usable draws)\n", M)); next }

  g   <- function(f) vapply(st, f, numeric(1))
  C17 <- g(function(z) z[["C_1917"]]); C85 <- g(function(z) z[["C_1985"]])
  RT  <- (C85 - C17) / N_PRE
  SI  <- g(function(z) z[["sigma_init"]]); SP <- g(function(z) z[["sigma_input"]])
  J17 <- g(function(z) z[["J_t0"]] * z[["sigma_init"]] * z[["sigma_input"]])
  J85 <- g(function(z) z[["J_t0"]] * z[["sigma_input"]])
  q <- function(v) unname(quantile(v, c(.05, .5, .95)))

  verdict <- paste(
    if (median(C17) < STOCK_FLOOR) sprintf("STOCK %.0f < floor", median(C17)) else "stock ok",
    if (median(RT) > OBS_RATE) sprintf("| rate %.1fx obs", median(RT) / OBS_RATE) else "| rate ok")

  cat(sprintf("%-8s %-16s %7.3f %8.2f %8.2f %8.1f %8.1f %+8.3f   %s\n",
              M, rid, median(SI), median(J17), median(J85),
              median(C17), median(C85), median(RT), verdict))
  cat(sprintf("%-8s %-16s %7s %8s %8s %8s %8s %8s   [%d draws | C_1917 %.1f-%.1f | rate %+.3f..%+.3f]\n",
              "", "", "", "", "", "", "", "", length(st),
              q(C17)[1], q(C17)[3], q(RT)[1], q(RT)[3]))

  rows[[M]] <- data.frame(model = M, run_id = rid, n_draw = length(st),
                          sigma_init = median(SI), sigma_input = median(SP),
                          J_1917 = median(J17), J_1985 = median(J85),
                          C_1917 = median(C17), C_1917_q05 = q(C17)[1], C_1917_q95 = q(C17)[3],
                          C_1985 = median(C85),
                          prerun_rate = median(RT), rate_q05 = q(RT)[1], rate_q95 = q(RT)[3])
}

if (length(rows)) {
  df <- do.call(rbind, rows)
  saveRDS(df, "doublechecks/init_state_plausibility.rds")
  cat("\nsaved -> doublechecks/init_state_plausibility.rds\n")
}
cat("\nNB sigma_init and sigma_input COUPLE through the 1917 flux; do not constrain them\n",
    "   independently. See memory sigma-input-physical-bounds.\n", sep = "")
