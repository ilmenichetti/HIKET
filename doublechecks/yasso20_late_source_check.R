# =============================================================================
# yasso20_late_source_check.R -- why does Yasso20's projection turn to a small
# source in the 2070s while the soil sits ~35% below its equilibrium? (2026-09-23)
# -----------------------------------------------------------------------------
# Hypothesis: the national litter input peaked around 2006 and declined, so in
# 2024 the fast pools (A, W, E, N) still hold carbon from the higher input of the
# preceding years and drain towards the lower 2024 input, while the slow humus pool
# gains only slowly. Three checks, on the production engines and draws:
#   (a) national mean litter input over the record (balanced plots)
#   (b) the baseline projection split into fast pools and humus
#   (c) the projection with litter held at its 2006-2024 MEAN instead of the 2024
#       value: if the hypothesis holds, the late decline shrinks or disappears
#   RESULT (2026-09-23): hypothesis FALSIFIED -- 2024 input is above the 2006-2024
#   mean and the fast pools do not drain. The late decline is the PHASE of the
#   recycled 20-year climate (decadal rates below); report rates per 20-yr cycle.
# Sources lines 1-176 of run_Yasso20_transient_predictive.R (config, loading,
# parameter assembly, engine wrappers; no file is written there).
#
# Run from repo root:  Rscript doublechecks/yasso20_late_source_check.R [N_DRAWS]
# =============================================================================
suppressMessages(source("manuscript/figures/run_ids.R"))
RUN_ID <- RID[["Yasso20"]]; cat("Yasso20 RUN_ID:", RUN_ID, "\n")
N_CHECK <- { a <- commandArgs(TRUE); if (length(a)) as.integer(a[1]) else 30L }
src <- readLines("Calibration_real_data_transient/run_Yasso20_transient_predictive.R", warn = FALSE)
cut <- grep("^# 4\\.  Posterior predictive simulation", src)[1] - 2L
e <- new.env()
local({ commandArgs <- function(trailingOnly = FALSE) if (trailingOnly) RUN_ID else c("R", RUN_ID)
        eval(parse(text = src[seq_len(cut)]), envir = environment()) }, envir = e)
suppressMessages(source("manuscript/figures/obs_basis.R"))
with(e, {
  set.seed(2025); draw_idx <- sample(nrow(posterior_phys), N_PP_DRAWS)   # as in Section 4
  plots <- as.character(intersect(plots_real, balanced_plots(obs_meta)))
  lit   <- function(df) rowSums(df[setdiff(names(df), c("plot_id", "year", "precip"))])

  # (a) national mean litter input, balanced plots
  J <- sapply(plots, function(p) lit(inputs_by_plot[[p]])); yrs <- inputs_by_plot[[plots[1]]]$year
  Jm <- rowMeans(J)
  cat(sprintf("\n(a) mean litter input (raw product, tC/ha/yr): 1990 %.2f | 2006 %.2f | 2015 %.2f | 2024 %.2f | 2006-2024 mean %.2f\n",
              Jm[yrs == 1990], Jm[yrs == 2006], Jm[yrs == 2015], Jm[yrs == 2024], mean(Jm[yrs >= 2006])))

  one <- function(d, pid, mean_input) {
    mp <- assemble_model_params(posterior_phys[draw_idx[d], ])
    clim <- climate_by_plot[[pid]]; inputs <- inputs_by_plot[[pid]]; lm <- litter_means[[pid]]
    n_ss <- min(STEADY_STATE_YEARS, nrow(clim))
    xs <- compute_xi_mean_yasso20_engine(clim[seq_len(n_ss), , drop = FALSE], mp)
    hist <- yasso20_run_engine(inputs, mp, steady_state_yasso20_engine(mp, lm, xs), clim)
    last <- max(hist$year); py <- seq(last + 1L, last + PROJ_YEARS)
    cr <- tail(clim, RECYCLE_YEARS * 12L)
    cp <- cr[((seq_len(PROJ_YEARS * 12L) - 1L) %% (RECYCLE_YEARS * 12L)) + 1L, , drop = FALSE]
    cp$year <- rep(py, each = 12L); cp$month <- rep(1:12, PROJ_YEARS); rownames(cp) <- NULL
    ip <- tail(inputs, 1L)[rep(1L, PROJ_YEARS), , drop = FALSE]; ip$year <- py; rownames(ip) <- NULL
    if (mean_input) {                                     # (c) litter held at the 2006-2024 mean
      lc <- setdiff(names(ip), c("plot_id", "year", "precip"))
      ip[lc] <- as.list(colMeans(inputs[inputs$year >= 2006, lc]))
    }
    pr <- yasso20_run_engine(ip, mp, attr(hist, "C_final"), cp)
    data.frame(year = c(last, pr$year),
               fast = c(sum(unlist(hist[nrow(hist), c("A","W","E","N")])), pr$A + pr$W + pr$E + pr$N),
               hum  = c(hist$H[nrow(hist)], pr$H))
  }
  agg <- function(mean_input) {
    r <- do.call(rbind, lapply(seq_len(N_CHECK), function(d)
      do.call(rbind, parallel::mclapply(plots, function(p) cbind(one(d, p, mean_input), draw = d),
                                        mc.cores = max(1L, parallel::detectCores() - 1L)))))
    aggregate(cbind(fast, hum) ~ year, r, mean)
  }
  rate <- function(a, y0, y1) with(a, ((fast + hum)[year == y1] - (fast + hum)[year == y0]) / (y1 - y0))
  b <- agg(FALSE)
  cat(sprintf("\n(b) baseline (2024 input), %d draws x %d plots, national mean tC/ha:\n", N_CHECK, length(plots)))
  for (y in c(2024, 2034, 2054, 2074, 2084)) with(b[b$year == y, ], cat(sprintf("    %d  fast %.1f  humus %.1f  total %.1f\n", y, fast, hum, fast + hum)))
  cat(sprintf("    change 2024-2084: fast %+.1f, humus %+.1f | total rate 2025-34 %+.3f, 2075-84 %+.3f\n",
              diff(b$fast[b$year %in% c(2024, 2084)]), diff(b$hum[b$year %in% c(2024, 2084)]), rate(b, 2024, 2034), rate(b, 2074, 2084)))
  cat("    decadal total rate (replayed climate years in brackets):\n")
  for (y0 in seq(2024, 2074, 10)) cat(sprintf("      %d-%d [%d-%d]  %+.3f\n", y0 + 1, y0 + 10,
      2005 + ((y0 - 2024) %% 20), 2014 + ((y0 - 2024) %% 20), rate(b, y0, y0 + 10)))
  cat(sprintf("    per 20-year climate cycle: 2025-44 %+.3f | 2045-64 %+.3f | 2065-84 %+.3f\n",
              rate(b, 2024, 2044), rate(b, 2044, 2064), rate(b, 2064, 2084)))
  c_ <- agg(TRUE)
  cat(sprintf("\n(c) litter held at its 2006-2024 mean: total rate 2025-34 %+.3f, 2075-84 %+.3f | fast change 2024-2084 %+.1f\n",
              rate(c_, 2024, 2034), rate(c_, 2074, 2084), diff(c_$fast[c_$year %in% c(2024, 2084)])))
})
