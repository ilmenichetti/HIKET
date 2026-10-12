# =============================================================================
# input_to_hold_sink.R -- how much litter input would hold a given sink (2026-10-12)
# -----------------------------------------------------------------------------
# Question (Lorenzo): given our assumptions, how much input would we need to keep the
# sink of the past for another 10 and 20 years?
#
# EXACT, no new runs: every model is linear in its input, and the predictive stage
# already projects each draw twice from 2025 -- litter held at 2024 (total_soc) and the
# same x MULT (total_soc_up, MULT = 1.2). Over any window the sink is then linear in the
# input multiplier m:
#     sink(m) = sink(1) + (m - 1)/(MULT - 1) * [sink(MULT) - sink(1)]
# so the m that holds a target sink is solved per draw.
# Targets: the observed 1985-2024 rate (+0.26), the observed 2006-2024 rate (+0.12), and
# each model's own modelled 1985-2024 rate. Climate: the recent 20 years recycled (no
# further warming), so these are LOWER bounds under warming. The 10-yr window replays
# the 2005-2014 climate only; the 20-yr window is one whole climate cycle.
# Basis: balanced plot set (n = 310), unweighted national mean, as every rate in the paper.
#
# Run from repo root:  Rscript doublechecks/input_to_hold_sink.R
# =============================================================================
suppressMessages({ source("manuscript/figures/run_ids.R"); source("manuscript/figures/obs_basis.R") })
cat("RUN_IDs:", paste(names(RID), RID, collapse = " | "), "\n")
MULT <- 1.2
TARGETS <- c(obs_1985_2024 = 0.26, obs_2006_2024 = 0.12)
q <- function(x) sprintf("%.2f [%.2f, %.2f]", median(x), quantile(x, .05), quantile(x, .95))

ir <- read.csv("Data/model_inputs/input_raw_monthly.csv")
lc <- grep("^C_(nwl|fwl|cwl)_", names(ir), value = TRUE)

res <- list()
for (m in FIG_MODELS) {
  pk <- readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_predictive_%s.rds", m, RID[[m]]))
  om <- readRDS(sprintf("Data/model_inputs/%s_inputs_%s.rds", m, RID[[m]]))$obs_meta
  bal <- balanced_plots(om)
  pp <- pk$projection_predictions; pp <- pp[as.integer(pp$plot_id) %in% bal, ]
  ph <- pk$posterior_predictions;  ph <- ph[as.integer(ph$plot_id) %in% bal, ]
  S <- tapply(pp$total_soc,    list(pp$draw, pp$year), mean)
  U <- tapply(pp$total_soc_up, list(pp$draw, pp$year), mean)
  H <- tapply(ph$total_soc,    list(ph$draw, ph$year), mean)
  yr <- as.integer(colnames(S)); hy <- as.integer(colnames(H))
  C24 <- H[, match(2024L, hy)]
  own <- (C24 - H[, match(1985L, hy)]) / (2024 - 1985)        # model's own 1985-2024 rate
  # effective litter input of 2024 (sigma_input x product), for absolute units
  post <- readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_%s.rds", m, RID[[m]]))
  si <- median(as.matrix(post)[, "sigma_input"])
  j24 <- mean(rowsum(rowSums(ir[ir$year == 2024 & ir$plot_id %in% bal, lc]),
                     ir$plot_id[ir$year == 2024 & ir$plot_id %in% bal]))
  for (w in c(10L, 20L)) {
    s1 <- (S[, match(2024L + w, yr)] - C24) / w
    su <- (U[, match(2024L + w, yr)] - C24) / w
    g  <- (su - s1) / (MULT - 1)                                 # sink per unit multiplier
    for (tg in c(names(TARGETS), "model_own_1985_2024")) {
      tv <- if (tg == "model_own_1985_2024") own else TARGETS[[tg]]
      need <- (tv - s1) / g                                      # m - 1
      res[[length(res) + 1]] <- data.frame(model = m, window = w, target = tg,
        target_rate = median(tv), baseline_sink = median(s1), gain_per_10pct = median(g) / 10,
        extra_input_pct = median(100 * need), lo = quantile(100 * need, .05), hi = quantile(100 * need, .95),
        extra_input_tC = median(need) * si * j24, row.names = NULL)
    }
  }
}
R <- do.call(rbind, res)
num <- sapply(R, is.numeric); R[num] <- lapply(R[num], round, 3)
for (tg in unique(R$target)) {
  cat(sprintf("\n=== target: %s ===\n", tg))
  print(R[R$target == tg, c("model", "window", "target_rate", "baseline_sink", "extra_input_pct", "lo", "hi", "extra_input_tC")],
        row.names = FALSE)
}
saveRDS(R, "doublechecks/figures/input_to_hold_sink.rds")
cat("\nextra_input_pct: % increase of the 2024 litter input, held from 2025, needed to keep the target sink\n",
    "as the mean over the window (median and 5-95% over draws); extra_input_tC: the same in tC/ha/yr of\n",
    "effective input (median sigma_input x the product's 2024 litter). No further warming assumed.\n")
