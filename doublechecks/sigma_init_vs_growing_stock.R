# =============================================================================
# sigma_init_vs_growing_stock.R   (2026-08-17)
#
# ⚠⚠ STALE AS OF 2026-08-18 -- DO NOT CITE ITS NUMBERS. It implements the PRE-P1
# engine, where the 1917 anchor was J_full_mean. P1 (2026-08-07) moved it to
# J_t0_mean, so the correct denominator is the 1985 flux, not the 1986-2024 mean.
# Consequences: the "volume-implied lower bound" of ~0.64 is wrong (the right
# comparison is 1400/1775 = 0.789, and 0.789^0.43 = 0.903), and the 0.818
# "inversion threshold" is obsolete -- with both ends on the same flux the pre-run
# builds iff sigma_init < 1, full stop. The whole "[0.64, 0.82] window" is wrong on
# both ends.
#
# ⚠⚠ RETIRED AS A BOUND, 2026-08-19 (decision, Lorenzo). Even corrected, this
# constrains the LITTER FLUX ratio, not the soil carbon stock. sigma_init also
# equilibrates the 1917 soil at that flux, and the NFI is silent on whether 1917
# soils had caught up -- after a century of slash-and-burn, raking and heavy
# cutting, plausibly they had not, and the one-parameter transient init cannot
# express "high flux, disequilibrated soil". Judging our own initial state by an
# EQUILIBRIUM-INIT benchmark would also re-import the convention the manuscript
# argues against (Lehtonen 2016, Palosuo 2008, Peltoniemi 2004).
#
# The bound is set aside so the correlated-likelihood run is judged without it;
# it may return. The replacement test is BIOLOGICAL PLAUSIBILITY of the derived
# state -- 1917 SOC stock and the implied pre-run rate:
#
#         >>>  doublechecks/init_state_plausibility.R  <<<
#
# Keep this script only as the record of where 0.90 came from. Do NOT quote its
# verdict column.
#
# THE POINT: sigma_init has an INDEPENDENT physical referent, and it is currently
# violated by every multi-pool model.
#
# sigma_init is, by construction, the ratio of the 1917 pre-run litter flux to the
# full-record mean flux:  J_1917 = J_full_mean * sigma_init * sigma_input, and the
# forward run uses J * sigma_input. sigma_input cancels. So
#
#        sigma_init  ==  (litter input in 1917) / (mean litter input 1986-2024)
#
# That is a quantity the NATIONAL FOREST INVENTORY can speak to, because litter
# input scales with the standing forest. NFI growing stock (Korhonen et al. 2024,
# Data/forest_history/nfi_growing_stock.csv) gives the 1917 : modern volume ratio.
#
# ⚠ THE COMPARISON IS CONSERVATIVE, and it matters which way. Litter is NOT
# proportional to stem volume -- foliage and fine-root litter track leaf area and
# turnover. In 1917 the forest was more heavily cut and younger, and young stands
# shed MORE litter per m3 than mature ones. Correcting for that raises the implied
# sigma_init ABOVE the volume ratio. So the volume ratio is a LOWER bound on the
# plausible value, and the discrepancy reported below is a floor, not a ceiling.
#
# The companion UPPER bound is in doublechecks/prerun_direction.R: the pre-run
# INVERTS (declines instead of building) once sigma_init > J_t0_mean/J_full_mean,
# median 0.818 over the plots. Together they bracket a narrow feasible window.
#
# Usage:  Rscript doublechecks/sigma_init_vs_growing_stock.R
# =============================================================================

suppressMessages(library(BayesianTools))

RUNS <- c(SP1     = "20260813_080305", TP2     = "20260813_080306",
          TP3     = "20260813_080306", Yasso07 = "20260814_105614",
          Yasso15 = "20260814_105757", Yasso20 = "20260814_105907")
LITTER_YEARS <- 1986:2024
PREINIT_YEAR <- 1917

# --- volume-implied lower bound ---------------------------------------------
gs     <- read.csv("Data/forest_history/nfi_growing_stock.csv")
v_mod  <- mean(approx(gs$midpoint_year, gs$total_stock_Mm3, LITTER_YEARS, rule = 2)$y)
v_init <- approx(gs$midpoint_year, gs$total_stock_Mm3, PREINIT_YEAR, rule = 2)$y
LOWER  <- v_init / v_mod

cat(sprintf("NFI growing stock %d          : %6.0f Mm3  (NFI1 %d = %.0f, curve flat to 1937)\n",
            PREINIT_YEAR, v_init, gs$midpoint_year[1], gs$total_stock_Mm3[1]))
cat(sprintf("time-weighted mean %d-%d : %6.0f Mm3\n",
            min(LITTER_YEARS), max(LITTER_YEARS), v_mod))
cat(sprintf("=> volume-implied sigma_init  : %6.2f   (a LOWER bound, see header)\n\n", LOWER))

# --- pre-run inversion upper bound, from the litter record itself ------------
UPPER <- NA
inp <- sprintf("Data/model_inputs/Yasso15_inputs_%s.rds", RUNS[["Yasso15"]])
if (file.exists(inp)) {
  lm_ <- readRDS(inp)$litter_means
  # Simple-model bundles carry J_{t0,full}_mean directly; Yasso bundles split litter
  # into nwl/fwl/cwl (each a vector over AWEN), so sum the size classes instead.
  tot <- function(z, tag) {
    k <- sprintf("J_%s_mean", tag)
    if (!is.null(z[[k]])) return(sum(z[[k]]))
    sum(unlist(z[sprintf(c("nwl_%s_mean","fwl_%s_mean","cwl_%s_mean"), tag)]))
  }
  r <- vapply(lm_, function(z) tot(z, "t0") / tot(z, "full"), numeric(1))
  UPPER <- median(r, na.rm = TRUE)
  cat(sprintf("pre-run inversion threshold   : %6.2f   (median J_t0/J_full over %d plots)\n\n",
              UPPER, length(r)))
}

# --- posteriors --------------------------------------------------------------
cat(sprintf("%-9s %8s %8s %8s   %-10s %s\n", "model", "q05", "median", "q95", "vs lower", "verdict"))
for (m in names(RUNS)) {
  f <- sprintf("Calibration_real_data_transient/runs/%s_posterior_%s.rds", m, RUNS[[m]])
  if (!file.exists(f)) { cat(sprintf("%-9s  (posterior not found)\n", m)); next }
  x <- getSample(readRDS(f)); j <- grep("sigma_init", colnames(x))[1]
  q <- unname(quantile(x[, j], c(.05, .5, .95)))
  cat(sprintf("%-9s %8.3f %8.3f %8.3f   %8.2fx   %s\n", m, q[1], q[2], q[3], LOWER / q[2],
              if (q[3] < LOWER) "ENTIRE posterior below the volume bound"
              else if (q[2] < LOWER) "median below the volume bound"
              else if (!is.na(UPPER) && q[2] > UPPER) "above the inversion threshold"
              else "inside the window"))
}
cat(sprintf("\nfeasible window ~ [%.2f, %.2f]\n", LOWER, UPPER))
cat("NB sigma_init and sigma_input COUPLE through the 1917 flux\n",
    "   (J_1917 = J_full * sigma_init * sigma_input); do not constrain them independently.\n", sep = "")
