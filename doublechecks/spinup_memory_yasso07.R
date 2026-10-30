# =============================================================================
# spinup_memory_yasso07.R   (2026-10-30)
#
# How long does PUBLISHED Yasso07 remember its starting state?
# The Finnish inventory (NID 2025, s. 6.4) starts Yasso07 in balance in 1971-76 and
# states that "approximately ten years of simulation since spin-up is enough to
# cancel out the effect of the spin-up level" (Peltoniemi et al. 2006).
#
# Test: under constant litter and climate Yasso07 is linear, so a starting error
# evolves on its own. A model started in balance (C_eq) is compared with a soil that
# really started at s * C_eq (s < 1: below balance, still recovering). The model stays
# flat; the soil gains. The difference in the ANNUAL CHANGE is what the inventory
# would misreport. Error placed (a) in all pools, (b) in the humus pool only.
#
# Setting: published Yasso07 parameters (YASSO07_DEFAULT_PARAMS), dataset-mean
# climate and litter of the current production bundle (tree litter, no sigma_input).
# Run from repo root: Rscript doublechecks/spinup_memory_yasso07.R
# =============================================================================
suppressMessages(source("manuscript/figures/run_ids.R"))
source("Model_functions_real_data_transient/Decomposition_functions/Yasso/yasso07_wrapper_transient.R")
dyn.load("Model_functions_real_data_transient/Decomposition_functions/Yasso/yasso07.so")

# === Reference forcing: dataset means ========================================
pkg <- readRDS(sprintf("Data/model_inputs/Yasso07_inputs_%s.rds", RID[["Yasso07"]]))
pl  <- as.character(pkg$plots_real)
clim <- data.frame(
  temp_mean      = mean(sapply(pl, function(p) mean(pkg$climate_by_plot[[p]]$temp_mean))),
  temp_amplitude = mean(sapply(pl, function(p) mean(pkg$climate_by_plot[[p]]$temp_amplitude))),
  precip         = mean(sapply(pl, function(p) mean(pkg$climate_by_plot[[p]]$precip))))
gm  <- function(f) rowMeans(sapply(pl, function(p) as.numeric(pkg$litter_means[[p]][[f]])))
nwl <- gm("nwl_mean"); fwl <- gm("fwl_mean"); cwl <- gm("cwl_mean")
P   <- YASSO07_DEFAULT_PARAMS
xi  <- compute_xi_mean_yasso07(clim, P[["beta1"]], P[["beta2"]], P[["gamma"]])
cat(sprintf("Reference: T %.2f C, amplitude %.2f, precip %.0f mm | litter %.2f tC/ha/yr | xi %.3f\n",
            clim$temp_mean, clim$temp_amplitude, clim$precip, sum(nwl, fwl, cwl), xi))

# === Balanced state and forward runs =========================================
Ceq <- yasso07_steady_state(P, nwl, fwl, cwl, xi)          # 15 = 3 size cohorts x AWENH
NY  <- 200L
idf <- data.frame(year = seq_len(NY), t(replicate(NY, c(nwl, fwl, cwl))))
names(idf) <- c("year", paste0(rep(c("nwl_", "fwl_", "cwl_"), each = 4), c("A","W","E","N")))
run <- function(C0) yasso07_run(idf, P, C0, rep(xi, NY))$total_soc
tot_eq <- sum(Ceq); H <- c(5, 10, 15)                       # humus index per cohort
cat(sprintf("Balanced stock %.1f tC/ha, of which humus %.0f%%\n\n", tot_eq, 100 * sum(Ceq[H]) / tot_eq))

# === Starting 10% / 20% below balance ========================================
case <- function(lab, C0) {
  soc  <- c(sum(C0), run(C0))                                # year 0..NY
  err  <- tot_eq - soc                                       # stock error of the balanced start
  rate <- diff(soc)                                          # soil's real gain = misreported change
  data.frame(case = lab, start_deficit = err[1],
             left_10 = err[11] / err[1], left_20 = err[21] / err[1], left_40 = err[41] / err[1],
             rate_yr1 = rate[1], rate_yr10 = rate[10], rate_yr20 = rate[20], rate_yr40 = rate[40])
}
mk_H <- function(s) { C <- Ceq; C[H] <- s * C[H]; C }
out <- rbind(case("all pools  -10%", 0.9 * Ceq), case("all pools  -20%", 0.8 * Ceq),
             case("humus only -10%", mk_H(0.9)), case("humus only -20%", mk_H(0.8)))
num <- sapply(out, is.numeric)
out[num] <- lapply(out[num], round, 3)
print(out, row.names = FALSE)
cat("\nleft_N = share of the starting stock error still present after N years;",
    "\nrate_yrN = the soil's real gain in year N (tC/ha/yr), which a balanced start reports as 0.",
    "\nCompare: observed 1985-2024 change on the balanced plots +0.26 tC/ha/yr.\n")
