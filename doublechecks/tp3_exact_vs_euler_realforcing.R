setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
suppressMessages({
  source("Calibration_real_data_transient/calibration_engine_transient.R")
suppressMessages(source("Model_functions_real_data_transient/Decomposition_functions/Yasso/yasso07_wrapper_transient.R"))
  source("Model_functions_real_data_transient/Decomposition_functions/SimpleModels/tp3_wrapper_transient.R")
})
RUN_ID <- "20260630_090644"
post <- readRDS(sprintf("Calibration_real_data_transient/runs/TP3_posterior_%s.rds", RUN_ID))
inp  <- readRDS(sprintf("Data/model_inputs/TP3_inputs_%s.rds", RUN_ID))
med  <- apply(BayesianTools::getSample(post), 2, median)

compute_xi_tp3 <- function(clim, mp) compute_xi_yasso07(clim$temp_mean, clim$temp_amplitude, clim$precip, mp["beta1"], mp["beta2"], mp["gamma"])

# Euler cascade step (the OLD integrator)
euler_run <- function(inputs, mp, C_init, xi) {
  n <- nrow(inputs); A<-S<-H<-numeric(n)
  cA<-C_init["A"];cS<-C_init["S"];cH<-C_init["H"]
  for (t in seq_len(n)) {
    J <- inputs$J_total[t]*mp["sigma_input"]; x<-xi[t]
    kA<-mp["alpha_A"]*x; kS<-mp["alpha_S"]*x; kH<-mp["alpha_H"]
    nA <- cA + J - kA*cA
    nS <- cS + mp["p_S"]*kA*cA - kS*cS
    nH <- cH + mp["p_H"]*kS*cS - kH*cH
    cA<-nA;cS<-nS;cH<-nH; A[t]<-cA;S[t]<-cS;H[t]<-cH
  }
  A+S+H
}
rough <- function(x){ x<-x[is.finite(x)]; mean(abs(diff(diff(x)))) }  # mean |2nd diff|

# run across a sample of plots
set.seed(1)
ids <- sample(inp$plots, 12)
cat(sprintf("%-8s %8s %8s %8s %8s %8s\n","plot","kA_med","kA_p95","frac>1","rough_ex","rough_eu"))
for (id in ids) {
  clim <- inp$climate_by_plot[[id]]; lm <- inp$litter_means[[id]]; ins <- inp$inputs_by_plot[[id]]
  xi <- compute_xi_tp3(clim, med)
  C0 <- tp3_transient_init(med, lm, mean(xi))
  ex <- tp3_run(ins, med, C0, xi)$total_soc
  eu <- euler_run(ins, med, C0, xi)
  kA <- med["alpha_A"]*xi
  cat(sprintf("%-8s %8.2f %8.2f %8.2f %8.4f %8.4f\n", id, median(kA), quantile(kA,.95),
              mean(kA>1), rough(ex), rough(eu)))
}
