# =============================================================================
# tp3_subannual_forcing_test.R -- does the ANNUAL time step matter?
# -----------------------------------------------------------------------------
# Re-runs the CURRENT TP3 posterior median with 12 exact monthly sub-steps per
# year (seasonal temperature reconstructed from temp_mean + temp_amplitude,
# litter flux constant within the year) and compares with the production annual
# scheme. TP3 is the test case because its fast pool turns over > once a year.
#
# 2026-09-23: repointed to the current run via manuscript/figures/run_ids.R
# (was hard-wired to 20260630_090644, pre-ICBM, not reproducible) and updated to
# the current parameterisation: alpha_A fixed (TP3_ALPHA_A_FIXED), xi on all
# three pools (C2). Quoted in the manuscript, sect. "Numerical integration".
#
# Run from repo root:  Rscript doublechecks/tp3_subannual_forcing_test.R
# =============================================================================

suppressMessages({
  source("Calibration_real_data_transient/calibration_engine_transient.R")
  source("Model_functions_real_data_transient/Decomposition_functions/Yasso/yasso07_wrapper_transient.R")
  source("Model_functions_real_data_transient/Decomposition_functions/SimpleModels/tp3_wrapper_transient.R")
  source("Prior_specs/TP3_priors.R")
  source("manuscript/figures/run_ids.R")
})
RUN_ID <- RID[["TP3"]]
cat("TP3 RUN_ID:", RUN_ID, "\n")
post <- readRDS(sprintf("Calibration_real_data_transient/runs/TP3_posterior_%s.rds", RUN_ID))
inp  <- readRDS(sprintf("Data/model_inputs/TP3_inputs_%s.rds", RUN_ID))
mp   <- c(apply(post, 2, median), alpha_A = TP3_ALPHA_A_FIXED)   # posterior is PHYSICAL

# === exact cascade step over a sub-step of length h (h = 1 == .tp3_step) ===
step_h <- function(cA,cS,cH,kA,kS,kH,pS,pH,J,h){
  if(abs(kA-kS)<1e-6)kS<-kS+1e-6; if(abs(kA-kH)<1e-6)kH<-kH+1e-6; if(abs(kS-kH)<1e-6)kH<-kH+2e-6
  Ass<-J/kA; Sss<-pS*J/kS; Hss<-pH*pS*J/kH
  l1<--kA;l2<--kS;l3<--kH; e1<-exp(l1*h);e2<-exp(l2*h);e3<-exp(l3*h)
  d21<-(e2-e1)/(l2-l1); d32<-(e3-e2)/(l3-l2)
  dd31<-e1/((l1-l2)*(l1-l3))+e2/((l2-l1)*(l2-l3))+e3/((l3-l1)*(l3-l2))
  a<-pS*kA; cc<-pH*kS; dA<-cA-Ass;dS<-cS-Sss;dH<-cH-Hss
  c(Ass+e1*dA, Sss+a*d21*dA+e2*dS, Hss+a*cc*dd31*dA+cc*d32*dS+e3*dH)
}
seas <- cos(2*pi*((1:12)-7)/12)                 # +1 July, -1 January

# === monthly run: 12 exact sub-steps a year, xi on all three pools (C2) ===
run_monthly <- function(ins, clim, mp, C0){
  n <- nrow(ins); tot <- numeric(n); cA <- C0[["A"]]; cS <- C0[["S"]]; cH <- C0[["H"]]
  pm <- 1 - exp(mp[["gamma"]]*clim$precip/1000)
  for (t in seq_len(n)) {
    J   <- ins$J_total[t]*mp[["sigma_input"]]
    Tm  <- clim$temp_mean[t] + clim$temp_amplitude[t]*seas
    xim <- exp(mp[["beta1"]]*Tm + mp[["beta2"]]*Tm^2) * pm[t]
    for (m in 1:12) {
      C <- step_h(cA,cS,cH, mp[["alpha_A"]]*xim[m], mp[["alpha_S"]]*xim[m], mp[["alpha_H"]]*xim[m],
                  mp[["p_S"]], mp[["p_H"]], J, 1/12)
      cA <- C[1]; cS <- C[2]; cH <- C[3]
    }
    tot[t] <- cA+cS+cH
  }
  tot
}

# === compare on all calibration plots, national mean trajectory ===
ids <- inp$plots_real
res <- lapply(ids, function(id){
  clim <- inp$climate_by_plot[[id]]; lm <- inp$litter_means[[id]]; ins <- inp$inputs_by_plot[[id]]
  xi   <- compute_xi_yasso07(clim$temp_mean, clim$temp_amplitude, clim$precip, mp["beta1"], mp["beta2"], mp["gamma"])
  n_ss <- min(inp$STEADY_STATE_YEARS, nrow(clim))
  xs   <- compute_xi_mean_yasso07(clim[seq_len(n_ss), , drop=FALSE], mp["beta1"], mp["beta2"], mp["gamma"])
  C0   <- tp3_transient_init(mp, lm, xs)
  data.frame(year = ins$year, ann = tp3_run(ins, mp, C0, xi)$total_soc, mon = run_monthly(ins, clim, mp, C0))
})
d  <- do.call(rbind, res)
nm <- aggregate(cbind(ann, mon) ~ year, d, mean)
cat(sprintf("plots %d | national mean, max |monthly - annual| over years: %.2f tC/ha (mean stock %.1f)\n",
            length(ids), max(abs(nm$mon - nm$ann)), mean(nm$ann)))
cat(sprintf("per plot, mean over years of (monthly - annual): median %.2f, 95th pct of |.| %.2f tC/ha\n",
            median(tapply(d$mon - d$ann, rep(seq_along(res), sapply(res, nrow)), mean)),
            quantile(abs(tapply(d$mon - d$ann, rep(seq_along(res), sapply(res, nrow)), mean)), .95)))
