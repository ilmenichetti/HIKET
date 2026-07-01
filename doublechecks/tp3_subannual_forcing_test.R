setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
suppressMessages({
  source("Calibration_real_data_transient/calibration_engine_transient.R")
  source("Model_functions_real_data_transient/Decomposition_functions/Yasso/yasso07_wrapper_transient.R")
  source("Model_functions_real_data_transient/Decomposition_functions/SimpleModels/tp3_wrapper_transient.R")
})
RUN_ID <- "20260630_090644"
post <- readRDS(sprintf("Calibration_real_data_transient/runs/TP3_posterior_%s.rds", RUN_ID))
inp  <- readRDS(sprintf("Data/model_inputs/TP3_inputs_%s.rds", RUN_ID))
med  <- apply(BayesianTools::getSample(post), 2, median)

# exact cascade step over arbitrary sub-step length h (h=1 == .tp3_step)
step_h <- function(cA,cS,cH,kA,kS,kH,pS,pH,J,h){
  if(abs(kA-kS)<1e-6)kS<-kS+1e-6; if(abs(kA-kH)<1e-6)kH<-kH+1e-6; if(abs(kS-kH)<1e-6)kH<-kH+2e-6
  Ass<-J/kA; Sss<-pS*J/kS; Hss<-pH*pS*J/kH
  l1<--kA;l2<--kS;l3<--kH; e1<-exp(l1*h);e2<-exp(l2*h);e3<-exp(l3*h)
  d21<-(e2-e1)/(l2-l1); d32<-(e3-e2)/(l3-l2)
  dd31<-e1/((l1-l2)*(l1-l3))+e2/((l2-l1)*(l2-l3))+e3/((l3-l1)*(l3-l2))
  a<-pS*kA; cc<-pH*kS; dA<-cA-Ass;dS<-cS-Sss;dH<-cH-Hss
  c(Ass+e1*dA, Sss+a*d21*dA+e2*dS, Hss+a*cc*dd31*dA+cc*d32*dS+e3*dH)
}
# temperature-only modifier at temperature T (matches compute_xi_yasso07 structure)
tmod <- function(T,b1,b2) exp(b1*T + b2*T^2)
seas <- cos(2*pi*((1:12)-7)/12)   # +1 in July, -1 in Jan; peak dev = temp_amp

# monthly sub-annual run: 12 exact sub-steps/yr, seasonal T reconstructed, litter flux constant
run_monthly <- function(ins, clim, mp, C0){
  n<-nrow(ins); tot<-numeric(n); cA<-C0["A"];cS<-C0["S"];cH<-C0["H"]
  pm <- 1 - exp(mp["gamma"]*clim$precip/1000)      # annual precip modifier (same as annual scheme)
  for(t in seq_len(n)){
    J<-ins$J_total[t]*mp["sigma_input"]
    Tm <- clim$temp_mean[t] + clim$temp_amplitude[t]*seas
    xim <- tmod(Tm, mp["beta1"], mp["beta2"]) * pm[t]   # 12 monthly xi
    for(m in 1:12){
      kA<-mp["alpha_A"]*xim[m]; kS<-mp["alpha_S"]*xim[m]; kH<-mp["alpha_H"]
      C<-step_h(cA,cS,cH,kA,kS,kH,mp["p_S"],mp["p_H"],J,1/12)
      cA<-C[1];cS<-C[2];cH<-C[3]
    }
    tot[t]<-cA+cS+cH
  }
  tot
}
compute_xi_tp3 <- function(clim, mp) compute_xi_yasso07(clim$temp_mean, clim$temp_amplitude, clim$precip, mp["beta1"], mp["beta2"], mp["gamma"])
rough <- function(x){x<-x[is.finite(x)]; mean(abs(diff(diff(x))))}

set.seed(1); ids<-sample(inp$plots,12)
cat(sprintf("%-8s %8s %8s %8s %8s %8s\n","plot","mean_ann","mean_mon","rgh_ann","rgh_mon","dMeanSOC"))
for(id in ids){
  clim<-inp$climate_by_plot[[id]]; lm<-inp$litter_means[[id]]; ins<-inp$inputs_by_plot[[id]]
  xi<-compute_xi_tp3(clim,med); C0<-tp3_transient_init(med,lm,mean(xi))
  ann<-tp3_run(ins,med,C0,xi)$total_soc
  mon<-run_monthly(ins,clim,med,C0)
  cat(sprintf("%-8s %8.1f %8.1f %8.3f %8.3f %8.2f\n",id,mean(ann),mean(mon),rough(ann),rough(mon),mean(mon)-mean(ann)))
}
