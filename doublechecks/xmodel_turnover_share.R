setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
suppressMessages({
  source("Calibration_real_data_transient/calibration_engine_transient.R")
  source("Model_functions_real_data_transient/Decomposition_functions/SimpleModels/sp1_wrapper_transient.R")
  source("Model_functions_real_data_transient/Decomposition_functions/SimpleModels/tp2_wrapper_transient.R")
  source("Model_functions_real_data_transient/Decomposition_functions/SimpleModels/tp3_wrapper_transient.R")
})
rid <- list(SP1="20260608_015026", TP2="20260608_020212", TP3="20260630_090644")
med <- list(); inp <- list()
for(m in names(rid)){
  p<-readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_%s.rds",m,rid[[m]]))
  med[[m]]<-apply(BayesianTools::getSample(p),2,median)
  inp[[m]]<-readRDS(sprintf("Data/model_inputs/%s_inputs_%s.rds",m,rid[[m]]))
}
cat("SP1 params:",paste(names(med$SP1),collapse=" "),"\n")
cat("TP2 params:",paste(names(med$TP2),collapse=" "),"\n")
cat("TP3 params:",paste(names(med$TP3),collapse=" "),"\n\n")
cxi <- function(clim,mp) compute_xi_yasso07(clim$temp_mean,clim$temp_amplitude,clim$precip,mp["beta1"],mp["beta2"],mp["gamma"])
rough <- function(x){x<-x[is.finite(x)];mean(abs(diff(diff(x))))}
set.seed(1)
# use common plot IDs present in all bundles
ids <- Reduce(intersect,lapply(inp,function(z)z$plots))
ids <- sample(ids, 30)
# --- per model: active-pool k*xi, active share, roughness ---
res <- data.frame()
for(m in names(rid)){
  mp<-med[[m]]; kax<-c(); ash<-c(); rg<-c()
  for(id in ids){
    clim<-inp[[m]]$climate_by_plot[[id]]; lm<-inp[[m]]$litter_means[[id]]; ins<-inp[[m]]$inputs_by_plot[[id]]
    xi<-cxi(clim,mp)
    if(m=="SP1"){ C0<-sp1_transient_init(mp,lm,mean(xi)); out<-sp1_run(ins,mp,C0,xi); tot<-out$total_soc; ash<-c(ash,1); ka<-mp["alpha"]*xi }
    if(m=="TP2"){ C0<-tp2_transient_init(mp,lm,mean(xi)); out<-tp2_run(ins,mp,C0,xi); tot<-out$total_soc; ash<-c(ash,mean(out$A/tot)); ka<-mp["alpha_A"]*xi }
    if(m=="TP3"){ C0<-tp3_transient_init(mp,lm,mean(xi)); out<-tp3_run(ins,mp,C0,xi); tot<-out$total_soc; ash<-c(ash,mean(out$A/tot)); ka<-mp["alpha_A"]*xi }
    kax<-c(kax,median(ka)); rg<-c(rg,rough(tot))
  }
  res<-rbind(res,data.frame(model=m, kAxi_med=median(kax), fastShare=mean(ash), roughness=mean(rg)))
}
print(res, row.names=FALSE)
