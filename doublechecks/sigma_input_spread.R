setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
suppressMessages(library(BayesianTools))
rid <- list(SP1="20260608_015026", TP2="20260608_020212", TP3="20260630_090644",
            Yasso07="20260611_032825", Yasso15="20260611_032825", Yasso20="20260611_033140")
cat(sprintf("%-8s %7s %7s %7s %7s | %7s %8s | %6s\n",
    "model","si2.5%","si_med","si97.5%","sinit_md","rawJ_md","effJ_md","sigmas"))
for(m in names(rid)){
  post <- readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_%s.rds", m, rid[[m]]))
  inp  <- readRDS(sprintf("Data/model_inputs/%s_inputs_%s.rds", m, rid[[m]]))
  s <- getSample(post)
  si <- s[,"sigma_input"]; sinit <- if("sigma_init"%in%colnames(s)) median(s[,"sigma_init"]) else NA
  q <- quantile(si, c(.025,.5,.975))
  Jt <- median(sapply(inp$plots, function(id) mean(inp$inputs_by_plot[[id]]$J_total)))
  cat(sprintf("%-8s %7.2f %7.2f %7.2f %7.2f | %7.2f %8.1f | %6.1f\n",
      m, q[1], q[2], q[3], sinit, Jt, Jt*q[2], log(q[2])/0.5))
}
cat("\nsi = sigma_input posterior; rawJ_md = median raw litter (tC/ha/yr); effJ_md = raw*sigma_input;\n")
cat("sigmas = log(si_med)/0.5 = distance of posterior median from lognormal(0,0.5) prior centre.\n")
