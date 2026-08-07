setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
suppressMessages(library(BayesianTools))
# Track the newest production posterior per model. Previously hardcoded to the
# 20260710 flux_pair run, which meant this check silently reported stale numbers
# after a re-calibration. Ablation posteriors are quarantined out of runs/ by
# doublechecks/quarantine_ablation_runs.R, so "newest" here is a production run.
rid <- sapply(c("SP1","TP2","TP3","Yasso07","Yasso15","Yasso20"), function(m) {
  fs <- list.files("Calibration_real_data_transient/runs",
                   pattern = sprintf("^%s_posterior_[0-9]{8}_[0-9]{6}\\.rds$", m))
  if (!length(fs)) stop("no posterior found for ", m)
  sub(sprintf("^%s_posterior_(.+)\\.rds$", m), "\\1", sort(fs, decreasing = TRUE)[1])
}, simplify = FALSE)
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
