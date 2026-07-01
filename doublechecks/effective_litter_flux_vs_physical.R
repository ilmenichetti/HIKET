setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
suppressMessages(library(BayesianTools))
rid <- list(SP1="20260608_015026", TP2="20260608_020212", TP3="20260630_090644",
            Yasso07="20260611_032825", Yasso15="20260611_032825", Yasso20="20260611_033140")
littercols <- function(df){ grep("^(nwl_|fwl_|cwl_)", names(df), value=TRUE) }
totlit <- function(df){
  if("J_total" %in% names(df)) return(mean(df$J_total))
  lc <- littercols(df); mean(rowSums(df[,lc,drop=FALSE]))
}
LO <- 0.5; HI <- 9   # physical bounds on total forest C input tC/ha/yr
cat(sprintf("%-8s %7s %7s %8s %8s   %s\n","model","rawJ","si_med","effJ_md","effJ_hi","verdict"))
for(m in names(rid)){
  post <- readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_%s.rds", m, rid[[m]]))
  inp  <- readRDS(sprintf("Data/model_inputs/%s_inputs_%s.rds", m, rid[[m]]))
  siq <- quantile(getSample(post)[,"sigma_input"], c(.5,.975))
  Jm  <- median(sapply(inp$plots, function(id) totlit(inp$inputs_by_plot[[id]])))
  eff <- Jm*siq[1]; effhi <- Jm*siq[2]
  v <- if(eff>HI) "ABOVE physical max" else if(eff<LO) "below min" else "physical OK"
  cat(sprintf("%-8s %7.2f %7.2f %8.1f %8.1f   %s\n", m, Jm, siq[1], eff, effhi, v))
}
cat(sprintf("\nPhysical window used: [%.1f, %.1f] tC/ha/yr total litter input.\n", LO, HI))
cat("Given shared raw litter ~2.5 (median), a cap effJ<=9 => sigma_input<=~3.6 (median plot),\n")
cat("and floor effJ>=0.5 => sigma_input>=~0.2. A homogeneous physical prior on sigma_input\n")
cat("~ [0.4, 3] (soft) would leave Yasso07/15 untouched, put Yasso20 at the edge, and bind SP1/TP2/TP3.\n")
