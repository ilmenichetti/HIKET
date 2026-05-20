# recover_tp2_posterior.R
# Reconstructs TP2 physical-space posterior from chains file.
# Run: apptainer_wrapper exec Rscript recover_tp2_posterior.R
# Check output before uncommenting the saveRDS line.

library(BayesianTools)
source("./Calibration_real_data/calibration_engine.R")

RUN_ID <- "20260512_193632"

param_spec <- list(
  list(names = "alpha_A",     type = "log"),
  list(names = "alpha_H",     type = "log"),
  list(names = "p_H",         type = "logit"),
  list(names = "beta1",       type = "log"),
  list(names = "beta2",       type = "unconstrained"),
  list(names = "gamma",       type = "unconstrained"),
  list(names = "sigma_init",  type = "log"),
  list(names = "sigma_input", type = "log")
)
transforms  <- build_transforms(param_spec)
to_original <- transforms$to_original
FREE_NAMES  <- transforms$param_names

chains  <- readRDS(file.path("./Calibration_real_data/runs",
                             sprintf("TP2_chains_%s.rds", RUN_ID)))
unc_mat <- getSample(chains, start = 2001)  # skip 2000 burnin
cat(sprintf("Unconstrained samples: %d x %d\n", nrow(unc_mat), ncol(unc_mat)))

phys_mat <- t(apply(unc_mat, 1, to_original))
colnames(phys_mat) <- FREE_NAMES
phys_df  <- as.data.frame(phys_mat)

# --- Sanity checks ---
cat("\nParameter ranges (physical space):\n")
for (nm in FREE_NAMES)
  cat(sprintf("  %-14s  [%9.5f, %9.5f]  median = %.5f\n",
              nm, min(phys_df[[nm]]), max(phys_df[[nm]]),
              median(phys_df[[nm]])))

# Expected ranges:
#   alpha_A    ~ 0.1 - 5.0   (positive rate)
#   alpha_H    ~ 0.0001 - 0.01
#   p_H        ~ 0.001 - 0.5  (fraction, 0-1)
#   beta1      ~ 0.01 - 0.5   (positive)
#   beta2      ~ any small negative
#   gamma      ~ -3 to -0.3   (negative)
#   sigma_init ~ 0.01 - 1.0
#   sigma_input~ 0.1 - 5.0

cat("\nAll alpha_A > 0:   ", all(phys_df$alpha_A > 0), "\n")
cat("All alpha_H > 0:   ", all(phys_df$alpha_H > 0), "\n")
cat("All p_H in (0,1):  ", all(phys_df$p_H > 0 & phys_df$p_H < 1), "\n")
cat("All beta1 > 0:     ", all(phys_df$beta1 > 0), "\n")
cat("All gamma < 0:     ", all(phys_df$gamma < 0), "\n")

# --- Uncomment to save once checks pass ---
saveRDS(phys_mat,
        file.path("./Calibration_real_data/runs",
                  sprintf("TP2_posterior_%s.rds", RUN_ID)))
cat("\nPosterior saved.\n")