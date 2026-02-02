# =============================================================================
# EXAMPLE: YASSO07 CALIBRATION RUN
# =============================================================================
# This script shows how to run Bayesian calibration for Yasso07
# 
# USAGE:
# 1. For testing on Mac: Use plot_subset with 4 plots, n_cores = NULL
# 2. For production on cluster: Use plot_subset = NULL (all plots), n_cores = 5-10
#
# Copy this code into your testing script after generating synthetic data
# =============================================================================

# Load calibration functions
source("./calibration_functions.R")

# Load model wrapper (must be loaded BEFORE calibration)
source("./Models_functions/wrapper_calibration_ready.R")

# Load synthetic data generation script
input_data <- read.csv("synthetic_input_data_template.csv")
site_data <- read.csv("synthetic_site_data_template.csv")


# =============================================================================
# SETUP: Define parameters to calibrate
# =============================================================================

# For Yasso07, we'll calibrate:
# - alpha_A: Decomposition rate of acid-soluble pool
# - beta1: Temperature sensitivity parameter
# - sigma: Observation error

params_yasso07 <- list(
  alpha_A = c(0.4, 1.0),
  beta1 = c(0.05, 0.15)
  # sigma is added automatically by the function
)

# Alternative: calibrate more parameters
# params_yasso07_extended <- list(
#   alpha_A = c(0.4, 1.0),
#   alpha_W = c(2.0, 6.0),
#   alpha_N = c(0.1, 0.4),
#   beta1 = c(0.05, 0.15),
#   beta2 = c(-0.003, 0.0),
#   sigma = c(1.0, 10.0)
# )

# =============================================================================
# EXTRACT OBSERVATIONS
# =============================================================================

# Extract SOC observations from synthetic data
obs_data <- extract_observations(input_data)

cat("Observations extracted:\n")
print(obs_data)

# =============================================================================
# RUN CALIBRATION
# =============================================================================

# Configuration for LOCAL TESTING (Mac)
test_config <- list(
  plot_subset = c("PLOT_001", "PLOT_002", "PLOT_003", "PLOT_004"),
  n_iter = 5000,        # Short run for testing
  n_chains = 4,         # Fewer chains
  n_cores = NULL,       # Sequential (easier to debug)
  spinup_years = 100,   # Short spinup (for Q initialization, hybrid approach)
  verbose = TRUE
)

# Configuration for CLUSTER PRODUCTION
# production_config <- list(
#   plot_subset = NULL,    # All plots
#   n_iter = 10000,       # Full run
#   n_chains = 3,         # More chains for convergence
#   n_cores = 5,          # Parallel processing
#   spinup_years = 100,
#   verbose = TRUE
# )

# Choose configuration (switch for Mac vs cluster)
config <- test_config  # Change to production_config on cluster

# =============================================================================
# RUN THE CALIBRATION
# =============================================================================

cat("\n")
cat("Starting Yasso07 calibration...\n")
cat("This will take a few minutes with 4 plots and 1000 iterations\n")
cat("\n")

# Run calibration
calibration_time <- system.time({
  yasso07_mcmc <- calibrate_model(
    model_name = "yasso07",
    input_df = input_data,
    site_df = site_data,
    obs_df = obs_data,
    params_to_calibrate = params_yasso07,
    sigma_prior = c(1.0, 10.0),
    n_iter = config$n_iter,
    n_chains = config$n_chains,
    n_cores = config$n_cores,
    plot_subset = config$plot_subset,
    spinup_years = config$spinup_years,
    sampler = "DEzs",
    verbose = config$verbose
  )
})

# Print timing
cat("\n=============================================================================\n")
cat("CALIBRATION TIMING\n")
cat("=============================================================================\n")
cat(sprintf("Elapsed time: %.2f seconds (%.2f minutes)\n", 
            calibration_time["elapsed"], 
            calibration_time["elapsed"]/60))
cat(sprintf("User time:    %.2f seconds\n", calibration_time["user"]))
cat(sprintf("System time:  %.2f seconds\n", calibration_time["sys"]))
cat("\n")


# =============================================================================
# EXAMINE RESULTS
# =============================================================================

# Print summary
print_calibration_summary(yasso07_mcmc, "Yasso07")

# Plot diagnostics
par(mfrow = c(2, 2))
plot(yasso07_mcmc)

# Get posterior samples
posterior_samples <- getSample(yasso07_mcmc)


# =============================================================================
# SAVE RESULTS
# =============================================================================
# 
# # Save MCMC output
# saveRDS(yasso07_mcmc, "./Calibration/tests/yasso07_calibration_test.rds")
# cat("\nResults saved to: yasso07_calibration_test.rds\n")

# =============================================================================
# POSTERIOR PREDICTIVE CHECK
# =============================================================================

# Get medians
post_median <- apply(posterior_samples, 2, median)
print(post_median)

# Get quantiles
post_quantiles <- apply(posterior_samples, 2, quantile, probs = c(0.025, 0.5, 0.975))
print(post_quantiles)

# Remove sigma (not a model parameter)
model_params_median <- post_median[names(post_median) != "sigma"]
print(model_params_median)

# Convert to Yasso07 format
yasso07_params_median <- vector_to_yasso07_params(model_params_median)

# Run predictions with calibrated parameters
predictions <- run_soc_models_multipplot_calibration(
  input_df = input_data,
  site_df = site_data,
  models = "yasso07",
  params_list = list(yasso07 = yasso07_params_median),
  spinup_years = 100
)

# Match predictions with observations
pred_at_obs <- merge(obs_data, predictions$summary, by = c("plot_id", "year"))

# Plot observed vs predicted
plot(pred_at_obs$soc_obs_tCha, pred_at_obs$yasso07,
     xlab = "Observed SOC (t C/ha)", 
     ylab = "Predicted SOC (t C/ha)",
     main = "Yasso07: Observed vs Predicted (Posterior Median)",
     pch = 19, cex = 1.5)
abline(0, 1, col = "red", lty = 2, lwd = 2)

# Calculate RMSE
rmse <- sqrt(mean((pred_at_obs$soc_obs_tCha - pred_at_obs$yasso07)^2))
cat(sprintf("\nRMSE: %.2f t C/ha\n", rmse))

# Calculate bias
bias <- mean(pred_at_obs$yasso07 - pred_at_obs$soc_obs_tCha)
cat(sprintf("Bias: %.2f t C/ha\n", bias))



# =============================================================================
# PLOT PARAMETER EVOLUTION
# =============================================================================
# Get priors from first chain
prior_lower <- yasso07_mcmc[[1]]$setup$prior$lower
prior_upper <- yasso07_mcmc[[1]]$setup$prior$upper

# Fix parameter names (sigma has empty name)
param_names <- names(prior_lower)
param_names[param_names == ""] <- "sigma"
names(prior_lower) <- param_names
names(prior_upper) <- param_names

# Extract posterior samples
posterior_samples <- getSample(yasso07_mcmc)

# Set up plot
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))

# Loop through parameters
for (i in 1:ncol(posterior_samples)) {
  param_name <- colnames(posterior_samples)[i]
  
  # Posterior samples
  post_samples <- posterior_samples[, i]
  
  # Prior bounds
  p_lower <- prior_lower[param_name]
  p_upper <- prior_upper[param_name]
  
  # Calculate posterior density
  dens <- density(post_samples)
  
  # Set up plot limits
  xlim <- c(min(p_lower, min(post_samples)), 
            max(p_upper, max(post_samples)))
  ylim <- c(0, max(dens$y, 1/(p_upper - p_lower)) * 1.1)
  
  # Create empty plot
  plot(1, type = "n",
       xlim = xlim,
       ylim = ylim,
       main = param_name,
       xlab = "Parameter value",
       ylab = "Density")
  
  # Add prior as filled polygon (red)
  prior_height <- 1 / (p_upper - p_lower)
  polygon(x = c(p_lower, p_upper, p_upper, p_lower),
          y = c(0, 0, prior_height, prior_height),
          col = rgb(1, 0, 0, 0.3),  # Red with 30% opacity
          border = "red",
          lwd = 2)
  
  # Add posterior as filled polygon (blue)
  polygon(x = c(dens$x, rev(dens$x)),
          y = c(dens$y, rep(0, length(dens$y))),
          col = rgb(0, 0, 1, 0.3),  # Blue with 30% opacity
          border = "blue",
          lwd = 2)
  
  # Add posterior median (solid line)
  post_median <- median(post_samples)
  abline(v = post_median, col = "darkblue", lwd = 2, lty = 1)
  
  # Add 95% CI (dashed lines)
  post_ci <- quantile(post_samples, c(0.025, 0.975))
  abline(v = post_ci, col = "darkblue", lwd = 1.5, lty = 2)
  
  # Legend (only on first plot)
  if (i == 1) {
    legend("topright", 
           legend = c("Prior", "Posterior", "Median", "95% CI"),
           col = c("red", "blue", "darkblue", "darkblue"),
           fill = c(rgb(1,0,0,0.3), rgb(0,0,1,0.3), NA, NA),
           border = c("red", "blue", NA, NA),
           lty = c(NA, NA, 1, 2),
           lwd = c(NA, NA, 2, 1.5),
           bty = "n",
           cex = 0.8)
  }
}

par(mfrow = c(1, 1))



# =============================================================================
# NOTES FOR SCALING TO CLUSTER
# =============================================================================

# When moving to cluster with 3700 plots:
# 
# 1. Change config to production_config
# 2. Set n_cores = 5-10 (depending on memory)
# 3. Use SLURM script:
#
# #!/bin/bash
# #SBATCH --nodes=1
# #SBATCH --ntasks-per-node=10
# #SBATCH --time=48:00:00
# #SBATCH --mem=32G
# 
# module load R/4.3.0
# Rscript calibrate_yasso07.R
#
# 4. Expect ~2-3 days runtime for full calibration
