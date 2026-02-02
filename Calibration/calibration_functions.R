# =============================================================================
# BAYESIAN CALIBRATION FUNCTIONS
# =============================================================================
# Flexible calibration framework for SOC models
# Works on Mac (local testing) and Linux cluster (production)
#
# Usage:
#   source("calibration_functions.R")
#   result <- calibrate_model(...)
# =============================================================================

if (!require(BayesianTools)) install.packages("BayesianTools")
library(BayesianTools)
library(parallel)

# =============================================================================
# PARAMETER CONVERSION FUNCTIONS
# =============================================================================
# Convert between parameter vectors (for MCMC) and model parameter lists

#' Convert parameter vector to Yasso07 parameter list
#' @param param_vec Named numeric vector from MCMC sampler
#' @return List compatible with yasso07_run()
vector_to_yasso07_params <- function(param_vec) {
  params <- YASSO07_DEFAULT_PARAMS  # Start with defaults
  
  # Update with calibrated parameters
  for (pname in names(param_vec)) {
    if (pname %in% names(params)) {
      params[[pname]] <- param_vec[pname]
    }
  }
  
  params
}

#' Convert parameter vector to Yasso15 parameter list
vector_to_yasso15_params <- function(param_vec) {
  params <- YASSO15_DEFAULT_PARAMS
  for (pname in names(param_vec)) {
    if (pname %in% names(params)) {
      params[[pname]] <- param_vec[pname]
    }
  }
  params
}

#' Convert parameter vector to Yasso20 parameter list
vector_to_yasso20_params <- function(param_vec) {
  params <- YASSO20_DEFAULT_PARAMS
  for (pname in names(param_vec)) {
    if (pname %in% names(params)) {
      params[[pname]] <- param_vec[pname]
    }
  }
  params
}

#' Convert parameter vector to RothC parameter list
vector_to_rothc_params <- function(param_vec) {
  params <- ROTHC_DEFAULT_PARAMS
  for (pname in names(param_vec)) {
    if (pname %in% names(params)) {
      params[[pname]] <- param_vec[pname]
    }
  }
  params
}

#' Convert parameter vector to Q-model parameter list
vector_to_q_params <- function(param_vec) {
  params <- Q_MODEL_DEFAULT_PARAMS
  for (pname in names(param_vec)) {
    if (pname %in% names(params)) {
      params[[pname]] <- param_vec[pname]
    }
  }
  params
}

# =============================================================================
# LIKELIHOOD FUNCTION FACTORY
# =============================================================================
# Creates likelihood functions for different models

#' Create likelihood function for a given model
#' 
#' @param model_name Character: "yasso07", "yasso15", "yasso20", "rothc", "q_model"
#' @param input_df Data frame with litter inputs and climate
#' @param site_df Data frame with site properties
#' @param obs_df Data frame with SOC observations (plot_id, year, soc_obs_tCha)
#' @param n_cores Number of cores for parallel plot processing (NULL = no parallel)
#' @param spinup_years Spinup length
#' @param verbose Print progress messages
#' @return Likelihood function for BayesianTools
create_likelihood <- function(model_name,
                              input_df,
                              site_df,
                              obs_df,
                              n_cores = NULL,
                              spinup_years = 100,
                              verbose = FALSE) {
  
  # Validate model name
  valid_models <- c("yasso07", "yasso15", "yasso20", "rothc", "q_model")
  if (!model_name %in% valid_models) {
    stop("model_name must be one of: ", paste(valid_models, collapse = ", "))
  }
  
  # Get plot IDs
  plot_ids <- unique(site_df$plot_id)
  n_plots <- length(plot_ids)
  
  if (verbose) {
    cat(sprintf("Creating likelihood for %s\n", model_name))
    cat(sprintf("  Plots: %d\n", n_plots))
    cat(sprintf("  Observations: %d\n", nrow(obs_df)))
    cat(sprintf("  Parallel cores: %s\n", ifelse(is.null(n_cores), "None", n_cores)))
  }
  
  # Pre-split data for efficiency
  input_split <- split(input_df, input_df$plot_id)
  site_split <- split(site_df, site_df$plot_id)
  
  # Select parameter conversion function
  param_converter <- switch(model_name,
    yasso07 = vector_to_yasso07_params,
    yasso15 = vector_to_yasso15_params,
    yasso20 = vector_to_yasso20_params,
    rothc = vector_to_rothc_params,
    q_model = vector_to_q_params
  )
  
  # Create likelihood function
  likelihood_fn <- function(param_vec) {
    
    # Check for valid parameters (all finite)
    if (any(!is.finite(param_vec))) {
      return(-Inf)
    }
    
    # Extract observation error (assumed to be last parameter)
    sigma <- param_vec[length(param_vec)]
    if (sigma <= 0) return(-Inf)  # Sigma must be positive
    
    # Convert parameter vector to model parameters
    model_params <- param_converter(param_vec)
    
    # Wrapper function for single plot
    run_single_plot <- function(pid) {
      tryCatch({
        result <- run_soc_models_oneplot_calibration(
          input_df = input_split[[pid]],
          site_row = site_split[[pid]],
          models = model_name,
          params_list = setNames(list(model_params), model_name),
          spinup_years = spinup_years,
          validate = FALSE
        )
        
        # Add plot_id
        result$plot_id <- pid
        result
        
      }, error = function(e) {
        # If model crashes, return NULL
        if (verbose) cat(sprintf("  Error in plot %s: %s\n", pid, e$message))
        return(NULL)
      })
    }
    
    # Run models (parallel or sequential)
    if (!is.null(n_cores) && n_cores > 1) {
      # Parallel execution
      cl <- makeCluster(n_cores)
      on.exit(stopCluster(cl))
      
      # Export necessary objects
      clusterExport(cl, 
                    c("input_split", "site_split", "model_params", "model_name",
                      "spinup_years", "run_soc_models_oneplot_calibration",
                      "YASSO07_DEFAULT_PARAMS", "YASSO15_DEFAULT_PARAMS",
                      "YASSO20_DEFAULT_PARAMS", "ROTHC_DEFAULT_PARAMS",
                      "Q_MODEL_DEFAULT_PARAMS"),
                    envir = environment())
      
      results_list <- parLapply(cl, plot_ids, run_single_plot)
      
    } else {
      # Sequential execution (good for debugging)
      results_list <- lapply(plot_ids, run_single_plot)
    }
    
    # Remove NULL results (crashed plots)
    results_list <- results_list[!sapply(results_list, is.null)]
    
    if (length(results_list) == 0) {
      return(-Inf)  # All plots failed
    }
    
    # Combine results
    pred_df <- dplyr::bind_rows(results_list)
    
    # Match predictions with observations
    matched <- merge(obs_df, pred_df, by = c("plot_id", "year"), all.x = TRUE)
    
    # Check for missing predictions
    if (any(is.na(matched[[model_name]]))) {
      n_missing <- sum(is.na(matched[[model_name]]))
      if (verbose) cat(sprintf("  Warning: %d observations have no predictions\n", n_missing))
      return(-Inf)
    }
    
    # Calculate log-likelihood (Gaussian)
    residuals <- matched$soc_obs_tCha - matched[[model_name]]
    ll <- sum(dnorm(residuals, mean = 0, sd = sigma, log = TRUE))
    
    return(ll)
  }
  
  return(likelihood_fn)
}

# =============================================================================
# MAIN CALIBRATION FUNCTION
# =============================================================================

#' Run Bayesian calibration for a SOC model
#' 
#' @param model_name Character: which model to calibrate
#' @param input_df Full input data (will be subset if plot_subset provided)
#' @param site_df Full site data (will be subset if plot_subset provided)
#' @param obs_df Full observation data (will be subset if plot_subset provided)
#' @param params_to_calibrate Named vector: c(param1 = c(lower, upper), ...)
#' @param sigma_prior Numeric vector: c(lower, upper) for observation error
#' @param n_iter Number of MCMC iterations
#' @param n_chains Number of MCMC chains
#' @param n_cores Number of cores for plot-level parallelization (NULL = sequential)
#' @param plot_subset Character vector of plot IDs to use (NULL = all plots)
#' @param spinup_years Spinup length
#' @param sampler MCMC sampler: "DEzs" (default), "DREAMzs", "Metropolis"
#' @param verbose Print progress
#' @return BayesianTools MCMC output object
calibrate_model <- function(model_name,
                            input_df,
                            site_df,
                            obs_df,
                            params_to_calibrate,
                            sigma_prior = c(1, 10),
                            n_iter = 10000,
                            n_chains = 3,
                            n_cores = NULL,
                            plot_subset = NULL,
                            spinup_years = 100,
                            sampler = "DEzs",
                            verbose = TRUE) {
  
  # -------------------------------------------------------------------------
  # Subset data if requested (for testing)
  # -------------------------------------------------------------------------
  if (!is.null(plot_subset)) {
    if (verbose) cat(sprintf("Subsetting to %d plots\n", length(plot_subset)))
    input_df <- input_df[input_df$plot_id %in% plot_subset, ]
    site_df <- site_df[site_df$plot_id %in% plot_subset, ]
    obs_df <- obs_df[obs_df$plot_id %in% plot_subset, ]
  }
  
  # -------------------------------------------------------------------------
  # Setup
  # -------------------------------------------------------------------------
  if (verbose) {
    cat("=============================================================================\n")
    cat(sprintf("BAYESIAN CALIBRATION: %s\n", toupper(model_name)))
    cat("=============================================================================\n")
    cat(sprintf("Model: %s\n", model_name))
    cat(sprintf("Plots: %d\n", length(unique(site_df$plot_id))))
    cat(sprintf("Observations: %d\n", nrow(obs_df)))
    cat(sprintf("Parameters: %s\n", paste(names(params_to_calibrate), collapse = ", ")))
    cat(sprintf("Iterations: %d\n", n_iter))
    cat(sprintf("Chains: %d\n", n_chains))
    cat(sprintf("Cores: %s\n", ifelse(is.null(n_cores), "Sequential", n_cores)))
    cat(sprintf("Sampler: %s\n\n", sampler))
  }
  
  # -------------------------------------------------------------------------
  # Create prior
  # -------------------------------------------------------------------------
  # Extract lower and upper bounds
  n_params <- length(params_to_calibrate)
  lower <- sapply(params_to_calibrate, function(x) x[1])
  upper <- sapply(params_to_calibrate, function(x) x[2])
  param_names <- names(params_to_calibrate)
  
  # Add sigma (observation error) as last parameter
  lower <- c(lower, sigma_prior[1])
  upper <- c(upper, sigma_prior[2])
  param_names <- c(param_names, "sigma")
  
  prior <- createUniformPrior(lower = lower, upper = upper, best = NULL)
  
  # -------------------------------------------------------------------------
  # Create likelihood
  # -------------------------------------------------------------------------
  if (verbose) cat("Creating likelihood function...\n")
  
  likelihood <- create_likelihood(
    model_name = model_name,
    input_df = input_df,
    site_df = site_df,
    obs_df = obs_df,
    n_cores = n_cores,
    spinup_years = spinup_years,
    verbose = verbose
  )
  
  # -------------------------------------------------------------------------
  # Create Bayesian setup
  # -------------------------------------------------------------------------
  if (verbose) cat("Creating Bayesian setup...\n")
  
  bayesian_setup <- createBayesianSetup(
    likelihood = likelihood,
    prior = prior,
    names = param_names
  )
  
  # -------------------------------------------------------------------------
  # Run MCMC
  # -------------------------------------------------------------------------
  if (verbose) {
    cat("\n")
    cat("Starting MCMC sampling...\n")
    cat(sprintf("Estimated time: %.1f minutes\n", 
                n_iter * length(unique(site_df$plot_id)) * 0.025 / 60))
    cat("\n")
  }
  
  start_time <- Sys.time()
  
  out <- runMCMC(
    bayesianSetup = bayesian_setup,
    sampler = sampler,
    settings = list(
      iterations = n_iter,
      nrChains = n_chains,
      message = verbose
    )
  )
  
  end_time <- Sys.time()
  elapsed <- as.numeric(difftime(end_time, start_time, units = "mins"))
  
  if (verbose) {
    cat("\n")
    cat("=============================================================================\n")
    cat("CALIBRATION COMPLETE\n")
    cat("=============================================================================\n")
    cat(sprintf("Elapsed time: %.1f minutes\n", elapsed))
    cat(sprintf("Time per iteration: %.2f seconds\n", elapsed * 60 / n_iter))
  }
  
  return(out)
}

# =============================================================================
# DIAGNOSTIC FUNCTIONS
# =============================================================================

#' Print calibration summary
#' @param mcmc_output Output from calibrate_model()
#' @param model_name Name of the model (for display)
print_calibration_summary <- function(mcmc_output, model_name = "") {
  
  cat("\n")
  cat("=============================================================================\n")
  cat(sprintf("CALIBRATION SUMMARY: %s\n", toupper(model_name)))
  cat("=============================================================================\n\n")
  
  # Get summary
  summ <- summary(mcmc_output)
  print(summ)
  
  cat("\n")
  
  # Convergence diagnostics - FIXED ERROR HANDLING
  n_chains <- tryCatch({
    mcmc_output$setup$numChains
  }, error = function(e) {
    NULL
  })
  
  if (!is.null(n_chains) && length(n_chains) > 0 && n_chains > 1) {
    cat("Convergence diagnostics (Gelman-Rubin):\n")
    gelm <- gelmanDiagnostics(mcmc_output)
    print(gelm)
    cat("\n")
  }
  
  # DIC
  cat("Model selection criterion:\n")
  dic_val <- DIC(mcmc_output)
  cat(sprintf("  DIC: %.2f\n", dic_val$DIC))
  cat("\n")
}

#' Plot calibration diagnostics
#' @param mcmc_output Output from calibrate_model()
#' @param model_name Name of the model (for plot title)
plot_calibration_diagnostics <- function(mcmc_output, model_name = "") {
  
  # Trace plots
  plot(mcmc_output, which = c("trace"))
  title(main = paste(model_name, "- Trace plots"), outer = TRUE, line = -1)
  
  # Marginal distributions
  plot(mcmc_output, which = c("marginal"))
  title(main = paste(model_name, "- Marginal posteriors"), outer = TRUE, line = -1)
  
  # Correlation plot
  correlationPlot(mcmc_output)
  title(main = paste(model_name, "- Parameter correlations"), outer = TRUE, line = -1)
}

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

#' Extract observations from input data
#' @param input_df Input data frame with soc_obs_tCha column
#' @return Data frame with plot_id, year, soc_obs_tCha
extract_observations <- function(input_df) {
  obs_df <- input_df %>%
    dplyr::filter(!is.na(soc_obs_tCha)) %>%
    dplyr::select(plot_id, year, soc_obs_tCha) %>%
    dplyr::distinct()
  
  return(obs_df)
}

cat("Calibration functions loaded successfully!\n")
cat("Main function: calibrate_model()\n")
cat("Available models: yasso07, yasso15, yasso20, rothc, q_model\n")
