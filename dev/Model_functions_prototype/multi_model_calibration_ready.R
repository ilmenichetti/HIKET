# =============================================================================
# CALIBRATION-READY FAST WRAPPER (SINGLE PLOT)
# =============================================================================
# Modified version with two critical fixes:
# 1. Robust years extraction (doesn't fail if all annual models excluded)
# 2. Custom parameter passing for Bayesian calibration
#
# DEPENDENCIES: This wrapper requires the following to be loaded:
# - Input_matching_functions.R (mapping functions)
# - Yasso07.R, Yasso15.R, Yasso20.R, RothC.R, Q_transient_T_hybrid_RCPP.R
# - dplyr package
# =============================================================================

# Load dependencies if not already loaded
if (!exists("map_input_to_yasso07")) {
  source("./Models_functions/Input_matching_functions.R")
}

if (!exists("YASSO07_DEFAULT_PARAMS")) {
  source("./Models_functions/Yasso07.R")
}

if (!exists("YASSO15_DEFAULT_PARAMS")) {
  source("./Models_functions/Yasso15.R")
}

if (!exists("YASSO20_DEFAULT_PARAMS")) {
  source("./Models_functions/Yasso20.R")
}

if (!exists("ROTHC_DEFAULT_PARAMS")) {
  source("./Models_functions/RothC.R")
}

if (!exists("Q_MODEL_DEFAULT_PARAMS")) {
  source("./Models_functions/Q_transient_T_hybrid_RCPP.R")
}

if (!require(dplyr)) {
  stop("dplyr package required. Install with: install.packages('dplyr')")
}

# Helper functions for RothC aggregation
aggregate_monthly_to_annual_mean <- function(monthly_soc, years) {
  unique_years <- unique(years)
  sapply(unique_years, function(y) mean(monthly_soc[years == y], na.rm = TRUE))
}

aggregate_monthly_to_annual_eoy <- function(monthly_soc, years, months) {
  unique_years <- unique(years)
  sapply(unique_years, function(y) {
    idx <- which(years == y & months == 12)
    if (length(idx) == 0) NA_real_ else monthly_soc[idx[1]]
  })
}

run_soc_models_oneplot_calibration <- function(
    input_df, 
    site_row = NULL,
    models = c("yasso07", "yasso15", "yasso20", "rothc", "q_model"),
    params_list = NULL,  # NEW: custom parameters per model
    spinup_years = NULL,
    rothc_annual_stat = c("eoy", "mean"),
    validate = FALSE) {
  
  rothc_annual_stat <- match.arg(rothc_annual_stat)
  if (validate) validate_input(input_df)
  
  # -------------------------------------------------------------------------
  # Extract parameters for each model (custom or default)
  # -------------------------------------------------------------------------
  get_params <- function(model_name, default_params) {
    if (!is.null(params_list) && model_name %in% names(params_list)) {
      return(params_list[[model_name]])
    }
    return(default_params)
  }
  
  params_y07   <- get_params("yasso07", YASSO07_DEFAULT_PARAMS)
  params_y15   <- get_params("yasso15", YASSO15_DEFAULT_PARAMS)
  params_y20   <- get_params("yasso20", YASSO20_DEFAULT_PARAMS)
  params_rothc <- get_params("rothc", ROTHC_DEFAULT_PARAMS)
  params_q     <- get_params("q_model", Q_MODEL_DEFAULT_PARAMS)
  
  # -------------------------------------------------------------------------
  # Map inputs (local only; don't store)
  # -------------------------------------------------------------------------
  y07_input <- NULL
  if (any(c("yasso07", "yasso15") %in% models)) {
    y07_input <- map_input_to_yasso07(input_df)
  }
  y20_input <- if ("yasso20" %in% models) map_input_to_yasso20(input_df) else NULL
  ro_input  <- if ("rothc"   %in% models) map_input_to_rothc(input_df, site_row) else NULL
  q_input   <- if ("q_model" %in% models) map_input_to_q_model(input_df) else NULL
  
  # -------------------------------------------------------------------------
  # FIX 1: Robust years extraction with fallback
  # -------------------------------------------------------------------------
  years <- NULL
  
  # Try to get years from mapped inputs
  if (!is.null(y07_input)) years <- y07_input$year
  if (is.null(years) && !is.null(y20_input)) years <- y20_input$year
  if (is.null(years) && !is.null(q_input))  years <- q_input$year
  if (is.null(years) && !is.null(ro_input)) years <- sort(unique(ro_input$year))
  
  # Fallback: extract from original input_df
  if (is.null(years)) {
    if ("year" %in% names(input_df)) {
      years <- sort(unique(input_df$year))
      warning("Years extracted from input_df (no models with mapped inputs selected)")
    } else {
      stop("Cannot determine years: no valid model inputs and input_df has no 'year' column",
           call. = FALSE)
    }
  }
  
  # Initialize summary data frame
  summary_df <- data.frame(year = years)
  
  # -------------------------------------------------------------------------
  # FIX 2: Run models with CUSTOM parameters (not hardcoded defaults)
  # -------------------------------------------------------------------------
  
  # YASSO07
  if ("yasso07" %in% models) {
    out <- yasso07_run(
      params = params_y07,  # <-- Custom params
      C0 = NULL, 
      input_df = y07_input,
      spinup_years = if (!is.null(spinup_years)) spinup_years else 5000
    )
    summary_df$yasso07 <- out$total_soc
  }
  
  # YASSO15
  if ("yasso15" %in% models) {
    out <- yasso15_run(
      params = params_y15,  # <-- Custom params
      C0 = NULL, 
      input_df = y07_input,
      spinup_years = if (!is.null(spinup_years)) spinup_years else 5000
    )
    summary_df$yasso15 <- out$total_soc
  }
  
  # YASSO20
  if ("yasso20" %in% models) {
    out <- yasso20_run(
      params = params_y20,  # <-- Custom params
      C0 = NULL, 
      input_df = y20_input, 
      spinup = TRUE
    )
    summary_df$yasso20 <- out$total_soc
  }
  
  # Q-MODEL
  if ("q_model" %in% models) {
    out <- q_model_run_hybrid_rcpp(
      params = params_q,  # <-- Custom params
      C0 = NULL, 
      input_df = q_input,
      site_row = site_row,
      spinup_years = if (!is.null(spinup_years)) spinup_years else 1500,
      verbose = FALSE
    )
    summary_df$q_model <- out$total_soc
  }
  
  # ROTHC
  if ("rothc" %in% models) {
    out <- rothc_run(
      params = params_rothc,  # <-- Custom params
      C0 = NULL, 
      input_df = ro_input,
      spinup_years = if (!is.null(spinup_years)) spinup_years else 10000
    )
    
    # Aggregate monthly to annual
    if (rothc_annual_stat == "mean") {
      summary_df$rothc <- aggregate_monthly_to_annual_mean(out$total_soc, ro_input$year)
    } else {
      summary_df$rothc <- aggregate_monthly_to_annual_eoy(out$total_soc, ro_input$year, ro_input$month)
    }
  }
  
  summary_df
}


# =============================================================================
# CALIBRATION-READY FAST WRAPPER (MULTIPLE PLOTS)
# =============================================================================

run_soc_models_multipplot_calibration <- function(
    input_df,
    site_df,
    models = c("yasso07", "yasso15", "yasso20", "rothc", "q_model"),
    params_list = NULL,  # NEW: custom parameters
    spinup_years = NULL,
    rothc_annual_stat = c("eoy", "mean"),
    validate = FALSE) {
  
  rothc_annual_stat <- match.arg(rothc_annual_stat)
  
  # Basic checks
  stopifnot(all(c("plot_id", "year") %in% names(input_df)))
  stopifnot("plot_id" %in% names(site_df))
  
  # Pre-split data once (major speedup for MCMC)
  input_split <- split(input_df, input_df$plot_id)
  site_split  <- split(site_df, site_df$plot_id)
  plot_ids <- intersect(names(input_split), names(site_split))
  
  # -------------------------------------------------------------------------
  # Run per plot with custom parameters
  # -------------------------------------------------------------------------
  res_by_plot <- lapply(plot_ids, function(pid) {
    
    dfp <- input_split[[pid]]
    
    # Ensure chronological ordering
    if ("month" %in% names(dfp)) {
      dfp <- dfp[order(dfp$year, dfp$month), ]
    } else {
      dfp <- dfp[order(dfp$year), ]
    }
    
    site_row <- site_split[[pid]]
    if (nrow(site_row) != 1) {
      stop("site_df must have 1 row per plot_id: ", pid, call. = FALSE)
    }
    
    # Call single-plot wrapper with custom parameters
    r <- run_soc_models_oneplot_calibration(
      input_df = dfp,
      site_row = site_row,
      models = models,
      params_list = params_list,  # <-- Pass through custom params
      spinup_years = spinup_years,
      rothc_annual_stat = rothc_annual_stat,
      validate = validate
    )
    
    # Add plot_id to output
    r$plot_id <- pid
    r <- r[, c("plot_id", setdiff(names(r), "plot_id")), drop = FALSE]
    
    r
  })
  
  names(res_by_plot) <- plot_ids
  
  # Optional: bind into single table for likelihood calculations
  summary_all <- dplyr::bind_rows(res_by_plot) %>%
    dplyr::arrange(plot_id, year)
  
  list(
    by_plot = res_by_plot,
    summary = summary_all  # Convenient for matching with observations
  )
}


# =============================================================================
# HELPER: Convert parameter vector to Q-model parameter list
# =============================================================================

vector_to_q_params <- function(param_vec) {
  # Assumes param_vec has named elements matching Q-model parameters
  # Example: c(eta11 = 0.36, beta = 7, e0 = 0.25, u00 = 0.0855, u01 = 0.0157)
  
  q_params <- Q_MODEL_DEFAULT_PARAMS  # Start with defaults
  
  # Update with calibrated parameters
  for (pname in names(param_vec)) {
    if (pname %in% names(q_params)) {
      q_params[[pname]] <- param_vec[pname]
    }
  }
  
  q_params
}


# =============================================================================
# EXAMPLE: BayesianTools likelihood function for Q-model
# =============================================================================

likelihood_q_model_example <- function(param_vec, 
                                       input_df, 
                                       site_df, 
                                       observations,
                                       spinup_years = 100,
                                       sigma = 5) {
  
  # Convert parameter vector to Q-model format
  params_q <- vector_to_q_params(param_vec)
  
  # Run model with custom parameters
  tryCatch({
    pred <- run_soc_models_multipplot_calibration(
      input_df = input_df,
      site_df = site_df,
      models = "q_model",
      params_list = list(q_model = params_q),  # <-- Custom parameters
      spinup_years = spinup_years
    )
    
    # Extract predictions at observation times/locations
    pred_summary <- pred$summary  # Has: plot_id, year, q_model
    
    # Match with observations
    matched <- merge(observations, pred_summary, 
                     by = c("plot_id", "year"),
                     all.x = TRUE)
    
    # Check for missing predictions
    if (any(is.na(matched$q_model))) {
      return(-Inf)  # Return very low likelihood if model failed
    }
    
    # Calculate log-likelihood (simple Gaussian)
    residuals <- matched$soc_obs_tCha - matched$q_model
    ll <- sum(dnorm(residuals, mean = 0, sd = sigma, log = TRUE))
    
    return(ll)
    
  }, error = function(e) {
    # If model crashes, return very low likelihood
    return(-Inf)
  })
}


# =============================================================================
# EXAMPLE USAGE
# =============================================================================

# Test with custom parameters (example: modify Q-model temperature sensitivity)
if (FALSE) {
  
  # Custom Q-model parameters
  custom_q_params <- Q_MODEL_DEFAULT_PARAMS
  custom_q_params$u01 <- 0.020  # Increase temperature sensitivity
  
  # Run with custom parameters
  results_custom <- run_soc_models_multipplot_calibration(
    input_df = input_data,
    site_df = site_data,
    models = "q_model",
    params_list = list(q_model = custom_q_params),
    spinup_years = 100
  )
  
  # Run with defaults for comparison
  results_default <- run_soc_models_multipplot_calibration(
    input_df = input_data,
    site_df = site_data,
    models = "q_model",
    params_list = NULL,  # Use defaults
    spinup_years = 100
  )
  
  # Compare
  plot(results_default$summary$year, results_default$summary$q_model, 
       type = "l", col = "black", 
       xlab = "Year", ylab = "SOC (t C/ha)",
       main = "Custom vs Default Parameters")
  lines(results_custom$summary$year, results_custom$summary$q_model, col = "red")
  legend("topright", c("Default", "Custom (u01=0.020)"), col = c("black", "red"), lty = 1)
}