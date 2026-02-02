# =============================================================================
# UNIFIED SOC MODEL COMPARISON WRAPPER 
# =============================================================================

source("./Models_functions/Input_matching_functions.R") # load the input matching functins to unify the input flow
source("./Models_functions/Yasso07.R") # load Yasso07 model functions
source("./Models_functions/Yasso15.R") # load Yasso15 model functions
source("./Models_functions/Yasso20.R") # load Yasso20 model functions
source("./Models_functions/RothC.R") # load RothC model functions
#source("./Models_functions/Q_transient_T_hybrid.R") # load Q-model functions and he wrapper with the hybrid initialization approach for Q-model
source("./Models_functions/Q_transient_T_hybrid_RCPP.R") # same but with compiled bits with rcpp



if(!require(Matrix)) install.packages("Matrix")
library(Matrix)

# -----------------------------------------------------------------------------
# MAIN WRAPPER FUNCTION
# -----------------------------------------------------------------------------

#' Run multiple SOC models on unified input data
#'
#' @param input_df Time series data (monthly or annual) with AWEN by cohort
#' @param site_row Single row from site_df with static site properties
#' @param models Character vector of models to run
#' @param spinup_years Optional override for models that use it (Yasso07/15/RothC/Q)
#' @param rothc_annual_stat How to annualise monthly RothC output for summary:
#'        "eoy" (December) or "mean" (annual mean of monthly states)
#' @return List with results from each model + results$summary
run_soc_models_oneplot <- function(input_df,
                           site_row = NULL,
                           models = c("yasso07", "yasso15", "yasso20", "rothc", "q_model"),
                           spinup_years = NULL,
                           rothc_annual_stat = c("eoy", "mean"),
                           verbose = FALSE,
                           validate = FALSE) {
  
  rothc_annual_stat <- match.arg(rothc_annual_stat)
  
  #validate inputs
  if(validate){validate_input(input_df)}
  
  results <- list()
  
  #if verbose, test the routing functions execution times
  if (verbose) {
    timings_map <- list()
  
  if ("yasso07" %in% models || "yasso15" %in% models) {
    tmap <- system.time({ yasso_input <- map_input_to_yasso07(input_df) })["elapsed"]
    timings_map$map_yasso07 <- tmap
    # reuse yasso_input for both models!
  }
  
  if ("yasso20" %in% models) {
    tmap <- system.time({ y20_input <- map_input_to_yasso20(input_df) })["elapsed"]
    timings_map$map_yasso20 <- tmap
  }
  
  if ("rothc" %in% models) {
    tmap <- system.time({ rothc_input <- map_input_to_rothc(input_df, site_row) })["elapsed"]
    timings_map$map_rothc <- tmap
  }
  
  if ("q_model" %in% models) {
    tmap <- system.time({ q_input <- map_input_to_q_model(input_df) })["elapsed"]
    timings_map$map_q <- tmap
  }
  
  print(timings_map)}
  
  
  # BENCHMARKING, open object to store execution times
  timings <- list()
  
  # ------------------------
  # YASSO07
  # ------------------------
  if ("yasso07" %in% models) {
    if(verbose){cat("Running Yasso07...\n")}
    
    yasso_input <- map_input_to_yasso07(input_df)
    yasso_spinup <- if (!is.null(spinup_years)) spinup_years else 5000
    
    t <- system.time({
      results$yasso07 <- yasso07_run(
        params = YASSO07_DEFAULT_PARAMS,
        C0 = NULL,
        input_df = yasso_input,
        spinup_years = yasso_spinup
      )
    })
    timings$yasso07 <- t["elapsed"]
    
    results$yasso07_input <- yasso_input
  }
  
  # ------------------------
  # YASSO15
  # ------------------------
  if ("yasso15" %in% models) {
    if(verbose){cat("Running Yasso15...\n")}
    
    y15_input <- map_input_to_yasso07(input_df) # alias of Yasso07 mapping
    y15_spinup <- if (!is.null(spinup_years)) spinup_years else 5000
    
    t <- system.time({
      results$yasso15 <- yasso15_run(
        params = YASSO15_DEFAULT_PARAMS,
        C0 = NULL,
        input_df = y15_input,
        spinup_years = y15_spinup
      )
    })
    timings$yasso15 <- t["elapsed"]
    
    results$yasso15_input <- y15_input
  }
  
  # ------------------------
  # YASSO20 (annual step, requires monthly forcing)
  # ------------------------
  if ("yasso20" %in% models) {
    if(verbose){cat("Running Yasso20...\n")}
    
    y20_input <- map_input_to_yasso20(input_df) # strict monthly-only

    t <- system.time({
      results$yasso20 <- yasso20_run(
        params = YASSO20_DEFAULT_PARAMS,
        C0 = NULL,
        input_df = y20_input,
        spinup = TRUE
      )
    })
    timings$yasso20 <- t["elapsed"]
    
    results$yasso20_input <- y20_input
  }
  
  # ------------------------
  # ROTHC
  # ------------------------
  if ("rothc" %in% models) {
    if(verbose){cat("Running RothC...\n")}
    
    rothc_input <- map_input_to_rothc(input_df, site_row)
    rothc_spinup <- if (!is.null(spinup_years)) spinup_years else 10000
    
    t <- system.time({
      results$rothc <- rothc_run(
        params = ROTHC_DEFAULT_PARAMS,
        C0 = NULL,
        input_df = rothc_input,
        spinup_years = rothc_spinup
      )
    })
    timings$rothc <- t["elapsed"]
    
    results$rothc_input <- rothc_input
  }
  
  # ------------------------
  # Q-MODEL
  # ------------------------
  if ("q_model" %in% models) {
    if(verbose){cat("Running Q-model...\n")}
    
    q_input <- map_input_to_q_model(input_df)
    q_spinup <- if (!is.null(spinup_years)) spinup_years else 1500
    
    t <- system.time({
      #results$q_model <- q_model_run_hybrid(
      results$q_model <- q_model_run_hybrid_rcpp(
        params = Q_MODEL_DEFAULT_PARAMS,
        C0 = NULL,
        input_df = q_input,
        site_row = site_row,
        spinup_years = q_spinup,
        verbose = verbose
      )
    })
    timings$q_model <- t["elapsed"]
    
    results$q_model_input <- q_input
  }
  
  # ------------------------
  # CREATE SUMMARY TABLE
  # ------------------------
  if (length(results) > 0) {
    
    # choose year vector from first annual model present, otherwise from RothC years
    years <- NULL
    if (!is.null(results$yasso07_input)) years <- years <- sort(unique(years)) #results$yasso07_input$year
    if (is.null(years) && !is.null(results$yasso15_input)) years <- results$yasso15_input$year
    if (is.null(years) && !is.null(results$yasso20_input)) years <- results$yasso20_input$year
    if (is.null(years) && !is.null(results$q_model_input))  years <- results$q_model_input$year
    if (is.null(years) && !is.null(results$rothc_input))    years <- unique(results$rothc_input$year)
    
    summary_df <- data.frame(year = years)
    
    if (!is.null(results$yasso07)) summary_df$yasso07 <- results$yasso07$total_soc
    if (!is.null(results$yasso15)) summary_df$yasso15 <- results$yasso15$total_soc
    if (!is.null(results$yasso20)) summary_df$yasso20 <- results$yasso20$total_soc
    if (!is.null(results$q_model))  summary_df$q_model  <- results$q_model$total_soc
    
    if (!is.null(results$rothc)) {
      ri <- results$rothc_input
      if (rothc_annual_stat == "mean") {
        summary_df$rothc <- aggregate_monthly_to_annual_mean(results$rothc$total_soc, ri$year)
      } else {
        if (!"month" %in% names(ri)) stop("RothC input missing `month`; cannot compute end-of-year.", call. = FALSE)
        summary_df$rothc <- aggregate_monthly_to_annual_eoy(results$rothc$total_soc, ri$year, ri$month)
      }
    }
    
    results$summary <- summary_df
  }
  
  # ------------------------
  # PRINT BENCHMARK SUMMARY
  # ------------------------
  if(verbose){  
  if (length(timings) > 0) {
    bench_df <- data.frame(
      model = names(timings),
      time_sec = round(unlist(timings), 3),
      row.names = NULL
    )
    
    cat("\nModel execution time (elapsed seconds):\n")
    print(bench_df, row.names = FALSE)}
  }
  
  
  results
}

# -----------------------------------------------------------------------------
# HELPERS (they aggregate monthly models to annual, currently just RothC)
# -----------------------------------------------------------------------------

# aggregate based on mean
# accessible via an option when running the wrapper,  rothc_annual_stat = c("eoy", "mean")
aggregate_monthly_to_annual_mean <- function(monthly_soc, years) {
  unique_years <- unique(years)
  sapply(unique_years, function(y) mean(monthly_soc[years == y], na.rm = TRUE))
}

# aggregate based on end-of-year (December), default
aggregate_monthly_to_annual_eoy <- function(monthly_soc, years, months) {
  unique_years <- unique(years)
  sapply(unique_years, function(y) {
    idx <- which(years == y & months == 12)
    if (length(idx) == 0) NA_real_ else monthly_soc[idx[1]]
  })
}






run_soc_models_multipplot <- function(input_df,
                                      site_df,
                                      models = c("yasso07","yasso15","yasso20","rothc","q_model"),
                                      spinup_years = NULL,
                                      rothc_annual_stat = c("eoy","mean"),
                                      verbose=FALSE,
                                      validate=FALSE) {
  
  rothc_annual_stat <- match.arg(rothc_annual_stat)
  
  # basic checks
  stopifnot(all(c("plot_id","year") %in% names(input_df)))
  stopifnot("plot_id" %in% names(site_df))
  
  plot_ids <- sort(unique(input_df$plot_id))
  
  # run per-plot
  res_by_plot <- lapply(plot_ids, function(pid) {
    
    dfp <- input_df %>%
      dplyr::filter(plot_id == pid) %>%
      dplyr::arrange(year, month)
    
    site_row <- site_df %>% dplyr::filter(plot_id == pid)
    
    if (nrow(site_row) != 1) {
      stop("site_df must have exactly 1 row per plot_id. plot_id=", pid,
           " has ", nrow(site_row), " rows.", call. = FALSE)
    }
    
    # IMPORTANT: pass the single-row site_row to Q + RothC mappings
    r <- run_soc_models_oneplot(
      input_df = dfp,
      site_row = site_row,
      models = models,
      spinup_years = spinup_years,
      rothc_annual_stat = rothc_annual_stat,
      verbose = verbose,
      validate = validate
    )
    
    # ensure summary exists and add plot_id (some of your models may drop it)
    if (!is.null(r$summary)) {
      r$summary$plot_id <- pid
      r$summary <- r$summary %>% dplyr::select(plot_id, dplyr::everything())
    }
    
    r
  })
  
  names(res_by_plot) <- plot_ids
  
  # bind summaries into one table: plot_id x year
  summary_all <- dplyr::bind_rows(lapply(res_by_plot, `[[`, "summary")) %>%
    dplyr::arrange(plot_id, year)
  
  list(
    by_plot = res_by_plot,
    summary = summary_all
  )
}





# =============================================================================
# FAST / LITE SOC WRAPPER – SINGLE PLOT
# =============================================================================
#
# Purpose:
#   Optimised wrapper intended for Bayesian calibration and fast model comparison.
#   It runs ONLY the selected model(s) and returns a minimal annual summary
#   (year + total SOC), avoiding all unnecessary allocations.
#
# Key differences vs standard wrapper:
#   - Does NOT return full model objects (pools, internal states, respiration).
#   - Does NOT store or return mapped inputs.
#   - Does NOT print timings or verbose messages.
#   - Immediately discards everything except total SOC.
#
# What it still guarantees:
#   - Uses the same mapping functions and model code as the standard wrapper.
#   - Identical numerical results for total SOC (given same parameters).
#   - Identical spinup logic for each model.
#
# Suitable for:
#   - BayesianTools likelihood evaluations
#   - Fast predictive model comparison (SOC trajectories, RMSE, log-lik)
#
# Not suitable for:
#   - Pool-level diagnostics
#   - Debugging mappings or internal dynamics
#   - Analyses requiring respiration or pool composition
# =============================================================================
run_soc_models_oneplot_fast <- function(input_df, site_row=NULL,
                                        models=c("yasso07","yasso15","yasso20","rothc","q_model"),
                                        spinup_years=NULL,
                                        rothc_annual_stat=c("eoy","mean"),
                                        validate=FALSE) {
  
  rothc_annual_stat <- match.arg(rothc_annual_stat)
  if (validate) validate_input(input_df)
  
  # --- mappings (local only; don't store) ---
  y07_input <- NULL
  if (any(c("yasso07","yasso15") %in% models)) {
    y07_input <- map_input_to_yasso07(input_df)
  }
  y20_input <- if ("yasso20" %in% models) map_input_to_yasso20(input_df) else NULL
  ro_input  <- if ("rothc"  %in% models) map_input_to_rothc(input_df, site_row) else NULL
  q_input   <- if ("q_model" %in% models) map_input_to_q_model(input_df) else NULL
  
  # --- years ---
  years <- NULL
  if (!is.null(y07_input)) years <- y07_input$year
  if (is.null(years) && !is.null(y20_input)) years <- y20_input$year
  if (is.null(years) && !is.null(q_input))  years <- q_input$year
  if (is.null(years) && !is.null(ro_input)) years <- sort(unique(ro_input$year))
  
  summary_df <- data.frame(year = years)
  
  # --- run models; immediately extract total_soc; drop the rest ---
  if ("yasso07" %in% models) {
    out <- yasso07_run(YASSO07_DEFAULT_PARAMS, NULL, y07_input,
                       spinup_years = if (!is.null(spinup_years)) spinup_years else 5000)
    summary_df$yasso07 <- out$total_soc
  }
  
  if ("yasso15" %in% models) {
    out <- yasso15_run(YASSO15_DEFAULT_PARAMS, NULL, y07_input,
                       spinup_years = if (!is.null(spinup_years)) spinup_years else 5000)
    summary_df$yasso15 <- out$total_soc
  }
  
  if ("yasso20" %in% models) {
    out <- yasso20_run(YASSO20_DEFAULT_PARAMS, NULL, y20_input, spinup = TRUE)
    summary_df$yasso20 <- out$total_soc
  }
  
  if ("q_model" %in% models) {
    out <- q_model_run_hybrid_rcpp(Q_MODEL_DEFAULT_PARAMS, NULL, q_input,
                                   site_row = site_row,
                                   spinup_years = if (!is.null(spinup_years)) spinup_years else 1500,
                                   verbose = FALSE)
    summary_df$q_model <- out$total_soc
  }
  
  if ("rothc" %in% models) {
    out <- rothc_run(ROTHC_DEFAULT_PARAMS, NULL, ro_input,
                     spinup_years = if (!is.null(spinup_years)) spinup_years else 10000)
    
    if (rothc_annual_stat == "mean") {
      summary_df$rothc <- aggregate_monthly_to_annual_mean(out$total_soc, ro_input$year)
    } else {
      summary_df$rothc <- aggregate_monthly_to_annual_eoy(out$total_soc, ro_input$year, ro_input$month)
    }
  }
  
  summary_df
}


# =============================================================================
# FAST / LITE SOC WRAPPER – MULTIPLE PLOTS
# =============================================================================
#
# Purpose:
#   High-throughput multi-plot driver built on run_soc_models_oneplot_fast().
#   Designed for MCMC and large Monte Carlo experiments.
#
# Key differences vs standard multi-plot wrapper:
#   - Pre-splits input and site data ONCE (no repeated filtering in likelihood).
#   - Runs per-plot models using the lite single-plot wrapper.
#   - Returns only compact annual SOC summaries per plot.
#   - No storage of intermediate model outputs or mapped inputs.
#
# What it still guarantees:
#   - Identical SOC predictions to the standard workflow.
#   - Correct plot-wise handling of site-level parameters.
#
# Output:
#   - List of per-plot annual SOC data frames
#   - (plot_id, year, SOC columns only)
#
# Suitable for:
#   - BayesianTools samplers (single-core or parallel)
#   - Ensemble simulations across thousands of plots
#
# Not suitable for:
#   - Pool diagnostics across plots
#   - Model debugging or validation of internal states
# =============================================================================
run_soc_models_multipplot_fast <- function(input_df,
                                           site_df,
                                           models = c("yasso07","yasso15","yasso20","rothc","q_model"),
                                           spinup_years = NULL,
                                           rothc_annual_stat = c("eoy","mean"),
                                           verbose=FALSE,
                                           validate=FALSE) {
  
  rothc_annual_stat <- match.arg(rothc_annual_stat)
  
  # basic checks
  stopifnot(all(c("plot_id","year") %in% names(input_df)))
  stopifnot("plot_id" %in% names(site_df))
  
  # --- Pre-split data once ---
  # This avoids repeatedly filtering the full data frame for each plot.
  input_split <- split(input_df, input_df$plot_id)
  site_split  <- split(site_df,  site_df$plot_id)
  
  # Use only plot_ids that are present in BOTH input and site tables
  plot_ids <- intersect(names(input_split), names(site_split))
  
  # --- Run per plot ---
  # Each iteration pulls the plot-specific data frame in O(1) and calls the fast single-plot wrapper.
  res_by_plot <- lapply(plot_ids, function(pid) {
    
    dfp <- input_split[[pid]]
    
    # Ensure chronological ordering (monthly or annual)
    if ("month" %in% names(dfp)) {
      dfp <- dfp[order(dfp$year, dfp$month), ]
    } else {
      dfp <- dfp[order(dfp$year), ]
    }
    
    site_row <- site_split[[pid]]
    if (nrow(site_row) != 1) stop("site_df must have 1 row per plot_id: ", pid, call.=FALSE)
    
    # Call fast per-plot wrapper (returns a small year x model table)
    r <- run_soc_models_oneplot_fast(
      input_df = dfp,
      site_row = site_row,
      models = models,
      spinup_years = spinup_years,
      rothc_annual_stat = rothc_annual_stat,
      validate = validate
    )
    
    # NOTE:
    # This function returns a data.frame already, not a list containing $summary.
    # If you want plot_id attached, do it here:
    r$plot_id <- pid
    r <- r[, c("plot_id", setdiff(names(r), "plot_id")), drop = FALSE]
    
    r
  })
  
  names(res_by_plot) <- plot_ids
  
  # Optional: bind into one table (often convenient for likelihood calculations)
  # summary_all <- dplyr::bind_rows(res_by_plot) |> dplyr::arrange(plot_id, year)
  
  list(
    by_plot = res_by_plot
    # , summary = summary_all
  )
}