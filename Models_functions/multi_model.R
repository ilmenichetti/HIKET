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
                           rothc_annual_stat = c("eoy", "mean")) {
  
  rothc_annual_stat <- match.arg(rothc_annual_stat)
  
  validate_input(input_df)
  
  results <- list()
  
  
  # BENCHMARKING, open object to store execution times
  timings <- list()
  
  # ------------------------
  # YASSO07
  # ------------------------
  if ("yasso07" %in% models) {
    cat("Running Yasso07...\n")
    
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
    cat("Running Yasso15...\n")
    
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
    cat("Running Yasso20...\n")
    
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
    cat("Running RothC...\n")
    
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
    cat("Running Q-model...\n")
    
    q_input <- map_input_to_q_model(input_df)
    q_spinup <- if (!is.null(spinup_years)) spinup_years else 1500
    
    t <- system.time({
      #results$q_model <- q_model_run_hybrid(
      results$q_model <- q_model_run_hybrid_rcpp(
        params = Q_MODEL_DEFAULT_PARAMS,
        C0 = NULL,
        input_df = q_input,
        site_row = site_row,
        spinup_years = q_spinup
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
  if (length(timings) > 0) {
    bench_df <- data.frame(
      model = names(timings),
      time_sec = round(unlist(timings), 3),
      row.names = NULL
    )
    
    cat("\nModel execution time (elapsed seconds):\n")
    print(bench_df, row.names = FALSE)
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
                                      rothc_annual_stat = c("eoy","mean")) {
  
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
      rothc_annual_stat = rothc_annual_stat
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
