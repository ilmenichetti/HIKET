# =============================================================================
# UNIFIED SOC MODEL COMPARISON WRAPPER (extended: Yasso07/Yasso15/Yasso20/RothC/Q)
# =============================================================================

source("./Models_functions/Input_matching_functions.R")
source("./Models_functions/Yasso07.R")
source("./Models_functions/Yasso15.R")
source("./Models_functions/Yasso20.R")
source("./Models_functions/RothC.R")
source("./Models_functions/Q_transient_T.R")
source("./Models_functions/Q_transient_T_hybrid.R")

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
run_soc_models <- function(input_df,
                           site_row = NULL,
                           models = c("yasso07", "yasso15", "yasso20", "rothc", "q_model"),
                           spinup_years = NULL,
                           rothc_annual_stat = c("eoy", "mean")) {
  
  rothc_annual_stat <- match.arg(rothc_annual_stat)
  
  validate_input(input_df)
  
  results <- list()
  
  # ------------------------
  # YASSO07
  # ------------------------
  if ("yasso07" %in% models) {
    cat("Running Yasso07...\n")
    
    yasso_input <- map_input_to_yasso07(input_df)
    yasso_spinup <- if (!is.null(spinup_years)) spinup_years else 5000
    
    results$yasso07 <- yasso07_run(
      params = YASSO07_DEFAULT_PARAMS,
      C0 = NULL,
      input_df = yasso_input,
      spinup_years = yasso_spinup
    )
    results$yasso07_input <- yasso_input
  }
  
  # ------------------------
  # YASSO15
  # ------------------------
  if ("yasso15" %in% models) {
    cat("Running Yasso15...\n")
    
    y15_input <- map_input_to_yasso07(input_df) # alias of Yasso07 mapping
    y15_spinup <- if (!is.null(spinup_years)) spinup_years else 5000
    
    results$yasso15 <- yasso15_run(
      params = YASSO15_DEFAULT_PARAMS,
      C0 = NULL,
      input_df = y15_input,
      spinup_years = y15_spinup
    )
    results$yasso15_input <- y15_input
  }
  
  # ------------------------
  # YASSO20 (annual step, requires monthly forcing)
  # ------------------------
  if ("yasso20" %in% models) {
    cat("Running Yasso20...\n")
    
    y20_input <- map_input_to_yasso20(input_df) # strict monthly-only
    results$yasso20 <- yasso20_run(
      params = YASSO20_DEFAULT_PARAMS,
      C0 = NULL,
      input_df = y20_input,
      spinup = TRUE
    )
    results$yasso20_input <- y20_input
  }
  
  # ------------------------
  # ROTHC
  # ------------------------
  if ("rothc" %in% models) {
    cat("Running RothC...\n")
    
    rothc_input <- map_input_to_rothc(input_df, site_row)
    rothc_spinup <- if (!is.null(spinup_years)) spinup_years else 10000
    
    results$rothc <- rothc_run(
      params = ROTHC_DEFAULT_PARAMS,
      C0 = NULL,
      input_df = rothc_input,
      spinup_years = rothc_spinup
    )
    results$rothc_input <- rothc_input
  }
  
  # ------------------------
  # Q-MODEL
  # ------------------------
  if ("q_model" %in% models) {
    cat("Running Q-model...\n")
    
    q_input <- map_input_to_q_model(input_df)
    q_spinup <- if (!is.null(spinup_years)) spinup_years else 1500
    
    results$q_model <- q_model_run_hybrid(
      params = Q_MODEL_DEFAULT_PARAMS,
      C0 = NULL,
      input_df = q_input,
      site_row = site_row,
      spinup_years = q_spinup
    )
    results$q_model_input <- q_input
  }
  
  # ------------------------
  # CREATE SUMMARY TABLE
  # ------------------------
  if (length(results) > 0) {
    
    # choose year vector from first annual model present, otherwise from RothC years
    years <- NULL
    if (!is.null(results$yasso07_input)) years <- results$yasso07_input$year
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
  
  results
}

# -----------------------------------------------------------------------------
# HELPERS
# -----------------------------------------------------------------------------

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

