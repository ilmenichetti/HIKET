# =============================================================================
# UNIFIED SOC MODEL COMPARISON WRAPPER
# =============================================================================
#
# Runs multiple SOC models with consistent inputs for structural comparison.
#
# Models supported:
#   - Yasso07: Chemical fractionation (AWENH pools)
#   - RothC:   DPM/RPM pools with clay-dependent stabilization
#   - Q-model: Continuous quality theory (requires species composition)
# =============================================================================

source("./Models_functions/Input_matching_functions.R")
source("./Models_functions/Yasso07.R")
source("./Models_functions/RothC.R")
source("./Models_functions/Q_transient_T.R")
source("./Models_functions/Q_transient_T_hybrid.R")



# Load required libraries
if(!require(Matrix)) install.packages("Matrix")
library(Matrix)

# -----------------------------------------------------------------------------
# MAIN WRAPPER FUNCTION
# -----------------------------------------------------------------------------

#' Run multiple SOC models on unified input data
#'
#' @param input_df Time series data (monthly or annual) with AWEN by cohort
#' @param site_row Single row from site_df with static site properties
#' @param models Character vector of models to run: "yasso07", "rothc", "q_model"
#' @param spinup_years Override default spinup duration
#' @return List with results from each model
run_soc_models <- function(input_df, 
                           site_row = NULL,
                           models = c("yasso07", "rothc", "q_model"),
                           spinup_years = NULL) {
  
  # Validate input
  validate_input(input_df)
  
  results <- list()
  
  # ------------------------
  # YASSO07
  # ------------------------
  if ("yasso07" %in% models) {
    cat("Running Yasso07...\n")
    
    # Map input to Yasso07 format (annual, 3 litter types)
    yasso_input <- map_input_to_yasso07(input_df)
    
    # Set spinup years
    yasso_spinup <- if (!is.null(spinup_years)) spinup_years else 5000
    
    # Run model
    results$yasso07 <- yasso07_run(
      params = YASSO07_DEFAULT_PARAMS,
      C0 = NULL,
      input_df = yasso_input,
      spinup_years = yasso_spinup
    )
  }
  
  # ------------------------
  # ROTHC
  # ------------------------
  if ("rothc" %in% models) {
    cat("Running RothC...\n")
    
    # Map input to RothC format (monthly, DPM/RPM split)
    rothc_input <- map_input_to_rothc(input_df, site_row)
    
    # Set spinup years
    rothc_spinup <- if (!is.null(spinup_years)) spinup_years else 10000
    
    # Run model
    results$rothc <- rothc_run(
      params = ROTHC_DEFAULT_PARAMS,
      C0 = NULL,
      input_df = rothc_input,
      spinup_years = rothc_spinup
    )
  }
  
  # ------------------------
  # Q-MODEL
  # ------------------------
  if ("q_model" %in% models) {
    cat("Running Q-model...\n")
    
    # Map input to Q-model format (annual, 5 litter types)
    q_input <- map_input_to_q_model(input_df)
    
    # Set spinup years
    q_spinup <- if (!is.null(spinup_years)) spinup_years else 1500
    
    # Run model
    #results$q_model <- q_model_run( #older model version without transient temperature
    #results$q_model <- q_model_run_transient(
    results$q_model <- q_model_run_hybrid(
        params = Q_MODEL_DEFAULT_PARAMS,
      C0 = NULL,
      input_df = q_input,
      site_row = site_row,
      spinup_years = q_spinup
    )
  }
  
  # ------------------------
  # CREATE SUMMARY TABLE
  # ------------------------
  if (length(results) > 0) {
    # Get years from first model result
    years <- if ("yasso07" %in% names(results)) {
      yasso_input$year
    } else if ("q_model" %in% names(results)) {
      q_input$year
    } else if ("rothc" %in% names(results)) {
      unique(rothc_input$year)
    }
    
    summary_df <- data.frame(year = years)
    
    if ("yasso07" %in% names(results)) {
      summary_df$yasso07 <- results$yasso07$total_soc
    }
    if ("rothc" %in% names(results)) {
      # Aggregate monthly to annual for comparison
      summary_df$rothc <- aggregate_monthly_to_annual(
        results$rothc$total_soc, 
        rothc_input$year
      )
    }
    if ("q_model" %in% names(results)) {
      summary_df$q_model <- results$q_model$total_soc
    }
    
    results$summary <- summary_df
  }
  
  results
}


# -----------------------------------------------------------------------------
# HELPER FUNCTIONS
# -----------------------------------------------------------------------------

#' Aggregate monthly SOC values to annual
aggregate_monthly_to_annual <- function(monthly_soc, years) {
  unique_years <- unique(years)
  annual_soc <- sapply(unique_years, function(y) {
    mean(monthly_soc[years == y])
  })
  annual_soc
}

#' Plot model comparison
plot_model_comparison <- function(results) {
  if (!"summary" %in% names(results)) {
    stop("No summary available. Run models first.")
  }
  
  if (!require(ggplot2)) {
    warning("ggplot2 not available, using base plot")
    
    summary_df <- results$summary
    years <- summary_df$year
    
    plot(years, summary_df[[2]], type = "l", 
         ylim = range(summary_df[, -1], na.rm = TRUE),
         xlab = "Year", ylab = "SOC (Mg C/ha)",
         main = "SOC Model Comparison",
         col = 1, lwd = 2)
    
    for (i in 3:ncol(summary_df)) {
      lines(years, summary_df[[i]], col = i-1, lwd = 2)
    }
    
    legend("topright", 
           legend = names(summary_df)[-1],
           col = 1:(ncol(summary_df)-1),
           lwd = 2)
    
  } else {
    library(ggplot2)
    library(tidyr)
    
    plot_df <- results$summary %>%
      pivot_longer(-year, names_to = "model", values_to = "soc")
    
    ggplot(plot_df, aes(x = year, y = soc, color = model)) +
      geom_line(linewidth = 1) +
      labs(
        title = "SOC Model Comparison",
        x = "Year",
        y = "SOC (Mg C/ha)",
        color = "Model"
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")
  }
}

#' Calculate model statistics
calculate_model_stats <- function(results) {
  if (!"summary" %in% names(results)) {
    stop("No summary available. Run models first.")
  }
  
  summary_df <- results$summary
  model_cols <- setdiff(names(summary_df), "year")
  
  stats <- data.frame(
    model = model_cols,
    mean_soc = sapply(model_cols, function(m) mean(summary_df[[m]], na.rm = TRUE)),
    min_soc = sapply(model_cols, function(m) min(summary_df[[m]], na.rm = TRUE)),
    max_soc = sapply(model_cols, function(m) max(summary_df[[m]], na.rm = TRUE)),
    sd_soc = sapply(model_cols, function(m) sd(summary_df[[m]], na.rm = TRUE)),
    cv_soc = sapply(model_cols, function(m) {
      sd(summary_df[[m]], na.rm = TRUE) / mean(summary_df[[m]], na.rm = TRUE) * 100
    })
  )
  
  rownames(stats) <- NULL
  stats
}

#' Export results to CSV
export_results <- function(results, output_dir = ".") {
  if (!"summary" %in% names(results)) {
    stop("No summary available. Run models first.")
  }
  
  # Summary table
  write.csv(results$summary, 
            file.path(output_dir, "soc_model_summary.csv"),
            row.names = FALSE)
  
  # Individual model outputs
  if ("yasso07" %in% names(results)) {
    yasso_out <- data.frame(
      year = results$summary$year,
      results$yasso07$pools,
      total_soc = results$yasso07$total_soc,
      respiration = results$yasso07$respiration
    )
    write.csv(yasso_out,
              file.path(output_dir, "yasso07_detailed.csv"),
              row.names = FALSE)
  }
  
  if ("rothc" %in% names(results)) {
    rothc_out <- data.frame(
      results$rothc$pools,
      total_soc = results$rothc$total_soc,
      respiration = results$rothc$respiration
    )
    write.csv(rothc_out,
              file.path(output_dir, "rothc_detailed.csv"),
              row.names = FALSE)
  }
  
  if ("q_model" %in% names(results)) {
    q_out <- data.frame(
      year = results$summary$year,
      results$q_model$pools,
      total_soc = results$q_model$total_soc,
      respiration = results$q_model$respiration
    )
    write.csv(q_out,
              file.path(output_dir, "q_model_detailed.csv"),
              row.names = FALSE)
  }
  
  cat(sprintf("Results exported to: %s\n", output_dir))
  invisible(TRUE)
}
