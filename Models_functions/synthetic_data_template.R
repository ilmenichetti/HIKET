# =============================================================================
# TESTING
# =============================================================================
#
# This script is for testing the whole SOC model comparison framework.
# It generates synthetic data, runs the models and compares results.
#
# =============================================================================

# Load required libraries
library(dplyr)
library(tidyr)

#load multi-model wrapper functions
source("./Models_functions/multi_model.R")


# =============================================================================
# STEP 1: Generate Synthetic Input Data
# =============================================================================

# Generate monthly input data for 2 plots, 2000-2020
generate_input_template <- function(plot_ids = c("PLOT_001", "PLOT_002"),
                                    years = 2000:2020) {
  
  grid <- expand.grid(
    plot_id = plot_ids,
    year = years,
    month = 1:12,
    stringsAsFactors = FALSE
  )
  grid <- grid[order(grid$plot_id, grid$year, grid$month), ]
  n <- nrow(grid)
  
  # AWEN proportions by cohort
  awen_props <- list(
    foliage      = c(A = 0.52, W = 0.10, E = 0.02, N = 0.36),
    branches     = c(A = 0.50, W = 0.05, E = 0.03, N = 0.42),
    stems        = c(A = 0.45, W = 0.02, E = 0.01, N = 0.52),
    fine_roots   = c(A = 0.44, W = 0.05, E = 0.01, N = 0.50),
    coarse_roots = c(A = 0.45, W = 0.03, E = 0.01, N = 0.51),
    understorey  = c(A = 0.50, W = 0.15, E = 0.03, N = 0.32)
  )
  
  # Annual litter totals (Mg C/ha/yr)
  litter_totals <- list(
    foliage = 1.2, branches = 0.4, stems = 0.1,
    fine_roots = 0.8, coarse_roots = 0.3, understorey = 0.3
  )
  
  cohorts <- names(litter_totals)
  
  # Generate litter columns (monthly = annual/12)
  for (coh in cohorts) {
    total <- litter_totals[[coh]] / 12
    props <- awen_props[[coh]]
    
    for (frac in c("A", "W", "E", "N")) {
      col_name <- paste0("C_", coh, "_", frac)
      base_value <- total * props[frac]
      
      # Variation by plot and year
      plot_effect <- ifelse(grid$plot_id == plot_ids[1], 1.0, 1.1)
      year_effect <- 1 + 0.1 * sin(2 * pi * (grid$year - min(years)) / 10)
      
      grid[[col_name]] <- round(base_value * plot_effect * year_effect, 6)
    }
  }
  
  # Climate
  base_temp <- ifelse(grid$plot_id == plot_ids[1], 3.5, 2.0)
  grid$temp_air <- round(base_temp + 12 * sin(2 * pi * (grid$month - 4) / 12) + rnorm(n, 0, 1), 1)
  
  monthly_precip <- c(40, 35, 35, 40, 45, 55, 70, 75, 60, 55, 50, 45)
  grid$precip <- round(monthly_precip[grid$month] * 
                         ifelse(grid$plot_id == plot_ids[1], 1.0, 1.1) * 
                         (1 + rnorm(n, 0, 0.2)), 1)
  
  grid$evap <- round(pmax(5, 20 + 40 * sin(2 * pi * (grid$month - 4) / 12) + rnorm(n, 0, 5)), 1)
  
  grid
}

# Generate site data
generate_site_template <- function(plot_ids = c("PLOT_001")) {  # Only one plot ID
  data.frame(
    plot_id = plot_ids,
    pine_prop = 0.3,      # Single value
    spruce_prop = 0.6,
    birch_prop = 0.1,
    clay = 18,
    soil_depth = 23,
    stringsAsFactors = FALSE
  )
}


# Generate data
input_data <- generate_input_template(plot_ids = "PLOT_001", years = 2000:2020)
site_data <- generate_site_template(plot_ids = "PLOT_001")

cat(sprintf("Generated input data: %d rows (months)\n", nrow(input_data)))
cat(sprintf("Time span: %d-%d\n", min(input_data$year), max(input_data$year)))

# =============================================================================
# STEP 2: Run Models
# =============================================================================

# Run all three models
results <- run_soc_models(
  input_df = input_data,
  site_row = site_data,
  models = c("yasso07", "rothc", "q_model"),
  spinup_years = 200  # Shorter spinup for demonstration
)

cat("\nModels completed successfully!\n")



# =============================================================================
# STEP 3: Examine Results
# =============================================================================

print(head(results$summary, 10))
cat("\n")
print(tail(results$summary, 10))

# Calculate statistics
stats <- calculate_model_stats(results)
print(stats)

# =============================================================================
# STEP 4: Visualize Results
# =============================================================================
plot_model_comparison(results)


# =============================================================================
# ADDITIONAL ANALYSES
# =============================================================================


# Turnover times
for(model_name in c("yasso07", "rothc", "q_model")) {
  if(model_name %in% names(results$summary)) {
    avg_soc <- mean(results$summary[[model_name]])
    
    # Get annual input
    if(model_name == "yasso07") {
      mapped <- map_input_to_yasso07(input_data)
      annual_input <- mean(mapped$C_nwl + mapped$C_fwl + mapped$C_cwl)
    } else if(model_name == "q_model") {
      mapped <- map_input_to_q_model(input_data)
      annual_input <- mean(mapped$C_needles + mapped$C_branches + mapped$C_stems + 
                             mapped$C_fine_roots + mapped$C_understorey)
    } else {
      mapped <- map_input_to_rothc(input_data, site_data)
      annual_input <- mean(aggregate(mapped$C_DPM + mapped$C_RPM, 
                                     by=list(mapped$year), FUN=sum)$x)
    }
    
    turnover <- avg_soc / annual_input
    cat(sprintf("  %s: %.1f years\n", model_name, turnover))
  }
}

# Model divergence
soc_range <- apply(results$summary[, -1], 1, function(x) {
  (max(x) - min(x)) / mean(x) * 100
})
cat(sprintf("  Mean CV across years: %.1f%%\n", mean(soc_range)))
cat(sprintf("  Range: %.1f%% - %.1f%%\n", min(soc_range), max(soc_range)))


# =============================================================================
# WRITE THE TEMPLATE WITH THE INPUT DATA
# =============================================================================
write.csv(input_data, "synthetic_input_data_template.csv", row.names = FALSE)
write.csv(site_data, "synthetic_site_data_template.csv", row.names = FALSE)
