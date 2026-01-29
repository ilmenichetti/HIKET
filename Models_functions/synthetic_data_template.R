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
  models = c("yasso07", "yasso15", "yasso20", "rothc", "q_model"),
  spinup_years = 200  # Shorter spinup for demonstration
)

cat("\nModels completed successfully!\n")



# =============================================================================
# STEP 3: Examine Results
# =============================================================================

print(head(results$summary, 10))
cat("\n")
print(tail(results$summary, 10))



# =============================================================================
# ADDITIONAL ANALYSES (ALL MODELS + 2 FIGURES)
# =============================================================================

# ---- Helper: annual input (Mg C/ha/yr) for each model ----
get_annual_input <- function(model_name, input_data, site_data) {
  
  if (model_name %in% c("yasso07", "yasso15")) {
    mapped <- map_input_to_yasso07(input_data)
    return(mean(mapped$C_nwl + mapped$C_fwl + mapped$C_cwl, na.rm = TRUE))
  }
  
  if (model_name == "yasso20") {
    mapped <- map_input_to_yasso20(input_data)  # strict monthly forcing required
    nwl <- mapped$nwl_A + mapped$nwl_W + mapped$nwl_E + mapped$nwl_N
    fwl <- mapped$fwl_A + mapped$fwl_W + mapped$fwl_E + mapped$fwl_N
    cwl <- mapped$cwl_A + mapped$cwl_W + mapped$cwl_E + mapped$cwl_N
    return(mean(nwl + fwl + cwl, na.rm = TRUE))
  }
  
  if (model_name == "q_model") {
    mapped <- map_input_to_q_model(input_data)
    return(mean(mapped$C_needles + mapped$C_branches + mapped$C_stems +
                  mapped$C_fine_roots + mapped$C_understorey, na.rm = TRUE))
  }
  
  if (model_name == "rothc") {
    mapped <- map_input_to_rothc(input_data, site_data)
    annual <- aggregate(mapped$C_DPM + mapped$C_RPM, by = list(year = mapped$year), FUN = sum)
    return(mean(annual$x, na.rm = TRUE))
  }
  
  stop(sprintf("Unknown model for annual input: %s", model_name))
}

# ---- 1) Turnover times for ALL models in the summary ----
model_cols <- setdiff(names(results$summary), "year")

turnover_df <- data.frame(
  model = model_cols,
  mean_soc = NA_real_,
  annual_input = NA_real_,
  turnover_years = NA_real_
)

for (m in model_cols) {
  avg_soc <- mean(results$summary[[m]], na.rm = TRUE)
  ann_in <- get_annual_input(m, input_data, site_data)
  
  turnover_df[turnover_df$model == m, "mean_soc"] <- avg_soc
  turnover_df[turnover_df$model == m, "annual_input"] <- ann_in
  turnover_df[turnover_df$model == m, "turnover_years"] <- avg_soc / ann_in
}

cat("\nTurnover times (mean SOC / mean annual input):\n")
print(turnover_df)

# ---- 2) Divergence across models per year ----
soc_divergence <- apply(results$summary[, model_cols, drop = FALSE], 1, function(x) {
  (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) / mean(x, na.rm = TRUE) * 100
})
cat(sprintf("\nModel divergence across years (range/mean):\n  Mean: %.1f%%\n  Range: %.1f%% - %.1f%%\n",
            mean(soc_divergence, na.rm = TRUE),
            min(soc_divergence, na.rm = TRUE),
            max(soc_divergence, na.rm = TRUE)))

# =============================================================================
# FIGURE 1: SOC time series comparison (all models)
# =============================================================================

png("./Figures/syntetic_data_soc_model_comparison_C_timeseries.png", width = 800, height = 600)

years <- results$summary$year

ylim_all <- range(results$summary[, model_cols, drop = FALSE], na.rm = TRUE)

plot(years, results$summary[[model_cols[1]]],
     type = "l", lwd = 2,
     xlab = "Year", ylab = "SOC (Mg C/ha)",
     ylim = ylim_all,
     main = "SOC model comparison (annual)")

if (length(model_cols) > 1) {
  for (i in 2:length(model_cols)) {
    lines(years, results$summary[[model_cols[i]]], lwd = 2, lty = i)
  }
}

legend("topleft",
       legend = model_cols,
       lty = seq_along(model_cols),
       lwd = 2,
       bty = "n")

dev.off()


# =============================================================================
# SET OF FIGUREs 2: Pairwise scatter vs a reference model (all models vs ref)
# =============================================================================

model_cols <- setdiff(names(results$summary), "year")

# Choose which models to cycle as references
ref_models <- c("yasso07", "yasso15", "yasso20", "rothc")
ref_models <- ref_models[ref_models %in% model_cols]

for (ref_model in ref_models) {
  
  # One PNG per reference model
  fn <- sprintf("./Figures/synthetic_data_soc_scatter_ref_%s.png", ref_model)
  png(fn, width = 1200, height = 800)
  
  others <- setdiff(model_cols, ref_model)
  x <- results$summary[[ref_model]]
  
  op <- par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1))
  
  for (m in head(others, 4)) {
    
    y <- results$summary[[m]]
    ok <- is.finite(x) & is.finite(y)
    
    plot(x[ok], y[ok],
         xlab = paste(ref_model, "SOC"),
         ylab = paste(m, "SOC"),
         main = paste(m, "vs", ref_model))
    
    abline(0, 1, lty = 2)
    
    rmse <- sqrt(mean((y[ok] - x[ok])^2))
    bias <- mean(y[ok] - x[ok])
    corv <- suppressWarnings(cor(x[ok], y[ok]))
    
    mtext(sprintf("r = %.2f   RMSE = %.2f   bias = %.2f",
                  corv, rmse, bias),
          side = 3, line = -1.2, cex = 0.8)
  }
  
  mtext(paste("SOC comparison – reference:", ref_model),
        outer = TRUE, line = -1, cex = 1.2)
  
  par(op)
  dev.off()
  
  cat("Saved:", fn, "\n")
}

# =============================================================================
# WRITE THE TEMPLATE WITH THE INPUT DATA
# =============================================================================
write.csv(input_data, "synthetic_input_data_template.csv", row.names = FALSE)
write.csv(site_data, "synthetic_site_data_template.csv", row.names = FALSE)
