# =============================================================================
# Compare Old vs New Q-model Implementation
# =============================================================================

library(dplyr)
library(tidyr)
library(zoo)

# -----------------------------------------------------------------------------
# Convert synthetic data to old Q-model format
# -----------------------------------------------------------------------------

convert_synthetic_to_old_format <- function(plot_data_mapped, site_data) {
  
  # plot_data_mapped should already be the output from map_input_to_q_model()
  # which has: year, temp_mean, C_needles, C_branches, C_stems, C_fine_roots, C_understorey
  
  # The old function expects:
  # litter_input: columns = region, year, litter_type (nwl/fwl/cwl), A, W, E, N, H, soil, litter_source
  # weather: columns = region, year, mean_T (and potentially ampli_T, sum_P for spinup)
  # species_proportions: data frame with columns Pine, Spruce, Birch (one row per region)
  
  # 1. Create litter_input in long format with AWENH
  # We need to partition total carbon into AWENH using typical ratios
  
  # Standard AWENH partitions (from Yasso literature)
  partition_needles <- c(A=0.52, W=0.17, E=0.12, N=0.08, H=0.11)
  partition_woody <- c(A=0.45, W=0.10, E=0.08, N=0.25, H=0.12)
  partition_roots <- c(A=0.50, W=0.15, E=0.13, N=0.10, H=0.12)
  
  # Create litter for each type
  litter_nwl <- plot_data_mapped %>%
    mutate(
      region = 1,
      litter_type = "nwl",
      litter_source = "tree",
      soil = 1,
      total_C = C_needles + C_understorey,
      A = total_C * partition_needles["A"],
      W = total_C * partition_needles["W"],
      E = total_C * partition_needles["E"],
      N = total_C * partition_needles["N"],
      H = total_C * partition_needles["H"]
    ) %>%
    select(region, year, litter_type, litter_source, soil, A, W, E, N, H)
  
  litter_fwl <- plot_data_mapped %>%
    mutate(
      region = 1,
      litter_type = "fwl",
      litter_source = "tree",
      soil = 1,
      total_C = C_branches,
      A = total_C * partition_woody["A"],
      W = total_C * partition_woody["W"],
      E = total_C * partition_woody["E"],
      N = total_C * partition_woody["N"],
      H = total_C * partition_woody["H"]
    ) %>%
    select(region, year, litter_type, litter_source, soil, A, W, E, N, H)
  
  litter_cwl <- plot_data_mapped %>%
    mutate(
      region = 1,
      litter_type = "cwl",
      litter_source = "tree",
      soil = 1,
      total_C = C_stems,
      A = total_C * partition_woody["A"],
      W = total_C * partition_woody["W"],
      E = total_C * partition_woody["E"],
      N = total_C * partition_woody["N"],
      H = total_C * partition_woody["H"]
    ) %>%
    select(region, year, litter_type, litter_source, soil, A, W, E, N, H)
  
  # Note: Fine roots are handled within the old Q-model as a fraction of nwl
  # The old model does: Lf = nwl * prop * 0.3
  # So we DON'T include C_fine_roots separately in litter_input
  
  litter_input <- bind_rows(litter_nwl, litter_fwl, litter_cwl)
  
  # 2. Create weather data
  weather <- plot_data_mapped %>%
    mutate(
      region = 1,
      mean_T = temp_mean,
      ampli_T = 10,  # Typical value
      sum_P = 600    # Typical annual precipitation
    ) %>%
    select(region, year, mean_T, ampli_T, sum_P)
  
  # Add years 1960-1990 for spinup (use average of simulation period)
  spinup_years_vec <- 1960:min(plot_data_mapped$year - 1)
  if(length(spinup_years_vec) > 0) {
    spinup_weather <- data.frame(
      region = 1,
      year = spinup_years_vec,
      mean_T = mean(weather$mean_T),
      ampli_T = 10,
      sum_P = 600
    )
    weather <- bind_rows(spinup_weather, weather)
  }
  
  # 3. Create species_proportions
  species_proportions <- data.frame(
    Pine = site_data$pine_prop,
    Spruce = site_data$spruce_prop,
    Birch = site_data$birch_prop
  )
  
  list(
    litter_input = litter_input,
    weather = weather,
    species_proportions = species_proportions
  )
}

# -----------------------------------------------------------------------------
# Run comparison
# -----------------------------------------------------------------------------

compare_q_models <- function(plot_data_mapped, site_data) {
  
  cat("=== COMPARING OLD VS NEW Q-MODEL IMPLEMENTATIONS ===\n\n")
  
  # Prepare data for old model
  cat("1. Converting data to old Q-model format...\n")
  old_format <- convert_synthetic_to_old_format(plot_data_mapped, site_data)
  
  cat("   Litter input rows:", nrow(old_format$litter_input), "\n")
  cat("   Weather rows:", nrow(old_format$weather), "\n")
  cat("   Species props:", paste(names(old_format$species_proportions), 
                                  unlist(old_format$species_proportions[1,]), 
                                  collapse=", "), "\n\n")
  
  # Run OLD model
  cat("2. Running OLD Q-model (FUNC_q_decomp)...\n")
  old_results <- FUNC_q_decomp(
    litter_input = old_format$litter_input,
    soil = 1,
    weather = old_format$weather,
    species_proportions = old_format$species_proportions
  )
  
  cat("   Old model output rows:", nrow(old_results), "\n")
  cat("   Old model columns:", paste(names(old_results), collapse=", "), "\n")
  cat("   Old model SOC range:", round(min(old_results$c_stock), 2), "-", 
      round(max(old_results$c_stock), 2), "Mg C/ha\n\n")
  
  # Run NEW model
  cat("3. Running NEW Q-model (q_model_run)...\n")
  new_results <- q_model_run(
    params = Q_MODEL_DEFAULT_PARAMS,
    input_df = plot_data_mapped,
    site_row = site_data,
    spinup_years = 1500
  )
  
  cat("   New model SOC length:", length(new_results$total_soc), "\n")
  cat("   New model SOC range:", round(min(new_results$total_soc), 2), "-", 
      round(max(new_results$total_soc), 2), "Mg C/ha\n\n")
  
  # Compare results
  cat("4. Comparing results...\n\n")
  
  # Match years (old model may have filtered some years)
  comparison <- data.frame(
    year = plot_data_mapped$year,
    new_soc = new_results$total_soc
  ) %>%
    left_join(
      old_results %>% select(year, old_soc = c_stock),
      by = "year"
    )
  
  # Calculate statistics
  comparison <- comparison %>%
    filter(!is.na(old_soc)) %>%
    mutate(
      diff = new_soc - old_soc,
      ratio = new_soc / old_soc,
      pct_diff = (new_soc - old_soc) / old_soc * 100
    )
  
  cat("   Years compared:", nrow(comparison), "\n")
  cat("   Mean old SOC:", round(mean(comparison$old_soc), 2), "Mg C/ha\n")
  cat("   Mean new SOC:", round(mean(comparison$new_soc), 2), "Mg C/ha\n")
  cat("   Mean difference:", round(mean(comparison$diff), 2), "Mg C/ha\n")
  cat("   Mean ratio (new/old):", round(mean(comparison$ratio), 3), "\n")
  cat("   Mean % difference:", round(mean(comparison$pct_diff), 1), "%\n\n")
  
  # Print first and last few years
  cat("   First 5 years:\n")
  print(head(comparison %>% select(year, old_soc, new_soc, diff, pct_diff), 5))
  cat("\n   Last 5 years:\n")
  print(tail(comparison %>% select(year, old_soc, new_soc, diff, pct_diff), 5))
  
  # Plot comparison
  cat("\n5. Generating comparison plot...\n")
  
  if(require(ggplot2, quietly = TRUE)) {
    p <- ggplot(comparison, aes(x = year)) +
      geom_line(aes(y = old_soc, color = "Old"), linewidth = 1) +
      geom_line(aes(y = new_soc, color = "New"), linewidth = 1) +
      labs(
        title = "Q-model Implementation Comparison",
        x = "Year",
        y = "SOC (Mg C/ha)",
        color = "Implementation"
      ) +
      theme_minimal()
    
    print(p)
  }
  
  invisible(list(
    comparison = comparison,
    old_results = old_results,
    new_results = new_results
  ))
}

# -----------------------------------------------------------------------------
# USAGE INSTRUCTIONS
# -----------------------------------------------------------------------------

