# Q-MODEL FUNCTIONS - COHORT-BASED DECOMPOSITION WITH TRANSIENT SPINUP
# ==============================================================================
# This implementation tracks individual litter cohorts through time, calculating
# how each cohort decomposes based on Q-model theory (Bosatta & Ågren 1991-2003).
# A 1500-year transient simulation is used to establish initial carbon stocks.


# Q-MODEL PARAMETERS ===========================================================

# Species-specific decomposition parameters
# q0n, q0w, q0f: initial quality parameters for needles, woody tissue, and fine roots
# eta11, beta, e0: microbial efficiency and quality transformation parameters
# tmaxb, tmaxs: maximum decomposition times for branches and stems
# u00, u01: base microbial growth rate and temperature modifier
q_model_params <- data.frame(
  Species = c("Pine", "Spruce", "Birch"),
  q0n = c(1.1, 1.01, 1.01),           
  q0w = c(1.06, 1.00, 1.00),          
  q0f = c(1.036, 1.036, 1.036),       
  eta11 = c(0.36, 0.36, 0.36),
  beta = c(7, 7, 7),
  e0 = c(0.25, 0.25, 0.25),
  tmaxb = 13,                          
  tmaxs = 60,                          
  u00 = 0.0855,                        
  u01 = 0.0157,                        
  stringsAsFactors = FALSE
)

# DECOMPOSITION FUNCTIONS ======================================================
# These calculate the remaining fraction of carbon after time t for different
# litter types, based on Q-model theory (Bosatta & Ågren 2003)

#' Needles and understory decomposition (simple power law)
calc_gn <- function(alpha_n, year_val, z) {
  return(1 / (1 + alpha_n * year_val)^z)
}

#' Fine roots decomposition (simple power law)
calc_gf <- function(alpha_f, year_val, z) {
  return(1 / (1 + alpha_f * year_val)^z)
}

#' Branch decomposition (accounts for gradual litter fall over tmaxb years)
calc_gb <- function(alpha_w, year_val, z, ib, tmaxb) {
  term1 <- 2 * ((1 + alpha_w * year_val)^(1 - z) - 
                  (1 + alpha_w * (year_val - ib))^(1 - z) * (1 - ib / tmaxb)) / 
    (alpha_w * (1 - z) * tmaxb)
  
  term2 <- 2 * ((1 + alpha_w * (year_val - ib))^(2 - z) - 
                  (1 + alpha_w * year_val)^(2 - z)) / 
    (alpha_w^2 * (1 - z) * (2 - z) * tmaxb^2)
  
  term3 <- (1 - ib / tmaxb)^2
  
  return(term1 + term2 + term3)
}

#' Stem decomposition (accounts for gradual litter fall over tmaxs years)
calc_gs <- function(alpha_w, year_val, z, is, tmaxs) {
  term1 <- 2 * ((1 + alpha_w * year_val)^(1 - z) - 
                  (1 + alpha_w * (year_val - is))^(1 - z) * (1 - is / tmaxs)) / 
    (alpha_w * (1 - z) * tmaxs)
  
  term2 <- 2 * ((1 + alpha_w * (year_val - is))^(2 - z) - 
                  (1 + alpha_w * year_val)^(2 - z)) / 
    (alpha_w^2 * (1 - z) * (2 - z) * tmaxs^2)
  
  term3 <- (1 - is / tmaxs)^2
  
  return(term1 + term2 + term3)
}



# MAIN Q-MODEL FUNCTION (REFACTORED) ===========================================

FUNC_q_decomp <- function(litter_input, soil, weather, species_proportions) {
  
  # STEP 1 & 2: PREPARE LITTER AND WEATHER DATA (Largely unchanged) =============
  
  if(soil == 1) {
    all_litter <- litter_input %>% 
      filter(soil == 1) %>% 
      ungroup() %>% 
      select(-soil)  
  } else {
    all_litter <- litter_input %>% 
      filter(soil == 2, litter_source %in% c("log", "nat")) %>% 
      ungroup() %>% 
      select(-soil)  
  }
  
  # Spinup weather: 1960-1990 averages for equilibrium calculation
  spinup_weather <- weather %>% 
    filter(year < 1991) %>% 
    group_by(region) %>% 
    summarize(mean_T = mean(mean_T)) %>% 
    ungroup()
  
  # Simulation weather: 30-year rolling mean to smooth interannual variation
  weather_roll <- weather %>%
    group_by(region) %>%
    mutate(roll_T = rollmean(mean_T, 30, align = "right", fill = NA)) %>%
    select(year, region, roll_T) %>%
    complete(year = min(year):max(litter_input$year)) %>%
    fill(roll_T, .direction = "downup")
  
  # STEP 3: PROCESS LITTER INPUTS (Unchanged) =================================
  
  if(soil == 1) {
    processed_litter <- all_litter %>% 
      filter(year > 1971) %>% 
      group_by(region, litter_type, year) %>% 
      summarise_at(c("A", "W", "E", "N", "H"), sum, na.rm = TRUE) %>%
      mutate(total_carbon = A + W + E + N + H) %>%
      ungroup()
  } else {
    processed_litter <- all_litter %>% 
      filter(!(litter_source == "nat" & !litter_type == "cwl"), year > 1971) %>% 
      group_by(region, litter_type, year) %>% 
      summarise_at(c("A", "W", "E", "N", "H"), sum, na.rm = TRUE) %>% 
      mutate(total_carbon = A + W + E + N + H) %>%
      ungroup()
  }
  
  # STEP 4: CALCULATE CONSTANT INPUTS FOR TRANSIENT SPINUP ====================
  # Use the average of the first 5 years of simulation data
  
  start_year_sim <- min(processed_litter$year)
  spinup_avg_years <- start_year_sim:(start_year_sim + 4)
  
  spinup_litter <- processed_litter %>%
    filter(year %in% spinup_avg_years) %>%
    group_by(region, litter_type) %>%
    summarise(total_carbon = mean(total_carbon, na.rm = TRUE)) %>%
    ungroup()
  
  # Convert to wide format
  spinup_wide <- spinup_litter %>%
    pivot_wider(names_from = litter_type, values_from = total_carbon, values_fill = 0) %>%
    left_join(spinup_weather, by = "region")
  
  # STEP 5: MAIN SIMULATION WITH INTEGRATED TRANSIENT SPINUP ==================
  # The model is run continuously for 1500 spinup years + N simulation years.
  
  spinup_duration <- 1500 # Years for the transient spinup
  
  # Prepare simulation data in wide format
  carbon_wide <- processed_litter %>%
    select(region, year, litter_type, total_carbon) %>%
    pivot_wider(names_from = litter_type, values_from = total_carbon, values_fill = 0) %>%
    left_join(weather_roll, by = c("region", "year"))
  
  all_results <- NULL
  
  # Loop over regions
  for(region_idx in 1:nrow(species_proportions)) {
    region_data_sim <- carbon_wide %>% filter(region == region_idx)
    n_years_sim <- nrow(region_data_sim)
    
    if(n_years_sim == 0) next
    
    cat(paste("Running transient spinup and simulation for region", region_idx, "...\n"))
    
    region_spinup_inputs <- spinup_wide %>% filter(region == region_idx)
    region_results <- NULL
    
    # Loop over species
    for(species_name in colnames(species_proportions)) {
      prop <- species_proportions[region_idx, species_name]
      if(prop <= 0) next
      
      params <- q_model_params[q_model_params$Species == species_name, ]
      
      # --- Create a combined dataframe for the entire run (spinup + simulation) ---
      
      # 1. Create data for the 1500-year spinup with constant inputs/climate
      spinup_df <- data.frame(
        year = (start_year_sim - spinup_duration):(start_year_sim - 1),
        region = region_idx,
        nwl = ifelse("nwl" %in% names(region_spinup_inputs), region_spinup_inputs$nwl, 0),
        fwl = ifelse("fwl" %in% names(region_spinup_inputs), region_spinup_inputs$fwl, 0),
        cwl = ifelse("cwl" %in% names(region_spinup_inputs), region_spinup_inputs$cwl, 0),
        roll_T = region_spinup_inputs$mean_T
      )
      
      # 2. Combine with the actual simulation data
      combined_data <- bind_rows(spinup_df, region_data_sim)
      total_years <- nrow(combined_data)
      
      # 3. Calculate species-specific litter inputs and decomp parameters for the whole period
      species_data_full <- combined_data %>%
        mutate(
          Ln = nwl * prop,
          Lf = nwl * prop * 0.3,
          Lb = fwl * prop,
          Ls = cwl * prop,
          Lu = nwl * prop * 0.1,
          fC = 0.5,
          mat = roll_T,
          u0 = params$u00 + params$u01 * mat,
          z = (1 - params$e0) / (params$beta * params$eta11 * params$e0),
          alfan = fC * params$beta * params$eta11 * u0 * params$q0n^params$beta,
          alfaf = fC * params$beta * params$eta11 * u0 * params$q0f^params$beta,
          alfaw = fC * params$beta * params$eta11 * u0 * params$q0w^params$beta
        )
      
      # --- Run the full simulation (spinup + main) ---
      
      # Initialize cohort tracking matrices for the full duration
      matrix_gn <- matrix(0, nrow = total_years, ncol = total_years)
      matrix_gf <- matrix(0, nrow = total_years, ncol = total_years)
      matrix_gb <- matrix(0, nrow = total_years, ncol = total_years)
      matrix_gs <- matrix(0, nrow = total_years, ncol = total_years)
      matrix_gu <- matrix(0, nrow = total_years, ncol = total_years)
      
      # Fill matrices: each input year's litter is tracked through all subsequent years
      for(input_year in 1:total_years) {
        
        # Decomp parameters are constant for spinup, variable for main simulation
        alfan <- species_data_full$alfan[input_year]
        alfaf <- species_data_full$alfaf[input_year]
        alfaw <- species_data_full$alfaw[input_year]
        z <- species_data_full$z[input_year]
        
        # Track this cohort through all future years
        for(sim_year in input_year:total_years) {
          decomp_time <- sim_year - input_year + 1
          
          ib <- min(decomp_time, params$tmaxb)
          is <- min(decomp_time, params$tmaxs)
          
          gn_val <- calc_gn(alfan, decomp_time, z)
          gf_val <- calc_gf(alfaf, decomp_time, z)
          gb_val <- calc_gb(alfaw, decomp_time, z, ib, params$tmaxb)
          gs_val <- calc_gs(alfaw, decomp_time, z, is, params$tmaxs)
          gu_val <- calc_gn(alfan, decomp_time, z)
          
          matrix_gn[sim_year, input_year] <- species_data_full$Ln[input_year] * gn_val
          matrix_gf[sim_year, input_year] <- species_data_full$Lf[input_year] * gf_val
          matrix_gb[sim_year, input_year] <- species_data_full$Lb[input_year] * gb_val
          matrix_gs[sim_year, input_year] <- species_data_full$Ls[input_year] * gs_val
          matrix_gu[sim_year, input_year] <- species_data_full$Lu[input_year] * gu_val
        }
      }
      
      # --- Extract results for the main simulation period ONLY ---
      
      species_results_full <- data.frame(
        region = region_idx,
        year = species_data_full$year,
        species = species_name
      ) %>%
        mutate(
          # Total stock is now the simple sum of all living cohorts
          c_total = rowSums(matrix_gn) + rowSums(matrix_gf) + rowSums(matrix_gb) +
            rowSums(matrix_gs) + rowSums(matrix_gu)
        )
      
      # Filter to keep only the years from the original input data
      species_results <- species_results_full %>%
        filter(year >= start_year_sim)
      
      region_results <- rbind(region_results, species_results)
    }
    
    all_results <- rbind(all_results, region_results)
  }
  
  # STEP 6: AGGREGATE AND CALCULATE EMISSION FACTORS (Unchanged logic) ==========
  
  q_final <- all_results %>%
    group_by(region, year) %>%
    summarise(c_total = sum(c_total), .groups = 'drop') %>%
    arrange(region, year) %>%
    group_by(region) %>%
    mutate(
      c_change = c_total - lag(c_total),
      EF_min = c_change
    ) %>%
    filter(!is.na(c_change)) %>%
    ungroup()
  
  # STEP 7: FORMAT OUTPUT (Unchanged logic) ====================================
  
  if(soil == 1) {
    q_final_processed <- q_final %>%
      arrange(region, year) %>%
      group_by(region) %>%
      mutate(
        EF_min = rollapply(EF_min, 5, mean, partial = TRUE),
        c_stock = c_total
      ) %>%
      select(region, year, EF_min, c_stock) %>%
      ungroup()
  } else {
    q_final_processed <- q_final %>%
      arrange(region, year) %>%
      group_by(region) %>%
      mutate(c_stock = c_total) %>%
      select(region, year, c_stock) %>%
      filter(year > 1989) %>%
      ungroup()
  }
  
  return(q_final_processed)
}