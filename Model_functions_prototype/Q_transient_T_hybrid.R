# =============================================================================
# Q-MODEL: HYBRID INITIALIZATION (STANDALONE, COMMENTED)
# =============================================================================
# Combines Swedish steady-state initialization with short spinup for speed.
# Uses transient (current-year) temperature for decomposition calculations.
#
# METHOD:
# 1. Calculate analytical steady-state SOC (C_ss)
# 2. Start with C_ss as bulk pool: C(t) = C_ss * (1 + alpha0*t)^(1-z)
# 3. Run short spinup (10-20 years) building explicit cohort tracking
# 4. During simulation: Total SOC = bulk_pool(t) + new_cohorts(t)
#
# SPEEDUP: 10-20 year spinup vs 1500 years = 75-150x faster
# =============================================================================

# =============================================================================
# SPECIES-SPECIFIC PARAMETERS
# =============================================================================
# Initial quality values (q0) for different tissue types by species
# Higher q0 = higher quality = faster decomposition

Q_SPECIES_PARAMS <- data.frame(
  species = c("pine", "spruce", "birch"),
  q0n = c(1.10, 1.01, 1.01),    # Needles/leaves quality
  q0w = c(1.06, 1.00, 1.00),    # Woody tissue quality (branches/stems)
  q0f = c(1.036, 1.036, 1.036), # Fine roots quality
  stringsAsFactors = FALSE
)

# =============================================================================
# MODEL PARAMETERS
# =============================================================================
# Default Q-model parameters controlling decomposition rates

Q_MODEL_DEFAULT_PARAMS <- list(
  eta11 = 0.36,        # Microbial efficiency
  beta = 7,            # Quality decline rate exponent
  e0 = 0.25,           # Initial decomposer efficiency
  u00 = 0.0855,        # Base microbial growth rate (1/year)
  u01 = 0.0157,        # Temperature sensitivity of growth (1/year/°C)
  tmaxb = 13,          # Branch fragmentation time (years)
  tmaxs = 60,          # Stem fragmentation time (years)
  fC = 0.5,            # Carbon fraction in dry matter
  spinup_years = 1500  # Default spinup (unused in hybrid method)
)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Calculate z parameter: controls how quality affects decomposition
calc_z <- function(e0, beta, eta11) {
  (1 - e0) / (beta * eta11 * e0)
}

# Calculate remaining fraction for simple tissues (transient temperature)
# Uses cumulative sum of alpha over cohort lifetime
calc_g_simple_transient <- function(alpha_vector, z) {
  cumulative_alpha_t <- sum(alpha_vector)  # Sum decomposition potential
  if(cumulative_alpha_t <= 0) return(1)
  1 / (1 + cumulative_alpha_t)^z           # Q-model decay formula
}

# Calculate remaining fraction for branches (with fragmentation)
# Uses average alpha as approximation for transient temperature
calc_g_branches_transient <- function(alpha_vector, decomp_time, z, ib, tmaxb) {
  if(decomp_time <= 0) return(1)
  
  alpha_avg <- mean(alpha_vector)  # Average decomposition rate
  if(alpha_avg <= 0) return(1)
  
  # Three-term fragmentation formula (Ågren & Bosatta 1998)
  term1 <- 2*((1 + alpha_avg*decomp_time)^(1-z) - 
              (1 + alpha_avg*(decomp_time-ib))^(1-z)*(1-ib/tmaxb)) / 
    (alpha_avg*(1-z)*tmaxb)
  
  term2 <- 2*((1 + alpha_avg*(decomp_time-ib))^(2-z) - 
              (1 + alpha_avg*decomp_time)^(2-z)) / 
    (alpha_avg^2*(1-z)*(2-z)*tmaxb^2)
  
  term3 <- (1 - ib/tmaxb)^2
  
  max(0, term1 + term2 + term3)
}

# Calculate remaining fraction for stems (similar to branches, longer time)
calc_g_stems_transient <- function(alpha_vector, decomp_time, z, is, tmaxs) {
  if(decomp_time <= 0) return(1)
  
  alpha_avg <- mean(alpha_vector)
  if(alpha_avg <= 0) return(1)
  
  # Same formula as branches but with stem fragmentation time
  term1 <- 2*((1 + alpha_avg*decomp_time)^(1-z) - 
              (1 + alpha_avg*(decomp_time-is))^(1-z)*(1-is/tmaxs)) / 
    (alpha_avg*(1-z)*tmaxs)
  
  term2 <- 2*((1 + alpha_avg*(decomp_time-is))^(2-z) - 
              (1 + alpha_avg*decomp_time)^(2-z)) / 
    (alpha_avg^2*(1-z)*(2-z)*tmaxs^2)
  
  term3 <- (1 - is/tmaxs)^2
  
  max(0, term1 + term2 + term3)
}

# =============================================================================
# STEADY-STATE CALCULATIONS
# =============================================================================
# Analytical equilibrium SOC for given inputs (used to initialize bulk pool)

# Simple tissues: needles, fine roots, understorey
calc_steady_state_simple <- function(L_input, alpha, z) {
  if(alpha <= 0 || L_input <= 0) return(0)
  if(abs(z - 1) < 0.001) {
    warning("z is very close to 1")
    return(L_input / alpha * 10)
  }
  L_input / (alpha * (z - 1))  # Analytical steady-state
}

# Branches: numerical integration over age distribution
calc_steady_state_branches <- function(L_input, alpha, z, tmaxb) {
  if(alpha <= 0 || L_input <= 0) return(0)
  
  max_age <- 500  # Integrate to 500 years
  ages <- 1:max_age
  
  # Calculate remaining fraction at each age
  g_values <- sapply(ages, function(t) {
    ib <- min(t, tmaxb)
    term1 <- 2*((1 + alpha*t)^(1-z) - (1 + alpha*(t-ib))^(1-z)*(1-ib/tmaxb)) / 
      (alpha*(1-z)*tmaxb)
    term2 <- 2*((1 + alpha*(t-ib))^(2-z) - (1 + alpha*t)^(2-z)) / 
      (alpha^2*(1-z)*(2-z)*tmaxb^2)
    term3 <- (1 - ib/tmaxb)^2
    max(0, term1 + term2 + term3)
  })
  
  L_input * sum(g_values)  # Steady-state = input × sum(g)
}

# Stems: same as branches but longer integration
calc_steady_state_stems <- function(L_input, alpha, z, tmaxs) {
  if(alpha <= 0 || L_input <= 0) return(0)
  
  max_age <- 1000  # Integrate to 1000 years (stems last longer)
  ages <- 1:max_age
  
  g_values <- sapply(ages, function(t) {
    is <- min(t, tmaxs)
    term1 <- 2*((1 + alpha*t)^(1-z) - (1 + alpha*(t-is))^(1-z)*(1-is/tmaxs)) / 
      (alpha*(1-z)*tmaxs)
    term2 <- 2*((1 + alpha*(t-is))^(2-z) - (1 + alpha*t)^(2-z)) / 
      (alpha^2*(1-z)*(2-z)*tmaxs^2)
    term3 <- (1 - is/tmaxs)^2
    max(0, term1 + term2 + term3)
  })
  
  L_input * sum(g_values)
}

# =============================================================================
# MAIN HYBRID Q-MODEL FUNCTION
# =============================================================================

q_model_run_hybrid <- function(params = Q_MODEL_DEFAULT_PARAMS, 
                               C0 = NULL, 
                               input_df, 
                               site_row = NULL, 
                               spinup_years = 20,
                               verbose = FALSE) {
  
  # ---------------------------------------------------------------------------
  # Setup dimensions
  # ---------------------------------------------------------------------------
  n_years_sim <- nrow(input_df)  # Simulation period length
  n_spinup <- spinup_years        # Short spinup period (10-20 years)
  
  if(verbose) {
    cat(sprintf("Q-model HYBRID initialization:\n"))
    cat(sprintf("  Spinup: %d years\n", n_spinup))
    cat(sprintf("  Simulation: %d years\n\n", n_years_sim))
  }
  
  # ---------------------------------------------------------------------------
  # Extract species proportions
  # ---------------------------------------------------------------------------
  if(is.null(site_row)) {
    # Default to pure spruce if no site data
    species_props <- c(pine = 0, spruce = 1, birch = 0)
  } else {
    # Get proportions from site data
    species_props <- c(
      pine = as.numeric(site_row$pine_prop),
      spruce = as.numeric(site_row$spruce_prop),
      birch = as.numeric(site_row$birch_prop)
    )
  }
  
  # ---------------------------------------------------------------------------
  # Extract litter inputs and temperature
  # ---------------------------------------------------------------------------
  required_cols <- c("C_needles", "C_branches", "C_fine_roots", "C_understorey")
  missing <- setdiff(required_cols, names(input_df))
  if (length(missing) > 0) {
    stop(sprintf("Missing columns: %s", paste(missing, collapse = ", ")))
  }
  
  C_needles <- as.numeric(input_df$C_needles)
  C_branches <- as.numeric(input_df$C_branches)
  C_stems <- if ("C_stems" %in% names(input_df)) {
    as.numeric(input_df$C_stems)
  } else {
    rep(0, n_years_sim)
  }
  C_fine_roots <- as.numeric(input_df$C_fine_roots)
  C_understorey <- as.numeric(input_df$C_understorey)
  temp_mean_vec <- as.numeric(input_df$temp_mean)
  
  # ---------------------------------------------------------------------------
  # Initialize output storage
  # ---------------------------------------------------------------------------
  total_soc <- numeric(n_years_sim)
  pools_mat <- matrix(0, n_years_sim, 5, 
                      dimnames = list(NULL, c("needles","fine_roots","branches","stems","understorey")))
  initial_pool_soc <- numeric(n_years_sim)  # Track bulk pool contribution
  
  # ---------------------------------------------------------------------------
  # Loop over species (pine, spruce, birch)
  # ---------------------------------------------------------------------------
  for(sp in names(species_props)) {
    prop <- as.numeric(species_props[sp])
    if(prop <= 0) next  # Skip species not present
    
    # Get species-specific parameters
    sp_params <- Q_SPECIES_PARAMS[Q_SPECIES_PARAMS$species == sp, ]
    if(nrow(sp_params) == 0) next
    
    if(verbose) cat(sprintf("%s (%.1f%%)\n", sp, prop * 100))
    
    # -------------------------------------------------------------------------
    # Calculate averages for steady-state initialization
    # -------------------------------------------------------------------------
    spinup_avg_years <- 1:min(5, n_years_sim)  # Use first 5 years
    avg_needles <- mean(C_needles[spinup_avg_years])
    avg_fine_roots <- mean(C_fine_roots[spinup_avg_years])
    avg_branches <- mean(C_branches[spinup_avg_years])
    avg_stems <- mean(C_stems[spinup_avg_years])
    avg_understorey <- mean(C_understorey[spinup_avg_years])
    avg_temp <- mean(temp_mean_vec)
    
    # -------------------------------------------------------------------------
    # Calculate decomposition parameters at average conditions
    # -------------------------------------------------------------------------
    u0_avg <- params$u00 + params$u01 * avg_temp  # Microbial growth rate
    z <- calc_z(params$e0, params$beta, params$eta11)  # Quality exponent
    
    # Extract species quality parameters
    q0n <- as.numeric(sp_params$q0n)  # Needles quality
    q0f <- as.numeric(sp_params$q0f)  # Fine roots quality
    q0w <- as.numeric(sp_params$q0w)  # Woody tissue quality
    
    # Calculate alpha (decomposition rate) for each tissue
    alfan_avg <- params$fC * params$beta * params$eta11 * u0_avg * q0n^params$beta
    alfaf_avg <- params$fC * params$beta * params$eta11 * u0_avg * q0f^params$beta
    alfaw_avg <- params$fC * params$beta * params$eta11 * u0_avg * q0w^params$beta
    
    # -------------------------------------------------------------------------
    # Calculate analytical steady-state SOC
    # -------------------------------------------------------------------------
    SOC_ss_needles <- calc_steady_state_simple(avg_needles * prop, alfan_avg, z)
    SOC_ss_fine_roots <- calc_steady_state_simple(avg_fine_roots * prop, alfaf_avg, z)
    SOC_ss_understorey <- calc_steady_state_simple(avg_understorey * prop, alfan_avg, z)
    SOC_ss_branches <- calc_steady_state_branches(avg_branches * prop, alfaw_avg, z, params$tmaxb)
    SOC_ss_stems <- calc_steady_state_stems(avg_stems * prop, alfaw_avg, z, params$tmaxs)
    
    SOC_ss_total <- SOC_ss_needles + SOC_ss_fine_roots + SOC_ss_branches + 
                    SOC_ss_stems + SOC_ss_understorey
    
    if(verbose) cat(sprintf("  Steady-state: %.1f gC/m²\n", SOC_ss_total))
    
    # -------------------------------------------------------------------------
    # Calculate alpha0 for bulk pool decomposition (Swedish formula)
    # -------------------------------------------------------------------------
    # This controls how fast the initial bulk pool decomposes
    i0 <- (avg_needles + avg_fine_roots + avg_branches + avg_stems + avg_understorey) * prop
    alpha0 <- -i0 / (SOC_ss_total * (1 - z))
    
    # -------------------------------------------------------------------------
    # Create combined timeline: spinup + simulation
    # -------------------------------------------------------------------------
    total_years <- n_spinup + n_years_sim
    
    # Litter inputs: constant during spinup, actual during simulation
    Ln <- c(rep(avg_needles * prop, n_spinup), C_needles * prop)
    Lf <- c(rep(avg_fine_roots * prop, n_spinup), C_fine_roots * prop)
    Lb <- c(rep(avg_branches * prop, n_spinup), C_branches * prop)
    Ls <- c(rep(avg_stems * prop, n_spinup), C_stems * prop)
    Lu <- c(rep(avg_understorey * prop, n_spinup), C_understorey * prop)
    
    # Temperature: constant during spinup, actual during simulation
    mat <- c(rep(avg_temp, n_spinup), temp_mean_vec)
    
    # -------------------------------------------------------------------------
    # Calculate time-varying decomposition rates
    # -------------------------------------------------------------------------
    u0 <- params$u00 + params$u01 * mat  # Temperature-dependent growth
    alfan <- params$fC * params$beta * params$eta11 * u0 * q0n^params$beta
    alfaf <- params$fC * params$beta * params$eta11 * u0 * q0f^params$beta
    alfaw <- params$fC * params$beta * params$eta11 * u0 * q0w^params$beta
    
    # -------------------------------------------------------------------------
    # Initialize cohort tracking matrices
    # -------------------------------------------------------------------------
    # Rows = simulation years, Columns = input years
    matrix_gn <- matrix(0, nrow = total_years, ncol = total_years)
    matrix_gf <- matrix(0, nrow = total_years, ncol = total_years)
    matrix_gb <- matrix(0, nrow = total_years, ncol = total_years)
    matrix_gs <- matrix(0, nrow = total_years, ncol = total_years)
    matrix_gu <- matrix(0, nrow = total_years, ncol = total_years)
    
    # -------------------------------------------------------------------------
    # Fill matrices: track each cohort's decomposition
    # -------------------------------------------------------------------------
    # For each input year, calculate remaining C at each future simulation year
    for(input_year in 1:total_years) {
      for(sim_year in input_year:total_years) {
        
        decomp_time <- sim_year - input_year + 1  # Years since input
        
        # Fragmentation progress
        ib <- min(decomp_time, params$tmaxb)  # Branch fragmentation
        is <- min(decomp_time, params$tmaxs)  # Stem fragmentation
        
        # Extract alpha vector: cohort experiences CURRENT temps, not input temp
        alpha_years <- input_year:sim_year
        alfan_vector <- alfan[alpha_years]
        alfaf_vector <- alfaf[alpha_years]
        alfaw_vector <- alfaw[alpha_years]
        
        # Calculate remaining fraction using transient temperature
        gn_val <- calc_g_simple_transient(alfan_vector, z)
        gf_val <- calc_g_simple_transient(alfaf_vector, z)
        gb_val <- calc_g_branches_transient(alfaw_vector, decomp_time, z, ib, params$tmaxb)
        gs_val <- calc_g_stems_transient(alfaw_vector, decomp_time, z, is, params$tmaxs)
        gu_val <- calc_g_simple_transient(alfan_vector, z)
        
        # Store remaining carbon
        matrix_gn[sim_year, input_year] <- Ln[input_year] * gn_val
        matrix_gf[sim_year, input_year] <- Lf[input_year] * gf_val
        matrix_gb[sim_year, input_year] <- Lb[input_year] * gb_val
        matrix_gs[sim_year, input_year] <- Ls[input_year] * gs_val
        matrix_gu[sim_year, input_year] <- Lu[input_year] * gu_val
      }
    }
    
    # -------------------------------------------------------------------------
    # Calculate bulk pool trajectory
    # -------------------------------------------------------------------------
    # Initial pool decomposes as single-quality material
    bulk_pool_trajectory <- numeric(total_years)
    for(year in 1:total_years) {
      t <- year - 1
      bulk_fraction <- (1 + alpha0 * t)^(1 - z)  # Swedish decomposition formula
      bulk_pool_trajectory[year] <- SOC_ss_total * bulk_fraction
    }
    
    if(verbose) cat(sprintf("  Bulk at spinup end: %.1f gC/m²\n", 
                           bulk_pool_trajectory[n_spinup]))
    
    # -------------------------------------------------------------------------
    # Extract simulation period results
    # -------------------------------------------------------------------------
    sim_start_idx <- n_spinup + 1
    sim_end_idx <- total_years
    
    # Sum across all input years (columns) for each simulation year (rows)
    sp_needles_new <- rowSums(matrix_gn[sim_start_idx:sim_end_idx, , drop = FALSE])
    sp_fine_roots_new <- rowSums(matrix_gf[sim_start_idx:sim_end_idx, , drop = FALSE])
    sp_branches_new <- rowSums(matrix_gb[sim_start_idx:sim_end_idx, , drop = FALSE])
    sp_stems_new <- rowSums(matrix_gs[sim_start_idx:sim_end_idx, , drop = FALSE])
    sp_understorey_new <- rowSums(matrix_gu[sim_start_idx:sim_end_idx, , drop = FALSE])
    
    bulk_pool_sim <- bulk_pool_trajectory[sim_start_idx:sim_end_idx]
    
    # -------------------------------------------------------------------------
    # Distribute bulk pool proportionally across tissues
    # -------------------------------------------------------------------------
    # Bulk represents old carbon not tracked explicitly
    bulk_prop_n <- SOC_ss_needles / SOC_ss_total
    bulk_prop_f <- SOC_ss_fine_roots / SOC_ss_total
    bulk_prop_b <- SOC_ss_branches / SOC_ss_total
    bulk_prop_s <- SOC_ss_stems / SOC_ss_total
    bulk_prop_u <- SOC_ss_understorey / SOC_ss_total
    
    # Total SOC = explicitly tracked cohorts + bulk pool
    sp_needles <- sp_needles_new + bulk_pool_sim * bulk_prop_n
    sp_fine_roots <- sp_fine_roots_new + bulk_pool_sim * bulk_prop_f
    sp_branches <- sp_branches_new + bulk_pool_sim * bulk_prop_b
    sp_stems <- sp_stems_new + bulk_pool_sim * bulk_prop_s
    sp_understorey <- sp_understorey_new + bulk_pool_sim * bulk_prop_u
    
    sp_total <- sp_needles + sp_fine_roots + sp_branches + sp_stems + sp_understorey
    
    if(verbose) cat(sprintf("  Start SOC: %.1f gC/m² (%.0f%% bulk)\n", 
                           sp_total[1], 100 * bulk_pool_sim[1] / sp_total[1]))
    
    # -------------------------------------------------------------------------
    # Accumulate across species
    # -------------------------------------------------------------------------
    total_soc <- total_soc + sp_total
    pools_mat[, "needles"] <- pools_mat[, "needles"] + sp_needles
    pools_mat[, "fine_roots"] <- pools_mat[, "fine_roots"] + sp_fine_roots
    pools_mat[, "branches"] <- pools_mat[, "branches"] + sp_branches
    pools_mat[, "stems"] <- pools_mat[, "stems"] + sp_stems
    pools_mat[, "understorey"] <- pools_mat[, "understorey"] + sp_understorey
    initial_pool_soc <- initial_pool_soc + bulk_pool_sim
  }
  
  # ---------------------------------------------------------------------------
  # Calculate respiration (mass balance)
  # ---------------------------------------------------------------------------
  C_input_total <- C_needles + C_branches + C_stems + C_fine_roots + C_understorey
  respiration <- numeric(n_years_sim)
  respiration[1] <- C_input_total[1]  # Year 1: all input respired (no previous SOC)
  for(t in 2:n_years_sim) {
    # Mass balance: respiration = input + SOC loss
    respiration[t] <- C_input_total[t] + total_soc[t-1] - total_soc[t]
  }
  
  if(verbose) cat(sprintf("\nFinal SOC: %.1f gC/m²\n\n", total_soc[1]))
  
  # ---------------------------------------------------------------------------
  # Return results
  # ---------------------------------------------------------------------------
  list(
    pools = pools_mat,               # Pool sizes by tissue type
    total_soc = total_soc,           # Total SOC trajectory
    respiration = respiration,       # Heterotrophic respiration
    initial_pool_soc = initial_pool_soc,  # Bulk pool contribution
    pool_names = c("needles", "fine_roots", "branches", "stems", "understorey"),
    species_props = species_props    # Species composition used
  )
}
