# =============================================================================
# Q-MODEL: HYBRID INITIALIZATION WITH SWEDISH APPROACH + SPINUP
# =============================================================================
#
# APPROACH (Swedish Heureka + age distribution refinement):
# 
# 1. Calculate analytical steady-state (C_ss) - correct total SOC amount
# 2. Start with C_ss as a BULK POOL that decomposes:
#    C_initial(t) = C_ss * (1 + alpha0 * t)^(1-z)  [Hyvönen et al. 1998, Eq. 7]
# 
# 3. Run SHORT SPINUP (10-20 years) where:
#    - Bulk pool gradually decomposes
#    - NEW cohorts enter each year and get tracked in matrices
#    - Proper age distribution builds up
# 
# 4. Start SIMULATION with:
#    - Much smaller bulk pool (mostly decomposed)
#    - 10-20 years of explicitly tracked cohorts
#    - Total SOC = C_initial(t) + C_new_cohorts(t)
#
# WHY THIS WORKS:
# - Start at RIGHT total SOC (no convergence drift)
# - Spinup transitions from "bulk single-quality" to "detailed multi-cohort"
# - By simulation start, most carbon is in explicit cohorts (not bulk)
# - Short spinup sufficient because we START at equilibrium
#
# SPEED:
# - 10-20 year spinup vs 1500 year spinup = 75-150x faster!
#
# =============================================================================

# =============================================================================
# LOAD BASE Q-MODEL FUNCTIONS
# =============================================================================

source("./Models_functions/Q_transient_T.R")

# =============================================================================
# STEADY-STATE CALCULATION FUNCTIONS (CORRECTED FOR FRAGMENTATION)
# =============================================================================

#' Calculate steady-state SOC for SIMPLE tissues
calc_steady_state_simple <- function(L_input, alpha, z) {
  if(alpha <= 0 || L_input <= 0) return(0)
  if(abs(z - 1) < 0.001) {
    warning("z is very close to 1")
    return(L_input / alpha * 10)
  }
  L_input / (alpha * (z - 1))
}

#' Calculate steady-state SOC for BRANCHES (numerical integration)
calc_steady_state_branches <- function(L_input, alpha, z, tmaxb) {
  if(alpha <= 0 || L_input <= 0) return(0)
  
  max_age <- 500
  ages <- 1:max_age
  
  g_values <- sapply(ages, function(t) {
    ib <- min(t, tmaxb)
    term1 <- 2*((1 + alpha*t)^(1-z) - (1 + alpha*(t-ib))^(1-z)*(1-ib/tmaxb)) / 
      (alpha*(1-z)*tmaxb)
    term2 <- 2*((1 + alpha*(t-ib))^(2-z) - (1 + alpha*t)^(2-z)) / 
      (alpha^2*(1-z)*(2-z)*tmaxb^2)
    term3 <- (1 - ib/tmaxb)^2
    max(0, term1 + term2 + term3)
  })
  
  L_input * sum(g_values)
}

#' Calculate steady-state SOC for STEMS (numerical integration)
calc_steady_state_stems <- function(L_input, alpha, z, tmaxs) {
  if(alpha <= 0 || L_input <= 0) return(0)
  
  max_age <- 1000
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
                               spinup_years = 20) {
  
  n_years_sim <- nrow(input_df)
  n_spinup <- spinup_years
  
  cat(sprintf("Q-model SWEDISH HYBRID initialization:\n"))
  cat(sprintf("  Spinup: %d years\n", n_spinup))
  cat(sprintf("  Simulation: %d years\n\n", n_years_sim))
  
  # Species proportions
  if(is.null(site_row)) {
    species_props <- c(pine = 0, spruce = 1, birch = 0)
  } else {
    species_props <- c(
      pine = as.numeric(site_row$pine_prop),
      spruce = as.numeric(site_row$spruce_prop),
      birch = as.numeric(site_row$birch_prop)
    )
  }
  
  # Extract inputs
  required_cols <- c("C_needles", "C_branches", "C_fine_roots", "C_understorey")
  missing <- setdiff(required_cols, names(input_df))
  if (length(missing) > 0) {
    stop(sprintf("Missing columns: %s", paste(missing, collapse = ", ")))
  }
  
  C_needles <- as.numeric(input_df$C_needles)
  C_branches <- as.numeric(input_df$C_branches)
  C_stems <- if ("C_stems" %in% names(input_df)) as.numeric(input_df$C_stems) else rep(0, n_years_sim)
  C_fine_roots <- as.numeric(input_df$C_fine_roots)
  C_understorey <- as.numeric(input_df$C_understorey)
  temp_mean_vec <- as.numeric(input_df$temp_mean)
  
  # Output storage
  total_soc <- numeric(n_years_sim)
  pools_mat <- matrix(0, n_years_sim, 5, 
                      dimnames = list(NULL, c("needles","fine_roots","branches","stems","understorey")))
  initial_pool_soc <- numeric(n_years_sim)
  
  # Loop over species
  for(sp in names(species_props)) {
    prop <- as.numeric(species_props[sp])
    if(prop <= 0) next
    
    sp_params <- Q_SPECIES_PARAMS[Q_SPECIES_PARAMS$species == sp, ]
    if(nrow(sp_params) == 0) next
    
    cat(sprintf("%s (%.1f%%)\n", sp, prop * 100))
    
    # Calculate averages for steady-state
    spinup_avg_years <- 1:min(5, n_years_sim)
    avg_needles <- mean(C_needles[spinup_avg_years])
    avg_fine_roots <- mean(C_fine_roots[spinup_avg_years])
    avg_branches <- mean(C_branches[spinup_avg_years])
    avg_stems <- mean(C_stems[spinup_avg_years])
    avg_understorey <- mean(C_understorey[spinup_avg_years])
    avg_temp <- mean(temp_mean_vec)
    
    # Calculate parameters
    u0_avg <- params$u00 + params$u01 * avg_temp
    z <- calc_z(params$e0, params$beta, params$eta11)
    q0n <- as.numeric(sp_params$q0n)
    q0f <- as.numeric(sp_params$q0f)
    q0w <- as.numeric(sp_params$q0w)
    
    alfan_avg <- params$fC * params$beta * params$eta11 * u0_avg * q0n^params$beta
    alfaf_avg <- params$fC * params$beta * params$eta11 * u0_avg * q0f^params$beta
    alfaw_avg <- params$fC * params$beta * params$eta11 * u0_avg * q0w^params$beta
    
    # Calculate analytical steady-state
    SOC_ss_needles <- calc_steady_state_simple(avg_needles * prop, alfan_avg, z)
    SOC_ss_fine_roots <- calc_steady_state_simple(avg_fine_roots * prop, alfaf_avg, z)
    SOC_ss_understorey <- calc_steady_state_simple(avg_understorey * prop, alfan_avg, z)
    SOC_ss_branches <- calc_steady_state_branches(avg_branches * prop, alfaw_avg, z, params$tmaxb)
    SOC_ss_stems <- calc_steady_state_stems(avg_stems * prop, alfaw_avg, z, params$tmaxs)
    
    SOC_ss_total <- SOC_ss_needles + SOC_ss_fine_roots + SOC_ss_branches + SOC_ss_stems + SOC_ss_understorey
    
    cat(sprintf("  Steady-state: %.1f gC/m²\n", SOC_ss_total))
    
    # Calculate alpha0 for bulk decomposition (Swedish formula)
    i0 <- (avg_needles + avg_fine_roots + avg_branches + avg_stems + avg_understorey) * prop
    alpha0 <- -i0 / (SOC_ss_total * (1 - z))
    
    # Create timeline
    total_years <- n_spinup + n_years_sim
    Ln <- c(rep(avg_needles * prop, n_spinup), C_needles * prop)
    Lf <- c(rep(avg_fine_roots * prop, n_spinup), C_fine_roots * prop)
    Lb <- c(rep(avg_branches * prop, n_spinup), C_branches * prop)
    Ls <- c(rep(avg_stems * prop, n_spinup), C_stems * prop)
    Lu <- c(rep(avg_understorey * prop, n_spinup), C_understorey * prop)
    mat <- c(rep(avg_temp, n_spinup), temp_mean_vec)
    
    u0 <- params$u00 + params$u01 * mat
    alfan <- params$fC * params$beta * params$eta11 * u0 * q0n^params$beta
    alfaf <- params$fC * params$beta * params$eta11 * u0 * q0f^params$beta
    alfaw <- params$fC * params$beta * params$eta11 * u0 * q0w^params$beta
    
    # Initialize matrices
    matrix_gn <- matrix(0, nrow = total_years, ncol = total_years)
    matrix_gf <- matrix(0, nrow = total_years, ncol = total_years)
    matrix_gb <- matrix(0, nrow = total_years, ncol = total_years)
    matrix_gs <- matrix(0, nrow = total_years, ncol = total_years)
    matrix_gu <- matrix(0, nrow = total_years, ncol = total_years)
    
    # Fill matrices with transient temperature
    for(input_year in 1:total_years) {
      for(sim_year in input_year:total_years) {
        decomp_time <- sim_year - input_year + 1
        ib <- min(decomp_time, params$tmaxb)
        is <- min(decomp_time, params$tmaxs)
        
        alpha_years <- input_year:sim_year
        alfan_vector <- alfan[alpha_years]
        alfaf_vector <- alfaf[alpha_years]
        alfaw_vector <- alfaw[alpha_years]
        
        gn_val <- calc_g_simple_transient(alfan_vector, z)
        gf_val <- calc_g_simple_transient(alfaf_vector, z)
        gb_val <- calc_g_branches_transient(alfaw_vector, decomp_time, z, ib, params$tmaxb)
        gs_val <- calc_g_stems_transient(alfaw_vector, decomp_time, z, is, params$tmaxs)
        gu_val <- calc_g_simple_transient(alfan_vector, z)
        
        matrix_gn[sim_year, input_year] <- Ln[input_year] * gn_val
        matrix_gf[sim_year, input_year] <- Lf[input_year] * gf_val
        matrix_gb[sim_year, input_year] <- Lb[input_year] * gb_val
        matrix_gs[sim_year, input_year] <- Ls[input_year] * gs_val
        matrix_gu[sim_year, input_year] <- Lu[input_year] * gu_val
      }
    }
    
    # Calculate bulk pool decomposition (SWEDISH APPROACH)
    bulk_pool_trajectory <- numeric(total_years)
    for(year in 1:total_years) {
      t <- year - 1
      bulk_fraction <- (1 + alpha0 * t)^(1 - z)
      bulk_pool_trajectory[year] <- SOC_ss_total * bulk_fraction
    }
    
    cat(sprintf("  Bulk at spinup end: %.1f gC/m²\n", bulk_pool_trajectory[n_spinup]))
    
    # Extract simulation period: BULK + NEW COHORTS
    sim_start_idx <- n_spinup + 1
    sim_end_idx <- total_years
    
    sp_needles_new <- rowSums(matrix_gn[sim_start_idx:sim_end_idx, , drop = FALSE])
    sp_fine_roots_new <- rowSums(matrix_gf[sim_start_idx:sim_end_idx, , drop = FALSE])
    sp_branches_new <- rowSums(matrix_gb[sim_start_idx:sim_end_idx, , drop = FALSE])
    sp_stems_new <- rowSums(matrix_gs[sim_start_idx:sim_end_idx, , drop = FALSE])
    sp_understorey_new <- rowSums(matrix_gu[sim_start_idx:sim_end_idx, , drop = FALSE])
    
    bulk_pool_sim <- bulk_pool_trajectory[sim_start_idx:sim_end_idx]
    
    # Distribute bulk proportionally
    bulk_prop_n <- SOC_ss_needles / SOC_ss_total
    bulk_prop_f <- SOC_ss_fine_roots / SOC_ss_total
    bulk_prop_b <- SOC_ss_branches / SOC_ss_total
    bulk_prop_s <- SOC_ss_stems / SOC_ss_total
    bulk_prop_u <- SOC_ss_understorey / SOC_ss_total
    
    sp_needles <- sp_needles_new + bulk_pool_sim * bulk_prop_n
    sp_fine_roots <- sp_fine_roots_new + bulk_pool_sim * bulk_prop_f
    sp_branches <- sp_branches_new + bulk_pool_sim * bulk_prop_b
    sp_stems <- sp_stems_new + bulk_pool_sim * bulk_prop_s
    sp_understorey <- sp_understorey_new + bulk_pool_sim * bulk_prop_u
    
    sp_total <- sp_needles + sp_fine_roots + sp_branches + sp_stems + sp_understorey
    
    cat(sprintf("  Start SOC: %.1f gC/m² (%.0f%% bulk)\n", 
                sp_total[1], 100 * bulk_pool_sim[1] / sp_total[1]))
    
    # Accumulate
    total_soc <- total_soc + sp_total
    pools_mat[, "needles"] <- pools_mat[, "needles"] + sp_needles
    pools_mat[, "fine_roots"] <- pools_mat[, "fine_roots"] + sp_fine_roots
    pools_mat[, "branches"] <- pools_mat[, "branches"] + sp_branches
    pools_mat[, "stems"] <- pools_mat[, "stems"] + sp_stems
    pools_mat[, "understorey"] <- pools_mat[, "understorey"] + sp_understorey
    initial_pool_soc <- initial_pool_soc + bulk_pool_sim
  }
  
  # Calculate respiration
  C_input_total <- C_needles + C_branches + C_stems + C_fine_roots + C_understorey
  respiration <- numeric(n_years_sim)
  respiration[1] <- C_input_total[1]
  for(t in 2:n_years_sim) {
    respiration[t] <- C_input_total[t] + total_soc[t-1] - total_soc[t]
  }
  
  cat(sprintf("\nFinal SOC: %.1f gC/m²\n\n", total_soc[1]))
  
  list(
    pools = pools_mat,
    total_soc = total_soc,
    respiration = respiration,
    initial_pool_soc = initial_pool_soc,
    pool_names = c("needles", "fine_roots", "branches", "stems", "understorey"),
    species_props = species_props
  )
}
