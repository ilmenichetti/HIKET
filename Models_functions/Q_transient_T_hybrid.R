# =============================================================================
# Q-MODEL WITH HYBRID INITIALIZATION AND TRANSIENT TEMPERATURE
# =============================================================================
#
# SPEEDUP: Uses analytical steady-state initialization followed by short spinup
# instead of full 1500-year spinup.
#
# Expected speedup: 10-20x faster!
# =============================================================================

# Load same parameters and functions as before
source("./Models_functions/Q_transient_T.R")

# =============================================================================
# ANALYTICAL STEADY-STATE CALCULATION
# =============================================================================

#' Calculate analytical steady-state SOC for Q-model
#'
#' For constant inputs L and temperature (constant alpha), the steady-state
#' is found by integrating the age distribution:
#'
#' SOC_ss = integral from 0 to infinity of L * g(t) dt
#'        = integral from 0 to infinity of L / (1 + alpha*t)^z dt
#'
#' This integral has an analytical solution!
#'
#' @param L_input Annual litter input (gC/m²/yr)
#' @param alpha Decomposition rate constant
#' @param z Quality decline exponent
#' @return Steady-state SOC (gC/m²)
calc_steady_state_analytical <- function(L_input, alpha, z) {
  
  if(alpha <= 0 || L_input <= 0) return(0)
  
  # For z != 1, the integral is:
  # SOC = L * (1/alpha) * (1/(z-1)) * [1 + alpha*t]^(1-z) evaluated from 0 to infinity
  #
  # At t=0: (1 + 0)^(1-z) = 1
  # At t→∞: (1 + alpha*t)^(1-z) → 0 (since z > 1 typically)
  #
  # So: SOC = L / (alpha * (z - 1))
  
  if(abs(z - 1) < 0.001) {
    # Special case z ≈ 1: use limiting form
    # SOC ≈ L / alpha * log(some large number)
    # This is approximate, but z is rarely exactly 1
    warning("z is very close to 1, steady-state calculation may be inaccurate")
    return(L_input / alpha * 10)  # Rough approximation
  }
  
  # General case
  SOC_ss <- L_input / (alpha * (z - 1))
  
  return(SOC_ss)
}

#' Calculate initial cohort distribution for analytical steady-state
#'
#' At steady-state with constant inputs, each cohort of age t contains:
#' C(t) = L / (1 + alpha*t)^z
#'
#' This creates the age distribution of SOC.
#'
#' @param L_input Annual litter input
#' @param alpha Decomposition rate constant  
#' @param z Quality decline exponent
#' @param max_age Maximum cohort age to track
#' @return Vector of cohort C amounts by age
calc_initial_cohort_distribution <- function(L_input, alpha, z, max_age = 2000) {
  
  if(alpha <= 0 || L_input <= 0) return(rep(0, max_age))
  
  ages <- 1:max_age
  
  # Each cohort: C(age) = L / (1 + alpha*age)^z
  C_by_age <- L_input / (1 + alpha * ages)^z
  
  return(C_by_age)
}

# =============================================================================
# MAIN Q-MODEL WITH HYBRID INITIALIZATION
# =============================================================================

q_model_run_hybrid <- function(params = Q_MODEL_DEFAULT_PARAMS, 
                               C0 = NULL, 
                               input_df, 
                               site_row = NULL, 
                               spinup_years = 200,  # Much shorter!
                               analytical_age = 2000) {  # Age of "initial" cohorts
  
  n_years_sim <- nrow(input_df)
  n_spinup <- spinup_years  # Use provided short spinup
  
  cat(sprintf("Running Q-model with hybrid initialization:\n"))
  cat(sprintf("  Analytical initialization: %d years back\n", analytical_age))
  cat(sprintf("  Active spinup: %d years\n", n_spinup))
  cat(sprintf("  Simulation: %d years\n", n_years_sim))
  
  # =========================================================================
  # SPECIES PROPORTIONS
  # =========================================================================
  
  if(is.null(site_row)) {
    species_props <- c(pine = 0, spruce = 1, birch = 0)
  } else {
    species_props <- c(
      pine = as.numeric(site_row$pine_prop),
      spruce = as.numeric(site_row$spruce_prop),
      birch = as.numeric(site_row$birch_prop)
    )
  }
  
  # =========================================================================
  # EXTRACT INPUT DATA
  # =========================================================================
  
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
  
  # =========================================================================
  # INITIALIZE OUTPUT STORAGE
  # =========================================================================
  
  total_soc <- numeric(n_years_sim)
  pools_mat <- matrix(0, n_years_sim, 5, 
                      dimnames = list(NULL, c("needles","fine_roots","branches","stems","understorey")))
  
  # =========================================================================
  # LOOP OVER SPECIES
  # =========================================================================
  
  for(sp in names(species_props)) {
    prop <- as.numeric(species_props[sp])
    if(prop <= 0) next
    
    sp_params <- Q_SPECIES_PARAMS[Q_SPECIES_PARAMS$species == sp, ]
    if(nrow(sp_params) == 0) {
      warning(sprintf("Unknown species: %s, skipping", sp))
      next
    }
    
    # =======================================================================
    # CALCULATE AVERAGE CONDITIONS FOR ANALYTICAL INITIALIZATION
    # =======================================================================
    
    start_year_sim <- 1
    spinup_avg_years <- start_year_sim:min(start_year_sim + 4, n_years_sim)
    
    avg_needles <- mean(C_needles[spinup_avg_years])
    avg_fine_roots <- mean(C_fine_roots[spinup_avg_years])
    avg_branches <- mean(C_branches[spinup_avg_years])
    avg_stems <- mean(C_stems[spinup_avg_years])
    avg_understorey <- mean(C_understorey[spinup_avg_years])
    avg_temp <- mean(temp_mean_vec)
    
    # =======================================================================
    # CALCULATE STEADY-STATE PARAMETERS
    # =======================================================================
    
    # Calculate u0 and alpha at average temperature
    u0_avg <- params$u00 + params$u01 * avg_temp
    z <- calc_z(params$e0, params$beta, params$eta11)
    
    # Species-specific parameters
    q0n <- as.numeric(sp_params$q0n)
    q0f <- as.numeric(sp_params$q0f)
    q0w <- as.numeric(sp_params$q0w)
    
    # Calculate alpha for each tissue type at average conditions
    alfan_avg <- params$fC * params$beta * params$eta11 * u0_avg * q0n^params$beta
    alfaf_avg <- params$fC * params$beta * params$eta11 * u0_avg * q0f^params$beta
    alfaw_avg <- params$fC * params$beta * params$eta11 * u0_avg * q0w^params$beta
    
    # =======================================================================
    # ANALYTICAL STEADY-STATE INITIALIZATION
    # =======================================================================
    # Instead of running 1500 years, we calculate what SOC would be
    # at steady-state and distribute it across cohorts
    
    cat(sprintf("\nInitializing %s (%.1f%% of stand):\n", sp, prop * 100))
    
    # Calculate steady-state SOC for each tissue type
    SOC_ss_needles <- calc_steady_state_analytical(avg_needles * prop, alfan_avg, z)
    SOC_ss_fine_roots <- calc_steady_state_analytical(avg_fine_roots * prop, alfaf_avg, z)
    SOC_ss_branches <- calc_steady_state_analytical(avg_branches * prop, alfaw_avg, z)
    SOC_ss_stems <- calc_steady_state_analytical(avg_stems * prop, alfaw_avg, z)
    SOC_ss_understorey <- calc_steady_state_analytical(avg_understorey * prop, alfan_avg, z)
    
    SOC_ss_total <- SOC_ss_needles + SOC_ss_fine_roots + SOC_ss_branches + 
      SOC_ss_stems + SOC_ss_understorey
    
    cat(sprintf("  Analytical steady-state SOC: %.1f gC/m²\n", SOC_ss_total))
    cat(sprintf("    Needles: %.1f, Fine roots: %.1f, Branches: %.1f\n", 
                SOC_ss_needles, SOC_ss_fine_roots, SOC_ss_branches))
    cat(sprintf("    Stems: %.1f, Understorey: %.1f\n", 
                SOC_ss_stems, SOC_ss_understorey))
    
    # =======================================================================
    # CREATE TIMELINE WITH SHORT ACTIVE SPINUP
    # =======================================================================
    # Timeline: [analytical_age years back] + [n_spinup active] + [n_sim]
    # But we only actively track the last (n_spinup + n_sim) years
    
    total_years <- n_spinup + n_years_sim
    
    # Create input vectors for active period
    Ln <- c(rep(avg_needles * prop, n_spinup), C_needles * prop)
    Lf <- c(rep(avg_fine_roots * prop, n_spinup), C_fine_roots * prop)
    Lb <- c(rep(avg_branches * prop, n_spinup), C_branches * prop)
    Ls <- c(rep(avg_stems * prop, n_spinup), C_stems * prop)
    Lu <- c(rep(avg_understorey * prop, n_spinup), C_understorey * prop)
    mat <- c(rep(avg_temp, n_spinup), temp_mean_vec)
    
    # Calculate alpha for each year in active period
    u0 <- params$u00 + params$u01 * mat
    alfan <- params$fC * params$beta * params$eta11 * u0 * q0n^params$beta
    alfaf <- params$fC * params$beta * params$eta11 * u0 * q0f^params$beta
    alfaw <- params$fC * params$beta * params$eta11 * u0 * q0w^params$beta
    
    # =======================================================================
    # INITIALIZE MATRICES WITH ANALYTICAL COHORTS
    # =======================================================================
    
    matrix_gn <- matrix(0, nrow = total_years, ncol = total_years)
    matrix_gf <- matrix(0, nrow = total_years, ncol = total_years)
    matrix_gb <- matrix(0, nrow = total_years, ncol = total_years)
    matrix_gs <- matrix(0, nrow = total_years, ncol = total_years)
    matrix_gu <- matrix(0, nrow = total_years, ncol = total_years)
    
    # Get initial cohort distributions (what would exist from past)
    init_dist_n <- calc_initial_cohort_distribution(avg_needles * prop, alfan_avg, z, analytical_age)
    init_dist_f <- calc_initial_cohort_distribution(avg_fine_roots * prop, alfaf_avg, z, analytical_age)
    init_dist_b <- calc_initial_cohort_distribution(avg_branches * prop, alfaw_avg, z, analytical_age)
    init_dist_s <- calc_initial_cohort_distribution(avg_stems * prop, alfaw_avg, z, analytical_age)
    init_dist_u <- calc_initial_cohort_distribution(avg_understorey * prop, alfan_avg, z, analytical_age)
    
    # Fill first row with initial SOC from "ancient" cohorts
    # These are cohorts that entered before our tracking window
    # We approximate their current amount based on age distribution
    
    # Sum of all cohorts older than our tracking window
    initial_old_SOC_n <- sum(init_dist_n[total_years:analytical_age])
    initial_old_SOC_f <- sum(init_dist_f[total_years:analytical_age])
    initial_old_SOC_b <- sum(init_dist_b[total_years:analytical_age])
    initial_old_SOC_s <- sum(init_dist_s[total_years:analytical_age])
    initial_old_SOC_u <- sum(init_dist_u[total_years:analytical_age])
    
    cat(sprintf("  Initial 'old' SOC (age>%d): %.1f gC/m²\n", 
                total_years, initial_old_SOC_n + initial_old_SOC_f + 
                  initial_old_SOC_b + initial_old_SOC_s + initial_old_SOC_u))
    
    # =======================================================================
    # FILL MATRICES: TRACK COHORTS WITH TRANSIENT TEMPERATURE
    # =======================================================================
    
    for(input_year in 1:total_years) {
      for(sim_year in input_year:total_years) {
        
        decomp_time <- sim_year - input_year + 1
        ib <- min(decomp_time, params$tmaxb)
        is <- min(decomp_time, params$tmaxs)
        
        # Extract alpha vector for transient temperature
        alpha_years <- input_year:sim_year
        alfan_vector <- alfan[alpha_years]
        alfaf_vector <- alfaf[alpha_years]
        alfaw_vector <- alfaw[alpha_years]
        
        # Calculate remaining fraction
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
    
    # =======================================================================
    # EXTRACT SIMULATION PERIOD AND ADD INITIAL OLD SOC
    # =======================================================================
    
    sim_start_idx <- n_spinup + 1
    sim_end_idx <- total_years
    
    sp_needles <- rowSums(matrix_gn[sim_start_idx:sim_end_idx, , drop = FALSE]) + initial_old_SOC_n
    sp_fine_roots <- rowSums(matrix_gf[sim_start_idx:sim_end_idx, , drop = FALSE]) + initial_old_SOC_f
    sp_branches <- rowSums(matrix_gb[sim_start_idx:sim_end_idx, , drop = FALSE]) + initial_old_SOC_b
    sp_stems <- rowSums(matrix_gs[sim_start_idx:sim_end_idx, , drop = FALSE]) + initial_old_SOC_s
    sp_understorey <- rowSums(matrix_gu[sim_start_idx:sim_end_idx, , drop = FALSE]) + initial_old_SOC_u
    
    sp_total <- sp_needles + sp_fine_roots + sp_branches + sp_stems + sp_understorey
    
    cat(sprintf("  Final SOC after spinup: %.1f gC/m²\n", sp_total[1]))
    
    # =======================================================================
    # ACCUMULATE ACROSS SPECIES
    # =======================================================================
    
    total_soc <- total_soc + sp_total
    pools_mat[, "needles"] <- pools_mat[, "needles"] + sp_needles
    pools_mat[, "fine_roots"] <- pools_mat[, "fine_roots"] + sp_fine_roots
    pools_mat[, "branches"] <- pools_mat[, "branches"] + sp_branches
    pools_mat[, "stems"] <- pools_mat[, "stems"] + sp_stems
    pools_mat[, "understorey"] <- pools_mat[, "understorey"] + sp_understorey
  }
  
  # =========================================================================
  # CALCULATE RESPIRATION
  # =========================================================================
  
  C_input_total <- C_needles + C_branches + C_stems + C_fine_roots + C_understorey
  
  respiration <- numeric(n_years_sim)
  respiration[1] <- C_input_total[1]
  for(t in 2:n_years_sim) {
    respiration[t] <- C_input_total[t] + total_soc[t-1] - total_soc[t]
  }
  
  # =========================================================================
  # RETURN RESULTS
  # =========================================================================
  
  cat(sprintf("\nFinal total SOC: %.1f gC/m²\n", total_soc[1]))
  
  list(
    pools = pools_mat,
    total_soc = total_soc,
    respiration = respiration,
    pool_names = c("needles", "fine_roots", "branches", "stems", "understorey"),
    species_props = species_props
  )
}