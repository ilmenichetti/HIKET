# =============================================================================
# Q-MODEL: Continuous Quality Theory with TRANSIENT TEMPERATURE
# =============================================================================
#
# MODIFICATION: This version uses transient (current-year) temperature instead
# of input-year temperature for decomposition calculations.
#
# KEY CHANGE: Instead of using analytical formula with constant alpha,
# we track cumulative alpha*t for each cohort as it decomposes.
#
# =============================================================================

# =============================================================================
# SPECIES-SPECIFIC PARAMETERS (unchanged)
# =============================================================================

Q_SPECIES_PARAMS <- data.frame(
  species = c("pine", "spruce", "birch"),
  q0n = c(1.10, 1.01, 1.01),    
  q0w = c(1.06, 1.00, 1.00),    
  q0f = c(1.036, 1.036, 1.036), 
  stringsAsFactors = FALSE
)

# =============================================================================
# COMMON MODEL PARAMETERS (unchanged)
# =============================================================================

Q_MODEL_DEFAULT_PARAMS <- list(
  eta11 = 0.36,   
  beta = 7,       
  e0 = 0.25,      
  u00 = 0.0855,   
  u01 = 0.0157,   
  tmaxb = 13,     
  tmaxs = 60,     
  fC = 0.5,       
  spinup_years = 1500
)

# =============================================================================
# DECOMPOSITION FUNCTIONS (unchanged)
# =============================================================================

calc_z <- function(e0, beta, eta11) {
  (1 - e0) / (beta * eta11 * e0)
}

calc_alpha <- function(fC, beta, eta11, u0, q0) {
  fC * beta * eta11 * u0 * q0^beta
}

calc_g_simple <- function(alpha, t, z) {
  if(t <= 0 || alpha <= 0) return(1)
  1 / (1 + alpha * t)^z
}

calc_g_branches <- function(alpha, t, z, ib, tmaxb) {
  if(t <= 0 || alpha <= 0) return(1)
  
  term1 <- 2*((1 + alpha*t)^(1-z) - (1 + alpha*(t-ib))^(1-z)*(1-ib/tmaxb)) / 
    (alpha*(1-z)*tmaxb)
  
  term2 <- 2*((1 + alpha*(t-ib))^(2-z) - (1 + alpha*t)^(2-z)) / 
    (alpha^2*(1-z)*(2-z)*tmaxb^2)
  
  term3 <- (1 - ib/tmaxb)^2
  
  max(0, term1 + term2 + term3)
}

calc_g_stems <- function(alpha, t, z, is, tmaxs) {
  if(t <= 0 || alpha <= 0) return(1)
  
  term1 <- 2*((1 + alpha*t)^(1-z) - (1 + alpha*(t-is))^(1-z)*(1-is/tmaxs)) / 
    (alpha*(1-z)*tmaxs)
  
  term2 <- 2*((1 + alpha*(t-is))^(2-z) - (1 + alpha*t)^(2-z)) / 
    (alpha^2*(1-z)*(2-z)*tmaxs^2)
  
  term3 <- (1 - is/tmaxs)^2
  
  max(0, term1 + term2 + term3)
}

# =============================================================================
# NEW FUNCTIONS FOR TRANSIENT TEMPERATURE
# =============================================================================

#' Calculate remaining fraction with transient temperature (simple tissues)
#' 
#' For needles, fine roots, understorey that decompose without fragmentation.
#' Uses cumulative alpha*t approach instead of analytical formula.
#'
#' @param alpha_vector Vector of alpha values from input_year to sim_year
#' @param z Quality decline exponent
#' @return Fraction of carbon remaining (0 to 1)
calc_g_simple_transient <- function(alpha_vector, z) {
  # Sum alpha over all years (cumulative decomposition potential)
  cumulative_alpha_t <- sum(alpha_vector)
  
  if(cumulative_alpha_t <= 0) return(1)
  
  # Apply Q-model formula with cumulative alpha*t
  1 / (1 + cumulative_alpha_t)^z
}

#' Calculate remaining fraction for branches with transient temperature
#' 
#' This is more complex because branches fragment gradually.
#' We need to account for both:
#' 1. Gradual colonization (fragmentation over tmaxb years)
#' 2. Transient temperature affecting decomposition
#'
#' APPROXIMATION: We use the average alpha over the decomposition period
#' with the analytical branch formula. This is simpler than full numerical
#' integration and should be reasonably accurate for slowly changing temperatures.
#'
#' @param alpha_vector Vector of alpha values from input_year to sim_year
#' @param decomp_time Total time since input
#' @param z Quality decline exponent
#' @param ib Effective time in soil
#' @param tmaxb Maximum colonization time
#' @return Fraction of carbon remaining
calc_g_branches_transient <- function(alpha_vector, decomp_time, z, ib, tmaxb) {
  if(decomp_time <= 0) return(1)
  
  # Use average alpha (approximation for transient temperature)
  alpha_avg <- mean(alpha_vector)
  
  if(alpha_avg <= 0) return(1)
  
  # Apply standard branch formula with average alpha
  calc_g_branches(alpha_avg, decomp_time, z, ib, tmaxb)
}

#' Calculate remaining fraction for stems with transient temperature
#' 
#' Same logic as branches but with longer fragmentation time.
calc_g_stems_transient <- function(alpha_vector, decomp_time, z, is, tmaxs) {
  if(decomp_time <= 0) return(1)
  
  # Use average alpha (approximation for transient temperature)
  alpha_avg <- mean(alpha_vector)
  
  if(alpha_avg <= 0) return(1)
  
  # Apply standard stem formula with average alpha
  calc_g_stems(alpha_avg, decomp_time, z, is, tmaxs)
}

# =============================================================================
# MAIN Q-MODEL FUNCTION WITH TRANSIENT TEMPERATURE
# =============================================================================

q_model_run_transient <- function(params = Q_MODEL_DEFAULT_PARAMS, C0 = NULL, 
                                  input_df, site_row = NULL, spinup_years = NULL) {
  
  # =========================================================================
  # SETUP (unchanged from original)
  # =========================================================================
  
  n_years_sim <- nrow(input_df)
  n_spinup <- if(is.null(spinup_years)) params$spinup_years else spinup_years
  
  # =========================================================================
  # SPECIES PROPORTIONS (unchanged)
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
  # EXTRACT INPUT DATA (unchanged)
  # =========================================================================
  
  required_cols <- c("C_needles", "C_branches", "C_fine_roots", "C_understorey")
  missing <- setdiff(required_cols, names(input_df))
  
  if (length(missing) > 0) {
    stop(
      sprintf(
        "Q-model requires columns: %s. Missing: %s",
        paste(required_cols, collapse = ", "),
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  C_needles <- as.numeric(input_df$C_needles)
  C_branches <- as.numeric(input_df$C_branches)
  C_stems <- if ("C_stems" %in% names(input_df)) {
    as.numeric(input_df$C_stems)
  } else {
    warning("Column 'C_stems' not found; setting to zero")
    rep(0, n_years_sim)
  }
  C_fine_roots <- as.numeric(input_df$C_fine_roots)
  C_understorey <- as.numeric(input_df$C_understorey)
  temp_mean_vec <- as.numeric(input_df$temp_mean)
  
  # =========================================================================
  # INITIALIZE OUTPUT STORAGE (unchanged)
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
    # SPINUP (unchanged)
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
    # CREATE COMBINED TIMELINE (unchanged)
    # =======================================================================
    
    total_years <- n_spinup + n_years_sim
    
    Ln <- c(rep(avg_needles * prop, n_spinup), C_needles * prop)
    Lf <- c(rep(avg_fine_roots * prop, n_spinup), C_fine_roots * prop)
    Lb <- c(rep(avg_branches * prop, n_spinup), C_branches * prop)
    Ls <- c(rep(avg_stems * prop, n_spinup), C_stems * prop)
    Lu <- c(rep(avg_understorey * prop, n_spinup), C_understorey * prop)
    mat <- c(rep(avg_temp, n_spinup), temp_mean_vec)
    
    # =======================================================================
    # CALCULATE DECOMPOSITION PARAMETERS (unchanged)
    # =======================================================================
    
    u0 <- params$u00 + params$u01 * mat
    z <- calc_z(params$e0, params$beta, params$eta11)
    
    q0n <- as.numeric(sp_params$q0n)
    q0f <- as.numeric(sp_params$q0f)
    q0w <- as.numeric(sp_params$q0w)
    
    alfan <- params$fC * params$beta * params$eta11 * u0 * q0n^params$beta
    alfaf <- params$fC * params$beta * params$eta11 * u0 * q0f^params$beta
    alfaw <- params$fC * params$beta * params$eta11 * u0 * q0w^params$beta
    
    # =======================================================================
    # MATRIX-BASED COHORT TRACKING - MODIFIED FOR TRANSIENT TEMPERATURE
    # =======================================================================
    # KEY CHANGE: Instead of using alpha[input_year] for entire trajectory,
    # we use alpha[input_year:sim_year] and calculate cumulative effect.
    
    matrix_gn <- matrix(0, nrow = total_years, ncol = total_years)
    matrix_gf <- matrix(0, nrow = total_years, ncol = total_years)
    matrix_gb <- matrix(0, nrow = total_years, ncol = total_years)
    matrix_gs <- matrix(0, nrow = total_years, ncol = total_years)
    matrix_gu <- matrix(0, nrow = total_years, ncol = total_years)
    
    # Fill matrices with transient temperature
    for(input_year in 1:total_years) {
      for(sim_year in input_year:total_years) {
        
        decomp_time <- sim_year - input_year + 1
        
        # For branches and stems
        ib <- min(decomp_time, params$tmaxb)
        is <- min(decomp_time, params$tmaxs)
        
        # =================================================================
        # KEY MODIFICATION: Extract alpha vector from input_year to sim_year
        # =================================================================
        # This is the vector of alpha values over the cohort's lifetime
        alpha_years <- input_year:sim_year
        
        alfan_vector <- alfan[alpha_years]
        alfaf_vector <- alfaf[alpha_years]
        alfaw_vector <- alfaw[alpha_years]
        
        # Calculate remaining fraction using TRANSIENT temperature functions
        gn_val <- calc_g_simple_transient(alfan_vector, z)
        gf_val <- calc_g_simple_transient(alfaf_vector, z)
        gb_val <- calc_g_branches_transient(alfaw_vector, decomp_time, z, ib, params$tmaxb)
        gs_val <- calc_g_stems_transient(alfaw_vector, decomp_time, z, is, params$tmaxs)
        gu_val <- calc_g_simple_transient(alfan_vector, z)  # Same as needles
        
        # Store remaining carbon
        matrix_gn[sim_year, input_year] <- Ln[input_year] * gn_val
        matrix_gf[sim_year, input_year] <- Lf[input_year] * gf_val
        matrix_gb[sim_year, input_year] <- Lb[input_year] * gb_val
        matrix_gs[sim_year, input_year] <- Ls[input_year] * gs_val
        matrix_gu[sim_year, input_year] <- Lu[input_year] * gu_val
      }
    }
    
    # =======================================================================
    # EXTRACT SIMULATION PERIOD RESULTS (unchanged)
    # =======================================================================
    
    sim_start_idx <- n_spinup + 1
    sim_end_idx <- total_years
    
    sp_needles <- rowSums(matrix_gn[sim_start_idx:sim_end_idx, , drop = FALSE])
    sp_fine_roots <- rowSums(matrix_gf[sim_start_idx:sim_end_idx, , drop = FALSE])
    sp_branches <- rowSums(matrix_gb[sim_start_idx:sim_end_idx, , drop = FALSE])
    sp_stems <- rowSums(matrix_gs[sim_start_idx:sim_end_idx, , drop = FALSE])
    sp_understorey <- rowSums(matrix_gu[sim_start_idx:sim_end_idx, , drop = FALSE])
    
    sp_total <- sp_needles + sp_fine_roots + sp_branches + sp_stems + sp_understorey
    
    # =======================================================================
    # ACCUMULATE ACROSS SPECIES (unchanged)
    # =======================================================================
    
    total_soc <- total_soc + sp_total
    pools_mat[, "needles"] <- pools_mat[, "needles"] + sp_needles
    pools_mat[, "fine_roots"] <- pools_mat[, "fine_roots"] + sp_fine_roots
    pools_mat[, "branches"] <- pools_mat[, "branches"] + sp_branches
    pools_mat[, "stems"] <- pools_mat[, "stems"] + sp_stems
    pools_mat[, "understorey"] <- pools_mat[, "understorey"] + sp_understorey
  }
  
  # =========================================================================
  # CALCULATE RESPIRATION (unchanged)
  # =========================================================================
  
  C_input_total <- C_needles + C_branches + C_stems + C_fine_roots + C_understorey
  
  respiration <- numeric(n_years_sim)
  respiration[1] <- C_input_total[1]
  for(t in 2:n_years_sim) {
    respiration[t] <- C_input_total[t] + total_soc[t-1] - total_soc[t]
  }
  
  # =========================================================================
  # RETURN RESULTS (unchanged)
  # =========================================================================
  
  list(
    pools = pools_mat,
    total_soc = total_soc,
    respiration = respiration,
    pool_names = c("needles", "fine_roots", "branches", "stems", "understorey"),
    species_props = species_props
  )
}