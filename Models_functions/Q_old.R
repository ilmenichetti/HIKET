# =============================================================================
# Q-MODEL: Continuous Quality Theory of Decomposition
# =============================================================================
#
# Cohort-based model with continuous quality decline (power-law decay).
# Annual timestep, transient spinup for initialization.
# Uses species-specific parameters weighted by composition.
#
# Reference: Bosatta & Ågren (1991, 2003)
# =============================================================================

# Species-specific decomposition parameters
Q_SPECIES_PARAMS <- data.frame(
  species = c("pine", "spruce", "birch"),
  q0n = c(1.10, 1.01, 1.01),    # needles/foliage initial quality
  q0w = c(1.06, 1.00, 1.00),    # woody initial quality
  q0f = c(1.036, 1.036, 1.036), # fine roots initial quality
  stringsAsFactors = FALSE
)

# Common parameters (not species-specific)
Q_MODEL_DEFAULT_PARAMS <- list(
  eta11 = 0.36, beta = 7, e0 = 0.25,
  u00 = 0.0855, u01 = 0.0157,
  tmaxb = 13, tmaxs = 60, fC = 0.5,
  spinup_years = 1500
)

# Power-law decay functions
calc_z <- function(e0, beta, eta11) (1 - e0) / (beta * eta11 * e0)
calc_alpha <- function(fC, beta, eta11, u0, q0) fC * beta * eta11 * u0 * q0^beta

calc_g_simple <- function(alpha, t, z) 1 / (1 + alpha * t)^z

calc_g_branches <- function(alpha, t, z, ib, tmaxb) {
  if(t <= 0 || alpha <= 0) return(1)
  term1 <- 2*((1 + alpha*t)^(1-z) - (1 + alpha*(t-ib))^(1-z)*(1-ib/tmaxb)) / (alpha*(1-z)*tmaxb)
  term2 <- 2*((1 + alpha*(t-ib))^(2-z) - (1 + alpha*t)^(2-z)) / (alpha^2*(1-z)*(2-z)*tmaxb^2)
  term3 <- (1 - ib/tmaxb)^2
  max(0, term1 + term2 + term3)
}

calc_g_stems <- function(alpha, t, z, is, tmaxs) {
  if(t <= 0 || alpha <= 0) return(1)
  term1 <- 2*((1 + alpha*t)^(1-z) - (1 + alpha*(t-is))^(1-z)*(1-is/tmaxs)) / (alpha*(1-z)*tmaxs)
  term2 <- 2*((1 + alpha*(t-is))^(2-z) - (1 + alpha*t)^(2-z)) / (alpha^2*(1-z)*(2-z)*tmaxs^2)
  term3 <- (1 - is/tmaxs)^2
  max(0, term1 + term2 + term3)
}

# Decompose a single cohort for ONE species
decompose_cohort_species <- function(C_input, t, temp, params, sp_params) {
  u0 <- params$u00 + params$u01 * temp
  z <- calc_z(params$e0, params$beta, params$eta11)
  
  # Extract numeric values from sp_params (may be data frame row)
  q0n <- as.numeric(sp_params$q0n)
  q0f <- as.numeric(sp_params$q0f)
  q0w <- as.numeric(sp_params$q0w)
  
  alpha_n <- calc_alpha(params$fC, params$beta, params$eta11, u0, q0n)
  alpha_f <- calc_alpha(params$fC, params$beta, params$eta11, u0, q0f)
  alpha_w <- calc_alpha(params$fC, params$beta, params$eta11, u0, q0w)
  
  ib <- min(t, params$tmaxb); is <- min(t, params$tmaxs)
  
  # Use as.numeric() to prevent name propagation issues
  c(needles = as.numeric(C_input["needles"]) * calc_g_simple(alpha_n, t, z),
    fine_roots = as.numeric(C_input["fine_roots"]) * calc_g_simple(alpha_f, t, z),
    branches = as.numeric(C_input["branches"]) * calc_g_branches(alpha_w, t, z, ib, params$tmaxb),
    stems = as.numeric(C_input["stems"]) * calc_g_stems(alpha_w, t, z, is, params$tmaxs),
    understorey = as.numeric(C_input["understorey"]) * calc_g_simple(alpha_n, t, z))
}

# Annual step for ONE species
q_model_annual_step_species <- function(cohorts, C_input, temp, params, sp_params) {
  new_cohort <- list(C_input = C_input, age = 0, temp_input = temp)
  cohorts <- c(cohorts, list(new_cohort))
  
  total_by_type <- c(needles=0, fine_roots=0, branches=0, stems=0, understorey=0)
  cohorts_keep <- list()
  
  for(cohort in cohorts) {
    cohort$age <- cohort$age + 1
    C_remain <- decompose_cohort_species(cohort$C_input, cohort$age, cohort$temp_input, params, sp_params)
    
    if(sum(C_remain, na.rm=TRUE) > 1e-6) {
      cohort$C_remaining <- C_remain
      cohorts_keep <- c(cohorts_keep, list(cohort))
      total_by_type <- total_by_type + C_remain
    }
  }
  
  list(cohorts = cohorts_keep, C_by_type = total_by_type, total_soc = sum(total_by_type, na.rm=TRUE))
}

# Spinup for ONE species
q_model_spinup_species <- function(C_input_annual, temp_mean, params, sp_params, n_years) {
  cohorts <- list()
  for(y in 1:n_years) {
    result <- q_model_annual_step_species(cohorts, C_input_annual, temp_mean, params, sp_params)
    cohorts <- result$cohorts
  }
  cohorts
}

# Main run function - loops over species weighted by proportions
#' @param params Common model parameters
#' @param C0 Initial cohorts (NULL for spinup)
#' @param input_df Annual input data from map_input_to_q_model()
#' @param site_row Single row from site_df with species proportions
#' @param spinup_years Override spinup duration
q_model_run <- function(params = Q_MODEL_DEFAULT_PARAMS, C0 = NULL, input_df, 
                        site_row = NULL, spinup_years = NULL) {
  
  n_years <- nrow(input_df)
  n_spinup <- if(is.null(spinup_years)) params$spinup_years else spinup_years
  
  # Get species proportions (default: 100% spruce if no site_row)
  if(is.null(site_row)) {
    species_props <- c(pine = 0, spruce = 1, birch = 0)
  } else {
    species_props <- c(
      pine = as.numeric(site_row$pine_prop),
      spruce = as.numeric(site_row$spruce_prop),
      birch = as.numeric(site_row$birch_prop)
    )
  }
  
  # Input columns
  if(!"C_needles" %in% names(input_df)) {
    stop("Q-model requires C_needles, C_branches, C_fine_roots, C_understorey columns")
  }
  
  C_needles <- as.numeric(input_df$C_needles)
  C_branches <- as.numeric(input_df$C_branches)
  C_stems <- if("C_stems" %in% names(input_df)) as.numeric(input_df$C_stems) else rep(0, n_years)
  C_fine_roots <- as.numeric(input_df$C_fine_roots)
  C_understorey <- as.numeric(input_df$C_understorey)
  temp_mean_vec <- as.numeric(input_df$temp_mean)
  
  # Storage for combined results
  total_soc <- respiration <- numeric(n_years)
  pools_mat <- matrix(0, n_years, 5, 
                      dimnames = list(NULL, c("needles","fine_roots","branches","stems","understorey")))
  
  # Run each species separately and combine
  species_results <- list()
  
  for(sp in names(species_props)) {
    prop <- as.numeric(species_props[sp])  # CRITICAL: Convert to numeric to avoid name issues
    if(prop <= 0) next
    
    sp_params <- Q_SPECIES_PARAMS[Q_SPECIES_PARAMS$species == sp, ]
    if(nrow(sp_params) == 0) {
      warning(sprintf("Unknown species: %s, skipping", sp))
      next
    }
    
    # Spinup for this species
    avg_input <- c(needles = mean(C_needles) * prop, 
                   fine_roots = mean(C_fine_roots) * prop,
                   branches = mean(C_branches) * prop, 
                   stems = mean(C_stems) * prop,
                   understorey = mean(C_understorey) * prop)
    
    cohorts <- q_model_spinup_species(avg_input, mean(temp_mean_vec), params, sp_params, n_spinup)
    
    # Run simulation for this species
    sp_soc <- numeric(n_years)
    sp_pools <- matrix(0, n_years, 5)
    
    for(t in 1:n_years) {
      C_input <- c(needles = C_needles[t] * prop, 
                   fine_roots = C_fine_roots[t] * prop,
                   branches = C_branches[t] * prop, 
                   stems = C_stems[t] * prop,
                   understorey = C_understorey[t] * prop)
      
      result <- q_model_annual_step_species(cohorts, C_input, temp_mean_vec[t], params, sp_params)
      cohorts <- result$cohorts
      sp_soc[t] <- result$total_soc
      sp_pools[t,] <- result$C_by_type
    }
    
    species_results[[sp]] <- list(soc = sp_soc, pools = sp_pools)
  }
  
  # Combine species results
  for(sp in names(species_results)) {
    total_soc <- total_soc + species_results[[sp]]$soc
    pools_mat <- pools_mat + species_results[[sp]]$pools
  }
  
  # Calculate respiration
  C_input_total <- C_needles + C_branches + C_stems + C_fine_roots + C_understorey
  respiration[1] <- C_input_total[1]
  for(t in 2:n_years) {
    respiration[t] <- C_input_total[t] + total_soc[t-1] - total_soc[t]
  }
  
  list(pools = pools_mat, total_soc = total_soc, respiration = respiration,
       pool_names = c("needles","fine_roots","branches","stems","understorey"),
       species_props = species_props)
}
