# =============================================================================
# Q-MODEL: Continuous Quality Theory of Decomposition (Matrix Architecture)
# =============================================================================
#
# IMPLEMENTATION APPROACH:
# This implementation uses a "cohort tracking" matrix approach:
# - Each year's litter input is tracked as a separate cohort
# - Each cohort decomposes according to its age (time since input)
# - Total SOC = sum of remaining carbon from all cohorts
#
# This version assumes inputs (C_needles, C_fine_roots, etc.) are provided
# as separate pools. Total C input = sum of all pools.
# =============================================================================

#TODO: check the following. It seems that there is a conceptual difficulty with transient temperatures and cohort tracking
# IMPORTANT: Use decomposition parameters from the INPUT year.
# This assumes the cohort's decomposition trajectory is set by
# conditions when it entered, not current conditions.
# (Alternative: use sim_year parameters for current conditions)
#TODO: update documentation for the temperature implementation in Q

# =============================================================================
# SPECIES-SPECIFIC PARAMETERS
# =============================================================================
# Different tree species have different initial litter quality (q0).
# Higher q0 = higher quality = faster initial decomposition.
#
# q0n: initial quality for needles/foliage
# q0w: initial quality for woody tissue (branches, stems)
# q0f: initial quality for fine roots

Q_SPECIES_PARAMS <- data.frame(
  species = c("pine", "spruce", "birch"),
  q0n = c(1.10, 1.01, 1.01),    # Pine needles decompose faster than spruce
  q0w = c(1.06, 1.00, 1.00),    # Pine wood slightly faster than spruce
  q0f = c(1.036, 1.036, 1.036), # Fine roots similar across species
  stringsAsFactors = FALSE
)

# =============================================================================
# COMMON MODEL PARAMETERS (not species-specific)
# =============================================================================

Q_MODEL_DEFAULT_PARAMS <- list(
  # Microbial efficiency parameters
  eta11 = 0.36,   # Carbon use efficiency parameter
  beta = 7,       # Quality-rate relationship exponent
  e0 = 0.25,      # Initial microbial growth efficiency
  
  # Temperature response
  u00 = 0.0855,   # Base microbial uptake rate at 0°C
  u01 = 0.0157,   # Temperature sensitivity (rate increase per °C)
  
  
  # Woody litter fragmentation times (years)
  # Branches and stems don't enter soil immediately - they are colonized by fungi gradually
  tmaxb = 13,     # Time for complete branch colonization
  tmaxs = 60,     # Time for complete stem colonization
  
  # Other parameters
  fC = 0.5,       # Carbon fraction of organic matter
  spinup_years = 1500  # Years to run for equilibrium initialization (default)
)

# =============================================================================
# DECOMPOSITION FUNCTIONS
# =============================================================================
# These calculate the fraction of carbon remaining after time t.

#' Calculate the z parameter (quality decline exponent)
#' 
#' z controls how quickly quality declines relative to mass loss.
#' Higher z = quality declines faster = decomposition slows more over time.
#'
#' @param e0 Initial efficiency
#' @param beta Quality-rate exponent  
#' @param eta11 Carbon use efficiency
#' @return z parameter (dimensionless)
calc_z <- function(e0, beta, eta11) {
  (1 - e0) / (beta * eta11 * e0)
}

#' Calculate alpha (decomposition rate constant)
#'
#' alpha combines all factors affecting decomposition rate:
#' - fC: carbon fraction
#' - beta, eta11: microbial parameters
#' - u0: temperature-dependent uptake rate
#' - q0: initial quality (species/tissue specific)
#'
#' @return alpha (yr^-1)
calc_alpha <- function(fC, beta, eta11, u0, q0) {
  fC * beta * eta11 * u0 * q0^beta
}

#' Power-law decay for needles, fine roots, understorey
#'
#' These tissues enter the soil immediately and decompose following
#' the basic Q-model equation: C(t)/C(0) = 1/(1 + alpha*t)^z
#'
#' @param alpha Decomposition rate constant
#' @param t Time since input (years)
#' @param z Quality decline exponent
#' @return Fraction of carbon remaining (0 to 1)
calc_g_simple <- function(alpha, t, z) {
  
  # Edge cases
  if(t <= 0 || alpha <= 0) return(1)
  
  # Power-law decay: slower than exponential at long times
  # At t=0: returns 1 (all carbon present)
  # As t→∞: approaches 0 (all carbon decomposed)
  1 / (1 + alpha * t)^z
}

#' Branch decomposition with gradual colonization
#'
#' Branches don't enter the soil all at once - they are gradually colonized
#' by fungi over tmaxb years. 
#'
#' The math accounts for:
#' - Material entering soil linearly over tmaxb years
#' - Each fragment decomposing according to its age in soil
#'
#' @param alpha Decomposition rate constant
#' @param t Time since branch death (years)
#' @param z Quality decline exponent
#' @param ib Effective time in soil = min(t, tmaxb)
#' @param tmaxb Maximum colonization time (13 years typically)
#' @return Fraction of carbon remaining
calc_g_branches <- function(alpha, t, z, ib, tmaxb) {
  if(t <= 0 || alpha <= 0) return(1)
  
  # This is the analytical solution for integrating decomposition
  # of continuously fragmenting material 
  term1 <- 2*((1 + alpha*t)^(1-z) - (1 + alpha*(t-ib))^(1-z)*(1-ib/tmaxb)) / 
    (alpha*(1-z)*tmaxb)
  
  term2 <- 2*((1 + alpha*(t-ib))^(2-z) - (1 + alpha*t)^(2-z)) / 
    (alpha^2*(1-z)*(2-z)*tmaxb^2)
  
  term3 <- (1 - ib/tmaxb)^2
  
  # Ensure non-negative (numerical precision issues)
  max(0, term1 + term2 + term3)
}

#' Stem decomposition with gradual fragmentation
#'
#' Same logic as branches but with longer colonization time (tmaxs = 60 years).
#' Large stems take decades to fully fragment and enter the soil.
#'
#' @param alpha Decomposition rate constant
#' @param t Time since stem death (years)
#' @param z Quality decline exponent
#' @param is Effective time in soil = min(t, tmaxs)
#' @param tmaxs Maximum fragmentation time (60 years typically)
#' @return Fraction of carbon remaining
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
# MAIN Q-MODEL FUNCTION
# =============================================================================
#'
#' Run Q-model simulation with matrix-based cohort tracking
#'
#' ALGORITHM OVERVIEW:
#' 1. Set up spinup period (1500 years) with constant average inputs
#' 2. Append actual simulation years with variable inputs
#' 3. For each input year, track that cohort through all future years
#' 4. Sum all cohorts to get total SOC for each year
#' 5. Return only the simulation period results
#'
#' @param params List of model parameters (default: Q_MODEL_DEFAULT_PARAMS)
#' @param C0 Initial pools - NOT USED, model does its own spinup
#' @param input_df Data frame with annual inputs from map_input_to_q_model()
#'        Required columns: C_needles, C_branches, C_stems, C_fine_roots, 
#'                          C_understorey, temp_mean
#' @param site_row Optional row with species proportions (pine_prop, spruce_prop, birch_prop)
#' @param spinup_years Override default spinup duration
#'
#' @return List with:
#'   - pools: matrix of pool values by year
#'   - total_soc: total SOC by year
#'   - respiration: CO2 release by year
#'   - pool_names: names of pools
#'   - species_props: species proportions used

q_model_run <- function(params = Q_MODEL_DEFAULT_PARAMS, C0 = NULL, input_df, 
                        site_row = NULL, spinup_years = NULL) {
  
  # =========================================================================
  # SETUP: Determine simulation parameters
  # =========================================================================
  
  n_years_sim <- nrow(input_df)  # Number of years to simulate
  n_spinup <- if(is.null(spinup_years)) params$spinup_years else spinup_years #override spinup years if specified, otherwise default from parameters
  
  
  # =========================================================================
  # SPECIES PROPORTIONS
  # =========================================================================
  # The model can handle mixed-species stands by running each species
  # separately and weighting by proportion.
  #
  # If no site_row provided, assumes pure spruce stand.
  
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
  # Input comes from map_input_to_q_model() which aggregates monthly data
  # to annual and sums AWEN fractions to total C per tissue type.
  #
  # These are INDEPENDENT pools, not derived from each other.
  # Total C input = C_needles + C_branches + C_stems + C_fine_roots + C_understorey

  
  # data integrity check
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
  # end of data integrity check
  
  C_needles <- as.numeric(input_df$C_needles)
  C_branches <- as.numeric(input_df$C_branches)
  
  # in some cases no stems inputs are present (pretty common case)
  C_stems <- if ("C_stems" %in% names(input_df)) {
    as.numeric(input_df$C_stems)
  } else {
    warning("Column 'C_stems' not found in input_df; setting now all stem inputs to zero")
    rep(0, n_years_sim)
  }

  C_fine_roots <- as.numeric(input_df$C_fine_roots)
  C_understorey <- as.numeric(input_df$C_understorey)
  temp_mean_vec <- as.numeric(input_df$temp_mean)
  
  # =========================================================================
  # INITIALIZE OUTPUT STORAGE
  # =========================================================================
  # Create the objects to accumulate results across species
  
  total_soc <- numeric(n_years_sim)
  pools_mat <- matrix(0, n_years_sim, 5, 
                      dimnames = list(NULL, c("needles","fine_roots","branches","stems","understorey")))
  
  # =========================================================================
  # LOOP OVER SPECIES
  # =========================================================================
  # Each species has different decomposition rates (q0 parameters).
  # We run the model separately for each species and weight by proportion.
  
  for(sp in names(species_props)) {
    prop <- as.numeric(species_props[sp])
    if(prop <= 0) next  # Skip species not present, restart loop from next species
    
    # Get species-specific parameters
    sp_params <- Q_SPECIES_PARAMS[Q_SPECIES_PARAMS$species == sp, ]
    if(nrow(sp_params) == 0) {
      warning(sprintf("Unknown species: %s, skipping", sp))
      next
    }
    
    # =======================================================================
    # SPINUP: Calculate average inputs for equilibrium initialization
    # =======================================================================
    # Use average of first 5 years of simulation data as constant spinup input.

    start_year_sim <- 1
    spinup_avg_years <- start_year_sim:min(start_year_sim + 4, n_years_sim)
    
    avg_needles <- mean(C_needles[spinup_avg_years])
    avg_fine_roots <- mean(C_fine_roots[spinup_avg_years])
    avg_branches <- mean(C_branches[spinup_avg_years])
    avg_stems <- mean(C_stems[spinup_avg_years])
    avg_understorey <- mean(C_understorey[spinup_avg_years])
    avg_temp <- mean(temp_mean_vec)  # Use mean of ALL years for spinup climate
    
    # =======================================================================
    # CREATE COMBINED TIMELINE (SPINUP + SIMULATION)
    # =======================================================================
    # The model runs continuously from year (-n_spinup+1) to year (n_years_sim).
    # Spinup years have constant inputs; simulation years have actual inputs.
    
    total_years <- n_spinup + n_years_sim
    
    # Create input vectors for entire timeline
    # [spinup_years with constant avg] + [simulation_years with actual data]
    #
    # Each input is weighted by species proportion (prop).
    # This handles mixed stands correctly.
    
    Ln <- c(rep(avg_needles * prop, n_spinup), C_needles * prop)
    Lf <- c(rep(avg_fine_roots * prop, n_spinup), C_fine_roots * prop)
    Lb <- c(rep(avg_branches * prop, n_spinup), C_branches * prop)
    Ls <- c(rep(avg_stems * prop, n_spinup), C_stems * prop)
    Lu <- c(rep(avg_understorey * prop, n_spinup), C_understorey * prop)
    mat <- c(rep(avg_temp, n_spinup), temp_mean_vec)  # Mean annual temperature
    
    # =======================================================================
    # CALCULATE DECOMPOSITION PARAMETERS FOR EACH YEAR
    # =======================================================================
    # Temperature affects decomposition rate through u0.
    # Higher temperature = higher u0 = faster decomposition.
    
    u0 <- params$u00 + params$u01 * mat  # Temperature-dependent uptake rate
    z <- calc_z(params$e0, params$beta, params$eta11)  # Quality decline exponent (constant)
    
    # Get species-specific initial quality parameters
    q0n <- as.numeric(sp_params$q0n)  # Needles quality
    q0f <- as.numeric(sp_params$q0f)  # Fine roots quality
    q0w <- as.numeric(sp_params$q0w)  # Wood quality
    
    # Calculate alpha (decomposition rate) for each tissue type and year
    # Alpha varies by year because u0 is already a vector that depends on temperature
    alfan <- params$fC * params$beta * params$eta11 * u0 * q0n^params$beta  # Needles
    alfaf <- params$fC * params$beta * params$eta11 * u0 * q0f^params$beta  # Fine roots
    alfaw <- params$fC * params$beta * params$eta11 * u0 * q0w^params$beta  # Wood (branches & stems)
    
    # =======================================================================
    # MATRIX-BASED COHORT TRACKING
    # =======================================================================
    # This is the core of the algorithm.
    #
    # We create a matrix where:
    #   - Rows = simulation years (when we measure SOC)
    #   - Columns = input years (when litter entered)
    #   - Value = remaining carbon from that cohort
    #
    # For example, matrix_gn[50, 10] = carbon remaining in year 50 from
    # needles that were input in year 10.
    #
    # Total SOC in year Y = sum of column Y across all input years that
    # have occurred by year Y.
    
    # Initialize tracking matrices for each tissue type
    matrix_gn <- matrix(0, nrow = total_years, ncol = total_years)  # Needles
    matrix_gf <- matrix(0, nrow = total_years, ncol = total_years)  # Fine roots
    matrix_gb <- matrix(0, nrow = total_years, ncol = total_years)  # Branches
    matrix_gs <- matrix(0, nrow = total_years, ncol = total_years)  # Stems
    matrix_gu <- matrix(0, nrow = total_years, ncol = total_years)  # Understorey
    
    # Fill matrices: track each cohort through time
    for(input_year in 1:total_years) {
      for(sim_year in input_year:total_years) {
        
        # How many years has this cohort been decomposing?
        decomp_time <- sim_year - input_year + 1
        
        # For branches and stems, calculate effective time in soil
        # (accounts for gradual colonization)
        ib <- min(decomp_time, params$tmaxb)  # Branches: max 13 years
        is <- min(decomp_time, params$tmaxs)  # Stems: max 60 years
        
        # IMPORTANT: Use decomposition parameters from the INPUT year.
        # This assumes the cohort's decomposition trajectory is set by
        # conditions when it entered, not current conditions.
        # (Alternative: use sim_year parameters for current conditions)
        
        # Calculate fraction remaining for each tissue type
        gn_val <- calc_g_simple(alfan[input_year], decomp_time, z)
        gf_val <- calc_g_simple(alfaf[input_year], decomp_time, z)
        gb_val <- calc_g_branches(alfaw[input_year], decomp_time, z, ib, params$tmaxb)
        gs_val <- calc_g_stems(alfaw[input_year], decomp_time, z, is, params$tmaxs)
        gu_val <- calc_g_simple(alfan[input_year], decomp_time, z)  # Same as needles
        
        # Store remaining carbon = input * fraction remaining
        matrix_gn[sim_year, input_year] <- Ln[input_year] * gn_val
        matrix_gf[sim_year, input_year] <- Lf[input_year] * gf_val
        matrix_gb[sim_year, input_year] <- Lb[input_year] * gb_val
        matrix_gs[sim_year, input_year] <- Ls[input_year] * gs_val
        matrix_gu[sim_year, input_year] <- Lu[input_year] * gu_val
      }
    }
    
    # =======================================================================
    # EXTRACT SIMULATION PERIOD RESULTS
    # =======================================================================
    # We only want results for the actual simulation years, not spinup.
    # Spinup was just to establish realistic initial SOC stocks.
    
    sim_start_idx <- n_spinup + 1  # First simulation year in the matrix
    sim_end_idx <- total_years     # Last simulation year
    
    # Sum across all cohorts (columns) for each simulation year (row)
    # rowSums gives total SOC from all historical inputs
    sp_needles <- rowSums(matrix_gn[sim_start_idx:sim_end_idx, , drop = FALSE])
    sp_fine_roots <- rowSums(matrix_gf[sim_start_idx:sim_end_idx, , drop = FALSE])
    sp_branches <- rowSums(matrix_gb[sim_start_idx:sim_end_idx, , drop = FALSE])
    sp_stems <- rowSums(matrix_gs[sim_start_idx:sim_end_idx, , drop = FALSE])
    sp_understorey <- rowSums(matrix_gu[sim_start_idx:sim_end_idx, , drop = FALSE])
    
    sp_total <- sp_needles + sp_fine_roots + sp_branches + sp_stems + sp_understorey
    
    # =======================================================================
    # ACCUMULATE ACROSS SPECIES
    # =======================================================================
    # Add this species' contribution to the totals.
    # (Already weighted by prop in the Ln, Lf, etc. vectors)
    
    total_soc <- total_soc + sp_total
    pools_mat[, "needles"] <- pools_mat[, "needles"] + sp_needles
    pools_mat[, "fine_roots"] <- pools_mat[, "fine_roots"] + sp_fine_roots
    pools_mat[, "branches"] <- pools_mat[, "branches"] + sp_branches
    pools_mat[, "stems"] <- pools_mat[, "stems"] + sp_stems
    pools_mat[, "understorey"] <- pools_mat[, "understorey"] + sp_understorey
  } # end of the species loop
  
  # =========================================================================
  # CALCULATE RESPIRATION (CO2 RELEASE)
  # =========================================================================
  # Respiration = Carbon in - Carbon stored = Input - ΔStock
  # For mass balance: Respiration[t] = Input[t] + Stock[t-1] - Stock[t]
  
  C_input_total <- C_needles + C_branches + C_stems + C_fine_roots + C_understorey
  
  respiration <- numeric(n_years_sim)
  respiration[1] <- C_input_total[1]  # First year: all input assumed respired (approx)
  for(t in 2:n_years_sim) {
    respiration[t] <- C_input_total[t] + total_soc[t-1] - total_soc[t]
  }
  
  # =========================================================================
  # RETURN RESULTS
  # =========================================================================
  
  list(
    pools = pools_mat,           # Matrix: year x pool
    total_soc = total_soc,       # Vector: total SOC by year
    respiration = respiration,   # Vector: CO2 release by year
    pool_names = c("needles", "fine_roots", "branches", "stems", "understorey"),
    species_props = species_props
  )
}