# =============================================================================
# Q-MODEL: Continuous Quality Theory of Decomposition (DEBUG VERSION)
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
  
  cat(sprintf("      [decompose_cohort] age=%d, temp=%.2f, u0=%.4f, z=%.4f\n", t, temp, u0, z))
  cat(sprintf("      [decompose_cohort] q0: needles=%.3f, fine_roots=%.3f, woody=%.3f\n", q0n, q0f, q0w))
  
  cat(sprintf("      [decompose_cohort] C_input names: %s\n", paste(names(C_input), collapse=", ")))
  cat(sprintf("      [decompose_cohort] C_input values: needles=%.6f, fine_roots=%.6f, branches=%.6f, stems=%.6f, understorey=%.6f\n",
              C_input["needles"], C_input["fine_roots"], C_input["branches"], C_input["stems"], C_input["understorey"]))
  
  alpha_n <- calc_alpha(params$fC, params$beta, params$eta11, u0, q0n)
  alpha_f <- calc_alpha(params$fC, params$beta, params$eta11, u0, q0f)
  alpha_w <- calc_alpha(params$fC, params$beta, params$eta11, u0, q0w)
  
  cat(sprintf("      [decompose_cohort] alpha: needles=%.4f, fine_roots=%.4f, woody=%.4f\n", alpha_n, alpha_f, alpha_w))
  
  ib <- min(t, params$tmaxb); is <- min(t, params$tmaxs)
  
  g_n <- calc_g_simple(alpha_n, t, z)
  g_f <- calc_g_simple(alpha_f, t, z)
  g_b <- calc_g_branches(alpha_w, t, z, ib, params$tmaxb)
  g_s <- calc_g_stems(alpha_w, t, z, is, params$tmaxs)
  
  cat(sprintf("      [decompose_cohort] g-values: needles=%.6f, fine_roots=%.6f, branches=%.6f, stems=%.6f\n",
              g_n, g_f, g_b, g_s))
  
  result <- c(needles = as.numeric(C_input["needles"]) * g_n,
              fine_roots = as.numeric(C_input["fine_roots"]) * g_f,
              branches = as.numeric(C_input["branches"]) * g_b,
              stems = as.numeric(C_input["stems"]) * g_s,
              understorey = as.numeric(C_input["understorey"]) * g_n)
  
  cat(sprintf("      [decompose_cohort] result values: needles=%.6f, fine_roots=%.6f, branches=%.6f, stems=%.6f, understorey=%.6f\n",
              result["needles"], result["fine_roots"], result["branches"], result["stems"], result["understorey"]))
  cat(sprintf("      [decompose_cohort] C_input total=%.6f, C_remain total=%.6f\n", 
              sum(C_input), sum(result, na.rm=TRUE)))
  
  result
}

# Annual step for ONE species
q_model_annual_step_species <- function(cohorts, C_input, temp, params, sp_params) {
  cat(sprintf("    [annual_step] Starting with %d cohorts, new input=%.6f\n", 
              length(cohorts), sum(C_input)))
  
  new_cohort <- list(C_input = C_input, age = 0, temp_input = temp)
  cohorts <- c(cohorts, list(new_cohort))
  
  cat(sprintf("    [annual_step] Now have %d cohorts (including new)\n", length(cohorts)))
  
  total_by_type <- c(needles=0, fine_roots=0, branches=0, stems=0, understorey=0)
  cohorts_keep <- list()
  
  n_filtered <- 0
  for(i in seq_along(cohorts)) {
    cohort <- cohorts[[i]]
    cohort$age <- cohort$age + 1
    
    cat(sprintf("    [annual_step] Processing cohort %d (age %d)\n", i, cohort$age))
    cat(sprintf("    [annual_step] Cohort %d C_input: names=%s\n", i, paste(names(cohort$C_input), collapse=", ")))
    cat(sprintf("    [annual_step] Cohort %d C_input: values=%s\n", i, paste(sprintf("%.6f", cohort$C_input), collapse=", ")))
    
    C_remain <- decompose_cohort_species(cohort$C_input, cohort$age, cohort$temp_input, params, sp_params)
    
    C_remain_sum <- sum(C_remain, na.rm=TRUE)
    cat(sprintf("    [annual_step] Cohort %d: C_remain_sum=%.9f\n", i, C_remain_sum))
    
    if(C_remain_sum > 1e-9) {
      cohort$C_remaining <- C_remain
      cohorts_keep <- c(cohorts_keep, list(cohort))
      total_by_type <- total_by_type + C_remain
      cat(sprintf("    [annual_step] Cohort %d KEPT\n", i))
    } else {
      n_filtered <- n_filtered + 1
      cat(sprintf("    [annual_step] Cohort %d FILTERED (below threshold)\n", i))
    }
  }
  
  cat(sprintf("    [annual_step] Kept %d cohorts, filtered %d cohorts\n", 
              length(cohorts_keep), n_filtered))
  cat(sprintf("    [annual_step] Total SOC = %.6f\n", sum(total_by_type, na.rm=TRUE)))
  
  list(cohorts = cohorts_keep, C_by_type = total_by_type, total_soc = sum(total_by_type, na.rm=TRUE))
}

# Spinup for ONE species
q_model_spinup_species <- function(C_input_annual, temp_mean, params, sp_params, n_years) {
  cat(sprintf("  [spinup] Starting spinup for %d years\n", n_years))
  cat(sprintf("  [spinup] Input: needles=%.6f, fine_roots=%.6f, branches=%.6f, stems=%.6f, understorey=%.6f\n",
              C_input_annual["needles"], C_input_annual["fine_roots"], 
              C_input_annual["branches"], C_input_annual["stems"], C_input_annual["understorey"]))
  cat(sprintf("  [spinup] Temperature: %.2f\n", temp_mean))
  
  cohorts <- list()
  for(y in 1:min(n_years, 5)) {  # Only show first 5 years
    cat(sprintf("  [spinup] === Year %d ===\n", y))
    result <- q_model_annual_step_species(cohorts, C_input_annual, temp_mean, params, sp_params)
    cohorts <- result$cohorts
    cat(sprintf("  [spinup] Year %d result: %d cohorts, SOC=%.6f\n", 
                y, length(cohorts), result$total_soc))
  }
  
  # Run remaining years silently
  if(n_years > 5) {
    cat(sprintf("  [spinup] Running years 6-%d silently...\n", n_years))
    for(y in 6:n_years) {
      result <- q_model_annual_step_species(cohorts, C_input_annual, temp_mean, params, sp_params)
      cohorts <- result$cohorts
    }
  }
  
  cat(sprintf("  [spinup] Final: %d cohorts\n", length(cohorts)))
  if(length(cohorts) > 0) {
    total_C <- sum(sapply(cohorts, function(c) sum(c$C_remaining, na.rm=TRUE)))
    cat(sprintf("  [spinup] Final SOC: %.6f\n", total_C))
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
  
  cat("[q_model_run] Starting Q-model run\n")
  
  n_years <- nrow(input_df)
  n_spinup <- if(is.null(spinup_years)) params$spinup_years else spinup_years
  
  cat(sprintf("[q_model_run] Simulation years: %d, Spinup years: %d\n", n_years, n_spinup))
  
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
  
  cat(sprintf("[q_model_run] Species proportions: pine=%.2f, spruce=%.2f, birch=%.2f\n",
              species_props["pine"], species_props["spruce"], species_props["birch"]))
  
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
  
  cat(sprintf("[q_model_run] First year inputs: needles=%.3f, fine_roots=%.3f, branches=%.3f, temp=%.2f\n",
              C_needles[1], C_fine_roots[1], C_branches[1], temp_mean_vec[1]))
  
  # Storage for combined results
  total_soc <- respiration <- numeric(n_years)
  pools_mat <- matrix(0, n_years, 5, 
                      dimnames = list(NULL, c("needles","fine_roots","branches","stems","understorey")))
  
  # Run each species separately and combine
  species_results <- list()
  
  for(sp in names(species_props)) {
    prop <- as.numeric(species_props[sp])  # Convert to numeric to avoid name issues
    if(prop <= 0) {
      cat(sprintf("[q_model_run] Skipping %s (proportion = 0)\n", sp))
      next
    }
    
    sp_params <- Q_SPECIES_PARAMS[Q_SPECIES_PARAMS$species == sp, ]
    if(nrow(sp_params) == 0) {
      warning(sprintf("Unknown species: %s, skipping", sp))
      next
    }
    
    cat(sprintf("\n[q_model_run] === Running Q-model for %s (%.0f%%) ===\n", sp, prop * 100))
    
    # Spinup for this species
    avg_input <- c(needles = mean(C_needles) * prop, 
                   fine_roots = mean(C_fine_roots) * prop,
                   branches = mean(C_branches) * prop, 
                   stems = mean(C_stems) * prop,
                   understorey = mean(C_understorey) * prop)
    
    cat(sprintf("[q_model_run] Species-weighted average inputs for spinup:\n"))
    cat(sprintf("  needles=%.6f, fine_roots=%.6f, branches=%.6f, stems=%.6f, understorey=%.6f\n",
                avg_input["needles"], avg_input["fine_roots"], avg_input["branches"],
                avg_input["stems"], avg_input["understorey"]))
    
    cohorts <- q_model_spinup_species(avg_input, mean(temp_mean_vec), params, sp_params, n_spinup)
    
    cat(sprintf("[q_model_run] %s spinup complete: %d cohorts\n", sp, length(cohorts)))
    
    # Run simulation for this species
    sp_soc <- numeric(n_years)
    sp_pools <- matrix(0, n_years, 5)
    
    cat(sprintf("[q_model_run] Starting main simulation for %s\n", sp))
    
    cat(sprintf("[q_model_run] Debug: C_needles[1] class=%s, names=%s\n", 
                class(C_needles[1]), paste(names(C_needles[1]), collapse=", ")))
    
    for(t in 1:min(n_years, 3)) {  # Only show first 3 years
      C_input <- c(needles = C_needles[t] * prop, 
                   fine_roots = C_fine_roots[t] * prop,
                   branches = C_branches[t] * prop, 
                   stems = C_stems[t] * prop,
                   understorey = C_understorey[t] * prop)
      
      cat(sprintf("[q_model_run] %s Year %d input: %.6f\n", sp, t, sum(C_input)))
      
      result <- q_model_annual_step_species(cohorts, C_input, temp_mean_vec[t], params, sp_params)
      cohorts <- result$cohorts
      sp_soc[t] <- result$total_soc
      sp_pools[t,] <- result$C_by_type
      
      cat(sprintf("[q_model_run] %s Year %d: %d cohorts, SOC=%.6f\n", 
                  sp, t, length(cohorts), result$total_soc))
    }
    
    # Run remaining years silently
    if(n_years > 3) {
      cat(sprintf("[q_model_run] Running years 4-%d silently...\n", n_years))
      for(t in 4:n_years) {
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
    }
    
    species_results[[sp]] <- list(soc = sp_soc, pools = sp_pools)
    cat(sprintf("[q_model_run] %s complete. Final SOC: %.6f\n", sp, sp_soc[n_years]))
  }
  
  # Combine species results
  cat("\n[q_model_run] Combining species results...\n")
  for(sp in names(species_results)) {
    cat(sprintf("  Adding %s: SOC[1]=%.6f\n", sp, species_results[[sp]]$soc[1]))
    total_soc <- total_soc + species_results[[sp]]$soc
    pools_mat <- pools_mat + species_results[[sp]]$pools
  }
  
  cat(sprintf("[q_model_run] Combined total SOC[1]=%.6f\n", total_soc[1]))
  
  # Calculate respiration
  C_input_total <- C_needles + C_branches + C_stems + C_fine_roots + C_understorey
  respiration[1] <- C_input_total[1]
  for(t in 2:n_years) {
    respiration[t] <- C_input_total[t] + total_soc[t-1] - total_soc[t]
  }
  
  cat("[q_model_run] Q-model run complete\n\n")
  
  list(pools = pools_mat, total_soc = total_soc, respiration = respiration,
       pool_names = c("needles","fine_roots","branches","stems","understorey"),
       species_props = species_props)
}
