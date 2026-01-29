# =============================================================================
# YASSO07 Soil Organic Carbon Model (3 litter types with wood size)
# =============================================================================
#
# Chemical fractionation model (AWENH pools), annual timestep.
# Handles 3 litter types: non-woody (nwl), fine woody (fwl), coarse woody (cwl)
# with different decomposition rates based on wood diameter.
#
# Reference: Tuomi et al. (2009) Ecol. Mod. 220:3362-3371
# =============================================================================

YASSO07_DEFAULT_PARAMS <- list(
  alpha_A=0.66, alpha_W=4.30, alpha_E=0.35, alpha_N=0.22,
  p_AW=0.32, p_AE=0.01, p_AN=0.93,
  p_WA=0.34, p_WE=0.00, p_WN=0.00,
  p_EA=0.00, p_EW=0.035, p_EN=0.005,
  p_NA=0.01, p_NW=0.0005, p_NE=0.03,
  beta1=0.096, beta2=-0.0014, gamma=-1.21,
  alpha_H=0.0033, p_H=0.04,
  delta1=-1.7, delta2=0.86, r=-0.306
)

yasso07_temp_mod <- function(temp_mean, temp_amp, beta1, beta2) {
  T1 <- temp_mean + 4*temp_amp/pi*(1/sqrt(2)-1)
  T2 <- temp_mean - 4*temp_amp/(sqrt(2)*pi)
  T3 <- temp_mean + 4*temp_amp/pi*(1-1/sqrt(2))
  T4 <- temp_mean + 4*temp_amp/(sqrt(2)*pi)
  mean(exp(beta1*c(T1,T2,T3,T4) + beta2*c(T1,T2,T3,T4)^2))
}

yasso07_precip_mod <- function(precip, gamma) 1 - exp(gamma * precip/1000)

yasso07_size_mod <- function(diameter, delta1, delta2, r) {
  if(is.null(diameter) || diameter <= 0) return(1)
  min(1, exp(-delta1 * diameter^delta2 * abs(r)))
}

yasso07_build_matrix <- function(params, temp_mod, precip_mod, size_mod = 1) {
  p <- params
  k <- -c(p$alpha_A, p$alpha_W, p$alpha_E, p$alpha_N) * temp_mod * precip_mod
  
  A <- matrix(0, 5, 5)
  A[1,1] <- k[1]; A[2,2] <- k[2]; A[3,3] <- k[3]; A[4,4] <- k[4] * size_mod
  
  # Transfers between AWEN
  A[2,1] <- -p$p_AW * k[1]; A[3,1] <- -p$p_AE * k[1]; A[4,1] <- -p$p_AN * k[1]
  A[1,2] <- -p$p_WA * k[2]; A[3,2] <- -p$p_WE * k[2]; A[4,2] <- -p$p_WN * k[2]
  A[1,3] <- -p$p_EA * k[3]; A[2,3] <- -p$p_EW * k[3]; A[4,3] <- -p$p_EN * k[3]
  A[1,4] <- -p$p_NA * k[4] * size_mod; A[2,4] <- -p$p_NW * k[4] * size_mod
  A[3,4] <- -p$p_NE * k[4] * size_mod
  
  # Humus
  A[5,5] <- -p$alpha_H * temp_mod * precip_mod
  A[5,1] <- -p$p_H * k[1]; A[5,2] <- -p$p_H * k[2]
  A[5,3] <- -p$p_H * k[3]; A[5,4] <- -p$p_H * k[4] * size_mod
  
  A
}

yasso07_annual_step <- function(pools, input_awen, params, temp_mean, temp_amp, precip, diameter = NULL) {
  temp_mod <- yasso07_temp_mod(temp_mean, temp_amp, params$beta1, params$beta2)
  precip_mod <- yasso07_precip_mod(precip, params$gamma)
  size_mod <- yasso07_size_mod(diameter, params$delta1, params$delta2, params$r)
  
  A <- yasso07_build_matrix(params, temp_mod, precip_mod, size_mod)
  
  input_vec <- c(input_awen, 0)  # No direct input to H
  
  # Matrix exponential solution
  expA <- as.matrix(Matrix::expm(A))
  
  if(max(abs(A)) < 1e-10) {
    new_pools <- pools + input_vec
  } else {
    Ainv <- solve(A)
    new_pools <- as.vector(expA %*% pools + (expA - diag(5)) %*% Ainv %*% input_vec)
  }
  
  resp <- sum(pools) + sum(input_vec) - sum(new_pools)
  
  list(pools = pmax(0, new_pools), respiration = max(0, resp))
}

# Steady state for one litter type
yasso07_steady_state <- function(params, input_awen, temp_mean, temp_amp, precip, diameter = NULL, n_years = 5000) {
  pools <- c(1, 1, 1, 1, 10)
  
  for(i in 1:n_years) {
    result <- yasso07_annual_step(pools, input_awen, params, temp_mean, temp_amp, precip, diameter)
    pools <- result$pools
  }
  pools
}

# Main run function - handles 3 litter types
#' @param params Model parameters
#' @param C0 Initial pools (list with nwl, fwl, cwl, each with 5-element vector)
#' @param input_df Annual input with columns for each litter type
yasso07_run <- function(params = YASSO07_DEFAULT_PARAMS, C0 = NULL, input_df, spinup_years = 5000) {
  
  if(!is.list(params)) params <- as.list(params)
  
  n_years <- nrow(input_df)
  
  # Wood diameters for each litter type (cm)
  wood_sizes <- list(
    nwl = 0,   # Non-woody (foliage, understorey, fine roots)
    fwl = 2,   # Fine woody (small branches)
    cwl = 15   # Coarse woody (large branches, stems)
  )
  
  # Storage for combined results
  pools_mat <- matrix(0, n_years, 5, dimnames = list(NULL, c("A","W","E","N","H")))
  total_soc <- respiration <- numeric(n_years)
  
  # Process each litter type separately
  litter_types <- c("nwl", "fwl", "cwl")
  
  for(lt in litter_types) {
    # Check if this litter type has input columns
    lt_a_col <- paste0(lt, "_A")
    if(!lt_a_col %in% names(input_df)) next
    
    input_A <- input_df[[paste0(lt, "_A")]]
    input_W <- input_df[[paste0(lt, "_W")]]
    input_E <- input_df[[paste0(lt, "_E")]]
    input_N <- input_df[[paste0(lt, "_N")]]
    
    # Initialize pools for this litter type
    if(is.null(C0) || !lt %in% names(C0)) {
      avg_input <- c(mean(input_A), mean(input_W), mean(input_E), mean(input_N))
      C0_lt <- yasso07_steady_state(
        params, avg_input,
        mean(input_df$temp_mean), 
        mean(input_df$temp_amplitude),
        mean(input_df$precip),
        diameter = wood_sizes[[lt]],
        n_years = spinup_years
      )
    } else {
      C0_lt <- C0[[lt]]
    }
    
    # Run simulation for this litter type
    pools_lt <- C0_lt
    lt_pools_mat <- matrix(0, n_years, 5)
    lt_pools_mat[1,] <- pools_lt
    
    for(t in 1:(n_years-1)) {
      input_awen <- c(input_A[t], input_W[t], input_E[t], input_N[t])
      
      result <- yasso07_annual_step(
        pools_lt, input_awen, params,
        input_df$temp_mean[t], 
        input_df$temp_amplitude[t],
        input_df$precip[t],
        diameter = wood_sizes[[lt]]
      )
      
      pools_lt <- result$pools
      lt_pools_mat[t+1,] <- pools_lt
      respiration[t+1] <- respiration[t+1] + result$respiration
    }
    
    # Add this litter type's pools to total
    pools_mat <- pools_mat + lt_pools_mat
  }
  
  # Calculate total SOC
  total_soc <- rowSums(pools_mat)
  
  list(
    pools = pools_mat, 
    total_soc = total_soc, 
    respiration = respiration,
    pool_names = c("A","W","E","N","H")
  )
}