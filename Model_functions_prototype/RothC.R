# =============================================================================
# ROTHC Soil Organic Carbon Model
# =============================================================================
#
# VERSION 2.2 - FIXED temperature modifier usage
# Last updated: 2026-01-27
#
# Five-pool model (DPM, RPM, BIO, HUM, IOM), monthly timestep.
#
# =============================================================================

ROTHC_DEFAULT_PARAMS <- list(
  k_DPM = 10.0,      # Decomposable plant material (yr^-1)
  k_RPM = 0.3,       # Resistant plant material (yr^-1)
  k_BIO = 0.66,      # Microbial biomass (yr^-1)
  k_HUM = 0.02,      # Humified organic matter (yr^-1)
  clay_effect = TRUE # Whether to use clay modifier
)

# Temperature rate modifier (returns value 0 to ~47.9)
# At 9.25°C this returns ~1.0 (reference conditions)
rothc_temp_mod <- function(temp) {
  ifelse(temp < -5, 0, 47.9 / (1 + exp(106 / (temp + 18.3))))
}

# Moisture rate modifier (topsoil moisture deficit)
rothc_moisture_mod <- function(precip, evap, clay, depth = 23, max_tsmd = NULL) {
  # Calculate max TSMD based on clay content and depth (from Fortran)
  # SMD15bar = -(20 + 1.3*clay - 0.01*clay^2)
  # SMD15barAdj = SMD15bar * depth / 23
  # SMD1bar = 0.444 * SMD15barAdj  (this is the "max TSMD" for bare soil)
  
  if (is.null(max_tsmd)) {
    SMD15bar <- -(20 + 1.3 * clay - 0.01 * clay^2)
    SMD15barAdj <- SMD15bar * depth / 23
    max_tsmd <- abs(0.444 * SMD15barAdj)  # Use absolute value
  }
  
  # Accumulated deficit
  tsmd <- pmax(0, evap - precip)
  tsmd <- pmin(tsmd, max_tsmd)
  
  # Rate modifier (0.2 to 1.0)
  pmax(0.2, 1.0 - 0.8 * tsmd / max_tsmd)
}

# Soil cover modifier
rothc_cover_mod <- function(soil_cover) {
  ifelse(soil_cover > 0.5, 0.6, 1.0)
}

# Clay modifier for BIO/HUM formation
rothc_clay_mod <- function(clay) {
  1.67 * (1.85 + 1.60 * exp(-0.0786 * clay))
}

# Monthly decomposition step
rothc_monthly_step <- function(pools, input_dpm, input_rpm, params, 
                               temp, precip, evap, clay, soil_cover = 1) {
  
  # Rate modifiers (NOT normalized - used directly as in Fortran)
  temp_mod <- rothc_temp_mod(temp)
  moist_mod <- rothc_moisture_mod(precip, evap, clay)
  cover_mod <- rothc_cover_mod(soil_cover)
  
  # Combined rate modifier (RateM in Fortran)
  rate_mod <- temp_mod * moist_mod * cover_mod
  
  # Monthly timestep
  tstep <- 1/12
  
  # Clay modifier for BIO/HUM partitioning
  x <- rothc_clay_mod(clay)
  f_CO2 <- x / (x + 1)
  f_BIO <- 0.46 / (x + 1)
  f_HUM <- 0.54 / (x + 1)
  
  if (!params$clay_effect) {
    f_CO2 <- 0.628
    f_BIO <- 0.172
    f_HUM <- 0.200
  }
  
  # Decomposition using exponential decay (as in Fortran)
  # C1 = C * exp(-RateM * k * tstep)
  C_DPM_remaining <- pools[1] * exp(-rate_mod * params$k_DPM * tstep)
  C_RPM_remaining <- pools[2] * exp(-rate_mod * params$k_RPM * tstep)
  C_BIO_remaining <- pools[3] * exp(-rate_mod * params$k_BIO * tstep)
  C_HUM_remaining <- pools[4] * exp(-rate_mod * params$k_HUM * tstep)
  
  # Amount decomposed
  C_DPM_loss <- pools[1] - C_DPM_remaining
  C_RPM_loss <- pools[2] - C_RPM_remaining
  C_BIO_loss <- pools[3] - C_BIO_remaining
  C_HUM_loss <- pools[4] - C_HUM_remaining
  
  # CO2 emissions
  resp <- (C_DPM_loss + C_RPM_loss + C_BIO_loss + C_HUM_loss) * f_CO2
  
  # Transfers to BIO and HUM
  to_BIO <- (C_DPM_loss + C_RPM_loss + C_BIO_loss + C_HUM_loss) * f_BIO
  to_HUM <- (C_DPM_loss + C_RPM_loss + C_BIO_loss + C_HUM_loss) * f_HUM
  
  # Update pools
  new_pools <- pools
  new_pools[1] <- C_DPM_remaining + input_dpm
  new_pools[2] <- C_RPM_remaining + input_rpm
  new_pools[3] <- C_BIO_remaining + to_BIO
  new_pools[4] <- C_HUM_remaining + to_HUM
  # new_pools[5] (IOM) unchanged
  
  list(pools = pmax(0, new_pools), respiration = max(0, resp))
}

# Steady state - ANALYTICAL SOLUTION
rothc_steady_state <- function(params, input_dpm, input_rpm, temp_mean, 
                               precip_mean, evap_mean, clay, n_years = 10000) {
  
  # Calculate IOM from clay (Falloon et al. equation)
  iom <- 0.049 * clay^1.139
  
  # Calculate average rate modifier (NOT normalized by 47.9!)
  temp_mod <- rothc_temp_mod(temp_mean)
  moist_mod <- rothc_moisture_mod(precip_mean / 12, evap_mean / 12, clay)
  cover_mod <- 0.6  # Assume vegetated
  
  # This is the annual average of the combined modifier
  # In reality you'd want to average monthly values, but for steady state
  # we use mean conditions
  rate_mod_annual <- temp_mod * moist_mod * cover_mod
  
  # Clay modifier for partitioning
  x <- rothc_clay_mod(clay)
  f_CO2 <- x / (x + 1)
  f_BIO <- 0.46 / (x + 1)
  f_HUM <- 0.54 / (x + 1)
  
  # Effective annual decomposition rates
  k_DPM <- params$k_DPM * rate_mod_annual
  k_RPM <- params$k_RPM * rate_mod_annual
  k_BIO <- params$k_BIO * rate_mod_annual
  k_HUM <- params$k_HUM * rate_mod_annual
  
  # Matrix A for steady-state solution
  A <- matrix(0, 4, 4)
  
  # Diagonal: net loss from each pool
  A[1,1] <- -k_DPM
  A[2,2] <- -k_RPM
  A[3,3] <- -k_BIO * (1 - f_BIO)  # Net loss (some returns to BIO)
  A[4,4] <- -k_HUM * (1 - f_HUM)  # Net loss (some returns to HUM)
  
  # Off-diagonal: transfers
  A[3,1] <- k_DPM * f_BIO   # DPM -> BIO
  A[3,2] <- k_RPM * f_BIO   # RPM -> BIO
  A[3,4] <- k_HUM * f_BIO   # HUM -> BIO
  A[4,1] <- k_DPM * f_HUM   # DPM -> HUM
  A[4,2] <- k_RPM * f_HUM   # RPM -> HUM
  A[4,3] <- k_BIO * f_HUM   # BIO -> HUM
  
  # Input vector (annual)
  b <- c(input_dpm, input_rpm, 0, 0)
  
  # Solve: A * C = -b
  C_active <- solve(A, -b)
  
  pools <- c(DPM = C_active[1], RPM = C_active[2], 
             BIO = C_active[3], HUM = C_active[4], IOM = iom)
  
  pmax(0, pools)
}

# Main run function
rothc_run <- function(params = ROTHC_DEFAULT_PARAMS, C0 = NULL, input_df, spinup_years = 10000) {
  
  if(!is.list(params)) params <- as.list(params)
  
  n_months <- nrow(input_df)
  
  pools_mat <- matrix(0, n_months, 5, 
                      dimnames = list(NULL, c("DPM","RPM","BIO","HUM","IOM")))
  total_soc <- respiration <- numeric(n_months)
  
  # Get input columns
  C_DPM <- input_df$C_DPM
  C_RPM <- input_df$C_RPM
  temp <- input_df$temp_air
  precip <- input_df$precip
  evap <- input_df$evap
  clay <- input_df$clay[1]
  soil_cover <- if("soil_cover" %in% names(input_df)) input_df$soil_cover else rep(1, n_months)
  
  # Initialize
  if(is.null(C0)) {
    C0 <- rothc_steady_state(
      params,
      mean(C_DPM) * 12,  # Annual input
      mean(C_RPM) * 12,
      mean(temp),
      mean(precip) * 12,
      mean(evap) * 12,
      clay,
      n_years = spinup_years
    )
  }
  
  pools <- C0
  iom_fixed <- pools[5]
  pools_mat[1,] <- pools
  total_soc[1] <- sum(pools)
  
  for(t in 1:(n_months-1)) {
    result <- rothc_monthly_step(
      pools, C_DPM[t], C_RPM[t], params,
      temp[t], precip[t], evap[t], clay, soil_cover[t]
    )
    
    pools <- result$pools
    pools[5] <- iom_fixed
    pools_mat[t+1,] <- pools
    total_soc[t+1] <- sum(pools)
    respiration[t+1] <- result$respiration
  }
  
  list(
    pools = pools_mat,
    total_soc = total_soc,
    respiration = respiration,
    pool_names = c("DPM","RPM","BIO","HUM","IOM")
  )
}


# =============================================================================
# Quick test
# =============================================================================
if (FALSE) {
  test_input <- data.frame(
    C_DPM = rep(0.1, 120),
    C_RPM = rep(0.2, 120),
    temp_air = rep(c(-8, -6, -2, 4, 10, 15, 17, 15, 10, 4, -2, -6), 10),
    precip = rep(40, 120),
    evap = rep(30, 120),
    clay = rep(15, 120)
  )
  
  result <- rothc_run(ROTHC_DEFAULT_PARAMS, input_df = test_input)
  
  cat("Steady-state pools (t C/ha):\n")
  print(round(result$pools[1,], 2))
  cat("\nTotal SOC:", round(result$total_soc[1], 2), "t C/ha\n")
}