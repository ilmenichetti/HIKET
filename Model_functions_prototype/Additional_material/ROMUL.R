# =============================================================================
# ROMULv Soil Organic Carbon Model
# =============================================================================
#
# Cohort-based decomposition with N cycling and volumetric soil moisture.
# Accepts monthly input, runs daily timestep internally.
#
# References:
#   Chertov et al. (2001), Linkosalo et al. (2013), Lehtonen et al. (2016)
# =============================================================================

ROMULV_DEFAULT_PARAMS <- list(
  kL = 0.12, DeposN = 2.0, ML = 0.1, gamma = 0.8, C_C = 0.5,
  DeltaB = 24.0, DeltaL = 12.8, k6 = 0.00006,
  ash_content = c(foliage=0.02, branches=0.02, stems=0.01,
                  fine_roots=0.02, coarse_roots=0.02, understorey=0.04),
  N_content = c(foliage=0.003, branches=0.004, stems=0.0014,
                fine_roots=0.0047, coarse_roots=0.0024, understorey=0.006),
  soil_temp_min = c(organic=-0.13, mineral=0.24),
  soil_temp_tau = c(organic=14.9, mineral=10.5),
  ET_ratio = 0.256, tau_soil = c(organic=0.89, mineral=9.42)
)

# Temperature response functions
romulv_f1 <- function(T) { r <- numeric(length(T)); r[T > -5 & T <= 1] <- 0.1595 + 0.0319*T[T > -5 & T <= 1]; r[T > 1 & T <= 35] <- 0.1754*exp(0.0871*T[T > 1 & T <= 35]); r[T > 35 & T <= 60] <- 8.791 - 0.1465*T[T > 35 & T <= 60]; r }
romulv_f2 <- function(T) { r <- numeric(length(T)); r[T > -5 & T <= 1] <- 0.1595 + 0.0319*T[T > -5 & T <= 1]; r[T > 1 & T <= 35] <- 0.1754*exp(0.0871*T[T > 1 & T <= 35]); r[T > 35 & T <= 60] <- 3.690 - 0.0615*T[T > 35 & T <= 60]; r }
romulv_f3 <- function(T) { r <- numeric(length(T)); r[T > -3 & T <= 7] <- 1.3; r[T > 7 & T <= 60] <- 1.472 - T[T > 7 & T <= 60]*0.0245; r }
romulv_f4 <- function(T) { r <- numeric(length(T)); r[T > -5 & T <= 1] <- 0.1595 + 0.0319*T[T > -5 & T <= 1]; r[T > 1 & T <= 20] <- 0.1754*exp(0.0871*T[T > 1 & T <= 20]); r[T > 20 & T <= 40] <- 1; r[T > 40 & T <= 80] <- 2.0 - 0.025*T[T > 40 & T <= 80]; r }
romulv_f5 <- function(T) { r <- numeric(length(T)); r[T > -5 & T <= 1] <- 0.078 + 0.0156*T[T > -5 & T <= 1]; r[T > 1 & T <= 13] <- 0.0675*exp(0.2088*T[T > 1 & T <= 13]); r[T > 13 & T <= 25] <- 1; r[T > 25 & T <= 50] <- 2.0 - 0.04*T[T > 25 & T <= 50]; r }
romulv_f6 <- function(T) { r <- numeric(length(T)); r[T > -5 & T <= 1] <- 0.1595 + 0.0319*T[T > -5 & T <= 1]; r[T > 1 & T <= 27.5] <- 0.1754*exp(0.0871*T[T > 1 & T <= 27.5]); r[T > 27.5 & T <= 35] <- 1.95; r[T > 35 & T <= 60] <- 4.68 - 0.078*T[T > 35 & T <= 60]; r }

romulv_moisture <- function(theta) { theta <- pmax(0, pmin(1, theta)); pmin(3.83*theta^1.25, 4.43*(1-theta)^0.8854, 1.0) }

# Rate coefficients
romulv_k1L <- function(ash, N) 0.0005 + 0.54*N
romulv_k1S <- function(ash, N) 0.0136 + 0.06*ash
romulv_k2L <- function(ash, N) 0.00060
romulv_k2S <- function(ash, N) 0.00126
romulv_k3L <- function(ash, N) 0.0089 + 0.78*N
romulv_k3S <- function(ash, N) if(ash < 0.18) 0.0394 - 0.21*ash else 0.0016
romulv_k4 <- function(ash, N) if(N <= 0.02) 0.05*N else 0.001
romulv_k5 <- function(ash, N) if(N <= 0.005) 0 else if(N >= 0.02) 0.007 else 0.007*(100*(2*N-0.01)/3)
romulv_k6 <- function() 0.00006

# Daily decomposition step
romulv_daily_step <- function(state, litter, litter_N, soil_temp, theta, params) {
  L <- state$L; F_ <- state$F; H <- state$H
  LN <- state$LN; FN <- state$FN; HN <- state$HN; Navail <- state$Navail
  ash <- params$ash_content; ML <- params$ML; gamma <- params$gamma
  DeltaB <- params$DeltaB; DeltaL <- params$DeltaL
  
  # Replace any NAs with 0
  L[is.na(L)] <- 0; F_[is.na(F_)] <- 0; LN[is.na(LN)] <- 0; FN[is.na(FN)] <- 0
  if(is.na(H)) H <- 0; if(is.na(HN)) HN <- 0; if(is.na(Navail)) Navail <- 0
  litter[is.na(litter)] <- 0; litter_N[is.na(litter_N)] <- 0
  
  org <- c(1,2,3,6); min_ <- c(4,5)
  F_org_sum <- sum(F_[org], na.rm=TRUE); FN_org_sum <- sum(FN[org], na.rm=TRUE)
  F_min_sum <- sum(F_[min_], na.rm=TRUE); FN_min_sum <- sum(FN[min_], na.rm=TRUE)
  NFconc_org <- if(F_org_sum > 1e-10) FN_org_sum/F_org_sum else 0.01
  NFconc_min <- if(F_min_sum > 1e-10) FN_min_sum/F_min_sum else 0.01
  
  moist_org <- romulv_moisture(theta[1]); moist_min <- romulv_moisture(theta[2])
  T_org <- soil_temp[1]; T_min <- soil_temp[2]
  
  fn <- litter_N / pmax(litter, 1e-10)
  kDL <- kDF <- kTL <- kTF <- kTF2 <- numeric(6)
  
  for(i in 1:6) {
    if(i %in% min_) {
      kDL[i] <- romulv_k1S(ash[i],fn[i])*romulv_f1(T_min)*moist_min
      kDF[i] <- romulv_k2S(ash[i],NFconc_min)*romulv_f2(T_min)*moist_min
      kTL[i] <- romulv_k3S(ash[i],fn[i])*romulv_f3(T_min)*moist_min
    } else {
      kDL[i] <- romulv_k1L(ash[i],fn[i])*romulv_f1(T_org)*moist_org
      kDF[i] <- romulv_k2L(ash[i],NFconc_org)*romulv_f2(T_org)*moist_org
      kTL[i] <- romulv_k3L(ash[i],fn[i])*romulv_f3(T_org)*moist_org
    }
    kTF[i] <- romulv_k4(ash[i],NFconc_min)*romulv_f4(T_min)*moist_min
    kTF2[i] <- romulv_k5(ash[i],NFconc_min)*romulv_f5(T_min)*moist_min
  }
  kDH <- romulv_k6()*romulv_f6(T_min)*moist_min
  
  dL <- litter - (kDL + kTL)*L
  dF <- kTL*L - (kDF + kTF + kTF2)*F_
  FH_fluxB <- sum(kTF*F_); FH_fluxL <- sum(kTF2*F_)
  FHN_fluxB <- sum(kTF*FN); FHN_fluxL <- sum(kTF2*FN)
  dH <- DeltaB*FHN_fluxB + DeltaL*FHN_fluxL - kDH*H
  
  N_diff <- 100*NFconc_min - 1.16*100*NFconc_org
  MF <- if(is.na(N_diff) || N_diff <= 0.44) 0.1 else if(N_diff <= 1.50) 0.5 else 1.0
  
  H_safe <- if(is.na(H) || is.na(HN)) FALSE else (HN > 1e-10 && H/HN/2 > 8)
  MH <- if(H_safe) 0.8 else 1.0
  
  dLN <- litter_N - (ML*kDL + kTL)*LN
  dFN <- kTL*LN - MF*kDF*FN - (kTF + kTF2)*FN
  dHN <- -kDH*MH*HN + gamma*(FHN_fluxB + FHN_fluxL)
  
  N_rel <- kDH*MH*HN + (1-gamma)*(FHN_fluxB+FHN_fluxL) + params$DeposN/365 +
    sum(ML*kDL*LN) + sum(MF*kDF*FN)
  dNavail <- N_rel - max(0, params$kL*Navail/365)
  
  resp <- params$C_C * (sum(kDL*L) + sum(kDF*F_) + kDH*H + FH_fluxB + FH_fluxL -
                          DeltaB*FHN_fluxB - DeltaL*FHN_fluxL)
  
  list(L=L+dL, F=F_+dF, H=H+dH, LN=LN+dLN, FN=FN+dFN, HN=HN+dHN,
       Navail=Navail+dNavail, respiration=resp)
}

# Monthly step (30 daily steps)
romulv_monthly_step <- function(state, litter, litter_N, temp_air, precip, SW_max, params) {
  # Ensure scalars
  temp_air <- as.numeric(temp_air)[1]
  precip <- as.numeric(precip)[1]
  SW_max <- as.numeric(SW_max)
  litter <- as.numeric(litter)
  litter_N <- as.numeric(litter_N)
  
  days <- 30
  daily_litter <- litter/days; daily_N <- litter_N/days
  daily_precip <- precip/days; ET <- 1.5
  
  if(is.null(state$soil_water)) state$soil_water <- SW_max * 0.6
  if(is.null(state$soil_temp)) state$soil_temp <- c(temp_air, temp_air)
  
  resp <- 0
  for(d in 1:days) {
    state$soil_temp <- state$soil_temp + (pmax(temp_air, params$soil_temp_min) - state$soil_temp)/params$soil_temp_tau
    
    # Soil water update
    sat <- SW_max/0.65; FC <- SW_max
    sw <- state$soil_water
    et <- c(params$ET_ratio*ET, (1-params$ET_ratio)*ET)
    overflow <- if(sw[1] > FC[1]) (sw[1]-FC[1])/params$tau_soil[1] else 0
    sw[1] <- pmin(sat[1], pmax(0, min(sw[1],FC[1]) + daily_precip - et[1]))
    if(sw[2] > FC[2]) sw[2] <- sw[2] - (sw[2]-FC[2])/params$tau_soil[2]
    sw[2] <- pmin(sat[2], pmax(0, sw[2] + overflow - et[2]))
    state$soil_water <- sw
    theta <- pmax(0, pmin(1, sw/sat))
    
    r <- romulv_daily_step(state, daily_litter, daily_N, state$soil_temp, theta, params)
    state$L <- r$L; state$F <- r$F; state$H <- r$H
    state$LN <- r$LN; state$FN <- r$FN; state$HN <- r$HN; state$Navail <- r$Navail
    resp <- resp + r$respiration
  }
  list(state=state, respiration=resp, total_soc=sum(state$L)+sum(state$F)+state$H)
}

# Main run function (expects monthly input)
romulv_run <- function(params = ROMULV_DEFAULT_PARAMS, C0 = NULL, input_df) {
  N_frac <- params$N_content
  
  if(!"month" %in% names(input_df)) stop("ROMULv requires monthly input data")
  
  years <- unique(input_df$year); n_years <- length(years)
  
  pools <- matrix(0, n_years, 3, dimnames=list(NULL, c("L","F","H")))
  total_soc <- respiration <- numeric(n_years)
  
  # Initialize
  state <- if(is.null(C0)) romulv_steady_state(params, input_df[1,], N_frac) else C0
  pools[1,] <- c(sum(state$L), sum(state$F), state$H)
  total_soc[1] <- sum(pools[1,])
  
  for(y in 1:(n_years-1)) {
    ydata <- input_df[input_df$year == years[y],]
    ann_resp <- 0
    
    for(m in 1:nrow(ydata)) {
      row <- ydata[m,]
      
      litter <- c(row$C_foliage, row$C_branches, row$C_stems,
                  row$C_fine_roots, row$C_coarse_roots, row$C_understorey)
      litter[is.na(litter)] <- 0
      
      litter_N <- if("N_foliage" %in% names(row)) {
        c(row$N_foliage, row$N_branches, row$N_stems,
          row$N_fine_roots, row$N_coarse_roots, row$N_understorey)
      } else litter * N_frac
      litter_N[is.na(litter_N)] <- 0
      
      SW_max <- c(row$SW_max_organic, row$SW_max_mineral)
      
      r <- romulv_monthly_step(state, litter, litter_N, row$temp_air, row$precip, SW_max, params)
      state <- r$state; ann_resp <- ann_resp + r$respiration
    }
    pools[y+1,] <- c(sum(state$L), sum(state$F), state$H)
    total_soc[y+1] <- r$total_soc; respiration[y+1] <- ann_resp
  }
  list(pools=pools, total_soc=total_soc, respiration=respiration, pool_names=c("L","F","H"))
}

# Steady state initialization
romulv_steady_state <- function(params, input_row, N_frac, n_iter=500) {
  state <- list(L=rep(1,6), F=rep(5,6), H=50, LN=rep(0.01,6), FN=rep(0.05,6),
                HN=2, Navail=10, soil_water=NULL, soil_temp=NULL)
  
  litter <- c(
    if("C_foliage" %in% names(input_row)) as.numeric(input_row$C_foliage) else 0.1,
    if("C_branches" %in% names(input_row)) as.numeric(input_row$C_branches) else 0.03,
    if("C_stems" %in% names(input_row)) as.numeric(input_row$C_stems) else 0,
    if("C_fine_roots" %in% names(input_row)) as.numeric(input_row$C_fine_roots) else 0.06,
    if("C_coarse_roots" %in% names(input_row)) as.numeric(input_row$C_coarse_roots) else 0.02,
    if("C_understorey" %in% names(input_row)) as.numeric(input_row$C_understorey) else 0.02
  )
  litter_N <- litter * N_frac
  temp <- if("temp_air" %in% names(input_row)) as.numeric(input_row$temp_air) else 4.5
  prec <- if("precip" %in% names(input_row)) as.numeric(input_row$precip) else 50
  SW_max <- c(
    if("SW_max_organic" %in% names(input_row)) as.numeric(input_row$SW_max_organic) else 20,
    if("SW_max_mineral" %in% names(input_row)) as.numeric(input_row$SW_max_mineral) else 120
  )
  
  for(i in 1:(n_iter*12)) {
    m <- ((i-1) %% 12) + 1
    temp_m <- temp + 12*sin(2*pi*(m-4)/12)
    prec_m <- prec * c(.07,.06,.07,.08,.08,.09,.10,.10,.09,.09,.09,.08)[m] * 12
    r <- romulv_monthly_step(state, litter, litter_N, temp_m, prec_m, SW_max, params)
    state <- r$state
  }
  state
}