# =============================================================================
# Yasso07 with optional intact wood pool (fragmentation delay)
# =============================================================================
#
# Modified from the original yasso07.light by Taru Palosuo (FMI, 2011)
# Based on Yasso07 description Tuomi & Liski 17.3.2008
#
# MODIFICATION:
# When UseIntactPool = TRUE, a 6th pool ("intact wood") is added.
# For woody litter (WoodySize > 0), a fraction f(d) of the AWEN input
# is routed to this pool instead of entering AWEN directly.
# The intact pool drains into AWEN at a first-order rate (FragRate),
# preserving the original chemical composition (AWEN fractions).
#
# This replaces the standard size-dependent rate reduction with a 
# fragmentation delay, producing a sigmoid-like mass loss curve for 
# coarse woody debris — conceptually similar to the Q-model's Tmax.
#
# When UseIntactPool = FALSE, the function behaves identically to
# the original yasso07.light.
#
# Copyright (C) <2017>  <Finnish Meteorological Institute>
# GNU General Public License v3
# =============================================================================

library(Matrix)

yasso07.intact = function(MeanTemperature, 
                          TemperatureAmplitude, 
                          Precipitation,
                          InitialCPool,       # length 5 (standard) or 6 (with intact pool)
                          LitterInput,        # AWENH, length 5
                          WoodySize, 
                          Yasso07Parameters,  # length 44
                          SimulationTime,
                          UseIntactPool = FALSE,
                          FragRate = 0.05) {  # first-order fragmentation rate [yr-1]
  
  # -----------------------------------------------------------
  # Parse inputs (same as original)
  # -----------------------------------------------------------
  MT = MeanTemperature
  TA = TemperatureAmplitude
  PR = Precipitation / 1000     # mm to meters
  LI = LitterInput
  PA = Yasso07Parameters
  WS = WoodySize
  TI = SimulationTime
  
  # -----------------------------------------------------------
  # Decomposition rates (alpha)
  # -----------------------------------------------------------
  alfa = c(-PA[1], -PA[2], -PA[3], -PA[4], -PA[35])
  
  # -----------------------------------------------------------
  # Transfer matrix p (5x5)
  # -----------------------------------------------------------
  p = matrix(c(
    -1,     PA[5],  PA[6],  PA[7],  0,
    PA[8],  -1,     PA[9],  PA[10], 0,
    PA[11], PA[12], -1,     PA[13], 0,
    PA[14], PA[15], PA[16], -1,     0,
    PA[36], PA[36], PA[36], PA[36], -1
  ), 5, 5, byrow = TRUE)
  
  # -----------------------------------------------------------
  # Climate modifiers (same as original)
  # -----------------------------------------------------------
  beta1 = PA[17]
  beta2 = PA[18]
  gamma = PA[26]
  
  T1 = MT + 4 * TA / pi * (1/sqrt(2) - 1)
  T2 = MT - 4 * TA / (sqrt(2) * pi)
  T3 = MT + 4 * TA / pi * (1 - 1/sqrt(2))
  T4 = MT + 4 * TA / (sqrt(2) * pi)
  
  temps = c(T1, T2, T3, T4)
  clim = mean(exp(beta1 * temps + beta2 * temps^2) * (1 - exp(gamma * PR)))
  k = alfa * clim
  
  # -----------------------------------------------------------
  # Size dependence parameters
  # -----------------------------------------------------------
  delta1 = PA[39]
  delta2 = PA[40]
  r      = PA[41]
  
  # Size modifier: 1.0 for d=0, decreasing for larger pieces
  size_dep = (1 + delta1 * WS + delta2 * WS^2)^r
  
  # ==========================================================
  # STANDARD YASSO07 (original behavior)
  # ==========================================================
  if (!UseIntactPool || WS == 0) {
    
    # Apply size dependence to AWEN rates (not H)
    k = c(k[1:4] * size_dep, k[5])
    
    A = p %*% diag(k)
    
    LC = as.array(solve(A) %*% (expm(A * TI) %*% (A %*% InitialCPool[1:5] + LI) - LI))
    return(LC)
  }
  
  # ==========================================================
  # MODIFIED YASSO07 WITH INTACT WOOD POOL
  # ==========================================================
  
  # Fraction routed to intact pool: Weibull CDF
  # Anchored at: d=0 -> 0, d=1 -> 0.05, d=40 -> ~0.99
  # d0 and p_shape are candidates for calibration later
  d0 = 11.4       # scale parameter [cm]
  p_shape = 1.22  # shape parameter
  f_intact = 1 - exp(-(WS / d0)^p_shape)
  
  # Keep size_dep on AWEN rates: fragmented material from large pieces
  # is still chunkier than fine litter. The intact pool adds a delay
  # ON TOP of the rate reduction, not instead of it.
  k_6 = c(k[1:4] * size_dep, k[5])
  
  # --- Build 6x6 system ---
  # Pools: A(1), W(2), E(3), N(4), H(5), I(6)
  
  A6 = matrix(0, 6, 6)
  
  # Top-left 5x5: standard AWEN-H dynamics (no size reduction on rates)
  A6[1:5, 1:5] = p %*% diag(k_6)
  
  # Column 6: transfers FROM intact pool TO AWEN pools
  # Released material has the same AWEN composition as the original litter
  awen_input = LI[1:4]
  awen_total = sum(awen_input)
  
  if (awen_total > 0) {
    awen_frac = awen_input / awen_total   # normalised AWEN fractions
  } else {
    awen_frac = c(0.25, 0.25, 0.25, 0.25)  # fallback (shouldn't happen)
  }
  
  # Fragmentation releases into AWEN at rate FragRate, partitioned by composition
  A6[1, 6] =  FragRate * awen_frac[1]    # I -> A
  A6[2, 6] =  FragRate * awen_frac[2]    # I -> W
  A6[3, 6] =  FragRate * awen_frac[3]    # I -> E
  A6[4, 6] =  FragRate * awen_frac[4]    # I -> N
  A6[5, 6] =  0                           # I -> H (none, fresh material)
  
  # Diagonal: intact pool loses mass at FragRate
  A6[6, 6] = -FragRate
  
  # Row 6, cols 1-5: no transfer FROM AWEN TO intact (one-way)
  # (already zero from initialization)
  
  # --- Split litter input ---
  # Direct to AWEN: (1 - f_intact) of original
  # To intact pool: f_intact of total AWEN carbon
  LI6 = c(LI[1:4] * (1 - f_intact),   # reduced direct AWEN input
          LI[5],                        # H input unchanged
          awen_total * f_intact)        # total AWEN carbon routed to intact
  
  # --- Initial conditions ---
  # Accept either 5 or 6 element initial pool
  if (length(InitialCPool) == 5) {
    IC6 = c(InitialCPool, 0)   # intact pool starts empty
  } else {
    IC6 = InitialCPool[1:6]
  }
  
  # --- Solve (same analytical solution, now 6x6) ---
  LC = as.array(solve(A6) %*% (expm(A6 * TI) %*% (A6 %*% IC6 + LI6) - LI6))
  
  # Return all 6 pools (user can sum 1:5 for comparable AWENH total)
  LC
}