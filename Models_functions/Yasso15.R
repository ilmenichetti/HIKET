# =============================================================================
# YASSO15 Soil Organic Carbon Model (R implementation)
# Structure compatible with YASSO07
#
# It also uses the same input mapping function
#
# =============================================================================

YASSO15_DEFAULT_PARAMS <- list(
  alpha_A = 0.66, alpha_W = 4.30, alpha_E = 0.35, alpha_N = 0.22,
  
  p_WA = 0.32, p_EA = 0.01, p_NA = 0.93,
  p_AW = 0.34, p_EW = 0.00, p_NW = 0.00,
  p_AE = 0.00, p_WE = 0.035, p_NE = 0.005,
  p_AN = 0.01, p_WN = 0.0005, p_EN = 0.03,
  
  beta1 = 0.096, beta2 = -0.0014,
  betaN1 = 0.096, betaN2 = -0.0014,
  betaH1 = 0.096, betaH2 = -0.0014,
  
  gamma = -1.21, gammaN = -1.21, gammaH = -1.21,
  
  p_H = 0.04, alpha_H = 0.0033,
  
  theta1 = -1.7, theta2 = 0.86, r = -0.306
)

yasso15_temp_mod <- function(temp_mean, temp_amp, beta1, beta2) {
  pi <- base::pi
  te <- c(
    temp_mean + 4*temp_amp*(1/sqrt(2)-1)/pi,
    temp_mean - 4*temp_amp/(sqrt(2)*pi),
    temp_mean + 4*temp_amp*(1-1/sqrt(2))/pi,
    temp_mean + 4*temp_amp/(sqrt(2)*pi)
  )
  mean(exp(beta1*te + beta2*te^2))
}

yasso15_precip_mod <- function(precip, gamma) {
  1 - exp(gamma * precip/1000)
}

yasso15_size_mod <- function(d, theta1, theta2, r) {
  if(is.null(d) || d == 0) return(1)
  min(1, (1 + theta1*d + theta2*d^2)^(-abs(r)))
}

yasso15_build_matrix <- function(params, climate, d, leac = 0) {
  
  temp_mean <- climate[1]
  precip <- climate[2]
  temp_amp <- climate[3]
  
  tem  <- yasso15_temp_mod(temp_mean, temp_amp, params$beta1,  params$beta2)
  temN <- yasso15_temp_mod(temp_mean, temp_amp, params$betaN1, params$betaN2)
  temH <- yasso15_temp_mod(temp_mean, temp_amp, params$betaH1, params$betaH2)
  
  tem  <- tem  * yasso15_precip_mod(precip, params$gamma)
  temN <- temN * yasso15_precip_mod(precip, params$gammaN)
  temH <- temH * yasso15_precip_mod(precip, params$gammaH)
  
  size_dep <- yasso15_size_mod(d, params$theta1, params$theta2, params$r)
  
  A <- matrix(0, 5, 5)
  
  A[1,1] <- -abs(params$alpha_A)*tem*size_dep
  A[2,2] <- -abs(params$alpha_W)*tem*size_dep
  A[3,3] <- -abs(params$alpha_E)*tem*size_dep
  A[4,4] <- -abs(params$alpha_N)*temN*size_dep
  A[5,5] <- -abs(params$alpha_H)*temH
  
  A[1,2] <- params$p_WA * abs(A[2,2])
  A[1,3] <- params$p_EA * abs(A[3,3])
  A[1,4] <- params$p_NA * abs(A[4,4])
  
  A[2,1] <- params$p_AW * abs(A[1,1])
  A[2,3] <- params$p_EW * abs(A[3,3])
  A[2,4] <- params$p_NW * abs(A[4,4])
  
  A[3,1] <- params$p_AE * abs(A[1,1])
  A[3,2] <- params$p_WE * abs(A[2,2])
  A[3,4] <- params$p_NE * abs(A[4,4])
  
  A[4,1] <- params$p_AN * abs(A[1,1])
  A[4,2] <- params$p_WN * abs(A[2,2])
  A[4,3] <- params$p_EN * abs(A[3,3])
  
  for(i in 1:4){
    A[5,i] <- params$p_H * abs(A[i,i])
  }
  
  for(i in 1:4){
    A[i,i] <- A[i,i] + leac * precip/1000
  }
  
  A
}

yasso15_annual_step <- function(pools, input_awen, params, climate, d = 0, leac = 0, steadystate = FALSE){
  
  A <- yasso15_build_matrix(params, climate, d, leac)
  b <- input_awen
  
  if(steadystate){
    xt <- -solve(A, b)
    resp <- sum(b) - sum(xt)
    return(list(pools = pmax(0, xt), respiration = max(0, resp)))
  }
  
  expA <- as.matrix(Matrix::expm(A))
  
  if(max(abs(A)) < 1e-12){
    new_pools <- pools + b
  } else {
    Ainv <- solve(A)
    new_pools <- as.vector(expA %*% pools + (expA - diag(5)) %*% Ainv %*% b)
  }
  
  resp <- sum(pools) + sum(b) - sum(new_pools)
  
  list(pools = pmax(0,new_pools), respiration = max(0,resp))
}

yasso15_steady_state <- function(params, input_awen, climate, d = 0){
  A <- yasso15_build_matrix(params, climate, d)
  -solve(A, input_awen)
}

yasso15_run <- function(params = YASSO15_DEFAULT_PARAMS, C0 = NULL, input_df, spinup_years = 5000){
  
  n_years <- nrow(input_df)
  
  wood_sizes <- list(nwl=0, fwl=2, cwl=15)
  litter_types <- c("nwl","fwl","cwl")
  
  pools_mat <- matrix(0, n_years, 5, dimnames=list(NULL,c("A","W","E","N","H")))
  respiration <- numeric(n_years)
  
  for(lt in litter_types){
    
    if(!paste0(lt,"_A") %in% names(input_df)) next
    
    input_A <- input_df[[paste0(lt,"_A")]]
    input_W <- input_df[[paste0(lt,"_W")]]
    input_E <- input_df[[paste0(lt,"_E")]]
    input_N <- input_df[[paste0(lt,"_N")]]
    
    if(is.null(C0) || !lt %in% names(C0)){
      avg_input <- c(mean(input_A),mean(input_W),mean(input_E),mean(input_N),0)
      climate <- c(mean(input_df$temp_mean),
                   mean(input_df$precip),
                   mean(input_df$temp_amplitude))
      C0_lt <- yasso15_steady_state(params, avg_input, climate, wood_sizes[[lt]])
    } else {
      C0_lt <- C0[[lt]]
    }
    
    pools_lt <- C0_lt
    lt_mat <- matrix(0,n_years,5)
    lt_mat[1,] <- pools_lt
    
    for(t in 1:(n_years-1)){
      input_awen <- c(input_A[t],input_W[t],input_E[t],input_N[t],0)
      climate <- c(input_df$temp_mean[t],
                   input_df$precip[t],
                   input_df$temp_amplitude[t])
      
      res <- yasso15_annual_step(pools_lt, input_awen, params, climate,
                                 d = wood_sizes[[lt]])
      pools_lt <- res$pools
      lt_mat[t+1,] <- pools_lt
      respiration[t+1] <- respiration[t+1] + res$respiration
    }
    
    pools_mat <- pools_mat + lt_mat
  }
  
  list(
    pools = pools_mat,
    total_soc = rowSums(pools_mat),
    respiration = respiration,
    pool_names = c("A","W","E","N","H")
  )
}
