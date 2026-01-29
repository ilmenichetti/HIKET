# -------------------------------
# Default parameters (theta[1:35]) using Global MAP values (DEzs)
# Source: Viskari et al. 2022 supplement, Table 4/5 (Global). :contentReference[oaicite:1]{index=1}
# -------------------------------
YASSO20_DEFAULT_PARAMS <- as.list(setNames(rep(0, 35), paste0("theta_", 1:35)))

# 1-4: alpha_A, alpha_W, alpha_E, alpha_N
YASSO20_DEFAULT_PARAMS$theta_1 <- 0.51
YASSO20_DEFAULT_PARAMS$theta_2 <- 5.19
YASSO20_DEFAULT_PARAMS$theta_3 <- 0.13
YASSO20_DEFAULT_PARAMS$theta_4 <- 0.10

# 5-16: p-terms in A_p (only a subset nonzero in the calibration note)
# Mapping per your A_p structure (column = source pool, row = destination pool):
# theta5 = pWA (W -> A)
# theta6 = pEA (E -> A)  [not listed in the supplement table; set to 0]
# theta7 = pNA (N -> A)  [set to 1 per note]
# theta8 = pAW (A -> W)  [set to 1 per note]
# theta9 = pEW (E -> W)
# theta10 = pNW (N -> W) [0]
# theta11 = pAE (A -> E) [0]
# theta12 = pWE (W -> E) [0]
# theta13 = pNE (N -> E) [0]
# theta14 = pAN (A -> N) [0]
# theta15 = pWN (W -> N)
# theta16 = pEN (E -> N) [0]
YASSO20_DEFAULT_PARAMS$theta_5  <- 0.50  # pWA
YASSO20_DEFAULT_PARAMS$theta_6  <- 0.00  # pEA (unknown in table)
YASSO20_DEFAULT_PARAMS$theta_7  <- 1.00  # pNA
YASSO20_DEFAULT_PARAMS$theta_8  <- 1.00  # pAW
YASSO20_DEFAULT_PARAMS$theta_9  <- 0.99  # pEW
YASSO20_DEFAULT_PARAMS$theta_10 <- 0.00
YASSO20_DEFAULT_PARAMS$theta_11 <- 0.00
YASSO20_DEFAULT_PARAMS$theta_12 <- 0.00
YASSO20_DEFAULT_PARAMS$theta_13 <- 0.00
YASSO20_DEFAULT_PARAMS$theta_14 <- 0.00
YASSO20_DEFAULT_PARAMS$theta_15 <- 0.16  # pWN
YASSO20_DEFAULT_PARAMS$theta_16 <- 0.00

# 17-21: leaching parameters (ignored here)
# keep as 0

# 22-23: beta_1, beta_2 (AWE)
YASSO20_DEFAULT_PARAMS$theta_22 <- 0.16
YASSO20_DEFAULT_PARAMS$theta_23 <- -0.002

# 24-25: beta_N1, beta_N2
YASSO20_DEFAULT_PARAMS$theta_24 <- 0.17
YASSO20_DEFAULT_PARAMS$theta_25 <- -0.005

# 26-27: beta_H1, beta_H2
YASSO20_DEFAULT_PARAMS$theta_26 <- 0.07
YASSO20_DEFAULT_PARAMS$theta_27 <- 0.00

# 28-30: gamma, gamma_N, gamma_H
YASSO20_DEFAULT_PARAMS$theta_28 <- -1.44
YASSO20_DEFAULT_PARAMS$theta_29 <- -2.00
YASSO20_DEFAULT_PARAMS$theta_30 <- -6.90

# 31-32: p_H, alpha_H (note: order in your code is theta[31] = p_H, theta[32] = alpha_H)
YASSO20_DEFAULT_PARAMS$theta_31 <- 0.004
YASSO20_DEFAULT_PARAMS$theta_32 <- 0.0015

# 33-35: woody size params phi1, phi2, r
YASSO20_DEFAULT_PARAMS$theta_33 <- -2.55
YASSO20_DEFAULT_PARAMS$theta_34 <- 1.24
YASSO20_DEFAULT_PARAMS$theta_35 <- 0.25
# -------------------------------
# Climate modifiers
# -------------------------------

yasso20_temp_mod <- function(monthly_temp, theta) {
  tem  <- mean(exp(theta[22]*monthly_temp + theta[23]*monthly_temp^2))
  temN <- mean(exp(theta[24]*monthly_temp + theta[25]*monthly_temp^2))
  temH <- mean(exp(theta[26]*monthly_temp + theta[27]*monthly_temp^2))
  list(tem = tem, temN = temN, temH = temH)
}

yasso20_precip_mod <- function(P, theta, tem, temN, temH) {
  tem  <- tem  * (1 - exp(theta[28]*P/1000))
  temN <- temN * (1 - exp(theta[29]*P/1000))
  temH <- temH * (1 - exp(theta[30]*P/1000))
  list(tem = tem, temN = temN, temH = temH)
}

yasso20_size_mod <- function(d, theta) {
  if(is.null(d) || d <= 0) return(1)
  min(1, (1 + theta[33]*d + theta[34]*d^2)^(-abs(theta[35])))
}

# -------------------------------
# Build system matrix
# -------------------------------

yasso20_build_matrix <- function(theta, tem, temN, temH, size_dep, leac = 0, precip = 0) {
  
  alpha <- abs(c(theta[1], theta[2], theta[3], theta[4], theta[32]))
  
  A_p <- matrix(c(
    -1,        theta[5],  theta[6],  theta[7],  0,
    theta[8],  -1,        theta[9],  theta[10], 0,
    theta[11], theta[12], -1,        theta[13], 0,
    theta[14], theta[15], theta[16], -1,        0,
    theta[31], theta[31], theta[31], theta[31], -1
  ), 5, 5, byrow = TRUE)
  
  k <- diag(c(
    tem  * alpha[1] * size_dep,
    tem  * alpha[2] * size_dep,
    tem  * alpha[3] * size_dep,
    temN * alpha[4] * size_dep,
    temH * alpha[5]               # no size effect in humus
  ))
  
  A <- A_p %*% k
  
  # Leaching (no humus)
  for(i in 1:4){
    A[i,i] <- A[i,i] + leac * precip / 1000
  }
  
  A
}

# -------------------------------
# Annual step
# -------------------------------

yasso20_annual_step <- function(pools, input_awen, params, monthly_temp, precip, diameter = 0, leac = 0, ss = FALSE) {
  
  theta <- unlist(params)
  
  tmods <- yasso20_temp_mod(monthly_temp, theta)
  pmods <- yasso20_precip_mod(precip, theta, tmods$tem, tmods$temN, tmods$temH)
  size_dep <- yasso20_size_mod(diameter, theta)
  
  A <- yasso20_build_matrix(theta, pmods$tem, pmods$temN, pmods$temH, size_dep, leac, precip)
  
  b <- input_awen
  
  if(pmods$tem <= 1e-16){
    new_pools <- pools + b
    resp <- 0
    return(list(pools = new_pools, respiration = resp))
  }
  
  if(ss){
    new_pools <- as.vector(solve(-A, b))
  } else {
    expA <- as.matrix(Matrix::expm(A))
    Ainv <- solve(A)
    new_pools <- as.vector(expA %*% pools + (expA - diag(5)) %*% Ainv %*% b)
  }
  
  resp <- sum(pools) + sum(b) - sum(new_pools)
  
  list(
    pools = pmax(0, new_pools),
    respiration = max(0, resp)
  )
}

# -------------------------------
# Steady state
# -------------------------------

yasso20_steady_state <- function(params, input_awen, monthly_temp, precip, diameter = 0, leac = 0) {
  theta <- unlist(params)
  
  tmods <- yasso20_temp_mod(monthly_temp, theta)
  pmods <- yasso20_precip_mod(precip, theta, tmods$tem, tmods$temN, tmods$temH)
  size_dep <- yasso20_size_mod(diameter, theta)
  
  A <- yasso20_build_matrix(theta, pmods$tem, pmods$temN, pmods$temH, size_dep, leac, precip)
  
  as.vector(solve(-A, input_awen))
}

# -------------------------------
# Main run function (multi-litter)
# -------------------------------

yasso20_run <- function(params = YASSO20_DEFAULT_PARAMS, C0 = NULL, input_df, spinup = TRUE) {
  
  if(!is.list(params)) params <- as.list(params)
  
  n_years <- nrow(input_df)
  
  wood_sizes <- list(
    nwl = 0,
    fwl = 2,
    cwl = 15
  )
  
  pools_mat <- matrix(0, n_years, 5, dimnames = list(NULL, c("A","W","E","N","H")))
  respiration <- numeric(n_years)
  
  litter_types <- c("nwl","fwl","cwl")
  
  for(lt in litter_types){
    
    lt_a_col <- paste0(lt,"_A")
    if(!lt_a_col %in% names(input_df)) next
    
    input_A <- input_df[[paste0(lt,"_A")]]
    input_W <- input_df[[paste0(lt,"_W")]]
    input_E <- input_df[[paste0(lt,"_E")]]
    input_N <- input_df[[paste0(lt,"_N")]]
    
    monthly_temp <- as.matrix(input_df[, paste0("T_",1:12)])
    precip <- input_df$precip
    
    if(is.null(C0) || !lt %in% names(C0)){
      avg_input <- c(mean(input_A), mean(input_W), mean(input_E), mean(input_N), 0)
      C0_lt <- yasso20_steady_state(
        params,
        avg_input,
        colMeans(monthly_temp),
        mean(precip),
        diameter = wood_sizes[[lt]],
        leac = 0
      )
    } else {
      C0_lt <- C0[[lt]]
    }
    
    pools_lt <- C0_lt
    lt_pools_mat <- matrix(0, n_years, 5)
    lt_pools_mat[1,] <- pools_lt
    
    for(t in 1:(n_years-1)){
      
      input_awen <- c(input_A[t], input_W[t], input_E[t], input_N[t], 0)
      
      res <- yasso20_annual_step(
        pools_lt,
        input_awen,
        params,
        monthly_temp[t,],
        precip[t],
        diameter = wood_sizes[[lt]],
        leac = 0,
        ss = FALSE
      )
      
      pools_lt <- res$pools
      lt_pools_mat[t+1,] <- pools_lt
      respiration[t+1] <- respiration[t+1] + res$respiration
    }
    
    pools_mat <- pools_mat + lt_pools_mat
  }
  
  total_soc <- rowSums(pools_mat)
  
  list(
    pools = pools_mat,
    total_soc = total_soc,
    respiration = respiration,
    pool_names = c("A","W","E","N","H")
  )
}
