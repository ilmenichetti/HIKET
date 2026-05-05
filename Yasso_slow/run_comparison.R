# =============================================================================
# Calibrate and compare: thinning + compressed rotation scenario
# =============================================================================

library(Matrix)
source("yasso07_intact.R")

dat = readRDS("yasso07_synthetic_data.rds")

# -----------------------------------------------------------------
# 1. Unpack
# -----------------------------------------------------------------
climate    = dat$climate
litter_nw  = dat$litter_nw
litter_fw  = dat$litter_fw
litter_cw  = dat$litter_cw
litter_st  = dat$litter_st
obs        = dat$observations
PA_base    = dat$parameters
awen_nw    = dat$awen_fractions$nonwoody
awen_fw    = dat$awen_fractions$fine_woody
awen_cw    = dat$awen_fractions$coarse_woody
awen_st    = dat$awen_fractions$stumps
d_nw       = dat$diameters["nonwoody"]
d_fw       = dat$diameters["fine_woody"]
d_cw       = dat$diameters["coarse_woody"]
d_st       = dat$diameters["stumps"]
n_total    = dat$scenario$n_total
cc1        = dat$scenario$cc1_year
cc2        = dat$scenario$cc2_year
cc3        = dat$scenario$cc3_year
obs_end    = dat$scenario$obs_end
FRAG_RATE  = dat$frag_rate

# Event years for plotting
events_cc   = c(cc1, cc2, cc3)
events_thin = c(dat$scenario$thin1_year, dat$scenario$thin2_year,
                dat$scenario$thin3_year, dat$scenario$thin4_year)

cat("Scenario:", n_total, "yr\n")
cat("CCs:", events_cc, "\n")
cat("Thinnings:", events_thin, "\n")
cat("FragRate:", round(FRAG_RATE, 4), "(half-life", round(log(2)/FRAG_RATE, 1), "yr)\n\n")

# -----------------------------------------------------------------
# 2. Forward run function
# -----------------------------------------------------------------
run_yasso = function(par_vec, PA, use_intact = FALSE) {
  
  k_mult = par_vec[1]; h_mult = par_vec[2]
  PA_mod = PA
  PA_mod[1:4] = PA[1:4] * k_mult
  PA_mod[35]  = PA[35]  * h_mult
  
  mean_MAT = mean(climate$MAT)
  mean_TA  = mean(climate$TA)
  mean_PR  = mean(climate$PR)
  
  rate_nw_eq = 2.0; rate_fw_eq = 0.8; rate_cw_eq = 0.5; rate_st_eq = 0.1
  
  n_state = if (use_intact) 6 else 5
  init = rep(0, n_state)
  
  p_nw = as.numeric(yasso07.intact(mean_MAT, mean_TA, mean_PR, init,
                                   c(rate_nw_eq / sum(awen_nw) * awen_nw, 0), d_nw, PA_mod, 3000,
                                   UseIntactPool = use_intact, FragRate = FRAG_RATE))
  p_fw = as.numeric(yasso07.intact(mean_MAT, mean_TA, mean_PR, init,
                                   c(rate_fw_eq / sum(awen_fw) * awen_fw, 0), d_fw, PA_mod, 3000,
                                   UseIntactPool = use_intact, FragRate = FRAG_RATE))
  p_cw = as.numeric(yasso07.intact(mean_MAT, mean_TA, mean_PR, init,
                                   c(rate_cw_eq / sum(awen_cw) * awen_cw, 0), d_cw, PA_mod, 3000,
                                   UseIntactPool = use_intact, FragRate = FRAG_RATE))
  p_st = as.numeric(yasso07.intact(mean_MAT, mean_TA, mean_PR, init,
                                   c(rate_st_eq / sum(awen_st) * awen_st, 0), d_st, PA_mod, 3000,
                                   UseIntactPool = use_intact, FragRate = FRAG_RATE))
  
  soc = numeric(n_total + 1)
  soc[1] = sum(p_nw) + sum(p_fw) + sum(p_cw) + sum(p_st)
  
  for (y in 1:n_total) {
    p_nw = as.numeric(yasso07.intact(climate$MAT[y], climate$TA[y], climate$PR[y],
                                     p_nw, litter_nw[y, ], d_nw, PA_mod, 1,
                                     UseIntactPool = use_intact, FragRate = FRAG_RATE))
    p_fw = as.numeric(yasso07.intact(climate$MAT[y], climate$TA[y], climate$PR[y],
                                     p_fw, litter_fw[y, ], d_fw, PA_mod, 1,
                                     UseIntactPool = use_intact, FragRate = FRAG_RATE))
    p_cw = as.numeric(yasso07.intact(climate$MAT[y], climate$TA[y], climate$PR[y],
                                     p_cw, litter_cw[y, ], d_cw, PA_mod, 1,
                                     UseIntactPool = use_intact, FragRate = FRAG_RATE))
    p_st = as.numeric(yasso07.intact(climate$MAT[y], climate$TA[y], climate$PR[y],
                                     p_st, litter_st[y, ], d_st, PA_mod, 1,
                                     UseIntactPool = use_intact, FragRate = FRAG_RATE))
    soc[y + 1] = sum(p_nw) + sum(p_fw) + sum(p_cw) + sum(p_st)
  }
  soc
}

# -----------------------------------------------------------------
# 3. Calibrate
# -----------------------------------------------------------------
obj_fun = function(par_vec, PA, use_intact, obs) {
  soc = tryCatch(run_yasso(par_vec, PA, use_intact),
                 error = function(e) rep(NA, n_total + 1))
  if (any(is.na(soc)) || any(soc < 0)) return(1e10)
  model_at_obs = soc[obs$year + 1]
  sum(((model_at_obs - obs$SOC_obs) / obs$SOC_obs_sd)^2)
}

cat("=== Calibrating standard Yasso07 ===\n")
fit_std = optim(par = c(1.0, 1.0), fn = obj_fun,
                PA = PA_base, use_intact = FALSE, obs = obs,
                method = "L-BFGS-B", lower = c(0.1, 0.1), upper = c(5.0, 5.0))
cat("k_mult:", round(fit_std$par[1], 4), "  h_mult:", round(fit_std$par[2], 4),
    "  obj:", round(fit_std$value, 2), "\n\n")

cat("=== Calibrating intact-pool Yasso07 ===\n")
fit_int = optim(par = c(1.0, 1.0), fn = obj_fun,
                PA = PA_base, use_intact = TRUE, obs = obs,
                method = "L-BFGS-B", lower = c(0.1, 0.1), upper = c(5.0, 5.0))
cat("k_mult:", round(fit_int$par[1], 4), "  h_mult:", round(fit_int$par[2], 4),
    "  obj:", round(fit_int$value, 2), "\n\n")

# -----------------------------------------------------------------
# 4. Predictions
# -----------------------------------------------------------------
soc_std = run_yasso(fit_std$par, PA_base, use_intact = FALSE)
soc_int = run_yasso(fit_int$par, PA_base, use_intact = TRUE)
soc_true = dat$soc_true_intact$total
years = 0:n_total

# -----------------------------------------------------------------
# 5. Plot (4 panels)
# -----------------------------------------------------------------
# Helper: add event markers
add_events = function(ylim) {
  abline(v = events_cc, lty = 2, col = "gray40")
  abline(v = events_thin, lty = 3, col = "gray60")
  text(cc1 + 0.5, ylim[2], "CC1", adj = c(0, 1), cex = 0.6, col = "gray40")
  text(cc2 + 0.5, ylim[2], "CC2", adj = c(0, 1), cex = 0.6, col = "gray40")
  text(cc3 + 0.5, ylim[2], "CC3", adj = c(0, 1), cex = 0.6, col = "gray40")
}

png("yasso07_intact_comparison.png", width = 11, height = 8, units = "in", res = 300)
par(mfrow = c(2, 2), mar = c(4, 4.5, 2, 1))

# Panel 1: SOC
ylim = range(c(soc_std, soc_int, soc_true,
               obs$SOC_obs + obs$SOC_obs_sd,
               obs$SOC_obs - obs$SOC_obs_sd), na.rm = TRUE)
plot(years, soc_std, type = "l", lwd = 2, col = "steelblue",
     ylim = ylim, xlab = "Year", ylab = "SOC [Mg C/ha]",
     main = "Standard vs intact-pool Yasso07 (truth = intact model)")
lines(years, soc_int, lwd = 2, col = "firebrick")
lines(years, soc_true, lwd = 1.5, lty = 3, col = "gray40")
points(obs$year, obs$SOC_obs, pch = 16, cex = 0.7)
arrows(obs$year, obs$SOC_obs - obs$SOC_obs_sd,
       obs$year, obs$SOC_obs + obs$SOC_obs_sd,
       angle = 90, code = 3, length = 0.02)
rect(obs_end, ylim[1], n_total, ylim[2], col = rgb(0,0,0,0.04), border = NA)
add_events(ylim)
legend("bottomleft",
       c("Standard (calibrated)", "Intact pool (calibrated)",
         "Truth (intact, uncalib.)", "Observations"),
       col = c("steelblue", "firebrick", "gray40", "black"),
       lty = c(1, 1, 3, NA), lwd = c(2, 2, 1.5, NA),
       pch = c(NA, NA, NA, 16), bty = "n", cex = 0.6)

# Panel 2: Residuals
res_std = soc_std - soc_true
res_int = soc_int - soc_true
ylim2 = range(c(res_std, res_int))
plot(years, res_std, type = "l", lwd = 2, col = "steelblue",
     ylim = ylim2, xlab = "Year", ylab = "Model - Truth [Mg C/ha]",
     main = "Residuals relative to true SOC")
lines(years, res_int, lwd = 2, col = "firebrick")
abline(h = 0, col = "gray70")
rect(obs_end, ylim2[1], n_total, ylim2[2], col = rgb(0,0,0,0.04), border = NA)
add_events(ylim2)
legend("bottomleft", c("Standard", "Intact pool"),
       col = c("steelblue", "firebrick"), lty = 1, lwd = 2, bty = "n", cex = 0.7)

# Panel 3: Difference
diff_soc = soc_int - soc_std
ylim3 = range(diff_soc)
plot(years, diff_soc, type = "l", lwd = 2, col = "darkorange",
     ylim = ylim3, xlab = "Year", ylab = "Intact - Standard [Mg C/ha]",
     main = "Difference between the two models")
abline(h = 0, col = "gray70")
rect(obs_end, ylim3[1], n_total, ylim3[2], col = rgb(0,0,0,0.04), border = NA)
add_events(ylim3)

# Panel 4: Litter inputs
cols = c("forestgreen", "steelblue", "darkorange", "firebrick")
ylim4 = c(0, max(cbind(dat$litter_total$rate_nw, dat$litter_total$rate_fw,
                       dat$litter_total$rate_cw, dat$litter_total$rate_st)) * 1.05)
matplot(dat$litter_total$year,
        cbind(dat$litter_total$rate_nw, dat$litter_total$rate_fw,
              dat$litter_total$rate_cw, dat$litter_total$rate_st),
        type = "l", lty = 1, lwd = 1.5, col = cols, ylim = ylim4,
        xlab = "Year", ylab = "Litter [Mg C/ha/yr]",
        main = "Litter inputs by diameter class")
add_events(ylim4)
legend("topright",
       paste0(c("Non-woody", "Fine woody", "Coarse woody", "Stumps"),
              " (d=", c(d_nw, d_fw, d_cw, d_st), ")"),
       col = cols, lty = 1, lwd = 1.5, bty = "n", cex = 0.55)
dev.off()

cat("Saved: yasso07_intact_comparison.png\n")

# -----------------------------------------------------------------
# 6. Single pulse
# -----------------------------------------------------------------
pulse_years = 60
mean_MAT = mean(climate$MAT)
mean_TA  = mean(climate$TA)
mean_PR  = mean(climate$PR)
pulse_awenh = c(awen_st, 0)
zero_input  = c(0, 0, 0, 0, 0)

p_std = pulse_awenh
pulse_std = numeric(pulse_years + 1); pulse_std[1] = sum(p_std)
for (y in 1:pulse_years) {
  p_std = as.numeric(yasso07.intact(mean_MAT, mean_TA, mean_PR,
                                    p_std, zero_input, d_st, PA_base, 1, UseIntactPool = FALSE))
  pulse_std[y + 1] = sum(p_std)
}

d0 = 11.4; p_shape = 1.22
f_intact = 1 - exp(-(d_st / d0)^p_shape)

p_int = c(pulse_awenh * (1 - f_intact), sum(pulse_awenh[1:4]) * f_intact)
pulse_int = numeric(pulse_years + 1); pulse_unfrag = numeric(pulse_years + 1)
pulse_int[1] = sum(p_int); pulse_unfrag[1] = p_int[6]
for (y in 1:pulse_years) {
  p_int = as.numeric(yasso07.intact(mean_MAT, mean_TA, mean_PR,
                                    p_int, zero_input, d_st, PA_base, 1, UseIntactPool = TRUE, FragRate = FRAG_RATE))
  pulse_int[y + 1] = sum(p_int); pulse_unfrag[y + 1] = p_int[6]
}

png("yasso07_intact_single_stump_pulse.png", width = 5.5, height = 3.5, units = "in", res = 300)
par(mar = c(4, 4.5, 2, 1))
plot(0:pulse_years, pulse_std, type = "l", lwd = 2.5, col = "steelblue",
     ylim = c(0, 1.05), xlab = "Years since pulse",
     ylab = "C remaining [Mg C/ha]",
     main = paste0("Single stump pulse (d = ", d_st, " cm)"))
lines(0:pulse_years, pulse_int, lwd = 2.5, col = "firebrick")
lines(0:pulse_years, pulse_unfrag, lwd = 1.5, col = "firebrick", lty = 2)
legend("topright",
       c("Standard Yasso07", "Intact pool (total C)", "Unfragmented wood"),
       col = c("steelblue", "firebrick", "firebrick"),
       lty = c(1, 1, 2), lwd = c(2.5, 2.5, 1.5), bty = "n", cex = 0.8)
dev.off()

cat("Saved: yasso07_intact_single_stump_pulse.png\n")

# -----------------------------------------------------------------
# 7. Summary
# -----------------------------------------------------------------
obs_idx  = obs$year + 1
pred_idx = (obs_end + 1):(n_total + 1)

cat("\n=== RMSE vs truth ===\n")
cat("Calibration period (yr 0-", obs_end, "):\n")
cat("  Standard:", round(sqrt(mean((soc_std[obs_idx] - soc_true[obs_idx])^2)), 2), "Mg C/ha\n")
cat("  Intact:  ", round(sqrt(mean((soc_int[obs_idx] - soc_true[obs_idx])^2)), 2), "Mg C/ha\n")
cat("Prediction period (yr", obs_end, "-", n_total, "):\n")
cat("  Standard:", round(sqrt(mean((soc_std[pred_idx] - soc_true[pred_idx])^2)), 2), "Mg C/ha\n")
cat("  Intact:  ", round(sqrt(mean((soc_int[pred_idx] - soc_true[pred_idx])^2)), 2), "Mg C/ha\n")