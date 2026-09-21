# =============================================================================
# basal_area_sign.R   (2026-09-04)
#
# THE QUESTION. Stand basal area both DRIVES the predictions (F6: the predicted
# cloud is cleanly ordered by BA quintile) and LEADS THE RESIDUALS in five of six
# models (F11). Those two facts together would mean the models do not merely lack
# a driver -- they MIS-SCALE one they already have. This resolves the direction.
#
# WHAT IS COMPUTED, per model, on the same basis as the RF (residual_log =
# log(obs) - log(pred), from each predictive bundle's residuals_df):
#   d log(pred) / d BA     how hard the model leans on basal area
#   d log(obs)  / d BA     how hard the DATA lean on it
#   d resid_log / d BA     the difference; < 0 => the model OVER-applies BA
# Reported for holdout, calibration and all observations, with 95% CIs.
#
# READING. The three slopes satisfy b_resid = b_obs - b_pred exactly (OLS is
# linear), so the test is really "does the model's BA response match the data's".
# =============================================================================
source("manuscript/figures/run_ids.R")
setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")

runs <- "Calibration_real_data_transient/runs"
sl <- function(y, x) {                     # slope, 95% CI, p, r
  ok <- is.finite(y) & is.finite(x); y <- y[ok]; x <- x[ok]
  f <- lm(y ~ x); ci <- confint(f)["x", ]
  c(b = unname(coef(f)["x"]), lo = ci[[1]], hi = ci[[2]],
    p = summary(f)$coefficients["x", 4], r = cor(x, y), n = length(y))
}

out <- list()
for (m in FIG_MODELS) {
  b  <- readRDS(file.path(runs, sprintf("%s_posterior_predictive_%s.rds", m, RID[[m]])))
  rd <- b$residuals_df
  rd <- rd[is.finite(rd$soc_obs_tCha) & is.finite(rd$soc_median) &
             is.finite(rd$basal_area_85) & rd$soc_obs_tCha > 0 & rd$soc_median > 0, ]
  rd$lo_ <- log(rd$soc_obs_tCha); rd$lp_ <- log(rd$soc_median)
  rd$rs_ <- rd$lo_ - rd$lp_                     # == residual_log
  for (set in c("holdout", "calib", "all")) {
    d <- switch(set, holdout = rd[rd$is_holdout %in% TRUE, ],
                     calib   = rd[!(rd$is_holdout %in% TRUE), ], rd)
    if (nrow(d) < 20) next
    out[[length(out) + 1]] <- data.frame(
      model = m, set = set, n = nrow(d),
      b_pred  = sl(d$lp_, d$basal_area_85)[["b"]],
      b_obs   = sl(d$lo_, d$basal_area_85)[["b"]],
      b_resid = sl(d$rs_, d$basal_area_85)[["b"]],
      lo      = sl(d$rs_, d$basal_area_85)[["lo"]],
      hi      = sl(d$rs_, d$basal_area_85)[["hi"]],
      p       = sl(d$rs_, d$basal_area_85)[["p"]],
      r_pred  = sl(d$lp_, d$basal_area_85)[["r"]],
      r_obs   = sl(d$lo_, d$basal_area_85)[["r"]])
  }
  rm(b, rd); gc()
}
res <- do.call(rbind, out)

cat("\n=== d log(SOC) / d basal_area  (per m2/ha) ==========================\n")
cat("b_pred = model's response;  b_obs = data's;  b_resid = b_obs - b_pred\n")
cat("b_resid < 0  =>  the model OVER-applies basal area\n\n")
for (s in c("holdout", "calib", "all")) {
  x <- res[res$set == s, ]
  if (!nrow(x)) next
  cat(sprintf("--- %s (n = %d) ---\n", toupper(s), x$n[1]))
  print(format(data.frame(model = x$model,
    b_pred = round(x$b_pred, 5), b_obs = round(x$b_obs, 5),
    b_resid = round(x$b_resid, 5),
    CI = sprintf("[%.4f, %.4f]", x$lo, x$hi),
    p = signif(x$p, 2), r_pred = round(x$r_pred, 3), r_obs = round(x$r_obs, 3)),
    justify = "right"), row.names = FALSE)
  cat("\n")
}
# --- CONFOUNDING CONTROL ------------------------------------------------------
# The marginal observed slope could be suppressed (or created) by geography: the
# north is sampled at 1/3 density and differs in both BA and SOC. Refit the OBSERVED
# relation with region and soil type as covariates; if the observed slope survives,
# the mismatch below is real, and if it shrinks the mismatch is larger still.
sa <- read.csv("Data/model_inputs/site_attributes.csv")
b  <- readRDS(file.path(runs, sprintf("Yasso15_posterior_predictive_%s.rds", RID[["Yasso15"]])))
rd <- b$residuals_df
rd$region <- factor(sa$region[match(rd$plot_id, sa$plot_id)])
rd$soil   <- factor(sa$soil_type[match(rd$plot_id, sa$plot_id)])
rd <- rd[is.finite(rd$soc_obs_tCha) & is.finite(rd$soc_median) & is.finite(rd$basal_area_85) &
           rd$soc_obs_tCha > 0 & rd$soc_median > 0 & !is.na(rd$region), ]
lo <- log(rd$soc_obs_tCha); lp <- log(rd$soc_median); BA <- rd$basal_area_85
D <- data.frame(lo = lo, lp = lp, BA = BA, region = rd$region, soil = rd$soil)
cat("=== observed BA slope, with controls (Yasso15 rows; obs are model-independent) ===\n")
for (lab in c("BA", "BA + region", "BA + region + soil")) {
  f <- lm(as.formula(paste("lo ~", lab)), data = D)
  cat(sprintf("  obs  ~ %-20s b = %8.5f   p = %.3g\n", lab,
              coef(f)[["BA"]], summary(f)$coefficients["BA", 4]))
}
for (lab in c("BA", "BA + region", "BA + region + soil")) {
  f <- lm(as.formula(paste("lp ~", lab)), data = D)
  cat(sprintf("  pred ~ %-20s b = %8.5f\n", lab, coef(f)[["BA"]]))
}
rng <- diff(range(BA))
cat(sprintf("\nAcross the observed BA range (%.0f m2/ha): models span x%.1f in SOC, data x%.2f\n",
            rng, exp(coef(lm(lp ~ BA, D))[["BA"]] * rng), exp(coef(lm(lo ~ BA, D))[["BA"]] * rng)))
cat(sprintf("Fully adjusted (region + soil): data span x%.2f  =>  model/data ratio %.0f\n",
            exp(coef(lm(lo ~ BA + region + soil, D))[["BA"]] * rng),
            coef(lm(lp ~ BA + region + soil, D))[["BA"]] / coef(lm(lo ~ BA + region + soil, D))[["BA"]]))

saveRDS(res, "doublechecks/basal_area_sign.rds")
cat("saved doublechecks/basal_area_sign.rds\n")
