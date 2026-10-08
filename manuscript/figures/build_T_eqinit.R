# =============================================================================
# build_T_eqinit.R   (2026-10-08)
#
# T_eqinit -- the two arms side by side, one row per model x arm:
#   fit       data log-likelihood (max and median over draws; Jacobian removed --
#             eqinit_draws.R) and dLL = eq - prod at the max. The arms are NESTED,
#             so dLL <= 0 up to sampling noise; one parameter fewer is cheap, but
#             the likelihood inherits the documented overconfidence (n_eff ~221),
#             so read dLL together with the rates, not alone.
#             Calibration / holdout RMSE and bias from the predictive stage.
#   params    sigma_input, sigma_init (= 1 in the eq arm), intrinsic transit time
#   rates     1985-2006, 2006-2024, 1985-2024 (median over draws), observed alongside
#   forecast  sink 2025-2044 and headroom (build_forward_scenarios.R definitions)
# Writes manuscript/figures/T_eqinit.csv and prints it.
#
# Usage:  Rscript manuscript/figures/build_T_eqinit.R   (after eqinit_draws.R)
# =============================================================================

setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
source("manuscript/figures/eqinit_common.R")
z <- load_eqinit_comparison()
D <- if (file.exists("manuscript/figures/eqinit_draws.rds")) readRDS("manuscript/figures/eqinit_draws.rds") else NULL
if (!is.null(D) && !identical(D$stamp, z$stamp)) {
  warning("eqinit_draws.rds is from other RUN_IDs -- fit and transit-time columns left empty"); D <- NULL }

rows <- do.call(rbind, lapply(z$models, function(m) do.call(rbind, lapply(c("prod", "eq"), function(a) {
  A <- z$arms[[m]][[a]]; dd <- if (!is.null(D)) D$d[[m]][[a]] else NULL
  r <- apply(A$rates, 2, median)
  data.frame(model = m, arm = a, n_free = if (!is.null(dd)) dd$n_free else NA,
             ll_max = if (!is.null(dd)) max(dd$ll) else NA, ll_median = if (!is.null(dd)) median(dd$ll) else NA,
             rmse_calib = A$metrics$calib$RMSE_mean, rmse_holdout = A$metrics$holdout$RMSE_mean,
             bias_calib = A$metrics$calib$bias_mean,
             sigma_input = A$pars["sigma_input", "med"], sigma_init = A$pars["sigma_init", "med"],
             mtt = if (!is.null(dd)) median(dd$thin$mtt, na.rm = TRUE) else NA,
             rate_85_06 = r[["1985-2006"]], rate_06_24 = r[["2006-2024"]], rate_85_24 = r[["1985-2024"]],
             sink_25_44 = median(A$fore$sink_first20), headroom_pct = 100 * median(A$fore$headroom))
}))))
rows$dLL <- ave(rows$ll_max, rows$model, FUN = function(x) x - x[1])
ob <- setNames(z$obs$rates$obs, z$obs$rates$interval)
rows <- rbind(rows, data.frame(model = "observed", arm = "", n_free = NA, ll_max = NA, ll_median = NA,
                               rmse_calib = NA, rmse_holdout = NA, bias_calib = NA, sigma_input = NA,
                               sigma_init = NA, mtt = NA, rate_85_06 = ob[["1985-2006"]],
                               rate_06_24 = ob[["2006-2024"]], rate_85_24 = ob[["1985-2024"]],
                               sink_25_44 = NA, headroom_pct = NA, dLL = NA))
write.csv(rows, "manuscript/figures/T_eqinit.csv", row.names = FALSE)
num <- vapply(rows, is.numeric, logical(1)); rows[num] <- lapply(rows[num], round, 3)
cat(sprintf("\nT_eqinit -- prod %s | eq %s%s\n\n", paste(unique(z$stamp$prod), collapse = "/"),
            paste(unique(z$stamp$eq), collapse = "/"),
            if (length(z$missing)) paste0(" | missing: ", paste(z$missing, collapse = ", ")) else ""))
print(rows, row.names = FALSE)
cat("\nwrote manuscript/figures/T_eqinit.csv\n")
