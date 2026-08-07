# =============================================================================
# production_fit_by_campaign.R   (2026-08-07)
#
# Per-campaign fit of the CURRENT production run, straight from the residuals
# CSVs written by run_*_transient_predictive.R. No forward re-runs needed.
#
# Reports BOTH R-squared definitions, because they disagree sharply here and
# quoting one against the other is a documented trap (NEXT_SESSION.md S2):
#
#   R2_cor = cor(obs, hat)^2                      -- ignores bias; what the
#                                                    predictive scripts report
#   R2_var = 1 - SS_res/SS_tot (variance explained) -- penalises bias; goes
#                                                    negative when a model is
#                                                    systematically off
#
# Usage:  Rscript doublechecks/production_fit_by_campaign.R
# =============================================================================

DIR_DIAG <- "Calibration_real_data_transient/diagnostics"
MODELS   <- c("SP1", "TP2", "TP3", "Yasso07", "Yasso15", "Yasso20")

# Campaign windows (VMI8 1985-86, Biosoil 2006, Komeetta 2024)
campaign_of <- function(y) ifelse(y <= 1990, "1985", ifelse(y <= 2015, "2006", "2024"))

newest <- function(model) {
  fs <- list.files(file.path(DIR_DIAG, model),
                   pattern = sprintf("^%s_residuals_[0-9]{8}_[0-9]{6}\\.csv$", model),
                   full.names = TRUE)
  if (!length(fs)) return(NA_character_)
  sort(fs, decreasing = TRUE)[1]
}

stat_block <- function(d) {
  obs <- d$soc_obs_tCha; hat <- d$soc_mean
  keep <- is.finite(obs) & is.finite(hat)
  obs <- obs[keep]; hat <- hat[keep]
  if (length(obs) < 3) return(NULL)
  data.frame(
    n        = length(obs),
    obs_med  = median(obs),
    hat_med  = median(hat),
    bias     = mean(hat - obs),
    RMSE     = sqrt(mean((hat - obs)^2)),
    R2_cor   = suppressWarnings(cor(obs, hat))^2,
    R2_var   = 1 - sum((obs - hat)^2) / sum((obs - mean(obs))^2))
}

out <- do.call(rbind, lapply(MODELS, function(m) {
  f <- newest(m)
  if (is.na(f)) return(NULL)
  d <- read.csv(f, stringsAsFactors = FALSE)
  d$campaign <- campaign_of(d$year)
  do.call(rbind, lapply(c("1985", "2006", "2024", "ALL"), function(cp) {
    sub <- if (cp == "ALL") d else d[d$campaign == cp, ]
    s <- stat_block(sub)
    if (is.null(s)) return(NULL)
    cbind(data.frame(model = m, campaign = cp), s)
  }))
}))

out[] <- lapply(out, function(x) if (is.numeric(x)) round(x, 3) else x)

cat("\n=== Per-campaign fit, current production run (all plots) ===\n\n")
for (cp in c("1985", "2006", "2024", "ALL")) {
  cat("---- campaign", cp, "----\n")
  print(out[out$campaign == cp, setdiff(names(out), "campaign")], row.names = FALSE)
  cat("\n")
}

cat("bias = mean(pred - obs) tC/ha.  R2_cor ignores bias; R2_var penalises it.\n\n")
