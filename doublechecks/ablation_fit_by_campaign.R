# =============================================================================
# ablation_fit_by_campaign.R   (2026-08-04)
#
# THE POINT: adjudicating C5 requires a metric the C5 weighting cannot touch.
#
# C5 multiplies the observation SD for the 1985 campaign. That changes the
# LIKELIHOOD ITSELF -- including its normalising constant -- so log-likelihoods
# from different C5 settings live on different scales and CANNOT be compared to
# decide which config "fits better". Neither can anything else read off the
# likelihood (DIC, WAIC, chain ll spread).
#
# What IS comparable is the model's UNWEIGHTED agreement with the observations,
# computed per campaign after the fact. If down-weighting 1985 is justified,
# then relaxing the pull of a suspect campaign should leave the model fitting the
# TRUSTED campaigns (2006, 2024) at least as well. If instead the 2006/2024 fit
# degrades as C5 strengthens, C5 is discarding signal rather than noise.
#
# So for each ablation config this script:
#   1. loads the posterior, takes the componentwise median parameter vector
#   2. runs the genuine forward model at those parameters (engine machinery,
#      sourced from the real calibration script -- see c3_preinit_shape_ablation.R)
#   3. reports RMSE / bias / R^2 per campaign year, UNWEIGHTED
#
# Read the 2006 and 2024 columns to judge C5; read 1985 to see what it did to the
# campaign it distrusts.
#
# Usage:  Rscript doublechecks/ablation_fit_by_campaign.R [MODEL]
# =============================================================================

args  <- commandArgs(trailingOnly = TRUE)
MODEL <- if (length(args) >= 1) args[[1]] else "TP2"

idx_f <- file.path("doublechecks", "ablation_logs",
                   sprintf("%s_ablation_index.csv", MODEL))
if (!file.exists(idx_f)) stop("No ablation index for ", MODEL)
idx <- read.csv(idx_f, stringsAsFactors = FALSE)

script <- sprintf("Calibration_real_data_transient/run_%s_transient_calibration.R", MODEL)
src    <- readLines(script, warn = FALSE)
cutix  <- grep("^t_run <- system.time\\(\\{", src)[1]
message(sprintf("[%s] sourcing setup ...", MODEL))
e <- new.env(parent = globalenv())
source(textConnection(paste(src[seq_len(cutix - 1L)], collapse = "\n")), local = e)

# --- engine bindings (args may be expressions, e.g. Yasso20's monthly x12) ----
start <- grep("ll_fn <- make_likelihood\\(", src)[1]
open <- 0L; end <- NA_integer_
for (i in seq(start, length(src))) {
  ch <- strsplit(src[i], "")[[1]]
  open <- open + sum(ch == "(") - sum(ch == ")")
  if (open == 0L) { end <- i; break }
}
ml <- as.list(str2lang(sub("^\\s*ll_fn\\s*<-\\s*", "",
                           paste(src[seq(start, end)], collapse = "\n"))))[-1]
argof <- function(nm, d = NULL) if (is.null(ml[[nm]])) d else eval(ml[[nm]], envir = e)

assemble_params <- argof("assemble_params")
compute_xi      <- argof("compute_xi")
compute_xi_mean <- argof("compute_xi_mean")
steady_state    <- argof("steady_state")
run_model       <- argof("run_model")
STEADY_N        <- argof("steady_state_n", NULL)

plots           <- get("plots",           e)
climate_by_plot <- get("climate_by_plot", e)
inputs_by_plot  <- get("inputs_by_plot",  e)
litter_means    <- get("litter_means",    e)
obs_meta        <- get("obs_meta",        e)
SOC_obs_all     <- get("SOC_obs_all",     e)
obs_years <- lapply(plots, function(pid)
  SOC_obs_all$year[as.character(SOC_obs_all$plot_id) == pid])
names(obs_years) <- plots

# --- forward at a PHYSICAL-scale parameter vector ----------------------------
# The stored posterior is already in original (physical) units, so it goes
# straight into assemble_params without to_original().
fwd <- function(pid, mp) {
  clim <- climate_by_plot[[pid]]; inputs <- inputs_by_plot[[pid]]
  lm <- litter_means[[pid]];      meta <- obs_meta[[pid]]
  if (any(is.na(meta$idx))) return(NULL)
  xa <- tryCatch(compute_xi(clim, mp), error = function(z) NULL); if (is.null(xa)) return(NULL)
  n_ss <- if (is.null(STEADY_N)) nrow(clim) else min(STEADY_N, nrow(clim))
  xs <- tryCatch(compute_xi_mean(clim[seq_len(n_ss), , drop = FALSE], mp),
                 error = function(z) NULL); if (is.null(xs)) return(NULL)
  ci <- tryCatch(steady_state(mp, lm, xs), error = function(z) NULL)
  if (is.null(ci) || any(!is.finite(ci)) || any(ci < 0)) return(NULL)
  ro <- tryCatch(run_model(inputs, mp, ci, xa), error = function(z) NULL); if (is.null(ro)) return(NULL)
  sh <- ro$total_soc[meta$idx]
  if (any(!is.finite(sh)) || any(sh <= 0)) return(NULL)
  data.frame(plot_id = pid, year = obs_years[[pid]], obs = meta$soc_obs, pred = sh,
             sigma_init = unname(mp["sigma_init"]))
}

load_post <- function(rid) {
  cand <- c(file.path("doublechecks","ablation_runs", sprintf("%s_posterior_%s.rds", MODEL, rid)),
            file.path("Calibration_real_data_transient","runs", sprintf("%s_posterior_%s.rds", MODEL, rid)))
  f <- cand[file.exists(cand)][1]; if (is.na(f)) NULL else readRDS(f)
}

# TWO R^2 definitions, deliberately both reported -- they answer different questions and
# quoting one against the other has already caused confusion here:
#   r2_cor  = cor(obs,pred)^2  -- what run_*_predictive.R reports (the "R2 0.05-0.11" in the
#             docs). Ignores bias and scale entirely: a model predicting 2x the truth with
#             perfect correlation scores 1.0. It measures pattern, not accuracy.
#   r2_var  = 1 - SSres/SStot  -- variance explained. Penalises the systematic
#             over-prediction hard, so it goes NEGATIVE whenever RMSE exceeds the
#             observations' own SD. Harsher, and computed per campaign here (within-campaign
#             variance is smaller than pooled, harsher again).
# Compare configs on RMSE and bias; use r2_cor only when comparing to older reported numbers.
metrics <- function(d) {
  r <- d$pred - d$obs
  c(n = nrow(d), rmse = sqrt(mean(r^2)), bias = mean(r),
    r2_var = 1 - sum(r^2) / sum((d$obs - mean(d$obs))^2),
    r2_cor = suppressWarnings(cor(d$obs, d$pred, use = "complete.obs")^2))
}

cat(sprintf("\n=====================================================================\n"))
cat(sprintf("UNWEIGHTED FIT BY CAMPAIGN  ::  %s\n", MODEL))
cat(sprintf("at each config's posterior-median parameters. These metrics do NOT\n"))
cat(sprintf("depend on the C5 weighting, so they are comparable ACROSS configs.\n"))
cat(sprintf("Judge C5 on 2006/2024: strengthening it must not degrade them.\n"))
cat(sprintf("=====================================================================\n"))

out <- list()
for (i in seq_len(nrow(idx))) {
  nm <- idx$config[i]; p <- load_post(idx$run_id[i])
  if (is.null(p)) { cat(sprintf("\n%-16s (posterior missing)\n", nm)); next }
  mp <- assemble_params(apply(p, 2, median))
  d  <- do.call(rbind, lapply(plots, fwd, mp = mp))
  cat(sprintf("\n%-16s  sigma_init = %.3f | sigma_input = %.3f\n",
              nm, median(p[, "sigma_init"]), median(p[, "sigma_input"])))
  cat(sprintf("  %-8s %6s %8s %8s %8s %8s\n", "campaign", "n", "RMSE", "bias", "R2var", "R2cor"))
  for (yr in c(1985L, 2006L, 2024L)) {
    dy <- d[d$year == yr, ]; if (!nrow(dy)) next
    m  <- metrics(dy)
    cat(sprintf("  %-8d %6d %8.2f %+8.2f %8.3f %8.3f\n",
                yr, m["n"], m["rmse"], m["bias"], m["r2_var"], m["r2_cor"]))
  }
  dt <- d[d$year %in% c(2006L, 2024L), ]
  mt <- metrics(dt)
  cat(sprintf("  %-8s %6d %8.2f %+8.2f %8.3f %8.3f   <- TRUSTED only\n",
              "2006+24", mt["n"], mt["rmse"], mt["bias"], mt["r2_var"], mt["r2_cor"]))
  # BIAS-CONTROLLED 1985 excess: the LOCO test must compare the held-out campaign against
  # the model's bias on campaigns it DID see. Without this control, any positively-biased
  # model would appear to "prove" that VMI8 reads low.
  d85 <- d[d$year == 1985L, ]
  if (nrow(d85)) {
    ex <- mean(d85$pred - d85$obs) - mt["bias"]
    cat(sprintf("  %-8s %6s %8s %+8.2f %8s %8s   <- 1985 excess OVER general bias\n",
                "excess", "", "", ex, "", ""))
  }
  out[[nm]] <- c(config = nm, trusted_rmse = unname(mt["rmse"]),
                 trusted_r2 = unname(mt["r2_cor"]),
                 excess85 = if (nrow(d85)) unname(mean(d85$pred - d85$obs) - mt["bias"]) else NA)
}

cat(sprintf("\n=====================================================================\n"))
cat("VERDICT AID: trusted-campaign (2006+2024) RMSE by config, lower = better\n")
for (nm in names(out))
  cat(sprintf("  %-16s RMSE %8.2f | R2cor %6.3f | 1985 excess %+7.2f\n", nm,
              as.numeric(out[[nm]]["trusted_rmse"]), as.numeric(out[[nm]]["trusted_r2"]),
              as.numeric(out[[nm]]["excess85"])))
cat("If RMSE is flat across C5 settings, C5 is not doing useful work.\n")
cat("For A5_LOCO: judge on the 1985 EXCESS (bias-controlled), not the raw 1985 bias.\n")
cat("=====================================================================\n")
