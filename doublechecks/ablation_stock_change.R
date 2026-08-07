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

# =============================================================================
# STOCK CHANGE by C5 setting   (2026-08-07)
#
# WHY: the C5 withdrawal was adjudicated on trusted-campaign RMSE, which is a
# STOCK metric. sigma_init does not move the stock much -- it moves the pre-run
# DIRECTION, and therefore the stock CHANGE, which is the inventory-relevant
# quantity and the one C5 turned out to control. This script re-runs the same
# forward machinery and reports paired per-plot stock change instead.
#
# Usage:  Rscript doublechecks/ablation_stock_change.R [MODEL]
# =============================================================================

bucket <- function(y) ifelse(y <= 1990, "c1985", ifelse(y <= 2015, "c2006", "c2024"))

cat(sprintf("\n=====================================================================\n"))
cat(sprintf("PAIRED STOCK CHANGE  ::  %s   (tC/ha/yr)\n", MODEL))
cat("Only plots observed in BOTH campaigns, so composition cannot drive it.\n")
cat(sprintf("Pre-run inverts (1917->1985 declines) when sigma_init > ~0.818.\n"))
cat(sprintf("=====================================================================\n"))
cat(sprintf("\n%-16s %8s %8s | %6s %8s %8s | %6s %8s %8s\n",
            "config","sig_init","sig_inp","n85_24","obs","pred","n06_24","obs","pred"))

for (i in seq_len(nrow(idx))) {
  nm <- idx$config[i]; p <- load_post(idx$run_id[i])
  if (is.null(p)) { cat(sprintf("%-16s (posterior missing)\n", nm)); next }
  mp <- assemble_params(apply(p, 2, median))
  d  <- do.call(rbind, lapply(plots, fwd, mp = mp))
  if (is.null(d) || !nrow(d)) { cat(sprintf("%-16s (forward failed)\n", nm)); next }
  d$cp <- bucket(d$year)

  pair <- function(a, b, yrs) {
    o <- reshape(d[, c("plot_id","cp","obs")],  idvar="plot_id", timevar="cp", direction="wide")
    q <- reshape(d[, c("plot_id","cp","pred")], idvar="plot_id", timevar="cp", direction="wide")
    m <- merge(o, q, by="plot_id")
    oa <- m[[paste0("obs.",a)]];  ob <- m[[paste0("obs.",b)]]
    pa <- m[[paste0("pred.",a)]]; pb <- m[[paste0("pred.",b)]]
    k  <- is.finite(oa)&is.finite(ob)&is.finite(pa)&is.finite(pb)
    if (sum(k) < 20) return(c(n=NA, obs=NA, pred=NA))
    c(n = sum(k), obs = mean(ob[k]-oa[k])/yrs, pred = mean(pb[k]-pa[k])/yrs)
  }
  r1 <- pair("c1985","c2024",39); r2 <- pair("c2006","c2024",18)
  cat(sprintf("%-16s %8.3f %8.3f | %6.0f %+8.3f %+8.3f | %6.0f %+8.3f %+8.3f\n",
      nm, median(p[,"sigma_init"]), median(p[,"sigma_input"]),
      r1["n"], r1["obs"], r1["pred"], r2["n"], r2["obs"], r2["pred"]))
}
cat("\nIf pred flips sign as C5 weakens, the C5 choice controls the reported sink.\n")
