# =============================================================================
# fraction_pin_trajectory.R   (2026-08-10)
#
# Question: if the transfer fractions were tightened, would F3 (mean SOC
# trajectory) and F4 (initialisation) still look reasonable? Specifically, does
# the ACCUMULATION -> PLATEAU shape survive, or does it only exist because the
# fractions ran to the p_AW ~ 0.99 corner?
#
# This runs the LIMIT case of tightening: fractions PINNED at their published
# prior centres, everything else at the posterior median, and sigma_input
# re-solved so the SOC level still matches observed. If the shape holds here it
# holds under any milder prior tightening.
#
# Reports the annual mean SOC trajectory for both configs, so the F3 shape can
# be read directly, plus the 1985 level (the F4 initialisation question) and the
# paired trends.
#
# Usage:  Rscript doublechecks/fraction_pin_trajectory.R [MODEL]
# =============================================================================

suppressWarnings(suppressMessages({
  a <- commandArgs(trailingOnly = TRUE)
  MODEL <- if (length(a) >= 1) a[[1]] else "Yasso15"
  library(BayesianTools)
}))
set.seed(2025)

src   <- readLines(file.path("Calibration_real_data_transient",
                             sprintf("run_%s_transient_calibration.R", MODEL)), warn = FALSE)
cutix <- grep("^t_run <- system.time\\(\\{", src)[1]
e <- new.env(parent = globalenv())
suppressMessages(source(textConnection(paste(src[seq_len(cutix-1L)], collapse="\n")), local = e))

start <- grep("ll_fn <- make_likelihood\\(", src)[1]
open <- 0L; end <- NA_integer_
for (i in seq(start, length(src))) {
  ch <- strsplit(src[i], "")[[1]]
  open <- open + sum(ch=="(") - sum(ch==")")
  if (open == 0L) { end <- i; break }
}
ml <- as.list(str2lang(sub("^\\s*ll_fn\\s*<-\\s*","",
                           paste(src[seq(start,end)], collapse="\n"))))[-1]
argof <- function(nm,d=NULL) if (is.null(ml[[nm]])) d else eval(ml[[nm]], envir=e)
to_original     <- argof("to_original")
assemble_params <- argof("assemble_params")
compute_xi      <- argof("compute_xi")
compute_xi_mean <- argof("compute_xi_mean")
steady_state    <- argof("steady_state")
run_model       <- argof("run_model")
STEADY_N        <- argof("steady_state_n", NULL)

plots <- get("plots", e); climate_by_plot <- get("climate_by_plot", e)
inputs_by_plot <- get("inputs_by_plot", e); litter_means <- get("litter_means", e)
obs_meta <- get("obs_meta", e); best_x <- get("best_x", e)
SOC_obs_all <- get("SOC_obs_all", e)
obs_years <- lapply(plots, function(pid) SOC_obs_all$year[as.character(SOC_obs_all$plot_id)==pid])
names(obs_years) <- plots

# full annual trajectory + SOC at observation years
fwd <- function(pid, mp) {
  clim <- climate_by_plot[[pid]]; inputs <- inputs_by_plot[[pid]]
  lm <- litter_means[[pid]]; meta <- obs_meta[[pid]]
  if (any(is.na(meta$idx))) return(NULL)
  xa <- tryCatch(compute_xi(clim, mp), error=function(z) NULL); if (is.null(xa)) return(NULL)
  n_ss <- if (is.null(STEADY_N)) nrow(clim) else min(STEADY_N, nrow(clim))
  xs <- tryCatch(compute_xi_mean(clim[seq_len(n_ss),,drop=FALSE], mp), error=function(z) NULL)
  if (is.null(xs)) return(NULL)
  C0 <- tryCatch(steady_state(mp, lm, xs), error=function(z) NULL)
  if (is.null(C0) || any(!is.finite(C0)) || any(C0 < 0)) return(NULL)
  ro <- tryCatch(run_model(inputs, mp, C0, xa), error=function(z) NULL); if (is.null(ro)) return(NULL)
  tr <- ro$total_soc
  yr <- if (!is.null(inputs$year) && length(inputs$year) == length(tr)) inputs$year else NULL
  hat <- tr[meta$idx]
  if (any(!is.finite(hat)) || any(hat <= 0)) return(NULL)
  list(traj = tr, years = yr, hat = hat, obs = meta$obs_years <- meta$soc_obs,
       oyr = obs_years[[pid]])
}

run_all <- function(p) { mp <- assemble_params(p)
  r <- lapply(plots, fwd, mp = mp); names(r) <- plots
  r[!vapply(r, is.null, logical(1))] }

soc_med <- function(r) median(unlist(lapply(r, `[[`, "hat")))
obs_med <- function(r) median(unlist(lapply(r, `[[`, "obs")))
paired  <- function(r, y0, y1, what="hat") {
  v <- vapply(r, function(z){ a<-z[[what]][z$oyr==y0]; b<-z[[what]][z$oyr==y1]
    if (length(a)!=1||length(b)!=1) NA_real_ else (b-a)/(y1-y0)}, numeric(1))
  mean(v, na.rm=TRUE) }
annual <- function(r) {
  yl <- lapply(r, `[[`, "years"); ok <- !vapply(yl, is.null, logical(1))
  if (!any(ok)) return(NULL)
  yrs <- yl[[which(ok)[1]]]
  M <- vapply(r[ok], function(z) if (length(z$traj)==length(yrs)) z$traj else rep(NA_real_,length(yrs)),
              numeric(length(yrs)))
  setNames(rowMeans(M, na.rm=TRUE), yrs)
}

# --- configs ------------------------------------------------------------------
p_def <- to_original(best_x)                      # published centres
rid <- sub(sprintf("^%s_posterior_(.+)\\.rds$", MODEL), "\\1",
           sort(list.files("Calibration_real_data_transient/runs",
                pattern=sprintf("^%s_posterior_[0-9]{8}_[0-9]{6}\\.rds$", MODEL)),
                decreasing=TRUE)[1])
smp <- getSample(readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_%s.rds",
                                 MODEL, rid)))
p_cal <- p_def
for (nm in intersect(names(p_cal), colnames(smp))) p_cal[nm] <- median(smp[, nm])

frac <- grep("^p_", names(p_cal), value = TRUE)
p_pin <- p_cal; p_pin[frac] <- p_def[frac]        # LIMIT case: fractions pinned

cat(sprintf("\n==============  %s  (run %s)  ==============\n", MODEL, rid))
cat("fractions pinned:", paste(frac, collapse=", "), "\n")
cat(sprintf("  calibrated p_AW = %.4f  ->  published p_AW = %.4f\n\n",
            p_cal[[frac[1]]], p_def[[frac[1]]]))

A <- run_all(p_cal)
tgt <- obs_med(A)
B0 <- run_all(p_pin)
p_pin["sigma_input"] <- p_pin["sigma_input"] * tgt / soc_med(B0)   # re-solve level
B  <- run_all(p_pin)

tr_A <- annual(A); tr_B <- annual(B)
cat(sprintf("%-22s %12s %12s\n", "", "CALIBRATED", "FRAC PINNED"))
cat(sprintf("%-22s %12.3f %12.3f\n", "sigma_input", p_cal["sigma_input"], p_pin["sigma_input"]))
cat(sprintf("%-22s %12.2f %12.2f\n", "median SOC (obs years)", soc_med(A), soc_med(B)))
cat(sprintf("%-22s %12.3f %12.3f\n", "paired 1985-2024", paired(A,1985,2024), paired(B,1985,2024)))
cat(sprintf("%-22s %12.3f %12.3f\n", "paired 2006-2024", paired(A,2006,2024), paired(B,2006,2024)))
cat(sprintf("%-22s %12.3f %12.3f\n", "  observed (same plots)",
            paired(A,1985,2024,"obs"), paired(A,2006,2024,"obs")))

if (!is.null(tr_A) && !is.null(tr_B)) {
  yrs <- as.integer(names(tr_A))
  key <- yrs[yrs %in% c(1985,1990,1995,2000,2006,2010,2015,2020,2024)]
  cat("\n--- annual MEAN SOC across plots (F3 shape) ---\n")
  cat(sprintf("%6s %12s %12s\n", "year", "calibrated", "frac pinned"))
  for (y in key) cat(sprintf("%6d %12.2f %12.2f\n", y, tr_A[as.character(y)], tr_B[as.character(y)]))
  cat(sprintf("\npeak year:  calibrated %d (%.2f) | pinned %d (%.2f)\n",
              yrs[which.max(tr_A)], max(tr_A), yrs[which.max(tr_B)], max(tr_B)))
  cat(sprintf("drop from peak to 2024: calibrated %+.2f | pinned %+.2f tC/ha\n",
              tr_A[length(tr_A)]-max(tr_A), tr_B[length(tr_B)]-max(tr_B)))
  dir.create("doublechecks/frac_pin", showWarnings = FALSE)
  write.csv(data.frame(year=yrs, calibrated=as.numeric(tr_A), pinned=as.numeric(tr_B)),
            sprintf("doublechecks/frac_pin/%s_trajectory.csv", MODEL), row.names=FALSE)
}
