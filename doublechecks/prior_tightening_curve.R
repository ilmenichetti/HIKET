# =============================================================================
# prior_tightening_curve.R   (2026-08-10)
#
# HOW TIGHT? The pinning test gave the destination (published values) but not
# the trade-off. This traces the path.
#
# lambda = fraction of the current posterior displacement from the published
# centre that is RETAINED, applied in the UNCONSTRAINED space the priors live
# in, to STRUCTURAL parameters only (fractions + climate + woody size):
#
#     x(lambda) = x_published + lambda * (x_calibrated - x_published)
#
#   lambda = 1  -> the run as calibrated
#   lambda = 0  -> pinned at published centres
#
# Auxiliaries (sigma_input, sigma_init) are NEVER pinned -- they carry the
# initialisation storyline. sigma_input is re-solved at every lambda so all
# configs hit the same observed median SOC; without that, a config that
# over-predicts the level also inflates its own trend by the same factor
# (SOC is exactly linear in sigma_input).
#
# MAPPING lambda BACK TO A PRIOR SD. For a weakly-identified parameter in the
# linear-Gaussian approximation the posterior displacement from the prior centre
# scales as sigma_prior^2 when the likelihood is weak, so
#
#     sigma_new  ~=  sigma_current * sqrt(lambda)
#
# i.e. lambda = 0.25 corresponds to HALVING the prior SD (fractions 0.4 -> 0.2).
# This is an approximation for choosing a target; the actual posterior needs a
# refit to confirm.
#
# Usage:  Rscript doublechecks/prior_tightening_curve.R [MODEL] [l1,l2,...]
# =============================================================================

suppressWarnings(suppressMessages({
  a <- commandArgs(trailingOnly = TRUE)
  MODEL <- if (length(a) >= 1) a[[1]] else "Yasso15"
  LAM   <- if (length(a) >= 2) as.numeric(strsplit(a[[2]], ",")[[1]])
           else c(1, 0.75, 0.5, 0.35, 0.25, 0.15, 0)
  library(BayesianTools)
}))
set.seed(2025)

src   <- readLines(file.path("Calibration_real_data_transient",
                             sprintf("run_%s_transient_calibration.R", MODEL)), warn = FALSE)
cutix <- grep("^t_run <- system.time", src)[1]
e <- new.env(parent = globalenv())
suppressMessages(source(textConnection(paste(src[seq_len(cutix-1L)], collapse="\n")), local = e))

st <- grep("ll_fn <- make_likelihood", src)[1]; op <- 0L; en <- NA_integer_
for (i in seq(st, length(src))) {
  ch <- strsplit(src[i], "")[[1]]; op <- op + sum(ch=="(") - sum(ch==")")
  if (op == 0L) { en <- i; break }
}
ml <- as.list(str2lang(sub("^\\s*ll_fn\\s*<-\\s*","",
                           paste(src[seq(st,en)], collapse="\n"))))[-1]
argof <- function(nm,d=NULL) if (is.null(ml[[nm]])) d else eval(ml[[nm]], envir=e)
to_original     <- argof("to_original")
assemble_params <- argof("assemble_params")
compute_xi      <- argof("compute_xi")
compute_xi_mean <- argof("compute_xi_mean")
steady_state    <- argof("steady_state")
run_model       <- argof("run_model")
STEADY_N        <- argof("steady_state_n", NULL)
to_unconstrained <- get("to_unconstrained", e)

plots <- get("plots", e); climate_by_plot <- get("climate_by_plot", e)
inputs_by_plot <- get("inputs_by_plot", e); litter_means <- get("litter_means", e)
obs_meta <- get("obs_meta", e); best_x <- get("best_x", e)
SOC_obs_all <- get("SOC_obs_all", e)
obs_years <- lapply(plots, function(p) SOC_obs_all$year[as.character(SOC_obs_all$plot_id)==p])
names(obs_years) <- plots

raw_J <- function(lm) if (!is.null(lm$J_total_mean)) unname(lm$J_total_mean) else
  sum(unlist(lm[intersect(c("nwl_mean","fwl_mean","cwl_mean"), names(lm))]))

fwd <- function(pid, mp) {
  clim <- climate_by_plot[[pid]]; inp <- inputs_by_plot[[pid]]
  lm <- litter_means[[pid]]; meta <- obs_meta[[pid]]
  if (any(is.na(meta$idx))) return(NULL)
  xa <- tryCatch(compute_xi(clim, mp), error=function(z) NULL); if (is.null(xa)) return(NULL)
  n_ss <- if (is.null(STEADY_N)) nrow(clim) else min(STEADY_N, nrow(clim))
  xs <- tryCatch(compute_xi_mean(clim[seq_len(n_ss),,drop=FALSE], mp), error=function(z) NULL)
  if (is.null(xs)) return(NULL)
  C0 <- tryCatch(steady_state(mp, lm, xs), error=function(z) NULL)
  if (is.null(C0) || any(!is.finite(C0)) || any(C0<0)) return(NULL)
  ro <- tryCatch(run_model(inp, mp, C0, xa), error=function(z) NULL); if (is.null(ro)) return(NULL)
  tr <- ro$total_soc; hat <- tr[meta$idx]
  if (any(!is.finite(hat)) || any(hat<=0)) return(NULL)
  list(traj=tr, years=if (!is.null(inp$year) && length(inp$year)==length(tr)) inp$year else NULL,
       hat=hat, obs=meta$soc_obs, oyr=obs_years[[pid]], J=raw_J(lm))
}
run_all <- function(p) { mp <- assemble_params(p); r <- lapply(plots, fwd, mp=mp)
  names(r) <- plots; r[!vapply(r, is.null, logical(1))] }
soc_med <- function(r) median(unlist(lapply(r, `[[`, "hat")))
obs_med <- function(r) median(unlist(lapply(r, `[[`, "obs")))
paired <- function(r,y0,y1,w="hat") mean(vapply(r, function(z){
  A<-z[[w]][z$oyr==y0]; B<-z[[w]][z$oyr==y1]
  if (length(A)!=1||length(B)!=1) NA_real_ else (B-A)/(y1-y0)}, numeric(1)), na.rm=TRUE)
annual <- function(r) { yl <- lapply(r,`[[`,"years"); ok <- !vapply(yl,is.null,logical(1))
  if (!any(ok)) return(NULL); yrs <- yl[[which(ok)[1]]]
  M <- vapply(r[ok], function(z) if (length(z$traj)==length(yrs)) z$traj else rep(NA_real_,length(yrs)),
              numeric(length(yrs)))
  setNames(rowMeans(M,na.rm=TRUE), yrs) }
mrt <- function(r, si) median(vapply(r, function(z){
  n<-length(z$traj); median(z$traj[seq(max(1,n-4),n)])/(si*z$J)}, numeric(1)), na.rm=TRUE)

p_def <- to_original(best_x)
rid <- sub(sprintf("^%s_posterior_(.+)\\.rds$", MODEL), "\\1",
           sort(list.files("Calibration_real_data_transient/runs",
                pattern=sprintf("^%s_posterior_[0-9]{8}_[0-9]{6}\\.rds$", MODEL)),
                decreasing=TRUE)[1])
smp <- getSample(readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_%s.rds",
                                 MODEL, rid)))
p_cal <- p_def; for (n in intersect(names(p_cal), colnames(smp))) p_cal[n] <- median(smp[,n])

AUX <- c("sigma_input","sigma_init")
x_def <- to_unconstrained(p_def); x_cal <- to_unconstrained(p_cal)
struct <- setdiff(names(p_def), AUX)

base <- run_all(p_cal); target <- obs_med(base)
o1 <- paired(base,1985,2024,"obs"); o2 <- paired(base,2006,2024,"obs")

cat(sprintf("\n===============  %s  (run %s)  ===============\n", MODEL, rid))
cat(sprintf("level-matched to observed median SOC = %.2f | sigma_init free at %.3f\n",
            target, p_cal[["sigma_init"]]))
cat(sprintf("OBSERVED trends: 1985-2024 %+.3f | 2006-2024 %+.3f\n\n", o1, o2))
cat(sprintf("%7s %9s %8s %9s %9s %9s %9s\n",
            "lambda","frac_SD","MRT","sig_in","t8524","t0624","pk->2024"))

out <- list()
for (lam in LAM) {
  x <- x_cal; x[struct] <- x_def[struct] + lam * (x_cal[struct] - x_def[struct])
  p <- to_original(x); p[AUX] <- p_cal[AUX]
  r0 <- run_all(p); if (!length(r0)) next
  p["sigma_input"] <- p["sigma_input"] * target / soc_med(r0)
  r <- run_all(p); tr <- annual(r)
  drop <- if (!is.null(tr)) tr[length(tr)] - max(tr) else NA_real_
  cat(sprintf("%7.2f %9.3f %8.2f %9.3f %+9.3f %+9.3f %+9.2f\n",
              lam, 0.4*sqrt(lam), mrt(r, p[["sigma_input"]]), p[["sigma_input"]],
              paired(r,1985,2024), paired(r,2006,2024), drop))
  out[[as.character(lam)]] <- tr
}
cat("\nfrac_SD = implied Tier-2 logit SD (current 0.400), via sigma*sqrt(lambda).\n")
cat("pk->2024 = drop from trajectory peak to 2024 (tC/ha); ~0 = flat plateau.\n")
