# =============================================================================
# prior_tightening_scenarios.R   (2026-08-10)
#
# Question: if we tightened ALL priors -- fractions AND the climate scaling --
# so each model stays near its published parameterisation, what happens to the
# F3 trajectory shape, the F4 initialisation level, and the bulk MRT?
#
# Three configs per model, all at the posterior median except where pinned:
#   CAL    as calibrated
#   FRAC   the 12 transfer fractions pinned at published centres
#   ALL    every STRUCTURAL parameter pinned at published centres (fractions,
#          climate scaling, woody size). Only the auxiliary sigma_input /
#          sigma_init remain free -- they are Tier-3 nuisance terms, not part of
#          any model's published parameterisation.
#
# LEVEL-MATCHING. sigma_input is re-solved in EVERY config so all three hit the
# same observed median SOC. Without this the comparison is rigged: SOC is exactly
# linear in sigma_input, so a config that over-predicts the level also inflates
# its own trend in tC/ha/yr by the same factor. Level-matching removes that and
# leaves shape, trend and MRT as the only differences.
#
# sigma_init is kept at its CALIBRATED value throughout (it is auxiliary, and it
# independently controls the stock-change sign -- pinning it would confound the
# question being asked).
#
# Usage:  Rscript doublechecks/prior_tightening_scenarios.R [MODEL]
# =============================================================================

suppressWarnings(suppressMessages({
  a <- commandArgs(trailingOnly = TRUE)
  MODEL <- if (length(a) >= 1) a[[1]] else "Yasso15"
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

AUX  <- c("sigma_input","sigma_init")
frac <- grep("^p_", names(p_cal), value=TRUE)
cfg <- list(CAL = p_cal,
            FRAC = { p <- p_cal; p[frac] <- p_def[frac]; p },
            ALL  = { p <- p_def; p[AUX] <- p_cal[AUX]; p })

base <- run_all(p_cal); target <- obs_med(base)
cat(sprintf("\n===============  %s  (run %s)  ===============\n", MODEL, rid))
cat(sprintf("all configs level-matched to observed median SOC = %.2f tC/ha\n", target))
cat(sprintf("sigma_init held at calibrated %.3f\n\n", p_cal[["sigma_init"]]))

res <- list()
for (k in names(cfg)) {
  p <- cfg[[k]]
  r0 <- run_all(p); if (!length(r0)) next
  p["sigma_input"] <- p["sigma_input"] * target / soc_med(r0)
  r <- run_all(p); res[[k]] <- list(p=p, r=r, tr=annual(r))
}

cat(sprintf("%-24s %10s %10s %10s\n", "", names(res)[1], names(res)[2], names(res)[3]))
row <- function(lab, f) cat(sprintf("%-24s %10s %10s %10s\n", lab,
  sprintf("%.3f", f(res[[1]])), sprintf("%.3f", f(res[[2]])), sprintf("%.3f", f(res[[3]]))))
row("sigma_input",        function(z) z$p[["sigma_input"]])
row("bulk MRT (yr)",      function(z) mrt(z$r, z$p[["sigma_input"]]))
row("paired 1985-2024",   function(z) paired(z$r,1985,2024))
row("paired 2006-2024",   function(z) paired(z$r,2006,2024))
cat(sprintf("%-24s %10.3f %10.3f\n", "  OBSERVED (same plots)",
            paired(base,1985,2024,"obs"), paired(base,2006,2024,"obs")))

if (!is.null(res[[1]]$tr)) {
  yrs <- as.integer(names(res[[1]]$tr))
  cat("\n--- annual MEAN SOC across plots (F3 shape) ---\n")
  cat(sprintf("%6s %10s %10s %10s\n", "year", names(res)[1], names(res)[2], names(res)[3]))
  for (y in intersect(c(1985,1995,2006,2010,2015,2024), yrs))
    cat(sprintf("%6d %10.2f %10.2f %10.2f\n", y,
        res[[1]]$tr[as.character(y)], res[[2]]$tr[as.character(y)], res[[3]]$tr[as.character(y)]))
  cat("\n")
  for (k in names(res)) {
    tr <- res[[k]]$tr
    cat(sprintf("%-6s peak %d (%.2f)  ->  2024 %.2f   drop %+.2f tC/ha\n",
        k, yrs[which.max(tr)], max(tr), tr[length(tr)], tr[length(tr)]-max(tr)))
  }
}
