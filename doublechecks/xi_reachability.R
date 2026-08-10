# =============================================================================
# xi_reachability.R   (2026-08-10)
#
# Can the calibration reach a SLOWER model at all? The MRT sweep scaled xi
# artificially; the sampler can only move xi through the climate parameters,
# which carry informative Tier-1 priors. So: draw from the model's OWN prior,
# push each draw through its OWN compute_xi, and ask where the xi needed to
# double MRT (= half the posterior xi) sits in that prior distribution.
#
# If half-xi is far in the prior tail, the slow-MRT solution is UNREACHABLE and
# the realised MRT is a consequence of fixed structure + informative climate
# priors, not something the SOC data chose.
#
# Usage:  Rscript doublechecks/xi_reachability.R [MODEL] [K]
# =============================================================================

suppressWarnings(suppressMessages({
  a     <- commandArgs(trailingOnly = TRUE)
  MODEL <- if (length(a) >= 1) a[[1]] else "Yasso15"
  K     <- if (length(a) >= 2) as.integer(a[[2]]) else 400L
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
argof <- function(nm, d=NULL) if (is.null(ml[[nm]])) d else eval(ml[[nm]], envir=e)
to_original     <- argof("to_original")
assemble_params <- argof("assemble_params")
compute_xi      <- argof("compute_xi")

plots <- get("plots", e); climate_by_plot <- get("climate_by_plot", e)
prior <- get("prior", e); best_x <- get("best_x", e)
sub   <- plots[seq_len(min(60, length(plots)))]        # subset: xi is plot-local

med_xi <- function(p_orig) {
  m <- assemble_params(p_orig)
  v <- unlist(lapply(sub, function(pid) {
    x <- tryCatch(compute_xi(climate_by_plot[[pid]], m), error=function(z) NULL)
    if (is.null(x)) return(NULL)
    if (is.list(x)) median(unlist(x)) else median(x)
  }))
  if (!length(v)) NA_real_ else median(v)
}

# posterior median (original space)
rid <- sub(sprintf("^%s_posterior_(.+)\\.rds$", MODEL), "\\1",
           sort(list.files("Calibration_real_data_transient/runs",
                pattern=sprintf("^%s_posterior_[0-9]{8}_[0-9]{6}\\.rds$", MODEL)),
                decreasing=TRUE)[1])
smp   <- getSample(readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_%s.rds",
                                   MODEL, rid)))
p_med <- to_original(best_x)
for (nm in intersect(names(p_med), colnames(smp))) p_med[nm] <- median(smp[, nm])

xi_post <- med_xi(p_med)
target  <- xi_post / 2

message(sprintf("[%s] drawing %d prior samples ...", MODEL, K))
draws <- prior$sampler(K)
if (is.null(dim(draws))) draws <- matrix(draws, nrow = K)
xi_prior <- vapply(seq_len(K), function(i)
  tryCatch(med_xi(to_original(draws[i, ])), error=function(z) NA_real_), numeric(1))
xi_prior <- xi_prior[is.finite(xi_prior)]

cat(sprintf("\n================  %s  (run %s)  ================\n", MODEL, rid))
cat(sprintf("posterior median xi          : %.4f\n", xi_post))
cat(sprintf("xi needed to DOUBLE MRT      : %.4f\n", target))
cat(sprintf("prior draws evaluated        : %d / %d\n", length(xi_prior), K))
cat(sprintf("prior xi  median             : %.4f\n", median(xi_prior)))
cat(sprintf("prior xi  90%% interval       : [%.4f, %.4f]\n",
            quantile(xi_prior,.05), quantile(xi_prior,.95)))
cat(sprintf("prior xi  full range         : [%.4f, %.4f]\n", min(xi_prior), max(xi_prior)))
frac <- mean(xi_prior <= target)
cat(sprintf("\nPRIOR MASS at or below the needed xi : %.2f%%  (%d of %d draws)\n",
            100*frac, sum(xi_prior <= target), length(xi_prior)))
cat(if (frac < 0.01)
      "=> effectively UNREACHABLE: the prior itself does not admit a model this slow.\n"
    else if (frac < 0.10)
      "=> reachable only in the far prior tail.\n"
    else "=> reachable within the prior.\n")
