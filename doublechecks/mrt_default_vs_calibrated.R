# =============================================================================
# mrt_default_vs_calibrated.R   (2026-08-10)
#
# Question: is the short bulk MRT something OUR calibration produced, or is it
# already there in each model's PUBLISHED parameterisation?
#
# If MRT at the published defaults ~= MRT after calibration, the short residence
# time is inherited from the model as published -- not a calibration artefact --
# and the finding transfers to the operational use of that model.
#
# NOTE on sigma_input: equilibrium MRT = C_ss / J is INVARIANT to sigma_input,
# because sigma_input scales the input and (in a linear system) the stock in
# exactly the same proportion. So this comparison isolates KINETICS + FRACTIONS
# + CLIMATE, and cannot be gamed by the input multiplier. The causal direction
# is therefore: fast kinetics FORCE a high sigma_input to reach the observed
# stock -- not the reverse.
#
# Two MRTs are reported:
#   MRT_eq   equilibrium, sum(C_ss) / J_ss  -- pure property of the parameters
#   MRT_fwd  end-of-run (2020-24) stock / effective flux -- what the run realises
#
# Usage:  Rscript doublechecks/mrt_default_vs_calibrated.R [MODEL]
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

raw_J_full <- function(lm) {
  if (!is.null(lm$J_full_mean)) return(unname(lm$J_full_mean))
  comp <- c("nwl_full_mean","fwl_full_mean","cwl_full_mean")
  hit  <- intersect(comp, names(lm))
  if (length(hit)) return(sum(unlist(lm[hit])))
  comp2 <- c("nwl_mean","fwl_mean","cwl_mean")
  sum(unlist(lm[intersect(comp2, names(lm))]))
}
raw_J_tot <- function(lm) {
  if (!is.null(lm$J_total_mean)) return(unname(lm$J_total_mean))
  sum(unlist(lm[intersect(c("nwl_mean","fwl_mean","cwl_mean"), names(lm))]))
}

mrt_at <- function(p) {
  mp <- assemble_params(p)
  si <- unname(p["sigma_input"]); s0 <- unname(p["sigma_init"])
  eq <- c(); fw <- c()
  for (pid in plots) {
    clim <- climate_by_plot[[pid]]; lm <- litter_means[[pid]]; meta <- obs_meta[[pid]]
    n_ss <- if (is.null(STEADY_N)) nrow(clim) else min(STEADY_N, nrow(clim))
    xs <- tryCatch(compute_xi_mean(clim[seq_len(n_ss),,drop=FALSE], mp), error=function(z) NULL)
    if (is.null(xs)) next
    Css <- tryCatch(steady_state(mp, lm, xs), error=function(z) NULL)
    if (is.null(Css) || any(!is.finite(Css))) next
    Jss <- raw_J_full(lm) * si * s0
    if (is.finite(Jss) && Jss > 0) eq <- c(eq, sum(Css)/Jss)

    xa <- tryCatch(compute_xi(clim, mp), error=function(z) NULL)
    if (is.null(xa)) next
    ro <- tryCatch(run_model(inputs_by_plot[[pid]], mp, Css, xa), error=function(z) NULL)
    if (is.null(ro)) next
    n <- length(ro$total_soc); tail_i <- seq(max(1, n-4), n)
    fw <- c(fw, median(ro$total_soc[tail_i]) / (si * raw_J_tot(lm)))
  }
  c(MRT_eq = median(eq, na.rm=TRUE), MRT_fwd = median(fw, na.rm=TRUE),
    n_eq = length(eq), sigma_input = si, sigma_init = s0)
}

p_def <- to_original(best_x)                       # published defaults / prior centres
rid <- sub(sprintf("^%s_posterior_(.+)\\.rds$", MODEL), "\\1",
           sort(list.files("Calibration_real_data_transient/runs",
                pattern=sprintf("^%s_posterior_[0-9]{8}_[0-9]{6}\\.rds$", MODEL)),
                decreasing=TRUE)[1])
smp <- getSample(readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_%s.rds",
                                 MODEL, rid)))
p_cal <- p_def
for (nm in intersect(names(p_cal), colnames(smp))) p_cal[nm] <- median(smp[, nm])

D <- mrt_at(p_def); C <- mrt_at(p_cal)

cat(sprintf("\n===============  %s  (calibrated run %s)  ===============\n", MODEL, rid))
cat(sprintf("%-26s %12s %12s\n", "", "PUBLISHED", "CALIBRATED"))
cat(sprintf("%-26s %12.2f %12.2f\n", "equilibrium MRT (yr)",  D["MRT_eq"],  C["MRT_eq"]))
cat(sprintf("%-26s %12.2f %12.2f\n", "forward-run MRT (yr)",  D["MRT_fwd"], C["MRT_fwd"]))
cat(sprintf("%-26s %12.3f %12.3f\n", "sigma_input",           D["sigma_input"], C["sigma_input"]))
cat(sprintf("%-26s %12.3f %12.3f\n", "sigma_init",            D["sigma_init"],  C["sigma_init"]))
cat(sprintf("%-26s %12d %12d\n",     "plots",       as.integer(D["n_eq"]), as.integer(C["n_eq"])))
cat(sprintf("\nratio calibrated/published (equilibrium MRT): %.2f\n", C["MRT_eq"]/D["MRT_eq"]))
cat("Equilibrium MRT is invariant to sigma_input -- this compares kinetics,\n")
cat("transfer fractions and climate response ONLY.\n")
