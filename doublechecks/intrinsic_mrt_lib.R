suppressWarnings(suppressMessages(library(BayesianTools)))
set.seed(2025)
DAT  <- "Model_functions_real_data_transient/Decomposition_functions/Yasso_original"
DCOL <- c("alpha_A","alpha_W","alpha_E","alpha_N","p_WA","p_EA","p_NA","p_AW","p_EW","p_NW",
          "p_AE","p_WE","p_NE","p_AN","p_WN","p_EN","w1","w2","w3","w4","w5",
          "beta1","beta2","betaN1","betaN2","betaH1","betaH2","gamma","gammaN","gammaH",
          "p_H","alpha_H","delta1","delta2","r")
# Auto-selected default (was a hard-coded 20260807 literal until 2026-08-20). Every
# consumer sources this lib and then overrides RID -- but only because they happen to
# source it BEFORE run_ids.R. A reorder would have silently reinstated a stale run, so
# the default is now the newest posterior rather than a literal.
RID <- vapply(c("Yasso07","Yasso15","Yasso20"), function(m) {
  fs <- list.files("Calibration_real_data_transient/runs",
                   pattern = sprintf("^%s_posterior_[0-9]{8}_[0-9]{6}\\.rds$", m))
  if (!length(fs)) stop("intrinsic_mrt_lib.R: no posterior found for ", m, call. = FALSE)
  sub(sprintf("^%s_posterior_(.+)\\.rds$", m), "\\1", sort(fs, decreasing = TRUE)[1])
}, character(1))

setup <- function(M) {
  src <- readLines(sprintf("Calibration_real_data_transient/run_%s_transient_calibration.R", M),
                   warn = FALSE)
  cut <- grep("^t_run <- system.time", src)[1]
  e <- new.env(parent = globalenv())
  suppressMessages(source(textConnection(paste(src[seq_len(cut-1)], collapse="\n")), local = e))
  st <- grep("ll_fn <- make_likelihood", src)[1]; op <- 0L; en <- NA_integer_
  for (i in seq(st, length(src))) {
    ch <- strsplit(src[i], "")[[1]]; op <- op + sum(ch=="(") - sum(ch==")")
    if (op == 0L) { en <- i; break }
  }
  ml <- as.list(str2lang(sub("^\\s*ll_fn\\s*<-\\s*","",
                             paste(src[seq(st,en)], collapse="\n"))))[-1]
  e$.to_original <- eval(ml[["to_original"]], envir = e)
  e$.assemble    <- eval(ml[["assemble_params"]], envir = e)
  e
}

# --- fixed reference condition, built once from the first model's data --------
ref <- local({
  e <- setup("Yasso15")
  plots <- get("plots", e); cbp <- get("climate_by_plot", e); lms <- get("litter_means", e)
  clim <- data.frame(
    temp_mean      = mean(vapply(plots, function(p) mean(cbp[[p]]$temp_mean),      numeric(1))),
    temp_amplitude = mean(vapply(plots, function(p) mean(cbp[[p]]$temp_amplitude), numeric(1))),
    precip         = mean(vapply(plots, function(p) mean(cbp[[p]]$precip),         numeric(1))))
  gm <- function(f) rowMeans(vapply(plots, function(p) as.numeric(lms[[p]][[f]]), numeric(4)))
  nwl <- gm("nwl_mean"); fwl <- gm("fwl_mean"); cwl <- gm("cwl_mean")
  tot <- sum(nwl) + sum(fwl) + sum(cwl)
  list(clim = clim, nwl = nwl/tot, fwl = fwl/tot, cwl = cwl/tot)   # sums to 1
})
cat(sprintf("\nReference condition (dataset means, unit input):\n"))
cat(sprintf("  T = %.2f C | amplitude = %.2f C | precip = %.0f mm\n",
            ref$clim$temp_mean, ref$clim$temp_amplitude, ref$clim$precip))
cat(sprintf("  litter split: nwl %.3f | fwl %.3f | cwl %.3f   (total %.3f)\n\n",
            sum(ref$nwl), sum(ref$fwl), sum(ref$cwl), sum(ref$nwl)+sum(ref$fwl)+sum(ref$cwl)))

mrt_fun <- function(M, e) {
  if (M == "Yasso07") {
    ss <- get("yasso07_steady_state", e); cxm <- get("compute_xi_mean_yasso07", e)
    function(p) {
      mp <- e$.assemble(p)
      xi <- cxm(ref$clim, mp[["beta1"]], mp[["beta2"]], mp[["gamma"]])
      sum(ss(mp, ref$nwl, ref$fwl, ref$cwl, xi))
    }
  } else {
    ss <- get("yasso15_steady_state", e); cxm <- get("compute_xi_mean_yasso15", e)
    .ypn <- sprintf("%s_PARAM_NAMES", toupper(M))
    YP <- if (exists(.ypn, envir = e, inherits = FALSE)) get(.ypn, envir = e) else NULL
    function(p) {
      mp <- e$.assemble(p)
      xi <- cxm(clim_ss = ref$clim, params = if (is.null(YP)) mp else mp[YP])
      sum(ss(mp, ref$nwl, ref$fwl, ref$cwl, xi, precip_mean = ref$clim$precip))
    }
  }
}

