# =============================================================================
# intrinsic_mrt.R   (2026-08-10)
#
# INTRINSIC model mean residence time: a property of the model's generator only.
#
#   MRT = 1' (-A)^-1 b    with b a UNIT input vector
#
# Implemented by calling each model's PURE steady-state routine with litter
# normalised to sum 1: MRT = sum(C_ss) in years.
#
# WHY THIS REPLACES THE EARLIER NUMBER. The engine binding named `steady_state`
# is NOT a steady state -- for Yasso it is `*_transient_init`, which equilibrates
# at the 1917 flux (J_full * sigma_init * sigma_input) and then RAMPS to 1985 on
# a flux carrying sigma_input only. So the earlier "equilibrium MRT" was
# C_1985 / (J * sigma_input * sigma_init): invariant to sigma_input (verified to
# 6 dp) but NOT to sigma_init (25.05 at 0.90 vs 33.84 at 0.35). Every MRT number
# computed that way is contaminated by an auxiliary parameter and is superseded.
#
# FIXED REFERENCE CONDITION -- identical for every model and every draw, so the
# comparison is of models, not of fits:
#   * unit total litter input, partitioned by the DATASET-MEAN AWEN x size
#     composition (nwl/fwl/cwl), so the size submodel is exercised realistically
#   * DATASET-MEAN climate (temp_mean, temp_amplitude, precip)
# Nothing from the SOC observations, sigma_input, sigma_init or the fit enters.
#
# MRT still depends on that reference (climate and litter chemistry); it is
# "model MRT under mean Finnish conditions", not a universal constant. The
# ours-vs-published comparison is unaffected by the choice, since both sides use
# the same reference.
#
# Usage:  Rscript doublechecks/intrinsic_mrt.R [N_DRAW]
# =============================================================================

suppressWarnings(suppressMessages({
  a <- commandArgs(trailingOnly = TRUE)
  N_DRAW <- if (length(a) >= 1) as.integer(a[[1]]) else 400L
  library(BayesianTools)
}))
set.seed(2025)

DAT  <- "Model_functions_real_data_transient/Decomposition_functions/Yasso_original"
DCOL <- c("alpha_A","alpha_W","alpha_E","alpha_N","p_WA","p_EA","p_NA","p_AW","p_EW","p_NW",
          "p_AE","p_WE","p_NE","p_AN","p_WN","p_EN","w1","w2","w3","w4","w5",
          "beta1","beta2","betaN1","betaN2","betaH1","betaH2","gamma","gammaN","gammaH",
          "p_H","alpha_H","delta1","delta2","r")
RID <- c(Yasso07="20260812_080941", Yasso15="20260812_080940", Yasso20="20260812_080940")  # 597031-33: sigma 0.80 + Student-t (2026-08-13)

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

out <- list()
for (M in names(RID)) {
  e <- setup(M); f <- mrt_fun(M, e)
  p_def <- e$.to_original(get("best_x", e))
  pub_pt <- tryCatch(f(p_def), error = function(z) NA_real_)

  smp <- getSample(readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_%s.rds",
                                   M, RID[[M]])))
  idx <- round(seq(1, nrow(smp), length.out = min(N_DRAW, nrow(smp))))
  ours <- vapply(idx, function(i) {
    p <- p_def; for (n in intersect(names(p), colnames(smp))) p[n] <- smp[i, n]
    tryCatch(f(p), error = function(z) NA_real_) }, numeric(1))

  pub <- NA
  df <- file.path(DAT, paste0(M, ".dat"))
  if (file.exists(df)) {
    X <- as.matrix(read.table(df)); colnames(X) <- DCOL
    j <- round(seq(1, nrow(X), length.out = min(N_DRAW, nrow(X))))
    pub <- vapply(j, function(i) {
      p <- p_def; nm <- intersect(names(p), DCOL); p[nm] <- X[i, nm]
      tryCatch(f(p), error = function(z) NA_real_) }, numeric(1))
    pub <- pub[is.finite(pub)]
  }
  ours <- ours[is.finite(ours)]
  out[[M]] <- list(ours = ours, pub = pub, pub_pt = pub_pt)

  cat(sprintf("=== %s ===\n", M))
  cat(sprintf("  published POINT     : %8.2f yr\n", pub_pt))
  if (length(pub) > 1)
    cat(sprintf("  published POSTERIOR : %8.2f  90%% [%.2f, %.2f]  n=%d\n",
                median(pub), quantile(pub,.05), quantile(pub,.95), length(pub)))
  cat(sprintf("  OUR posterior       : %8.2f  90%% [%.2f, %.2f]  n=%d\n\n",
              median(ours), quantile(ours,.05), quantile(ours,.95), length(ours)))
}
saveRDS(list(out = out, ref = ref), "doublechecks/intrinsic_mrt.rds")
cat("wrote doublechecks/intrinsic_mrt.rds\n")
