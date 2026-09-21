# =============================================================================
# woody_compartment_test.R   (2026-09-04)
#
# THE QUESTION (discussion memo D4). Yasso's stock is dead organic matter AND
# soil organic matter, and woody litter enters it. Our target is a SOIL
# measurement. If woody-derived carbon is in the model total but not in a soil
# core, the mismatch grows with standing biomass -- and a too-steep predicted
# gradient against basal area is exactly what we observe.
#
# THE TEST. Recompute each plot's modelled stock with the woody litter inputs
# zeroed, and see whether the between-plot gradient against basal area collapses
# toward the observed one.
#
#   full      as calibrated
#   no_cwl    coarse woody litter zeroed  (stumps, dead wood, logging residues --
#             the class that cannot be in a soil core)
#   no_woody  coarse AND fine woody zeroed (branches, woody roots as well)
#
# ⚠ The modelled stock used here is the transient-init ENDPOINT (the 1985 state),
# obtained from the same routine the calibration uses. It is a per-plot stock
# driven by that plot's own litter, which is what the gradient question is about.
# It is NOT the 2024 prediction; the forward run adds 39 more years of the same
# inputs and does not change the between-plot structure materially.
#
# ⚠ Zeroing is applied to ALL THREE litter variants a plot carries (_mean,
# _full_mean, _t0_mean) so the pre-run and the endpoint stay consistent.
#
# Usage:  Rscript doublechecks/woody_compartment_test.R [N_DRAW]
# =============================================================================
suppressWarnings(suppressMessages({
  a <- commandArgs(trailingOnly = TRUE)
  N_DRAW <- if (length(a) >= 1) as.integer(a[[1]]) else 5L
  library(BayesianTools)
}))
set.seed(2025)
MODELS   <- strsplit(Sys.getenv("HIKET_MODELS", "SP1,TP2,TP3,Yasso07,Yasso15,Yasso20"), ",")[[1]]
DIR_RUNS <- "Calibration_real_data_transient/runs"
NCORE    <- max(1L, parallel::detectCores() - 1L)

latest_run_id <- function(m) {
  fs <- list.files(DIR_RUNS, pattern = sprintf("^%s_posterior_[0-9]{8}_[0-9]{6}\\.rds$", m))
  if (!length(fs)) return(NA_character_)
  sub(sprintf("^%s_posterior_(.+)\\.rds$", m), "\\1", sort(fs, decreasing = TRUE)[1])
}
# same harness as init_state_plausibility.R
setup <- function(M) {
  src <- readLines(sprintf("Calibration_real_data_transient/run_%s_transient_calibration.R", M), warn = FALSE)
  cut <- grep("^t_run <- system.time", src)[1]
  e <- new.env(parent = globalenv())
  invisible(capture.output(suppressMessages(
    source(textConnection(paste(src[seq_len(cut - 1)], collapse = "\n")), local = e))))
  st <- grep("ll_fn <- make_likelihood", src)[1]; op <- 0L; en <- NA_integer_
  for (i in seq(st, length(src))) {
    ch <- strsplit(src[i], "")[[1]]; op <- op + sum(ch == "(") - sum(ch == ")")
    if (op == 0L) { en <- i; break }
  }
  ml <- as.list(str2lang(sub("^\\s*ll_fn\\s*<-\\s*", "", paste(src[seq(st, en)], collapse = "\n"))))[-1]
  for (k in c("to_original", "assemble_params", "compute_xi_mean", "steady_state"))
    assign(paste0(".", k), eval(ml[[k]], envir = e), envir = e)
  e
}
zero_woody <- function(lm, what) {
  for (v in c("_mean", "_full_mean", "_t0_mean")) {
    for (cls in what) {
      nm <- paste0(cls, v)
      if (!is.null(lm[[nm]])) lm[[nm]][] <- 0
    }
  }
  lm
}
stock_by_plot <- function(e, pf, plots, cbp, lms, what = character(0)) {
  mp <- e$.assemble_params(pf)
  v <- unlist(parallel::mclapply(plots, function(pid) {
    xi <- tryCatch(e$.compute_xi_mean(cbp[[pid]], mp), error = function(z) NULL)
    if (is.null(xi)) return(NA_real_)
    l <- if (length(what)) zero_woody(lms[[pid]], what) else lms[[pid]]
    tryCatch(sum(e$.steady_state(mp, l, xi)), error = function(z) NA_real_)
  }, mc.cores = NCORE))
  setNames(v, plots)
}
slope <- function(C, ba) {
  ok <- is.finite(C) & is.finite(ba) & C > 0
  unname(coef(lm(log(C[ok]) ~ ba[ok]))[2])
}

# observed reference, from any predictive bundle (observations are model-independent)
rid0 <- latest_run_id("Yasso15")
rd   <- readRDS(file.path(DIR_RUNS, sprintf("Yasso15_posterior_predictive_%s.rds", rid0)))$residuals_df
BA_of <- function(plots) rd$basal_area_85[match(as.integer(plots), rd$plot_id)]
oo <- !duplicated(rd$plot_id)
b_obs <- slope(rd$soc_obs_tCha[oo], rd$basal_area_85[oo])
RNG <- diff(range(rd$basal_area_85, na.rm = TRUE))

cat("\n=== Does the between-plot gradient come from woody litter? ===\n")
cat(sprintf("observed  d log(SOC)/dBA = %+.4f   (span x%.2f over %.0f m2/ha)\n\n", b_obs,
            exp(b_obs * RNG), RNG))
cat(sprintf("%-8s %10s %10s %10s | %8s %8s %8s | %s\n", "model",
            "full", "no_cwl", "no_woody", "x full", "x noCWL", "x noWdy", "woody share of input"))
res <- list()
for (M in MODELS) {
  rid <- latest_run_id(M); if (is.na(rid)) next
  e <- tryCatch(setup(M), error = function(z) NULL); if (is.null(e)) next
  ch <- file.path(DIR_RUNS, sprintf("%s_chains_%s.rds", M, rid)); if (!file.exists(ch)) next
  s_ <- do.call(rbind, lapply(readRDS(ch), function(z) getSample(z, parametersOnly = FALSE, start = 2)))
  free <- names(get("best_x", e)); plots <- get("plots", e)
  cbp <- get("climate_by_plot", e); lms <- get("litter_means", e)
  ba  <- BA_of(plots)
  ws  <- mean(vapply(plots, function(p) { l <- lms[[p]]
           w <- sum(l$fwl_mean, l$cwl_mean); w / (w + sum(l$nwl_mean)) }, numeric(1)), na.rm = TRUE)
  idx <- sample(nrow(s_), min(N_DRAW, nrow(s_)))
  B <- matrix(NA_real_, length(idx), 3, dimnames = list(NULL, c("full","no_cwl","no_woody")))
  for (i in seq_along(idx)) {
    pf <- tryCatch(e$.to_original(s_[idx[i], free]), error = function(z) NULL); if (is.null(pf)) next
    B[i,"full"]     <- slope(stock_by_plot(e, pf, plots, cbp, lms), ba)
    B[i,"no_cwl"]   <- slope(stock_by_plot(e, pf, plots, cbp, lms, "cwl"), ba)
    B[i,"no_woody"] <- slope(stock_by_plot(e, pf, plots, cbp, lms, c("fwl","cwl")), ba)
  }
  m <- colMeans(B, na.rm = TRUE)
  cat(sprintf("%-8s %10.4f %10.4f %10.4f | %8.2f %8.2f %8.2f | %.3f\n", M,
              m[1], m[2], m[3], exp(m[1]*RNG), exp(m[2]*RNG), exp(m[3]*RNG), ws))
  res[[M]] <- data.frame(model = M, run_id = rid, n_draw = sum(is.finite(B[,1])),
                         b_full = m[1], b_no_cwl = m[2], b_no_woody = m[3],
                         woody_share = ws, b_obs = b_obs, range_ba = RNG)
}
if (length(res)) {
  df <- do.call(rbind, res); saveRDS(df, "doublechecks/woody_compartment_test.rds")
  cat(sprintf("\nobserved slope %+.4f (span x%.2f). saved doublechecks/woody_compartment_test.rds\n",
              b_obs, exp(b_obs*RNG)))
}
