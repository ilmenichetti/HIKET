# =============================================================================
# published_kinetics_input_refit.R   (2026-09-02)
#
# THE QUESTION (Lorenzo). For the projection experiment we need a "published
# Yasso" arm that AGREES WITH THE DATA over the observed window, so that any
# difference in EXTRAPOLATION is attributable to the parameters and not to a
# level mismatch. So: holding the published kinetics fixed, what litter input
# would Yasso need to hit the three observed SOC campaign means?
#
# WHAT IS HELD FIXED. Everything except the auxiliaries: rates, fractions,
# climate parameters and woody size are the PUBLISHED point values
# (to_original(best_x), the same vector whose MRT is 33.47/30.38/19.03 in
# intrinsic_mrt.R). Only sigma_input -- and optionally sigma_init -- is free.
#
# WHAT IS MEASURED. Cross-plot mean SOC at the three campaigns on the BALANCED
# plot set (n = 310), against observed 64.5 / 71.5 / 73.6 tC/ha -- the decided
# reporting basis. This is a LEVEL-and-SHAPE match on the aggregate, which is
# what the projection experiment needs; it is NOT a calibration.
#
# ⚠ This is a refit of ONE auxiliary parameter, not of the model. It puts the
# published kinetics at the point on the MRT x sigma_input ridge that agrees
# with our data, which is exactly what makes the forward comparison fair.
#
# ⚠⚠ sigma_init IS PINNED TO OUR POSTERIOR MEDIAN IN BOTH ARMS (Lorenzo,
# 2026-09-02). The published default is 0.900 against our 0.53-0.58, so leaving
# it free would make the arms differ in TWO auxiliaries and put an
# initialisation difference inside a comparison that is supposed to be about
# inputs vs turnover. Pinned, the contrast is exactly two coordinates on the
# ridge: sigma_input against MRT. The sigma_init = 0.900 variant is still
# computed and reported, but only as a sensitivity.
#
# Usage:  Rscript doublechecks/published_kinetics_input_refit.R
# =============================================================================

suppressWarnings(suppressMessages(library(BayesianTools)))
set.seed(2025)
MODELS <- c("Yasso07", "Yasso15", "Yasso20")
NCORE  <- max(1L, parallel::detectCores() - 1L)

# --- environment with the calibration machinery, and the make_likelihood args -
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
  for (nm in c("to_original","assemble_params","compute_xi","compute_xi_mean",
               "steady_state","run_model"))
    assign(paste0(".", nm), eval(ml[[nm]], envir = e), envir = e)
  e
}

# --- predicted cross-plot mean SOC at the three campaigns ---------------------
campaign_of <- function(y) ifelse(y <= 2000L, 1L, ifelse(y <= 2015L, 2L, 3L))

pred_means <- function(e, p_free, keep) {
  mp <- e$.assemble_params(p_free)
  plots <- get("plots", e); cbp <- get("climate_by_plot", e)
  ibp <- get("inputs_by_plot", e); lms <- get("litter_means", e); om <- get("obs_meta", e)
  plots <- plots[plots %in% keep]
  res <- parallel::mclapply(plots, function(pid) {
    clim <- cbp[[pid]]; inputs <- ibp[[pid]]; lm <- lms[[pid]]; meta <- om[[pid]]
    if (any(is.na(meta$idx))) return(NULL)
    xa <- tryCatch(e$.compute_xi(clim, mp), error = function(z) NULL); if (is.null(xa)) return(NULL)
    xs <- tryCatch(e$.compute_xi_mean(clim, mp), error = function(z) NULL); if (is.null(xs)) return(NULL)
    C0 <- tryCatch(e$.steady_state(mp, lm, xs), error = function(z) NULL)
    if (is.null(C0) || any(!is.finite(C0)) || any(C0 < 0)) return(NULL)
    ro <- tryCatch(e$.run_model(inputs, mp, C0, xa), error = function(z) NULL); if (is.null(ro)) return(NULL)
    s <- ro$total_soc[meta$idx]; if (any(!is.finite(s)) || any(s <= 0)) return(NULL)
    data.frame(camp = campaign_of(1984L + meta$idx), pred = s, obs = meta$soc_obs)
  }, mc.cores = NCORE)
  d <- do.call(rbind, res[!vapply(res, is.null, logical(1))])
  if (is.null(d)) return(NULL)
  a <- aggregate(cbind(pred, obs) ~ camp, d, mean)
  list(pred = a$pred, obs = a$obs, n = nrow(d) / 3)
}

cat(sprintf("cores: %d\n\n", NCORE))
OUT <- list()

# OUR posterior medians -- needed BEFORE the published arm, because sigma_init is
# now pinned to them.
rid <- vapply(MODELS, function(m) {
  fs <- list.files("Calibration_real_data_transient/runs",
                   pattern = sprintf("^%s_posterior_[0-9]{8}_[0-9]{6}\\.rds$", m))
  sub(sprintf("^%s_posterior_(.+)\\.rds$", m), "\\1", sort(fs, decreasing = TRUE)[1])
}, character(1))
cat("our RUN_IDs:", paste(names(rid), rid, sep="=", collapse=" | "), "\n\n")

fit_input <- function(e, p_base, keep) {
  obj <- function(s) {
    p <- p_base; p["sigma_input"] <- s
    r <- pred_means(e, p, keep)
    if (is.null(r)) return(1e6)
    sum((log(r$pred) - log(r$obs))^2)
  }
  op <- optimize(obj, interval = c(0.3, 6), tol = 1e-3)
  p <- p_base; p["sigma_input"] <- op$minimum
  list(s = op$minimum, r = pred_means(e, p, keep), rms = sqrt(op$objective/3))
}

for (M in MODELS) {
  e <- setup(M)
  p_def <- e$.to_original(get("best_x", e))              # PUBLISHED point, physical
  om    <- get("obs_meta", e)
  keep  <- names(om)[vapply(om, function(z) length(z$soc_obs) >= 3L, logical(1))]

  smp <- getSample(readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_%s.rds",
                                   M, rid[[M]])))
  p_ours <- p_def
  for (n in intersect(names(p_ours), colnames(smp))) p_ours[n] <- median(smp[, n])
  si_ours <- unname(p_ours["sigma_init"])

  # PRIMARY: published kinetics, sigma_init PINNED to ours, sigma_input free
  pb <- p_def; pb["sigma_init"] <- si_ours
  A  <- fit_input(e, pb, keep)

  # SENSITIVITY: published kinetics with the published sigma_init 0.900
  B  <- fit_input(e, p_def, keep)

  # OURS, same metric
  r_ours <- pred_means(e, p_ours, keep)
  rms_ours <- sqrt(sum((log(r_ours$pred) - log(r_ours$obs))^2)/3)

  # published kinetics at sigma_input = 1, sigma_init pinned to ours
  p1 <- pb; p1["sigma_input"] <- 1; r1 <- pred_means(e, p1, keep)

  cat(sprintf("=== %s ===  (sigma_init PINNED at our %.3f in both arms)\n", M, si_ours))
  cat(sprintf("  observed                       : %.1f  %.1f  %.1f\n", A$r$obs[1], A$r$obs[2], A$r$obs[3]))
  cat(sprintf("  published kin., sigma_input 1.000: %.1f  %.1f  %.1f\n", r1$pred[1], r1$pred[2], r1$pred[3]))
  cat(sprintf("  published kin., BEST %.3f       : %.1f  %.1f  %.1f   rms log %.4f   effJ %.2f\n",
              A$s, A$r$pred[1], A$r$pred[2], A$r$pred[3], A$rms, A$s*2.511))
  cat(sprintf("  OURS,            sigma_input %.3f: %.1f  %.1f  %.1f   rms log %.4f   effJ %.2f\n",
              p_ours["sigma_input"], r_ours$pred[1], r_ours$pred[2], r_ours$pred[3],
              rms_ours, p_ours["sigma_input"]*2.511))
  cat(sprintf("  [sensitivity: sigma_init 0.900 -> sigma_input %.3f, rms log %.4f]\n\n",
              B$s, B$rms))

  OUT[[M]] <- list(
    pinned = list(sigma_init = si_ours, s_hat = A$s, pred = A$r$pred, rms_log = A$rms,
                  pred_s1 = r1$pred),
    sens_090 = list(sigma_init = unname(p_def["sigma_init"]), s_hat = B$s, rms_log = B$rms),
    ours = list(s_input = unname(p_ours["sigma_input"]), s_init = si_ours,
                pred = r_ours$pred, rms_log = rms_ours),
    obs = A$r$obs, n = A$r$n)
}

saveRDS(OUT, "doublechecks/published_kinetics_input_refit.rds")
cat("\nwrote doublechecks/published_kinetics_input_refit.rds\n")
