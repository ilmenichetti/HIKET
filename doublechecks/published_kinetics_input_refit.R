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
# ⚠ This is a refit of two auxiliary parameters, not of the model. It puts the
# published kinetics at the point on the MRT x sigma_input ridge that agrees
# with our data, which is exactly what makes the forward comparison fair.
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

for (M in MODELS) {
  e <- setup(M)
  p_def <- e$.to_original(get("best_x", e))              # PUBLISHED point, physical
  om    <- get("obs_meta", e)
  keep  <- names(om)[vapply(om, function(z) length(z$soc_obs) >= 3L, logical(1))]

  # --- 1-D: sigma_input only, sigma_init at the published/default value -------
  s_init0 <- unname(p_def["sigma_init"])
  obj <- function(s) {
    p <- p_def; p["sigma_input"] <- s
    r <- pred_means(e, p, keep)
    if (is.null(r)) return(1e6)
    sum((log(r$pred) - log(r$obs))^2)                    # match all three campaigns
  }
  op <- optimize(obj, interval = c(0.3, 6), tol = 1e-3)
  s_hat <- op$minimum
  p <- p_def; p["sigma_input"] <- s_hat; r <- pred_means(e, p, keep)

  # baseline: sigma_input = 1 (tree litter as given), published kinetics
  p1 <- p_def; p1["sigma_input"] <- 1; r1 <- pred_means(e, p1, keep)

  cat(sprintf("=== %s ===  (published kinetics fixed; sigma_init = %.3f)\n", M, s_init0))
  cat(sprintf("  observed campaign means      : %.1f  %.1f  %.1f\n", r$obs[1], r$obs[2], r$obs[3]))
  cat(sprintf("  sigma_input = 1.000          : %.1f  %.1f  %.1f\n", r1$pred[1], r1$pred[2], r1$pred[3]))
  cat(sprintf("  BEST sigma_input = %.3f      : %.1f  %.1f  %.1f   (rms log %.4f, n=%d)\n",
              s_hat, r$pred[1], r$pred[2], r$pred[3], sqrt(op$objective/3), r$n))
  cat(sprintf("  => effective flux = %.2f tC/ha/yr  (our arm B for this model differs; see below)\n\n",
              s_hat * 2.511))
  OUT[[M]] <- list(s_hat = s_hat, sigma_init = s_init0, pred = r$pred, obs = r$obs,
                   pred_s1 = r1$pred, rms_log = sqrt(op$objective/3), n = r$n)
}

# --- the same metric for OUR calibration, so the comparison is like-for-like --
# Without this the published arm's fit has nothing to be judged against, and the
# whole point is whether the observed window can TELL THEM APART.
cat("\n================ OUR calibration, same metric ================\n")
rid <- vapply(MODELS, function(m) {
  fs <- list.files("Calibration_real_data_transient/runs",
                   pattern = sprintf("^%s_posterior_[0-9]{8}_[0-9]{6}\\.rds$", m))
  sub(sprintf("^%s_posterior_(.+)\\.rds$", m), "\\1", sort(fs, decreasing = TRUE)[1])
}, character(1))
cat("our RUN_IDs:", paste(names(rid), rid, sep="=", collapse=" | "), "\n\n")

for (M in MODELS) {
  e <- setup(M)
  p_def <- e$.to_original(get("best_x", e))
  om    <- get("obs_meta", e)
  keep  <- names(om)[vapply(om, function(z) length(z$soc_obs) >= 3L, logical(1))]
  smp   <- getSample(readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_%s.rds",
                                     M, rid[[M]])))
  p <- p_def
  for (n in intersect(names(p), colnames(smp))) p[n] <- median(smp[, n])
  r <- pred_means(e, p, keep)
  rms <- sqrt(sum((log(r$pred) - log(r$obs))^2) / 3)
  cat(sprintf("%-8s sigma_input %.3f  sigma_init %.3f : %.1f %.1f %.1f  (obs %.1f %.1f %.1f)  rms log %.4f\n",
              M, p["sigma_input"], p["sigma_init"], r$pred[1], r$pred[2], r$pred[3],
              r$obs[1], r$obs[2], r$obs[3], rms))
  OUT[[M]]$ours <- list(s_input = unname(p["sigma_input"]), s_init = unname(p["sigma_init"]),
                        pred = r$pred, rms_log = rms)
}

saveRDS(OUT, "doublechecks/published_kinetics_input_refit.rds")
cat("\nwrote doublechecks/published_kinetics_input_refit.rds\n")
