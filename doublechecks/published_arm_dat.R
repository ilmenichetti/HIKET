# =============================================================================
# published_arm_dat.R   (2026-09-02)
#
# Builds the PUBLISHED arm for the forward experiment from the FMI .dat posterior
# samples -- NOT from to_original(best_x).
#
# ⚠⚠ WHY THIS FILE EXISTS. published_kinetics_input_refit.R used
# to_original(best_x) as "the published parameterisation". For Yasso07 and
# Yasso15 that is fine (their point MRT 33.47 / 30.38 matches the published
# posterior 30.44 for Y15). For YASSO20 IT IS NOT: Prior_specs/Yasso20_priors.R
# takes the 12 transfer-fraction centres from YASSO15, because Yasso20's
# fractions were freed. So best_x for Yasso20 is a HYBRID -- Yasso20 climate and
# size with Yasso15 fractions -- whose MRT is 19.03, against the genuine
# published Yasso20 posterior median of 25.04. A 32% error, and it inverted the
# sign of the published-vs-ours comparison for that model.
#
# F12_mrt_yasso.png plots the .dat POSTERIOR where it exists, which is why the
# figure disagreed with the refit. The figure was right.
#
# Yasso07 has no .dat file, so it keeps the point vector (its own published one).
#
# For each model: take a representative published draw (the one whose intrinsic
# MRT is closest to the published posterior median), then refit sigma_input alone
# with sigma_init PINNED at our posterior median, scoring the three campaign
# means on the balanced set.
#
# Usage:  Rscript doublechecks/published_arm_dat.R
# =============================================================================

suppressWarnings(suppressMessages(library(BayesianTools)))
set.seed(2025)
MODELS <- c("Yasso07", "Yasso15", "Yasso20")
NCORE  <- max(1L, parallel::detectCores() - 1L)
DAT    <- "Model_functions_real_data_transient/Decomposition_functions/Yasso_original"
DCOL   <- c("alpha_A","alpha_W","alpha_E","alpha_N","p_WA","p_EA","p_NA","p_AW","p_EW","p_NW",
            "p_AE","p_WE","p_NE","p_AN","p_WN","p_EN","w1","w2","w3","w4","w5",
            "beta1","beta2","betaN1","betaN2","betaH1","betaH2","gamma","gammaN","gammaH",
            "p_H","alpha_H","delta1","delta2","r")


# --- shared machinery (same as published_kinetics_input_refit.R) --------------
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
campaign_of <- function(y) ifelse(y <= 2000L, 1L, ifelse(y <= 2015L, 2L, 3L))
pred_means <- function(e, p_free, keep, mp = NULL) {
  if (is.null(mp)) mp <- e$.assemble_params(p_free)
  plots <- get("plots", e); cbp <- get("climate_by_plot", e)
  ibp <- get("inputs_by_plot", e); lms <- get("litter_means", e); om <- get("obs_meta", e)
  plots <- plots[plots %in% keep]
  res <- parallel::mclapply(plots, function(pid) {
    clim <- cbp[[pid]]; inputs <- ibp[[pid]]; lm <- lms[[pid]]; meta <- om[[pid]]
    if (any(is.na(meta$idx))) return(NULL)
    xa <- tryCatch(e$.compute_xi(clim, mp), error=function(z) NULL); if (is.null(xa)) return(NULL)
    xs <- tryCatch(e$.compute_xi_mean(clim, mp), error=function(z) NULL); if (is.null(xs)) return(NULL)
    C0 <- tryCatch(e$.steady_state(mp, lm, xs), error=function(z) NULL)
    if (is.null(C0) || any(!is.finite(C0)) || any(C0 < 0)) return(NULL)
    ro <- tryCatch(e$.run_model(inputs, mp, C0, xa), error=function(z) NULL); if (is.null(ro)) return(NULL)
    sv <- ro$total_soc[meta$idx]; if (any(!is.finite(sv)) || any(sv <= 0)) return(NULL)
    data.frame(camp = campaign_of(1984L + meta$idx), pred = sv, obs = meta$soc_obs)
  }, mc.cores = NCORE)
  d <- do.call(rbind, res[!vapply(res, is.null, logical(1))])
  if (is.null(d)) return(NULL)
  a <- aggregate(cbind(pred, obs) ~ camp, d, mean)
  list(pred = a$pred, obs = a$obs, n = nrow(d)/3)
}
fit_input <- function(e, p_base, keep, mp_base = NULL) {
  mk <- function(s) { if (is.null(mp_base)) { p <- p_base; p["sigma_input"] <- s
                        e$.assemble_params(p) } else { m <- mp_base; m["sigma_input"] <- s; m } }
  obj <- function(s) { r <- pred_means(e, NULL, keep, mp = mk(s))
    if (is.null(r)) return(1e6); sum((log(r$pred) - log(r$obs))^2) }
  op <- optimize(obj, interval = c(0.3, 6), tol = 1e-3)
  list(s = op$minimum, r = pred_means(e, NULL, keep, mp = mk(op$minimum)),
       rms = sqrt(op$objective/3))
}

# --- BUDGET CLAMP (2026-09-02) ----------------------------------------------
# Published Yasso20 is calibrated with the AWEN transfer fractions PLUS p_H
# exhausting the pool budget EXACTLY (median max = 1.000000; 100% of draws above
# 0.9999, 60.8% marginally above 1.0 on floating point). Yasso15 sits at 0.99939
# and only 10.5% of its draws reach the threshold.
#
# model_step (yasso07.f90 / yasso15.f90) declares a matrix "not well conditioned"
# when a column's off-diagonal sum exceeds 0.9999 x |diagonal| -- which is exactly
# "outgoing fractions + p_H > 0.9999" -- and falls back to an EXPLICIT EULER step.
# With alpha_W ~ 4.4 and xi up to 3.4 that step has factor ~ -14, so it oscillates
# and diverges. The fallback, not the parameters, is the problem: the exact path is
# stable AND fast (1 ms) at 0.9998, and SOC converges smoothly as the budget rises
#   0.9900 -> 51.27 | 0.9950 -> 52.01 | 0.9980 -> 52.47 | 0.9998 -> 52.75
# so the limit at 1.0 is ~52.76 and clamping to 0.9998 costs ~0.02%.
#
# ⬜ The PROPER fix is in the Fortran (raise the threshold; fall back only when the
# budget genuinely exceeds 1, where the system really does create mass). That
# touches a .so shared by Yasso15 and Yasso20 and needs its own verification, so it
# is deliberately NOT done here. This clamp is the minimal, documented alternative.
BUDGET_MAX <- 0.9998
OUTG <- list(A=c("p_AW","p_AE","p_AN"), W=c("p_WA","p_WE","p_WN"),
             E=c("p_EA","p_EW","p_EN"), N=c("p_NA","p_NW","p_NE"))
clamp_budget <- function(mp) {
  pH <- mp[["p_H"]]
  for (g in OUTG) {
    g <- intersect(g, names(mp)); if (!length(g)) next
    sm <- sum(mp[g]); tgt <- BUDGET_MAX - pH
    if (is.finite(sm) && sm > tgt) mp[g] <- mp[g] * (tgt / sm)
  }
  mp
}

MR <- readRDS("doublechecks/intrinsic_mrt.rds")
rid <- vapply(MODELS, function(m) {
  fs <- list.files("Calibration_real_data_transient/runs",
                   pattern = sprintf("^%s_posterior_[0-9]{8}_[0-9]{6}\\.rds$", m))
  sub(sprintf("^%s_posterior_(.+)\\.rds$", m), "\\1", sort(fs, decreasing=TRUE)[1])
}, character(1))
cat("cores:", NCORE, "| our RUN_IDs:", paste(names(rid), rid, sep="=", collapse=" | "), "\n\n")

OUT <- list()
for (M in MODELS) {
  e <- setup(M)
  p_def <- e$.to_original(get("best_x", e))
  om <- get("obs_meta", e)
  keep <- names(om)[vapply(om, function(z) length(z$soc_obs) >= 3L, logical(1))]
  smp <- getSample(readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_%s.rds", M, rid[[M]])))
  p_ours <- p_def; for (n in intersect(names(p_ours), colnames(smp))) p_ours[n] <- median(smp[, n])
  si_ours <- unname(p_ours["sigma_init"])

  # ---- PUBLISHED vector: full override (incl. the fixed_rates block) + clamp ----
  df <- file.path(DAT, paste0(M, ".dat"))
  if (file.exists(df)) {
    X <- as.matrix(read.table(df)); colnames(X) <- DCOL
    pub <- MR$out[[M]]$pub
    j2  <- round(seq(1, nrow(X), length.out = min(400, nrow(X))))
    k   <- j2[which.min(abs(pub - median(pub)))]
    p <- p_def; nm <- intersect(names(p), DCOL); p[nm] <- X[k, nm]
    mp_pub <- e$.assemble_params(p)
    # ⚠ assemble keeps fixed_rates -- for Yasso20 those are YASSO15's alphas and
    # w1-w5. Override the whole DCOL block so this really is the published model.
    for (v in intersect(names(mp_pub), DCOL)) mp_pub[v] <- X[k, v]
    basis <- sprintf(".dat draw %d (full override + budget clamp)", k)
  } else {
    mp_pub <- e$.assemble_params(p_def)          # Yasso07: its own published point
    basis  <- "published POINT (no .dat file)"
  }
  mp_pub <- clamp_budget(mp_pub)
  mp_pub["sigma_init"] <- si_ours                # PINNED to ours in both arms

  # intrinsic MRT of the arm actually used
  ss <- if (M == "Yasso07") get("yasso07_steady_state", e) else get("yasso15_steady_state", e)
  YP <- get(sprintf("%s_PARAM_NAMES", toupper(M)), envir = e)
  mrt_of <- function(mp) {
    if (M == "Yasso07") {
      xi <- get("compute_xi_mean_yasso07", e)(MR$ref$clim, mp[["beta1"]], mp[["beta2"]], mp[["gamma"]])
      sum(ss(mp, MR$ref$nwl, MR$ref$fwl, MR$ref$cwl, xi))
    } else {
      xi <- get("compute_xi_mean_yasso15", e)(clim_ss = MR$ref$clim, params = mp[YP])
      sum(ss(mp, MR$ref$nwl, MR$ref$fwl, MR$ref$cwl, xi, precip_mean = MR$ref$clim$precip))
    }
  }
  mrt_pub <- tryCatch(mrt_of(mp_pub), error = function(z) NA_real_)

  A <- fit_input(e, NULL, keep, mp_base = mp_pub)
  r_ours <- pred_means(e, p_ours, keep)
  rms_ours <- sqrt(sum((log(r_ours$pred) - log(r_ours$obs))^2)/3)
  mrt_ours <- median(MR$out[[M]]$ours)

  cat(sprintf("=== %s ===  basis: %s\n", M, basis))
  cat(sprintf("  observed                     : %.1f %.1f %.1f\n", A$r$obs[1],A$r$obs[2],A$r$obs[3]))
  cat(sprintf("  PUBLISHED MRT %5.2f  s_inp %.3f : %.1f %.1f %.1f   rms %.4f  effJ %.2f\n",
              mrt_pub, A$s, A$r$pred[1],A$r$pred[2],A$r$pred[3], A$rms, A$s*2.511))
  cat(sprintf("  OURS      MRT %5.2f  s_inp %.3f : %.1f %.1f %.1f   rms %.4f  effJ %.2f\n",
              mrt_ours, p_ours["sigma_input"], r_ours$pred[1],r_ours$pred[2],r_ours$pred[3],
              rms_ours, p_ours["sigma_input"]*2.511))
  cat(sprintf("  contrast: MRT x%.2f  input x%.2f   ridge product %.1f vs %.1f (%+.1f%%)\n\n",
              mrt_ours/mrt_pub, p_ours["sigma_input"]/A$s, mrt_pub*A$s, mrt_ours*p_ours["sigma_input"],
              100*(mrt_ours*p_ours["sigma_input"]/(mrt_pub*A$s)-1)))
  # ⚠ store mp_pub WITH the fitted sigma_input. It was saved without it once, and
  # every downstream run then used the prior centre and sat ~25 tC/ha low.
  mp_pub["sigma_input"] <- A$s
  OUT[[M]] <- list(basis=basis, mrt_pub=mrt_pub, s_pub=A$s, pred_pub=A$r$pred, rms_pub=A$rms,
                   mrt_ours=mrt_ours, s_ours=unname(p_ours["sigma_input"]),
                   pred_ours=r_ours$pred, rms_ours=rms_ours, sigma_init=si_ours,
                   obs=A$r$obs, mp_pub=mp_pub)
}

saveRDS(OUT, "doublechecks/published_arm_dat.rds")
cat("wrote doublechecks/published_arm_dat.rds\n")
