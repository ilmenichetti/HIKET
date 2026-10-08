# =============================================================================
# eqinit_common.R   (2026-10-08)
#
# Shared extraction for the equilibrium-init counterfactual figures (F16, F17, F19,
# T_eqinit). Reads both arms' predictive bundles and returns, per model x arm:
#   traj      draws x years (1985-2084): national mean SOC, history + projection
#   rates     draws x 3: mean annual change between campaigns, each plot on its
#             OWN observation years (as F5), so the denominators are the real spans
#   fore      per draw: C_2024, C_eq, headroom, sink over the first / last 20-yr
#             projection cycle (as build_forward_scenarios.R)
#   flux      effective litter input sigma_input x J(year), posterior median + 90%
#   pars      sigma_input / sigma_init posterior summaries
#   metrics   calibration / holdout metrics from the predictive stage
# plus the shared observations (campaign markers, observed rates).
#
# BASIS: balanced plot set (all three campaigns), unweighted, whole profile, true
# observation years -- obs_basis.R. Both arms use the SAME plot set (checked).
# Uncertainty: 90% interval over posterior draws of the national mean (parameter
# uncertainty of the mean, not plot scatter).
#
# Cached in manuscript/figures/eqinit_comparison.rds, stamped with both arms'
# RUN_IDs; rebuilt automatically when any RUN_ID changes.
# =============================================================================

source("manuscript/figures/run_ids_eqinit.R")
source("manuscript/figures/obs_basis.R")
source("manuscript/figures/model_palette.R")

EQ_CACHE  <- "manuscript/figures/eqinit_comparison.rds"
INTERVALS <- list("1985-2006" = c(1L, 2L), "2006-2024" = c(2L, 3L), "1985-2024" = c(1L, 3L))
ARM_LAB   <- c(prod = "transient start", eq = "equilibrium start")
ARM_LTY   <- c(prod = 1, eq = 2)

q90 <- function(x) c(med = median(x, na.rm = TRUE), lo = unname(quantile(x, .05, na.rm = TRUE)),
                     hi = unname(quantile(x, .95, na.rm = TRUE)))

# Per balanced plot: the true observation year of each campaign (NA if absent).
plot_campaign_years <- function(om, plots) {
  do.call(rbind, lapply(as.character(plots), function(p) {
    yr <- 1984L + om[[p]]$idx; cp <- campaign_of(yr)
    data.frame(plot_id = p, c1 = yr[cp == 1][1], c2 = yr[cp == 2][1], c3 = yr[cp == 3][1])
  }))
}

# Observed rates on the same basis: mean of per-plot (obs_b - obs_a)/dt, 95% CI.
observed_rates <- function(om, plots) {
  do.call(rbind, lapply(names(INTERVALS), function(iv) {
    ab <- INTERVALS[[iv]]
    r <- unlist(lapply(as.character(plots), function(p) {
      yr <- 1984L + om[[p]]$idx; cp <- campaign_of(yr); s <- om[[p]]$soc_obs
      ia <- which(cp == ab[1])[1]; ib <- which(cp == ab[2])[1]
      if (is.na(ia) || is.na(ib)) return(NULL)
      (s[ib] - s[ia]) / (yr[ib] - yr[ia])
    }))
    se <- sd(r) / sqrt(length(r))
    data.frame(interval = iv, obs = mean(r), lo = mean(r) - 1.96 * se, hi = mean(r) + 1.96 * se,
               n = length(r))
  }))
}

# Total annual litter per plot-year from an input bundle (simple: J_total; Yasso: AWEN x size).
litter_total <- function(inp) {
  if ("J_total" %in% names(inp)) return(inp$J_total)
  rowSums(inp[, grep("^(nwl|fwl|cwl)_", names(inp)), drop = FALSE])
}

extract_arm <- function(m, arm, bal) {
  f   <- arm_files(m, arm)
  b   <- readRDS(f$pred); ib <- readRDS(f$inp)
  bal_c <- as.character(bal)
  # --- trajectories: draw x year national mean ---------------------------------
  pp <- b$posterior_predictions;  pp <- pp[as.character(pp$plot_id) %in% bal_c, ]
  pj <- b$projection_predictions; pj <- pj[as.character(pj$plot_id) %in% bal_c, ]
  H  <- tapply(pp$total_soc, list(pp$draw, pp$year), mean)
  P  <- tapply(pj$total_soc, list(pj$draw, pj$year), mean)
  stopifnot(identical(rownames(H), rownames(P)))
  traj <- cbind(H, P)
  # --- rates, per plot on its own campaign years ---------------------------------
  cy <- plot_campaign_years(ib$obs_meta, bal)
  key <- paste(pp$plot_id, pp$year)
  at_c <- function(k) {                                  # draw x plot matrix at campaign k
    yk <- cy[[paste0("c", k)]]
    sub <- pp[key %in% paste(cy$plot_id, yk), c("plot_id", "draw", "total_soc")]
    M <- tapply(sub$total_soc, list(sub$draw, sub$plot_id), mean)
    M[, cy$plot_id, drop = FALSE]
  }
  C <- lapply(1:3, at_c)
  rates <- sapply(names(INTERVALS), function(iv) {
    ab <- INTERVALS[[iv]]
    dt <- cy[[paste0("c", ab[2])]] - cy[[paste0("c", ab[1])]]
    rowMeans(sweep(C[[ab[2]]] - C[[ab[1]]], 2, dt, "/"), na.rm = TRUE)
  })
  # --- forecast quantities (build_forward_scenarios.R definitions) ---------------
  first <- !duplicated(pj[c("draw", "plot_id")])
  C0 <- tapply(pj$C_last[first], pj$draw[first], mean)
  CE <- tapply(pj$C_eq[first],   pj$draw[first], mean)
  yr <- as.integer(colnames(P)); y0 <- min(yr) - 1L; y1 <- max(yr)
  atP <- function(y) P[, match(y, yr)]
  fore <- data.frame(C_2024 = C0, C_eq = CE, headroom = CE / C0 - 1,
                     sink_first20 = (atP(y0 + 20L) - C0) / 20,
                     sink_last20  = (atP(y1) - atP(y1 - 20L)) / 20)
  # --- parameters and effective litter flux --------------------------------------
  post <- readRDS(f$post)
  si   <- post[, "sigma_input"]
  sinit <- if ("sigma_init" %in% colnames(post)) post[, "sigma_init"] else rep(1, nrow(post))
  J <- do.call(rbind, lapply(bal_c, function(p) {
    x <- ib$inputs_by_plot[[p]]; data.frame(year = x$year, J = litter_total(x)) }))
  Jy <- tapply(J$J, J$year, mean)
  flux <- data.frame(year = as.integer(names(Jy)), J = as.numeric(Jy),
                     med = as.numeric(Jy) * median(si),
                     lo  = as.numeric(Jy) * quantile(si, .05),
                     hi  = as.numeric(Jy) * quantile(si, .95))
  list(traj = traj, rates = rates, fore = fore, flux = flux,
       pars = rbind(sigma_input = q90(si), sigma_init = q90(sinit)),
       metrics = list(calib = b$metrics_calib, holdout = b$metrics_holdout),
       rid = f$rid, n_plots = length(unique(pp$plot_id)))
}

load_eqinit_comparison <- function(force = FALSE) {
  stamp <- list(prod = RID[EQ_MODELS], eq = RID_EQ[EQ_MODELS])
  if (!force && file.exists(EQ_CACHE)) {
    z <- readRDS(EQ_CACHE)
    if (identical(z$stamp, stamp)) { message("eqinit cache: up to date"); return(z) }
    message("eqinit cache: RUN_IDs changed -- rebuilding")
  }
  out <- list(stamp = stamp, models = EQ_MODELS, missing = EQ_MISSING, arms = list(), obs = list())
  for (m in EQ_MODELS) {
    om  <- readRDS(arm_files(m, "prod")$inp)$obs_meta
    bal <- balanced_plots(om)
    bal_eq <- balanced_plots(readRDS(arm_files(m, "eq")$inp)$obs_meta)
    if (!setequal(bal, bal_eq)) stop(m, ": the two arms have different balanced plot sets")
    message(sprintf("  %s: extracting both arms (%d balanced plots)", m, length(bal)))
    out$arms[[m]] <- list(prod = extract_arm(m, "prod", bal), eq = extract_arm(m, "eq", bal))
    if (!length(out$obs)) out$obs <- list(camp = obs_campaigns(om, bal),
                                          rates = observed_rates(om, bal), n = length(bal))
  }
  saveRDS(out, EQ_CACHE)
  out
}
