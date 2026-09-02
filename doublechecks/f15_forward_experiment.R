# =============================================================================
# f15_forward_experiment.R   (2026-09-02)
#
# THE FORWARD EXPERIMENT. Two parameterisations of the same model that agree on
# everything observable, asked what the soil does under warming.
#
#   ARM "published" : published kinetics, sigma_input refitted so the arm matches
#                     the three observed campaign means, sigma_init PINNED to ours
#   ARM "ours"      : arm B posterior (RUN_IDs 20260831_1624*)
#
# Forcing: observed 1985-2024, then 60 projection years recycling the last 20
# observed years; warming applied to the PROJECTION ONLY. Response = warm - control.
#
# BOOTSTRAP. Our arm: posterior draws. Published Yasso15/20: FMI .dat draws.
# Yasso07 has NO .dat, so its published arm is a POINT -- a line against a band.
#
# Each published draw is placed on its own model's ridge at the level that matches
# the observed window: sigma_input = (s_ref * MRT_ref) / MRT_draw. Legitimate
# because SOC is EXACTLY linear in sigma_input (verified to machine precision), so
# this fixes the level without touching the kinetics -- and the RELATIVE response
# is invariant to it anyway.
#
# ⚠ Inputs are frozen under warming (no productivity response), so BOTH arms are
# biased toward loss. This is a sensitivity, not a projection.
#
# Usage:  Rscript doublechecks/f15_forward_experiment.R [N_DRAW] [N_PLOT]
# =============================================================================

suppressWarnings(suppressMessages(library(BayesianTools)))
ar <- commandArgs(trailingOnly = TRUE)
N_DRAW <- if (length(ar) >= 1) as.integer(ar[[1]]) else 40L
N_PLOT <- if (length(ar) >= 2) as.integer(ar[[2]]) else NA_integer_
set.seed(2025)

src <- readLines("doublechecks/published_arm_dat.R")
cut <- grep("^MR <- readRDS", src)[1]
eval(parse(text = paste(src[seq_len(cut-1)], collapse = "\n")), envir = globalenv())
MR  <- readRDS("doublechecks/intrinsic_mrt.rds")
ARM <- readRDS("doublechecks/published_arm_dat.rds")
rid <- vapply(c("Yasso07","Yasso15","Yasso20"), function(m) {
  fs <- list.files("Calibration_real_data_transient/runs",
                   pattern = sprintf("^%s_posterior_[0-9]{8}_[0-9]{6}\\.rds$", m))
  sub(sprintf("^%s_posterior_(.+)\\.rds$", m), "\\1", sort(fs, decreasing = TRUE)[1])
}, character(1))
N_PROJ <- 60L; N_RECYCLE <- 20L; DTS <- c(0, 2, 5)

extend <- function(inp, clim) {
  ny <- nrow(inp); rpy <- nrow(clim)/ny
  iy <- inp[(ny-N_RECYCLE+1):ny, , drop=FALSE]
  cy <- clim[(nrow(clim)-N_RECYCLE*rpy+1):nrow(clim), , drop=FALSE]
  reps <- ceiling(N_PROJ/N_RECYCLE)
  ie <- do.call(rbind, replicate(reps, iy, simplify=FALSE))[seq_len(N_PROJ), , drop=FALSE]
  ce <- do.call(rbind, replicate(reps, cy, simplify=FALSE))[seq_len(N_PROJ*rpy), , drop=FALSE]
  ie$year <- max(inp$year) + seq_len(N_PROJ)
  ce$year <- rep(max(inp$year) + seq_len(N_PROJ), each = rpy)
  list(inp = rbind(inp, ie), clim = rbind(clim, ce), n_obs = ny, rpy = rpy)
}
TCOL <- function(cl) if ("temp_air" %in% names(cl)) "temp_air" else "temp_mean"

# cross-plot MEAN TRAJECTORY (length n_obs + N_PROJ) for one parameter vector
traj_mean <- function(e, mp, keep, dT) {
  ibp <- get("inputs_by_plot", e); cbp <- get("climate_by_plot", e)
  lms <- get("litter_means", e); pl <- get("plots", e); pl <- pl[pl %in% keep]
  v <- parallel::mclapply(pl, function(pid) {
    x <- extend(ibp[[pid]], cbp[[pid]]); cl <- x$clim
    tc <- TCOL(cl); j <- (x$n_obs*x$rpy + 1):nrow(cl)
    if (dT != 0) cl[[tc]][j] <- cl[[tc]][j] + dT
    xa <- tryCatch(e$.compute_xi(cl, mp), error=function(z) NULL); if (is.null(xa)) return(NULL)
    xs <- tryCatch(e$.compute_xi_mean(cl[seq_len(x$n_obs*x$rpy), , drop=FALSE], mp),
                   error=function(z) NULL); if (is.null(xs)) return(NULL)
    C0 <- tryCatch(e$.steady_state(mp, lms[[pid]], xs), error=function(z) NULL)
    if (is.null(C0) || any(!is.finite(C0)) || any(C0 < 0)) return(NULL)
    ro <- tryCatch(e$.run_model(x$inp, mp, C0, xa), error=function(z) NULL); if (is.null(ro)) return(NULL)
    t <- ro$total_soc; if (any(!is.finite(t)) || any(t <= 0)) return(NULL)
    t
  }, mc.cores = NCORE)
  v <- v[!vapply(v, is.null, logical(1))]
  if (!length(v)) return(NULL)
  rowMeans(do.call(cbind, v))
}

OUT <- list()
for (M in c("Yasso07","Yasso15","Yasso20")) {
  e <- setup(M); p_def <- e$.to_original(get("best_x", e))
  om <- get("obs_meta", e)
  keep <- names(om)[vapply(om, function(z) length(z$soc_obs) >= 3L, logical(1))]
  if (!is.na(N_PLOT)) keep <- keep[seq_len(min(N_PLOT, length(keep)))]
  YP <- get(sprintf("%s_PARAM_NAMES", toupper(M)), envir = e)
  ss <- if (M=="Yasso07") get("yasso07_steady_state", e) else get("yasso15_steady_state", e)
  mrt_of <- function(mp) tryCatch({
    if (M=="Yasso07") { xi <- get("compute_xi_mean_yasso07", e)(MR$ref$clim, mp[["beta1"]], mp[["beta2"]], mp[["gamma"]])
      sum(ss(mp, MR$ref$nwl, MR$ref$fwl, MR$ref$cwl, xi))
    } else { xi <- get("compute_xi_mean_yasso15", e)(clim_ss=MR$ref$clim, params=mp[YP])
      sum(ss(mp, MR$ref$nwl, MR$ref$fwl, MR$ref$cwl, xi, precip_mean=MR$ref$clim$precip)) }
  }, error=function(z) NA_real_)

  # ---- assemble the draw sets --------------------------------------------
  smp <- getSample(readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_%s.rds", M, rid[[M]])))
  idx <- round(seq(1, nrow(smp), length.out = N_DRAW))
  ours_mp <- lapply(idx, function(i) { p <- p_def
    for (n in intersect(names(p), colnames(smp))) p[n] <- smp[i, n]; e$.assemble_params(p) })

  df <- file.path(DAT, paste0(M, ".dat"))
  if (file.exists(df)) {
    X <- as.matrix(read.table(df)); colnames(X) <- DCOL
    jj <- round(seq(1, nrow(X), length.out = N_DRAW))
    ref_prod <- ARM[[M]]$mrt_pub * ARM[[M]]$s_pub
    pub_mp <- lapply(jj, function(i) { p <- p_def; nm <- intersect(names(p), DCOL); p[nm] <- X[i, nm]
      mp <- e$.assemble_params(p); for (v in intersect(names(mp), DCOL)) mp[v] <- X[i, v]
      mp <- clamp_budget(mp); mp["sigma_init"] <- ARM[[M]]$sigma_init
      m <- mrt_of(mp); mp["sigma_input"] <- if (is.finite(m) && m > 0) ref_prod/m else ARM[[M]]$s_pub
      mp })
  } else {
    # Yasso07: a POINT, no band. ⚠ Force the fitted sigma_input at read time --
    # ARM$mp_pub has been stale before, and for Yasso07 there is no ridge-scaling
    # step to overwrite it, so the arm silently ran ~15% low.
    mp <- ARM[[M]]$mp_pub; mp["sigma_input"] <- ARM[[M]]$s_pub; pub_mp <- list(mp)
  }

  for (arm in c("published","ours")) {
    mps <- if (arm=="published") pub_mp else ours_mp
    t0 <- Sys.time()
    res <- lapply(seq_along(mps), function(i) {
      mp <- mps[[i]]
      tr <- lapply(DTS, function(d) traj_mean(e, mp, keep, d)); names(tr) <- paste0("dT", DTS)
      if (any(vapply(tr, is.null, logical(1)))) return(NULL)
      list(traj = tr, mrt = mrt_of(mp), s_input = unname(mp[["sigma_input"]]))
    })
    res <- res[!vapply(res, is.null, logical(1))]
    OUT[[M]][[arm]] <- res
    cat(sprintf("%-8s %-10s %3d/%3d draws  %.1f min\n", M, arm, length(res), length(mps),
                as.numeric(difftime(Sys.time(), t0, units="mins"))))
  }
}
saveRDS(list(out=OUT, n_draw=N_DRAW, n_proj=N_PROJ, dts=DTS, arm=ARM, rid=rid),
        "doublechecks/f15_forward_experiment.rds")
cat("\nwrote doublechecks/f15_forward_experiment.rds\n")
