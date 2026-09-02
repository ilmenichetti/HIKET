# =============================================================================
# climate_response_arms.R   (2026-09-02)
#
# THE CHECK BEFORE THE FIGURE. F15's panel (c) rests entirely on Yasso20 landing
# at ZERO: it has essentially no displacement between the published arm and ours
# (MRT x1.02, input x0.97), so if the climate response really tracks position on
# the MRT x sigma_input ridge, Yasso20 MUST NOT separate. Yasso07 (MRT x0.73,
# input x1.27) is the positive control -- it must.
#
#   null holds  -> build F15 as designed
#   null fails  -> the effect is NOT ridge position, and the figure is wrong
#
# SETUP. Forward from the transient init through 1985-2024 on observed forcing,
# then 60 projection years recycling the last 20 observed years (the same
# stationary-climate convention the predictive stage uses). Warming is applied to
# the PROJECTION segment only. Response = SOC(2084, warm) - SOC(2084, control),
# reported absolute (tC/ha) and relative (%).
#
# ⚠ Relative response is EXACTLY invariant to sigma_input (SOC is linear in it,
# verified to machine precision), so the % response isolates KINETICS. The
# absolute response carries both. Both are reported.
# ⚠ Inputs are frozen under warming -- no productivity response -- so BOTH arms
# are biased toward loss. This is a sensitivity, not a projection.
#
# Usage:  Rscript doublechecks/climate_response_arms.R [N_DRAW]
# =============================================================================

suppressWarnings(suppressMessages(library(BayesianTools)))
a <- commandArgs(trailingOnly = TRUE)
N_DRAW <- if (length(a) >= 1) as.integer(a[[1]]) else 0L    # 0 = medians only
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
N_PROJ <- 60L; N_RECYCLE <- 20L

# ⚠ Yasso07 carries ANNUAL climate (temp_mean/temp_amplitude, 1 row/yr); Yasso15
# and Yasso20 carry MONTHLY (temp_air, 12 rows/yr). Detect, do not assume.
extend <- function(inp, clim) {
  ny  <- nrow(inp); rpy <- nrow(clim) / ny            # rows of climate per year
  iy  <- inp[(ny-N_RECYCLE+1):ny, , drop = FALSE]
  nm  <- nrow(clim); cy <- clim[(nm-N_RECYCLE*rpy+1):nm, , drop = FALSE]
  reps <- ceiling(N_PROJ / N_RECYCLE)
  ie <- do.call(rbind, replicate(reps, iy, simplify = FALSE))[seq_len(N_PROJ), , drop = FALSE]
  ce <- do.call(rbind, replicate(reps, cy, simplify = FALSE))[seq_len(N_PROJ*rpy), , drop = FALSE]
  ie$year <- max(inp$year) + seq_len(N_PROJ)
  ce$year <- rep(max(inp$year) + seq_len(N_PROJ), each = rpy)
  list(inp = rbind(inp, ie), clim = rbind(clim, ce), n_obs = ny, rpy = rpy)
}
TCOL <- function(cl) if ("temp_air" %in% names(cl)) "temp_air" else "temp_mean"

soc_2084 <- function(e, mp, keep, dT) {
  ibp <- get("inputs_by_plot", e); cbp <- get("climate_by_plot", e)
  lms <- get("litter_means", e); pl <- get("plots", e); pl <- pl[pl %in% keep]
  v <- parallel::mclapply(pl, function(pid) {
    x <- extend(ibp[[pid]], cbp[[pid]])
    cl <- x$clim; tc <- TCOL(cl); j <- (x$n_obs*x$rpy + 1):nrow(cl)
    if (dT != 0) cl[[tc]][j] <- cl[[tc]][j] + dT
    xa <- tryCatch(e$.compute_xi(cl, mp),      error = function(z) NULL); if (is.null(xa)) return(NA_real_)
    xs <- tryCatch(e$.compute_xi_mean(cl[seq_len(x$n_obs*x$rpy), , drop=FALSE], mp),
                                                error = function(z) NULL); if (is.null(xs)) return(NA_real_)
    C0 <- tryCatch(e$.steady_state(mp, lms[[pid]], xs), error = function(z) NULL)
    if (is.null(C0) || any(!is.finite(C0)) || any(C0 < 0)) return(NA_real_)
    ro <- tryCatch(e$.run_model(x$inp, mp, C0, xa), error = function(z) NULL); if (is.null(ro)) return(NA_real_)
    tail(ro$total_soc, 1)
  }, mc.cores = NCORE)
  mean(unlist(v), na.rm = TRUE)
}

cat(sprintf("projection: %d yr recycling the last %d observed years; warming on the projection only\n\n",
            N_PROJ, N_RECYCLE))
for (M in c("Yasso07", "Yasso20")) {
  e <- setup(M); p_def <- e$.to_original(get("best_x", e))
  om <- get("obs_meta", e)
  keep <- names(om)[vapply(om, function(z) length(z$soc_obs) >= 3L, logical(1))]
  rid_m <- rid[[M]]
  smp <- getSample(readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_%s.rds", M, rid_m)))
  p_ours <- p_def; for (n in intersect(names(p_ours), colnames(smp))) p_ours[n] <- median(smp[, n])
  mp_ours <- e$.assemble_params(p_ours)
  # ⚠ ARM$mp_pub was saved BEFORE fit_input applied the fitted sigma_input, so the
  # published arm would otherwise run at the prior centre and sit ~25 tC/ha low.
  mp_pub  <- ARM[[M]]$mp_pub; mp_pub["sigma_input"] <- ARM[[M]]$s_pub

  cat(sprintf("=== %s ===  published MRT %.2f / sigma_input %.3f | ours MRT %.2f / %.3f\n",
              M, ARM[[M]]$mrt_pub, ARM[[M]]$s_pub, ARM[[M]]$mrt_ours, ARM[[M]]$s_ours))
  for (dT in c(2, 5)) {
    rows <- lapply(list(published = mp_pub, ours = mp_ours), function(mp) {
      c0 <- soc_2084(e, mp, keep, 0); cw <- soc_2084(e, mp, keep, dT)
      c(ctrl = c0, warm = cw, abs = cw - c0, rel = 100*(cw/c0 - 1))
    })
    p <- rows$published; o <- rows$ours
    cat(sprintf("  +%d C : published %6.2f -> %6.2f  (%+6.2f tC/ha, %+6.2f%%) | ours %6.2f -> %6.2f  (%+6.2f tC/ha, %+6.2f%%)\n",
                dT, p["ctrl"], p["warm"], p["abs"], p["rel"], o["ctrl"], o["warm"], o["abs"], o["rel"]))
    cat(sprintf("         DIFFERENCE ours - published : %+6.2f tC/ha   %+6.2f percentage points\n",
                o["abs"] - p["abs"], o["rel"] - p["rel"]))
  }
  cat("\n")
}
