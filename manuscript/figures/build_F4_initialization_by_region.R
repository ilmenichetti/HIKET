source("manuscript/figures/run_ids.R")   # auto-selects current RUN_IDs
source("Calibration_real_data_transient/prerun_climate.R")   # hiket_ss_clim(): the run's pre-run climate
setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
# F4 BY REGION (2026-10-30) -- build_F4_initialization.R split into South and North,
# same panels and layout: per region a column of SOC (2/3) over input flux (1/3),
# 1917-2084, initialization | calibration | forecast. Decided with Samuli: regional
# split for measured and modelled trends; this is the candidate replacement for F4.
# Reconstruction code (wrapper sources + spin-up helpers) is copied VERBATIM from
# build_F4_initialization.R -- keep the two in step. Differences from F4:
#   * per-plot spin-up and input are kept, then averaged within region;
#   * the spin-up is averaged over the region's BALANCED plots (F4 averages over all
#     calibrated plots and then shifts the curve onto the balanced t0 stock);
#   * both columns share y limits, so the regional levels compare by eye.
# "North" = design weight 3 (soc_homogenized_plot.csv), as F_trends_by_region.
# Output: manuscript/figures/F4_initialization_by_region.png

rid <- as.list(RID)
source("manuscript/figures/model_palette.R")
source("manuscript/figures/obs_basis.R")
col <- MODEL_COL

.om_basis <- readRDS(sprintf("Data/model_inputs/Yasso20_inputs_%s.rds", RID[["Yasso20"]]))$obs_meta
BAL <- balanced_plots(.om_basis)
P   <- read.csv("Data/SOC_homogeneized/soc_homogenized_plot.csv")
WT  <- tapply(P$weight, P$plot_id, function(w) w[1])
REG <- list(South = BAL[WT[as.character(BAL)] == 1], North = BAL[WT[as.character(BAL)] == 3])
stopifnot(!anyNA(unlist(REG)), sum(lengths(REG)) == length(BAL))
message("F4 by region ", basis_note(BAL), "; South ", length(REG$South), ", North ", length(REG$North))
CACHE <- sprintf("manuscript/figures/F4reg_cache_bal_c3_%s.rds",
                 substr(paste(RID[FIG_MODELS], collapse = "-"), 1, 120))

if (!file.exists(CACHE)) {
  suppressMessages(library(BayesianTools))
  source("./Model_functions_real_data/Decomposition_functions/Yasso/yasso07_wrapper.R")
  source("./Model_functions_real_data/input_compatibility_layer.R")
  source("./Model_functions_real_data_transient/Decomposition_functions/SimpleModels/sp1_wrapper_transient.R")
  source("./Model_functions_real_data_transient/Decomposition_functions/SimpleModels/tp2_wrapper_transient.R")
  source("./Model_functions_real_data_transient/Decomposition_functions/SimpleModels/tp3_wrapper_transient.R")
  source("./Model_functions_real_data_transient/Decomposition_functions/Yasso/yasso07_wrapper_transient.R")
  source("./Model_functions_real_data_transient/Decomposition_functions/Yasso/yasso15_wrapper_transient.R")
  source("./Model_functions_real_data_transient/Decomposition_functions/Yasso/yasso20_wrapper_transient.R")
  source("./Prior_specs/Yasso07_priors.R"); source("./Prior_specs/Yasso15_priors.R"); source("./Prior_specs/Yasso20_priors.R")
  dyn.load("./Model_functions_real_data_transient/Decomposition_functions/Yasso/yasso07.so")
  dyn.load("./Model_functions_real_data_transient/Decomposition_functions/Yasso/yasso15.so")

  pre_shape <- function(lm, n = 68L)
    if (!is.null(lm$preinit_shape) && length(lm$preinit_shape) == n)
      lm$preinit_shape else (seq_len(n) - 1L) / (n - 1L)

  interp_awen <- function(v1917, v1985, n, f){
    m<-outer(1-f, v1917)+outer(f, v1985); colnames(m)<-names(v1917); m }
  spinup_yasso <- function(pm, lm, xm, PN, steady_fn, run_fn, xi_list){
    n<-68L; params<-pm[PN]; shp<-pre_shape(lm, n)
    n17<-lm$nwl_t0_mean*pm["sigma_init"]*pm["sigma_input"]; n85<-lm$nwl_t0_mean*pm["sigma_input"]
    f17<-lm$fwl_t0_mean*pm["sigma_init"]*pm["sigma_input"]; f85<-lm$fwl_t0_mean*pm["sigma_input"]
    c17<-lm$cwl_t0_mean*pm["sigma_init"]*pm["sigma_input"]; c85<-lm$cwl_t0_mean*pm["sigma_input"]
    if(xi_list) C0 <- steady_fn(params=params, nwl_mean=n17, fwl_mean=f17, cwl_mean=c17, xi_ss=xm, precip_mean=lm$precip_mean)
    else        C0 <- steady_fn(params=params, nwl_mean=n17, fwl_mean=f17, cwl_mean=c17, xi_mean=xm)
    idf <- data.frame(year=seq_len(n), interp_awen(n17,n85,n,shp), interp_awen(f17,f85,n,shp), interp_awen(c17,c85,n,shp))
    names(idf) <- c("year","nwl_A","nwl_W","nwl_E","nwl_N","fwl_A","fwl_W","fwl_E","fwl_N","cwl_A","cwl_W","cwl_E","cwl_N")
    if(xi_list){ xa<-list(xi_awe=rep(xm$xi_awe,n),xi_n=rep(xm$xi_n,n),xi_h=rep(xm$xi_h,n))
      out<-run_fn(input_df=idf, params=params, C_init=C0, xi_arrays=xa, precip=rep(lm$precip_mean,n))
    } else out<-run_fn(input_df=idf, params=params, C_init=C0, xi_array=rep(xm,n))
    c(sum(C0), out$total_soc[1:67])   # years 1917..1984
  }
  spinup_simple <- function(pm, lm, xm, step){
    si<-unname(pm["sigma_init"]); sinp<-unname(pm["sigma_input"])
    J17<-lm$J_t0_mean*si*sinp; J85<-lm$J_t0_mean*sinp; tr<-numeric(68); shp<-pre_shape(lm)
    if(step=="sp1"){ k<-unname(pm["alpha"])*xm; C<-J17/k
      for(i in 1:68){ Css<-(J17+(J85-J17)*shp[i])/k; C<-Css+(C-Css)*exp(-k); tr[i]<-C }; C0<-J17/k }
    if(step=="tp2"){ kA<-unname(pm["alpha_A"])*xm; kH<-unname(pm["alpha_H"])*xm; pH<-unname(pm["p_H"])
      A<-J17/kA; H<-pH*J17/kH; C0<-A+H
      for(i in 1:68){ J<-J17+(J85-J17)*shp[i]; s<-tp2_step(A,H,J,kA,kH,pH); A<-unname(s["A"]); H<-unname(s["H"]); tr[i]<-A+H } }
    if(step=="tp3"){ kA<-unname(pm["alpha_A"])*xm; kS<-unname(pm["alpha_S"])*xm; kH<-unname(pm["alpha_H"]); pS<-unname(pm["p_S"]); pH<-unname(pm["p_H"])
      C<-c(A=J17/kA, S=pS*J17/kS, H=pH*pS*J17/kH); C0<-sum(C)
      for(i in 1:68){ J<-J17+(J85-J17)*shp[i]; C<-.tp3_step(C["A"],C["S"],C["H"],kA,kS,kH,pS,pH,J); tr[i]<-sum(C) } }
    c(C0, tr[1:67])
  }


  # --- per-region aggregates of what F4 draws ---
  agg_stored <- function(m) {                      # 1985-2084 posterior band, per region
    b <- readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_predictive_%s.rds", m, rid[[m]]))
    pp <- rbind(b$posterior_predictions[, c("plot_id","year","draw","total_soc")],
                b$projection_predictions[, c("plot_id","year","draw","total_soc")])
    lapply(REG, function(pl) {
      x <- pp[pp$plot_id %in% as.character(pl), c("year","draw","total_soc")]
      ydr <- interaction(x$year, x$draw, drop = TRUE, lex.order = TRUE)
      mbar <- as.vector(rowsum(x$total_soc, ydr)) / as.vector(rowsum(rep(1, nrow(x)), ydr))
      yr <- as.integer(do.call(rbind, strsplit(levels(ydr), ".", fixed = TRUE))[, 1])
      data.frame(year = sort(unique(yr)), m = tapply(mbar, yr, mean),
                 lo = tapply(mbar, yr, quantile, 0.025, names = FALSE),
                 hi = tapply(mbar, yr, quantile, 0.975, names = FALSE)) })
  }
  .Jt0 <- function(L) if (!is.null(L$J_t0_mean)) sum(L$J_t0_mean) else
                      sum(L$nwl_t0_mean) + sum(L$fwl_t0_mean) + sum(L$cwl_t0_mean)
  .Jyr <- function(d) { k <- grep("^(nwl|fwl|cwl)_", names(d))
                        if (length(k)) rowSums(d[, k, drop = FALSE]) else d$J_total }

  spin <- list(); stored <- list(); flux <- list()
  for (m in names(rid)) {
    message("  reconstructing ", m)
    pkg  <- readRDS(sprintf("Data/model_inputs/%s_inputs_%s.rds", m, rid[[m]]))
    post <- readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_%s.rds", m, rid[[m]]))
    pm <- apply(getSample(post), 2, median); SSY <- pkg$STEADY_STATE_YEARS
    plots <- intersect(as.character(pkg$plots_real), as.character(BAL))
    hiket_set_prerun_clim(pkg$climate_ss_by_plot, plots)
    if (m == "TP2") { source("Prior_specs/TP2_priors.R"); pm["alpha_A"] <- TP2_ALPHA_A_FIXED }
    if (m == "TP3") { source("Prior_specs/TP3_priors.R"); pm["alpha_A"] <- TP3_ALPHA_A_FIXED }
    stored[[m]] <- agg_stored(m)
    if (m %in% c("SP1","TP2","TP3")) {
      step <- tolower(m)
      Mtx <- sapply(plots, function(pid){ clim<-pkg$climate_by_plot[[pid]]; n<-min(SSY,nrow(clim))
        xm<-unname(compute_xi_mean_yasso07(hiket_ss_clim(pid, clim, n), pm["beta1"], pm["beta2"], pm["gamma"]))
        spinup_simple(pm, pkg$litter_means[[pid]], xm, step) })
    } else {
      if(m=="Yasso07"){ PN<-YASSO07_PARAM_NAMES; DEF<-YASSO07_DEFAULT_PARAMS; steady<-yasso07_steady_state; runf<-yasso07_run; xil<-FALSE; xifn<-function(cl) compute_xi_mean_yasso07(cl, pm["beta1"],pm["beta2"],pm["gamma"]) }
      if(m=="Yasso15"){ PN<-YASSO15_PARAM_NAMES; DEF<-YASSO15_DEFAULT_PARAMS; steady<-yasso15_steady_state; runf<-yasso15_run; xil<-TRUE; xifn<-function(cl) compute_xi_mean_yasso15(cl, pm[PN]) }
      if(m=="Yasso20"){ PN<-YASSO20_PARAM_NAMES; DEF<-YASSO20_DEFAULT_PARAMS; steady<-yasso15_steady_state; runf<-yasso15_run; xil<-TRUE; xifn<-function(cl) compute_xi_mean_yasso20(cl, pm[PN]) }
      free <- setdiff(names(getSample(post)[1,]), c("sigma_init","sigma_input"))
      fixnm <- setdiff(PN, free); pmf <- c(DEF[fixnm], pm[free])[PN]; pmf <- c(pmf, sigma_init=unname(pm["sigma_init"]), sigma_input=unname(pm["sigma_input"]))
      mo <- if(m=="Yasso07") SSY else SSY*12L
      Mtx <- sapply(plots, function(pid){ clim<-pkg$climate_by_plot[[pid]]; n<-min(mo,nrow(clim))
        xm<-xifn(hiket_ss_clim(pid, clim, n)); spinup_yasso(pmf, pkg$litter_means[[pid]], xm, PN, steady, runf, xil) })
    }
    spin[[m]] <- lapply(REG, function(pl)
      data.frame(year = 1917:1984, m = rowMeans(Mtx[, intersect(colnames(Mtx), as.character(pl)), drop = FALSE], na.rm = TRUE)))

    # input flux: same three regimes as F4, per region
    S   <- getSample(post); nd <- min(2000L, nrow(S))
    Sdr <- S[sample.int(nrow(S), nd), c("sigma_init","sigma_input"), drop = FALSE]
    shp <- pre_shape(pkg$litter_means[[plots[1]]])
    flux[[m]] <- lapply(REG, function(pl) {
      pr   <- intersect(plots, as.character(pl))
      Jt0  <- mean(vapply(pr, function(p) .Jt0(pkg$litter_means[[p]]), numeric(1)), na.rm = TRUE)
      ann  <- do.call(rbind, lapply(pr, function(p) { d <- pkg$inputs_by_plot[[p]]
                data.frame(year = d$year, J = .Jyr(d)) }))
      Jann <- tapply(ann$J, ann$year, mean, na.rm = TRUE); yv <- as.integer(names(Jann))
      pre  <- outer(Sdr[,"sigma_input"], rep(1, length(shp))) * Jt0 *
              (outer(Sdr[,"sigma_init"], 1 - shp) + outer(rep(1, nd), shp))
      cal  <- outer(Sdr[,"sigma_input"], as.numeric(Jann))
      prj  <- outer(Sdr[,"sigma_input"], rep(as.numeric(Jann)[length(Jann)], 2084 - max(yv)))
      allJ <- cbind(pre, cal, prj)
      data.frame(year = c(1917:1984, yv, (max(yv)+1):2084), m = apply(allJ, 2, median),
                 lo = apply(allJ, 2, quantile, 0.025, names = FALSE),
                 hi = apply(allJ, 2, quantile, 0.975, names = FALSE)) })
  }
  saveRDS(list(spin = spin, stored = stored, flux = flux), CACHE)
}

cache <- readRDS(CACHE)
# per region: anchor each spin-up to the model's own 1985 stock (as F4), print a check
for (r in names(REG)) for (m in names(rid)) {
  s <- cache$spin[[m]][[r]]; t0 <- with(cache$stored[[m]][[r]], m[year == 1985])
  cat(sprintf("%-5s %-8s spin 1917=%.1f 1984=%.1f | stored 1985=%.1f (shift %+.2f)\n",
              r, m, s$m[1], tail(s$m, 1), t0, t0 - tail(s$m, 1)))
  cache$spin[[m]][[r]]$m <- s$m + (t0 - tail(s$m, 1))
}

ribbon <- function(d, keep, c0){ i<-d$year %in% keep
  polygon(c(d$year[i], rev(d$year[i])), c(d$lo[i], rev(d$hi[i])), col=adjustcolor(c0,0.12), border=NA) }
CM <- lapply(REG, function(pl) obs_campaigns(.om_basis, pl))
SY <- lapply(REG, function(pl) obs_subyears(.om_basis, pl, camp = 1L))
# shared limits across the two regions
yl <- range(unlist(lapply(names(REG), function(r) c(
        unlist(lapply(names(rid), function(m) c(cache$spin[[m]][[r]]$m, cache$stored[[m]][[r]]$lo,
                                                  cache$stored[[m]][[r]]$hi))),
        CM[[r]]$lo, CM[[r]]$hi, SY[[r]]$lo, SY[[r]]$hi))), finite = TRUE)
yl <- yl + c(-1, 1) * 0.04 * diff(yl)
fl <- range(unlist(lapply(cache$flux, function(z) lapply(z, function(d) c(d$lo, d$hi)))), finite = TRUE)
fl <- fl + c(-1, 1) * 0.06 * diff(fl)

# Print size: full text width (6.53 in), 9 pt. Columns = regions; in each column the
# SOC panel (2/3) sits on the input panel (1/3) and they share the x axis, as in F4.
png("manuscript/figures/F4_initialization_by_region.png", width = 6.53, height = 4.9, units = "in",
    res = 300, pointsize = 9)
layout(matrix(1:4, 2), heights = c(2, 1))
for (r in names(REG)) {
  first <- r == names(REG)[1]
  stored <- lapply(cache$stored, `[[`, r); spin <- lapply(cache$spin, `[[`, r)
  flux <- lapply(cache$flux, `[[`, r); cm <- CM[[r]]; sy <- SY[[r]]; yc <- round(cm$year)

  par(mar = c(0, if (first) 3.6 else 1.2, 1.4, 3.6), mgp = c(2.2, 0.5, 0), las = 1, tcl = -0.3)
  par(cex = 1)
  plot(NA, xlim = c(1917, 2084), ylim = yl, xlab = "", xaxt = "n", yaxt = if (first) "s" else "n",
       ylab = if (first) expression("Mean SOC  (tC ha"^-1*")") else "", main = "")
  if (!first) axis(2, labels = FALSE)
  mtext(sprintf("%s (%d plots)", r, length(REG[[r]])), side = 3, line = 0.2, adj = 0, font = 2)
  rect(1917, yl[1], 1985, yl[2], col = adjustcolor("grey85",  0.40), border = NA)
  rect(2024, yl[1], 2084, yl[2], col = adjustcolor("#9ec7e8", 0.22), border = NA)
  abline(v = c(1985, 2024), col = "grey60", lty = 2)
  ytop <- yl[2] - 0.03 * diff(yl)
  text(1951,   ytop, "initialization", cex = 0.8, col = "grey45", font = 3)
  text(2004.5, ytop, "calibration",    cex = 0.8, col = "grey45", font = 3)
  text(2054,   ytop, "forecast",       cex = 0.8, col = "grey45", font = 3)
  for (m in names(rid)) {
    ribbon(stored[[m]], 1985:2084, col[m])
    lines(spin[[m]]$year, spin[[m]]$m, col = col[m], lwd = 1.8)
    lines(stored[[m]]$year, stored[[m]]$m, col = col[m], lwd = 1.8)
  }
  arrows(sy$year, sy$lo, sy$year, sy$hi, angle = 90, code = 3, length = 0.02, col = SUBYEAR_COL, lwd = 1.2)
  points(sy$year, sy$m, pch = 21, bg = "white", col = SUBYEAR_COL, cex = 0.8, lwd = 1.4)
  lines(cm$year, cm$m, col = "black", lwd = 2.6)
  arrows(cm$year, cm$lo, cm$year, cm$hi, angle = 90, code = 3, length = 0.03, col = "firebrick", lwd = 1.8)
  points(cm$year, cm$m, pch = 21, bg = "firebrick", col = "black", cex = 1.15)
  # half-width columns: labels staggered in two rows (Biosoil lower), or they collide
  ylab_c <- yl[2] - c(0.11, 0.21, 0.11) * diff(yl)
  segments(cm$year, cm$hi, cm$year, ylab_c - 0.05 * diff(yl), col = adjustcolor("firebrick", 0.35), lty = 3)
  text(cm$year, ylab_c, sprintf("%s\n%d", CAMPAIGN_LABELS, yc), cex = 0.75, col = "firebrick")
  if (!first) legend("bottomright", bty = "n", cex = 0.8, pch = c(21, 21), pt.bg = c("firebrick", "white"),
                     col = c("black", SUBYEAR_COL), pt.cex = c(1.15, 0.8), pt.lwd = c(1, 1.4),
                     legend = c("campaign mean", "VMI8 sampling year (plot subset)"))
  lab  <- sort(sapply(stored, function(d) d$m[d$year == 2084]), decreasing = TRUE)
  laby <- lab; mind <- 0.055 * diff(yl)
  for (k in 2:length(laby)) if (laby[k-1] - laby[k] < mind) laby[k] <- laby[k-1] - mind
  text(2084, laby, names(lab), pos = 4, cex = 0.75, col = col[names(lab)], font = 2, xpd = NA)

  par(mar = c(3.2, if (first) 3.6 else 1.2, 0, 3.6), mgp = c(2.0, 0.5, 0))
  plot(NA, xlim = c(1917, 2084), ylim = fl, xlab = "Year", yaxt = if (first) "s" else "n",
       ylab = if (first) expression("Input  (tC ha"^-1*" yr"^-1*")") else "", main = "")
  if (!first) axis(2, labels = FALSE)
  rect(1917, fl[1], 1985, fl[2], col = adjustcolor("grey85",  0.40), border = NA)
  rect(2024, fl[1], 2084, fl[2], col = adjustcolor("#9ec7e8", 0.22), border = NA)
  abline(v = c(1985, 2024), col = "grey60", lty = 2)
  for (m in names(rid)) { d <- flux[[m]]
    polygon(c(d$year, rev(d$year)), c(d$lo, rev(d$hi)), col = adjustcolor(col[m], 0.12), border = NA) }
  for (m in names(rid)) lines(flux[[m]]$year, flux[[m]]$m, col = col[m], lwd = 1.8)
  if (first) text(1951, fl[2] - 0.17 * diff(fl), "shape after Liski\net al. (2006), assumed",
                  cex = 0.7, col = "grey40", font = 3)
}
dev.off(); cat("Wrote F4_initialization_by_region.png\n")
