source("manuscript/figures/run_ids.R")   # auto-selects current RUN_IDs
setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
# F4 -- THE INITIALIZATION PROBLEM (protagonist). TWO STACKED PANELS, shared x axis
# (layout, heights 2:1; x ticks and label on the LOWER panel only).
#   (a) TOP (2/3) : FULL SOC trajectory incl. the transient-init spin-up 1917->1985
#               (median line; modelled init device: constant climate + litter ramped
#               along the GROWING-STOCK shape from a below-equilibrium 1917 anchor)
#               + observed window 1985-2024 (mean + ribbon) + the observed campaign
#               path + the model-observation gap at campaign 1.
#   (b) BOTTOM (1/3) : the MEAN LITTER INPUT FLUX driving the panel above, 1917-2084,
#               with 95% posterior bands. This is the figure's own driver made visible:
#               the spin-up's flat-then-steep 1917-1985 ramp is an ASSUMPTION (see
#               M&M "the pre-run input shape"), and the reader should see it.
#   ⚠ 2026-08-20: the spin-up here hard-coded a LINEAR ramp while the wrappers have
#     used the growing-stock shape since C3 (2026-07-16). Fixed; caches re-keyed "c3".
#   ⚠ The old header quoted "63->102->105" and "+26 tC/ha": both come from the SOC
#     series RETIRED 2026-08-04 (no stoniness correction). On the corrected target the
#     campaign path is ~61/65/67 and the campaign-1 gap is NEGATIVE (~-9, an UNDER-
#     prediction). Do not reintroduce those numbers.
# Ribbons: posterior 2.5-97.5% of the cross-plot mean; same per-model colour, alpha 0.12.
# Heavy reconstruction is cached to F4_cache.rds (delete it to recompute).

rid <- as.list(RID)
source("manuscript/figures/model_palette.R")   # shared per-model palette (Temperature Diverging)
source("manuscript/figures/obs_basis.R")      # shared observed-SOC basis (see that file)
col <- MODEL_COL

# Balanced plot set: model curves and observed markers must describe ONE population
# (2026-08-12). A model curve is a single line and cannot change population by year.
.om_basis <- readRDS(sprintf("Data/model_inputs/Yasso20_inputs_%s.rds", RID[["Yasso20"]]))$obs_meta
BAL <- balanced_plots(.om_basis)
message("F4/F3 ", basis_note(BAL))
# Cache key includes the RUN_IDs, so a re-calibration invalidates it automatically.
# It used to be a fixed filename guarded by file.exists(), which meant that after a
# re-calibration this script "rebuilt" the figure from the PREVIOUS run's cache and
# reported success -- F4 and F3 (which reads this cache) were both silently stale.
# "bal" in the key: the cache now holds BALANCED-plot-set aggregates, so a cache
# written before 2026-08-12 must not be reused.
# "c3" in the key: caches written before 2026-08-20 hold a spin-up built on a LINEAR
# ramp, which the calibration has not used since C3. They must not be reused.
CACHE <- sprintf("manuscript/figures/F4_cache_bal_c3_%s.rds",
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

  # --- (i) stored trajectory bands 1985-2084 (posterior draws) ---
  agg_stored <- function(m) {
    b <- readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_predictive_%s.rds", m, rid[[m]]))
    pp <- rbind(b$posterior_predictions[, c("plot_id","year","draw","total_soc")],
                b$projection_predictions[, c("plot_id","year","draw","total_soc")])
    pp <- pp[pp$plot_id %in% BAL, c("year","draw","total_soc")]   # basis: see obs_basis.R
    ydr <- interaction(pp$year, pp$draw, drop = TRUE, lex.order = TRUE)
    mbar <- as.vector(rowsum(pp$total_soc, ydr)) / as.vector(rowsum(rep(1, nrow(pp)), ydr))
    yr <- as.integer(do.call(rbind, strsplit(levels(ydr), ".", fixed = TRUE))[, 1])
    data.frame(year = sort(unique(yr)), m = tapply(mbar, yr, mean),
               lo = tapply(mbar, yr, quantile, 0.025, names=FALSE),
               hi = tapply(mbar, yr, quantile, 0.975, names=FALSE))
  }

  # --- (ii) spin-up 1917-1984 cross-plot-mean median line ---
  # FIXED 2026-08-17: the 1917 anchor uses J_t0_mean, matching P1 in the wrappers
  # (tp2_wrapper_transient.R ~354, yasso15_wrapper_transient.R ~282). This script
  # still used J_full_mean, so the drawn initialization sat on a DIFFERENT flux
  # basis from the calibration that produced the posterior -- and sigma_init means
  # J_1917/J_1985, not J_1917/J_full.
  # SHAPE OF THE 1917->1985 RAMP. Until 2026-08-20 this script hard-coded a LINEAR
  # ramp, but the wrappers have followed the growing-stock shape since C3
  # (2026-07-16): tp2_wrapper_transient.R ~376, yasso15_wrapper_transient.R ~308 both
  # read lm$preinit_shape. So the figure named after the initialization problem was
  # drawing an initialization the calibration never ran -- same class of defect as the
  # J_full_mean/J_t0_mean bug fixed above. Now read from the bundle, linear fallback,
  # exactly as the wrappers do.
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

  spin <- list(); stored <- list(); flux <- list()
  for(m in names(rid)){
    message("  reconstructing ", m)
    pkg <- readRDS(sprintf("Data/model_inputs/%s_inputs_%s.rds", m, rid[[m]]))
    post <- readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_%s.rds", m, rid[[m]]))
    pm <- apply(getSample(post), 2, median); SSY <- pkg$STEADY_STATE_YEARS; plots <- pkg$plots_real
    # TP2/TP3: alpha_A is externally anchored to ICBM (C1) and held FIXED, so it is
    # NOT a posterior column and must be re-injected exactly as the calibration and
    # predictive scripts do. Without it the cascade coefficients are NA and the
    # spin-up fails for every plot.
    if (m == "TP2") { source("Prior_specs/TP2_priors.R"); pm["alpha_A"] <- TP2_ALPHA_A_FIXED }
    if (m == "TP3") { source("Prior_specs/TP3_priors.R"); pm["alpha_A"] <- TP3_ALPHA_A_FIXED }
    stored[[m]] <- agg_stored(m)
    if(m %in% c("SP1","TP2","TP3")){
      step <- tolower(m); if(step=="sp1"){}
      Mtx <- sapply(plots, function(pid){ clim<-pkg$climate_by_plot[[pid]]; n<-min(SSY,nrow(clim))
        xm<-unname(compute_xi_mean_yasso07(clim[seq_len(n),,drop=FALSE], pm["beta1"], pm["beta2"], pm["gamma"]))
        spinup_simple(pm, pkg$litter_means[[pid]], xm, step) })
    } else {
      if(m=="Yasso07"){ PN<-YASSO07_PARAM_NAMES; DEF<-YASSO07_DEFAULT_PARAMS; steady<-yasso07_steady_state; runf<-yasso07_run; xil<-FALSE; xifn<-function(cl) compute_xi_mean_yasso07(cl, pm["beta1"],pm["beta2"],pm["gamma"]); SS<-min(SSY,NA) }
      if(m=="Yasso15"){ PN<-YASSO15_PARAM_NAMES; DEF<-YASSO15_DEFAULT_PARAMS; steady<-yasso15_steady_state; runf<-yasso15_run; xil<-TRUE; xifn<-function(cl) compute_xi_mean_yasso15(cl, pm[PN]) }
      if(m=="Yasso20"){ PN<-YASSO20_PARAM_NAMES; DEF<-YASSO20_DEFAULT_PARAMS; steady<-yasso15_steady_state; runf<-yasso15_run; xil<-TRUE; xifn<-function(cl) compute_xi_mean_yasso20(cl, pm[PN]) }
      free <- setdiff(names(getSample(post)[1,]), c("sigma_init","sigma_input"))
      fixnm <- setdiff(PN, free); pmf <- c(DEF[fixnm], pm[free])[PN]; pmf <- c(pmf, sigma_init=unname(pm["sigma_init"]), sigma_input=unname(pm["sigma_input"]))
      mo <- if(m=="Yasso07") SSY else SSY*12L
      Mtx <- sapply(plots, function(pid){ clim<-pkg$climate_by_plot[[pid]]; n<-min(mo,nrow(clim))
        xm<-xifn(clim[seq_len(n),,drop=FALSE]); spinup_yasso(pmf, pkg$litter_means[[pid]], xm, PN, steady, runf, xil) })
    }
    spin[[m]] <- data.frame(year=1917:1984, m=rowMeans(Mtx, na.rm=TRUE))

    # --- (iii) MEAN LITTER INPUT FLUX 1917-2084, with posterior uncertainty ---
    # Same three regimes the models actually see:
    #   1917-1984  J(t) = J_t0 * sigma_input * [sigma_init + shape(t)*(1-sigma_init)]
    #   1985-2024  observed annual litter * sigma_input
    #   2025-2084  held at the 2024 value * sigma_input (what the predictive stage does:
    #              run_*_transient_predictive.R ~305 repeats the last input row)
    # Only sigma_init/sigma_input carry posterior uncertainty; the litter series itself
    # is data. Bands are therefore a faithful read of what the MCMC says about the flux.
    S    <- getSample(post); nd <- min(2000L, nrow(S))
    Sdr  <- S[sample.int(nrow(S), nd), c("sigma_init","sigma_input"), drop = FALSE]
    shp  <- pre_shape(pkg$litter_means[[plots[1]]])
    # Bundles differ by family: SP1/TP2/TP3 carry the aggregated J_t0_mean / J_total,
    # the Yassos carry AWEN components that must be summed. Handle both.
    # NB the Yasso *_t0_mean are length-4 AWEN vectors, not scalars -- sum them.
    .Jt0 <- function(L) if (!is.null(L$J_t0_mean)) sum(L$J_t0_mean) else
                        sum(L$nwl_t0_mean) + sum(L$fwl_t0_mean) + sum(L$cwl_t0_mean)
    .Jyr <- function(d) { k <- grep("^(nwl|fwl|cwl)_", names(d))
                          if (length(k)) rowSums(d[, k, drop = FALSE]) else d$J_total }
    Jt0  <- mean(vapply(plots, function(p) .Jt0(pkg$litter_means[[p]]), numeric(1)), na.rm = TRUE)
    ann  <- do.call(rbind, lapply(plots, function(p) {
              d <- pkg$inputs_by_plot[[p]]
              data.frame(year = d$year, J = .Jyr(d)) }))
    Jann <- tapply(ann$J, ann$year, mean, na.rm = TRUE)
    yv   <- as.integer(names(Jann))
    pre  <- outer(Sdr[,"sigma_input"], rep(1, length(shp))) * Jt0 *
            (outer(Sdr[,"sigma_init"], 1 - shp) + outer(rep(1, nd), shp))
    cal  <- outer(Sdr[,"sigma_input"], as.numeric(Jann))
    prj  <- outer(Sdr[,"sigma_input"], rep(as.numeric(Jann)[length(Jann)], 2084 - max(yv)))
    allJ <- cbind(pre, cal, prj)
    flux[[m]] <- data.frame(year = c(1917:1984, yv, (max(yv)+1):2084),
                            m  = apply(allJ, 2, median),
                            lo = apply(allJ, 2, quantile, 0.025, names = FALSE),
                            hi = apply(allJ, 2, quantile, 0.975, names = FALSE))
  }
  saveRDS(list(spin=spin, stored=stored, flux=flux), CACHE)
}

cache <- readRDS(CACHE); spin <- cache$spin; stored <- cache$stored; flux <- cache$flux

# observed campaign means +/- 95% CI
om <- readRDS(sprintf("Data/model_inputs/Yasso20_inputs_%s.rds", RID[["Yasso20"]]))$obs_meta
# ONE marker per campaign at its TRUE mean year (VMI8 ~1989, not 1985), on the
# BALANCED set -- the basis this figure declares. Addressed by campaign INDEX:
# a literal `year == 1985` no longer matches anything (see obs_basis.R).
cm  <- obs_campaigns(om, BAL)
sy  <- obs_subyears(om, BAL, camp = 1L)   # VMI8 per-sampling-year subsets
yc  <- round(cm$year)                       # modelled years to read the models at
m_c1   <- sapply(stored, function(d) d$m[d$year == yc[1]])
mstart <- mean(m_c1)
gap    <- mstart - cm$m[1]                  # SIGNED: >0 over-prediction, <0 under
obs_rate <- (cm$m[2] - cm$m[1]) / (cm$year[2] - cm$year[1])
mod_rate <- (mean(sapply(stored, function(d) d$m[d$year == yc[2]])) - mstart) /
            (cm$year[2] - cm$year[1])
# validation print
for(m in names(rid)) cat(sprintf("%-8s spin 1917=%.1f 1984=%.1f | stored%d=%.1f\n",
                                 m, spin[[m]]$m[1], tail(spin[[m]]$m,1), yc[1], m_c1[m]))
# Anchor each spin-up to the model's own t0 (1985) stock, so the spin-up and the
# stored trajectory JOIN. Must be the t0 year, NOT campaign 1's mean year (~1989):
# anchoring on the campaign left a visible step at 1985 in the merged panel.
m_t0 <- sapply(stored, function(d) d$m[d$year == 1985])
for(m in names(rid)) spin[[m]]$m <- spin[[m]]$m + (m_t0[m] - tail(spin[[m]]$m, 1))

ribbon <- function(d, keep, c0){ i<-d$year %in% keep
  polygon(c(d$year[i], rev(d$year[i])), c(d$lo[i], rev(d$hi[i])), col=adjustcolor(c0,0.12), border=NA) }

## ONE panel, 1917-2084: initialization | calibration | forecast, separated by shaded
## bands and rules. Was two panels; merged 2026-08-17.
## ylim is COMPUTED, never hardcoded -- the old ylim=c(55,114) sat above the anchored
## spin-up (which starts near 20 tC/ha), so the initialization was drawn off-scale and
## the panel silently showed nothing of the thing the figure is named after.
png("manuscript/figures/F4_initialization.png", width = 12.6, height = 8.6, units = "in", res = 200)
# Two stacked panels sharing one x axis: SOC (2/3) over input flux (1/3). Ticks and the
# "Year" label live on the LOWER panel only; the panels touch (top has bottom mar 0,
# bottom has top mar 0), so the shared axis reads as one.
layout(matrix(1:2, ncol = 1), heights = c(2, 1))
par(mar = c(0, 4.8, 2.6, 4.4), mgp = c(2.8, 0.7, 0), las = 1)

yl <- range(unlist(lapply(names(rid), function(m)
        c(spin[[m]]$m, stored[[m]]$lo, stored[[m]]$hi))), cm$lo, cm$hi, sy$lo, sy$hi, finite = TRUE)
yl <- yl + c(-1, 1) * 0.04 * diff(yl)

plot(NA, xlim = c(1917, 2084), ylim = yl, xlab = "", xaxt = "n",
     ylab = "Mean SOC across plots (tC/ha)", main = "")
rect(1917, yl[1], 1985, yl[2], col = adjustcolor("grey85",  0.40), border = NA)
rect(2024, yl[1], 2084, yl[2], col = adjustcolor("#9ec7e8", 0.22), border = NA)
abline(v = c(1985, 2024), col = "grey60", lty = 2)
ytop <- yl[2] - 0.03 * diff(yl)
text(1951,   ytop, "initialization", cex = 0.8, col = "grey45", font = 3)
text(2004.5, ytop, "calibration",    cex = 0.8, col = "grey45", font = 3)
text(2054,   ytop, "forecast",       cex = 0.8, col = "grey45", font = 3)

for (m in names(rid)) {
  ribbon(stored[[m]], 1985:2084, col[m])
  lines(spin[[m]]$year, spin[[m]]$m, col = col[m], lwd = 2)
  s <- stored[[m]]; lines(s$year, s$m, col = col[m], lwd = 2)
}
arrows(sy$year, sy$lo, sy$year, sy$hi, angle = 90, code = 3, length = 0.03, col = SUBYEAR_COL, lwd = 1.4)
points(sy$year, sy$m, pch = 21, bg = "white", col = SUBYEAR_COL, cex = 0.95, lwd = 1.6)
lines(cm$year, cm$m, col = "black", lwd = 3)
arrows(cm$year, cm$lo, cm$year, cm$hi, angle = 90, code = 3, length = 0.04, col = "firebrick", lwd = 2)
points(cm$year, cm$m, pch = 21, bg = "firebrick", col = "black", cex = 1.5)
text(cm$year, cm$lo, sprintf("%s\n%d", CAMPAIGN_LABELS, yc), pos = c(2,1,1), offset = 0.9,
     cex = 0.7, col = "firebrick")
legend("bottomright", bty = "n", cex = 0.72, pch = c(21, 21), pt.bg = c("firebrick", "white"),
       col = c("black", SUBYEAR_COL), pt.cex = c(1.5, 0.95), pt.lwd = c(1, 1.6),
       legend = c("campaign mean", "VMI8 sampling year (plot subset)"))

# per-model labels at the right edge, nudged apart where they coincide
lab  <- sort(sapply(stored, function(d) d$m[d$year == 2084]), decreasing = TRUE)
laby <- lab; mind <- 0.032 * diff(yl)
for (k in 2:length(laby)) if (laby[k-1] - laby[k] < mind) laby[k] <- laby[k-1] - mind
text(2084, laby, names(lab), pos = 4, cex = 0.7, col = col[names(lab)], font = 2, xpd = NA)

## ---- LOWER PANEL: the litter input flux driving the panel above -------------
par(mar = c(4.2, 4.8, 0, 4.4))
fl <- range(unlist(lapply(flux, function(d) c(d$lo, d$hi))), finite = TRUE)
fl <- fl + c(-1, 1) * 0.06 * diff(fl)
plot(NA, xlim = c(1917, 2084), ylim = fl, xlab = "Year",
     ylab = expression("Mean litter input (tC ha"^-1*" yr"^-1*")"), main = "")
rect(1917, fl[1], 1985, fl[2], col = adjustcolor("grey85",  0.40), border = NA)
rect(2024, fl[1], 2084, fl[2], col = adjustcolor("#9ec7e8", 0.22), border = NA)
abline(v = c(1985, 2024), col = "grey60", lty = 2)
for (m in names(rid)) {
  d <- flux[[m]]
  polygon(c(d$year, rev(d$year)), c(d$lo, rev(d$hi)),
          col = adjustcolor(col[m], 0.12), border = NA)
}
for (m in names(rid)) lines(flux[[m]]$year, flux[[m]]$m, col = col[m], lwd = 2)
# Name the assumption on the figure: the 1917-1985 segment is not data.
text(1951, fl[1] + 0.06 * diff(fl), "growing-stock ramp (assumed)",
     cex = 0.68, col = "grey40", font = 3)
dev.off()
cat(sprintf("\nstart mean=%.1f; gap=%+.1f; obs rate=%.2f mod rate=%.2f\nWrote F4_initialization.png\n", mstart, gap, obs_rate, mod_rate))
