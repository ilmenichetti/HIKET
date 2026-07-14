setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
# F4 (redesign, 2-panel) -- THE INITIALIZATION PROBLEM (protagonist).
#   (a) LEFT  : FULL trajectory incl. the transient-init spin-up 1917->1985 (median line;
#               modelled init device: constant climate + linearly ramped litter from a
#               below-equilibrium 1917 anchor) + observed window 1985-2024 (mean + ribbon)
#               + observed accumulation path (63->102->105) + the +26 tC/ha 1985 gap.
#   (b) RIGHT : FUTURE forecast 2025-2084 (projection mean + ribbon) -- structural divergence.
# Ribbons: posterior 2.5-97.5% of the cross-plot mean; same per-model colour, alpha 0.12.
# Heavy reconstruction is cached to F4_cache.rds (delete it to recompute).

rid <- list(SP1="20260710_104903", TP2="20260710_104904", TP3="20260710_104904",
            Yasso07="20260710_104902", Yasso15="20260710_104902", Yasso20="20260710_102431")
col <- c(SP1="#5d4037", TP2="#c9922b", TP3="#f2c200",
         Yasso07="#1f6fb4", Yasso15="#d1495b", Yasso20="#2e8b57")
CACHE <- "manuscript/figures/F4_cache.rds"

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
    pp <- rbind(b$posterior_predictions[, c("year","draw","total_soc")],
                b$projection_predictions[, c("year","draw","total_soc")])
    ydr <- interaction(pp$year, pp$draw, drop = TRUE, lex.order = TRUE)
    mbar <- as.vector(rowsum(pp$total_soc, ydr)) / as.vector(rowsum(rep(1, nrow(pp)), ydr))
    yr <- as.integer(do.call(rbind, strsplit(levels(ydr), ".", fixed = TRUE))[, 1])
    data.frame(year = sort(unique(yr)), m = tapply(mbar, yr, mean),
               lo = tapply(mbar, yr, quantile, 0.025, names=FALSE),
               hi = tapply(mbar, yr, quantile, 0.975, names=FALSE))
  }

  # --- (ii) spin-up 1917-1984 cross-plot-mean median line ---
  interp_awen <- function(v1917, v1985, n){ f<-(seq_len(n)-1)/(n-1)
    m<-outer(1-f, v1917)+outer(f, v1985); colnames(m)<-names(v1917); m }
  spinup_yasso <- function(pm, lm, xm, PN, steady_fn, run_fn, xi_list){
    n<-68L; params<-pm[PN]
    n17<-lm$nwl_full_mean*pm["sigma_init"]*pm["sigma_input"]; n85<-lm$nwl_t0_mean*pm["sigma_input"]
    f17<-lm$fwl_full_mean*pm["sigma_init"]*pm["sigma_input"]; f85<-lm$fwl_t0_mean*pm["sigma_input"]
    c17<-lm$cwl_full_mean*pm["sigma_init"]*pm["sigma_input"]; c85<-lm$cwl_t0_mean*pm["sigma_input"]
    if(xi_list) C0 <- steady_fn(params=params, nwl_mean=n17, fwl_mean=f17, cwl_mean=c17, xi_ss=xm, precip_mean=lm$precip_mean)
    else        C0 <- steady_fn(params=params, nwl_mean=n17, fwl_mean=f17, cwl_mean=c17, xi_mean=xm)
    idf <- data.frame(year=seq_len(n), interp_awen(n17,n85,n), interp_awen(f17,f85,n), interp_awen(c17,c85,n))
    names(idf) <- c("year","nwl_A","nwl_W","nwl_E","nwl_N","fwl_A","fwl_W","fwl_E","fwl_N","cwl_A","cwl_W","cwl_E","cwl_N")
    if(xi_list){ xa<-list(xi_awe=rep(xm$xi_awe,n),xi_n=rep(xm$xi_n,n),xi_h=rep(xm$xi_h,n))
      out<-run_fn(input_df=idf, params=params, C_init=C0, xi_arrays=xa, precip=rep(lm$precip_mean,n))
    } else out<-run_fn(input_df=idf, params=params, C_init=C0, xi_array=rep(xm,n))
    c(sum(C0), out$total_soc[1:67])   # years 1917..1984
  }
  spinup_simple <- function(pm, lm, xm, step){
    si<-unname(pm["sigma_init"]); sinp<-unname(pm["sigma_input"])
    J17<-lm$J_full_mean*si*sinp; J85<-lm$J_t0_mean*sinp; tr<-numeric(68)
    if(step=="sp1"){ k<-unname(pm["alpha"])*xm; C<-J17/k
      for(i in 1:68){ Css<-(J17+(J85-J17)*(i-1)/67)/k; C<-Css+(C-Css)*exp(-k); tr[i]<-C }; C0<-J17/k }
    if(step=="tp2"){ kA<-unname(pm["alpha_A"])*xm; kH<-unname(pm["alpha_H"])*xm; pH<-unname(pm["p_H"])
      A<-J17/kA; H<-pH*J17/kH; C0<-A+H
      for(i in 1:68){ J<-J17+(J85-J17)*(i-1)/67; s<-tp2_step(A,H,J,kA,kH,pH); A<-unname(s["A"]); H<-unname(s["H"]); tr[i]<-A+H } }
    if(step=="tp3"){ kA<-unname(pm["alpha_A"])*xm; kS<-unname(pm["alpha_S"])*xm; kH<-unname(pm["alpha_H"]); pS<-unname(pm["p_S"]); pH<-unname(pm["p_H"])
      C<-c(A=J17/kA, S=pS*J17/kS, H=pH*pS*J17/kH); C0<-sum(C)
      for(i in 1:68){ J<-J17+(J85-J17)*(i-1)/67; C<-.tp3_step(C["A"],C["S"],C["H"],kA,kS,kH,pS,pH,J); tr[i]<-sum(C) } }
    c(C0, tr[1:67])
  }

  spin <- list(); stored <- list()
  for(m in names(rid)){
    message("  reconstructing ", m)
    pkg <- readRDS(sprintf("Data/model_inputs/%s_inputs_%s.rds", m, rid[[m]]))
    post <- readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_%s.rds", m, rid[[m]]))
    pm <- apply(getSample(post), 2, median); SSY <- pkg$STEADY_STATE_YEARS; plots <- pkg$plots_real
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
  }
  saveRDS(list(spin=spin, stored=stored), CACHE)
}

cache <- readRDS(CACHE); spin <- cache$spin; stored <- cache$stored

# observed campaign means +/- 95% CI
om <- readRDS("Data/model_inputs/Yasso20_inputs_20260710_102431.rds")$obs_meta
obs <- do.call(rbind, lapply(names(om), function(p){ z<-om[[p]]; if(!length(z$soc_obs)) return(NULL)
  data.frame(year=1984L+z$idx, soc=z$soc_obs) }))
cm <- aggregate(soc~year, obs, function(x) c(m=mean(x), lo=mean(x)-1.96*sd(x)/sqrt(length(x)), hi=mean(x)+1.96*sd(x)/sqrt(length(x))))
cm <- data.frame(year=cm$year, m=cm$soc[,"m"], lo=cm$soc[,"lo"], hi=cm$soc[,"hi"])
m1985 <- sapply(stored, function(d) d$m[d$year==1985]); mstart <- mean(m1985); gap <- mstart - cm$m[cm$year==1985]
obs_rate <- (cm$m[cm$year==2006]-cm$m[cm$year==1985])/(2006-1985)
mod_rate <- (mean(sapply(stored,function(d)d$m[d$year==2006]))-mstart)/(2006-1985)
# validation print
for(m in names(rid)) cat(sprintf("%-8s spin 1917=%.1f 1984=%.1f | stored1985=%.1f\n", m, spin[[m]]$m[1], tail(spin[[m]]$m,1), m1985[m]))
# anchor each spin-up to the calibrated 1985 state (removes the median-vs-mean kink;
# the spin-up shows the recovery SHAPE, its level is pinned to the posterior 1985 stock)
for(m in names(rid)) spin[[m]]$m <- spin[[m]]$m + (m1985[m] - tail(spin[[m]]$m, 1))

ribbon <- function(d, keep, c0){ i<-d$year %in% keep
  polygon(c(d$year[i], rev(d$year[i])), c(d$lo[i], rev(d$hi[i])), col=adjustcolor(c0,0.12), border=NA) }

png("manuscript/figures/F4_initialization.png", width = 11.8, height = 5.6, units = "in", res = 200)
layout(matrix(1:2, nrow = 1), widths = c(1.35, 1))
par(mar = c(4.0, 4.6, 3.2, 1.0), mgp = c(2.6, 0.7, 0), las = 1)

## (a) full trajectory: spin-up 1917 -> observed 2024
plot(NA, xlim=c(1917,2024), ylim=c(55,114), xlab="Year", ylab="Mean SOC across plots (tC/ha)",
     main="(a)  Full trajectory: transient-init spin-up + observed window")
rect(1917,54,1985,116, col=adjustcolor("grey85",0.35), border=NA)
text(1951, 57, "modelled transient initialization (1917->1985)", cex=0.72, col="grey45", font=3)
abline(v=1985, col="grey70", lty=2)
for(m in names(rid)){
  ribbon(stored[[m]], 1985:2024, col[m])
  lines(spin[[m]]$year, spin[[m]]$m, col=col[m], lwd=2, lty=1)
  s<-stored[[m]]; i<-s$year<=2024 & s$year>=1985; lines(s$year[i], s$m[i], col=col[m], lwd=2)
}
lines(cm$year, cm$m, col="black", lwd=3)
arrows(cm$year, cm$lo, cm$year, cm$hi, angle=90, code=3, length=0.04, col="firebrick", lwd=2)
points(cm$year, cm$m, pch=21, bg="firebrick", col="black", cex=1.5)
text(cm$year, cm$lo, c("VMI8\n1985","Biosoil\n2006","Komeetta\n2024"), pos=c(4,1,2), offset=0.8, cex=0.68, col="firebrick")
arrows(1986.4, cm$m[cm$year==1985], 1986.4, mstart, angle=90, code=3, length=0.04, col="grey20", lwd=1.6)
text(1990, 72, sprintf("+%.0f tC/ha\n1985 over-prediction", gap), pos=4, cex=0.74, col="grey15", font=2)
legend("bottomright", bty="n", cex=0.76, lwd=2, col=col, legend=names(col), title="model means", ncol=2)

## (b) future forecast 2024-2084
par(mar = c(4.0, 4.2, 3.2, 3.6))
yl2 <- range(sapply(stored, function(d){ i<-d$year>=2024; c(d$lo[i], d$hi[i]) }))
plot(NA, xlim=c(2024,2084), ylim=yl2, xlab="Year", ylab="", main="(b)  Future forecast: structural divergence")
for(m in names(rid)){ ribbon(stored[[m]], 2024:2084, col[m]); s<-stored[[m]]; i<-s$year>=2024; lines(s$year[i], s$m[i], col=col[m], lwd=2.2) }
lab <- sort(sapply(stored, function(d) d$m[d$year==2084]), decreasing=TRUE)
laby <- lab; mind <- 0.03*diff(yl2)              # spread near-coincident labels
for(k in 2:length(laby)) if(laby[k-1]-laby[k] < mind) laby[k] <- laby[k-1]-mind
text(2084, laby, names(lab), pos=4, cex=0.66, col=col[names(lab)], font=2, xpd=NA)
dev.off()
cat(sprintf("\nstart mean=%.1f; gap=+%.1f; obs rate=%.2f mod rate=%.2f\nWrote F4_initialization.png\n", mstart, gap, obs_rate, mod_rate))
