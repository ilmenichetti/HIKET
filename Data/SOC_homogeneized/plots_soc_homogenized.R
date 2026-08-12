# =============================================================================
# plots_soc_homogenized.R  — diagnostic figures for the homogenized SOC baseline
# Sourced by build_soc_homogenized.R (expects `layers`, `plot_tab`, `PLOTDIR`).
# Base R, PNG only. Checks: overall means, means by region, layers, depth, coverage.
# =============================================================================

CAMPS <- c("VMI8","Biosoil","Komeetta"); YRS <- c(1985,2006,2024)
COL_040 <- "#2C6FB5"; COL_1M <- "#c9922b"
wmean <- function(x,w) weighted.mean(x,w,na.rm=TRUE)
wse   <- function(x,w){ x<-x[!is.na(x)]; w<-w[seq_along(x)]; m<-weighted.mean(x,w)
  sqrt(sum(w^2)*weighted.mean((x-m)^2,w)/(sum(w)^2/length(w)))/sqrt(length(x)) } # design-ish SE

# --- 1. Overall trajectory: weighted mean 0-40 vs 1m per campaign -------------
# Balanced comparison: plots present in ALL three campaigns (avoids composition shifts).
plot_tab <- as.data.frame(plot_tab)                  # drop tibble quirks for base-R plotting
bal <- names(which(table(plot_tab$plot_id) == 3))
pt  <- plot_tab[plot_tab$plot_id %in% bal, ]
png(file.path(PLOTDIR,"01_trajectory_overall.png"), 900, 620, res=110)
agg <- do.call(rbind, lapply(YRS, function(y){ d<-pt[pt$year==y,]
  data.frame(year=y, m040=wmean(d$soc_0_40_Mgha,d$weight), m1m=wmean(d$soc_profile_Mgha,d$weight),
             s040=sd(d$soc_0_40_Mgha)/sqrt(nrow(d)), s1m=sd(d$soc_profile_Mgha)/sqrt(nrow(d)), n=nrow(d))}))
yl <- range(0, agg$m1m*1.15)
plot(agg$year, agg$m1m, type="o", pch=19, col=COL_1M, lwd=2, ylim=yl, xaxt="n",
     xlab="", ylab="SOC (Mg C / ha)", main="Homogenized SOC baseline — weighted mean trajectory")
arrows(agg$year, agg$m1m-1.96*agg$s1m, agg$year, agg$m1m+1.96*agg$s1m, angle=90, code=3, length=.04, col=COL_1M)
lines(agg$year, agg$m040, type="o", pch=17, col=COL_040, lwd=2)
arrows(agg$year, agg$m040-1.96*agg$s040, agg$year, agg$m040+1.96*agg$s040, angle=90, code=3, length=.04, col=COL_040)
axis(1, YRS, paste0(CAMPS,"\n",YRS))
text(agg$year, agg$m1m,  sprintf("%.1f",agg$m1m),  pos=3, col=COL_1M, cex=.9, font=2)
text(agg$year, agg$m040, sprintf("%.1f",agg$m040), pos=1, col=COL_040, cex=.9, font=2)
legend("bottomright", c("to 1 m (extrapolated)","0-40 cm (measured)"), col=c(COL_1M,COL_040),
       pch=c(19,17), lwd=2, bty="n")
mtext(sprintf("balanced subset present in all 3 campaigns: n = %d plots", agg$n[1]), 3, cex=.8)
dev.off()

# --- 2. Trajectory by region (North / South) --------------------------------
# Balanced subset (same plots all 3 campaigns) so the split is not confounded by
# year-to-year changes in which plots were sampled.
png(file.path(PLOTDIR,"02_trajectory_by_region.png"), 900, 620, res=110)
regs <- c("South","North"); rc <- c(South="#4a9c5d", North="#8c4a9c")
agg2 <- expand.grid(year=YRS, reg=regs)
agg2$m <- mapply(function(y,r){ d<-pt[pt$year==y & pt$region_label==r,]
  if(!nrow(d)) NA else mean(d$soc_profile_Mgha)}, agg2$year, agg2$reg)   # equal weights within a region
agg2$s <- mapply(function(y,r){ d<-pt[pt$year==y & pt$region_label==r,]
  if(!nrow(d)) NA else sd(d$soc_profile_Mgha)/sqrt(nrow(d))}, agg2$year, agg2$reg)
agg2$n <- mapply(function(y,r) sum(pt$year==y & pt$region_label==r), agg2$year, agg2$reg)
plot(NA, xlim=range(YRS)+c(-1,1), ylim=c(0, max(agg2$m,na.rm=TRUE)*1.15), xaxt="n",
     xlab="", ylab="SOC to 1 m (Mg C / ha)", main="By region (balanced subset)")
for(r in regs){ a<-agg2[agg2$reg==r,]
  lines(a$year, a$m, type="o", pch=19, col=rc[r], lwd=2)
  arrows(a$year, a$m-1.96*a$s, a$year, a$m+1.96*a$s, angle=90, code=3, length=.04, col=rc[r]) }
axis(1, YRS, paste0(CAMPS,"\n",YRS))
legend("topleft", sprintf("%s (n=%d)", regs, agg2$n[match(regs,agg2$reg)]), col=rc[regs], pch=19, lwd=2, bty="n")
dev.off()

# --- 3. C stock by layer, per campaign ---------------------------------------
png(file.path(PLOTDIR,"03_layer_boxplots.png"), 1050, 620, res=110)
lord <- c("organic","0-5cm","0-10cm","5-20cm","10-20cm","20-40cm")
L <- layers[layers$layer %in% lord,]; L$layer <- factor(L$layer, levels=lord)
L$campaign <- factor(L$campaign, levels=CAMPS)
boxplot(C_Mgha ~ campaign + layer, data=L, las=2, cex.axis=.6, outline=FALSE,
        col=rep(c("#9ecae1","#a1d99b","#fdae6b"),length(lord)),
        main="C stock by layer and campaign", ylab="C (Mg/ha)", xlab="")
legend("topright", CAMPS, fill=c("#9ecae1","#a1d99b","#fdae6b"), bty="n", cex=.8)
dev.off()

# --- 4. Depth: measured 0-40 vs depth-capped whole profile (per plot, 2024) ---
png(file.path(PLOTDIR,"04_depth_0_40_vs_profile.png"), 760, 720, res=110)
d24 <- plot_tab[plot_tab$year==2024,]
plot(d24$soc_0_40_Mgha, d24$soc_profile_Mgha, pch=21, bg="#c9922b55", col="#00000040",
     xlab="0-40 cm measured (Mg/ha)", ylab="whole profile, capped at soil depth (Mg/ha)",
     main="Depth extrapolation contribution (Komeetta 2024)", asp=1)
abline(0,1,lty=2,col="grey40")
mtext(sprintf("median deep addition = %.1f Mg/ha (%.0f%%); thin-soil plots sit on the 1:1 line (no extrapolation)",
      median(d24$soc_deep_Mgha,na.rm=TRUE),
      100*median(d24$soc_deep_Mgha/d24$soc_0_40_Mgha,na.rm=TRUE)), 3, cex=.7)
dev.off()

# --- 7. Extrapolation validation: predicted vs MEASURED 40-80 cm (2006, Krs 204) ---
png(file.path(PLOTDIR,"07_extrap_validation_40_80.png"), 760, 720, res=110)
v <- plot_tab[plot_tab$year==2006 & !is.na(plot_tab$soc_40_80_meas_Mgha) &
              !is.na(plot_tab$soc_40_80_pred_Mgha) & plot_tab$soc_40_80_meas_Mgha>0, ]
lim <- c(0, quantile(c(v$soc_40_80_meas_Mgha, v$soc_40_80_pred_Mgha), .98, na.rm=TRUE))
plot(v$soc_40_80_meas_Mgha, v$soc_40_80_pred_Mgha, pch=21, bg="#4a9c5d55", col="#00000040",
     xlim=lim, ylim=lim, asp=1,
     xlab="MEASURED 40-80 cm (Krs 204, Mg/ha)", ylab="extrapolated 40-80 cm (Mg/ha)",
     main="Extrapolation validated against measured deep layer")
abline(0,1,lty=2,col="grey40")
fit <- lm(soc_40_80_pred_Mgha ~ soc_40_80_meas_Mgha, v); abline(fit, col="firebrick", lwd=2)
mtext(sprintf("n=%d  median pred/obs=%.2f  bias=%.1f Mg/ha  (slight under = conservative)",
      nrow(v), median(v$soc_40_80_pred_Mgha/v$soc_40_80_meas_Mgha),
      mean(v$soc_40_80_pred_Mgha - v$soc_40_80_meas_Mgha)), 3, cex=.75)
dev.off()

# --- 5. Coverage: plots per campaign, all-three, region ----------------------
png(file.path(PLOTDIR,"05_coverage.png"), 900, 560, res=110)
par(mfrow=c(1,2))
np <- sapply(YRS, function(y) length(unique(plot_tab$plot_id[plot_tab$year==y])))
all3 <- sum(table(plot_tab$plot_id)==3)
bp <- barplot(c(np, `all 3`=all3), col=c("#9ecae1","#a1d99b","#fdae6b","grey60"),
        names.arg=c(CAMPS,"all 3"), main="Plot coverage", ylab="n plots", las=2, cex.names=.8)
text(bp, c(np,all3), c(np,all3), pos=3, cex=.8, xpd=NA)
rs <- table(plot_tab$region_source[!duplicated(plot_tab$plot_id)])
nice <- c(biosoil_design   = "from Biosoil\ndesign file\n(official N/S)",
          site_raw_region  = "from site_raw\nregion field",
          latitude_ETRS    = "inferred from\nplot northing\n(site_raw)",
          latitude_sitekey = "inferred from\nplot northing\n(site key)",
          default_South    = "no geodata:\ndefaulted\nto South")
names(rs) <- ifelse(names(rs) %in% names(nice), nice[names(rs)], names(rs))
par(mar=c(6,4,4,1))
bp2 <- barplot(rs, col="#2C6FB5", main="How each plot's North/South was assigned",
               ylab="n plots", las=1, cex.names=.68, cex.main=.95)
text(bp2, rs, rs, pos=3, cex=.8, xpd=NA)
mtext("North weight 3, South weight 1", side=1, line=4.4, cex=.7, col="grey40")
dev.off()

# --- 6. Distribution of total SOC (0-40) per campaign ------------------------
png(file.path(PLOTDIR,"06_hist_total.png"), 1000, 380, res=110)
par(mfrow=c(1,3))
for(y in YRS){ d<-plot_tab$soc_0_40_Mgha[plot_tab$year==y]
  hist(d, breaks=30, col="#9ecae1", border="white", main=paste0(CAMPS[YRS==y]," ",y),
       xlab="SOC 0-40 cm (Mg/ha)", xlim=c(0,200))
  abline(v=median(d,na.rm=TRUE), col="tomato", lwd=2, lty=2) }
dev.off()

# --- 8. New baseline vs OLD (inflated) calibration target, per plot x campaign ---
oldcands <- list.files(file.path(ROOT,"Data/model_inputs"), pattern="Yasso20_inputs_.*\\.rds", full.names=TRUE)
if (length(oldcands)) {
  om  <- readRDS(tail(sort(oldcands),1))$obs_meta
  old <- do.call(rbind, lapply(names(om), function(p){ z<-om[[p]]; if(!length(z$soc_obs)) return(NULL)
    data.frame(plot_id=as.integer(p), year=1984L+z$idx, old=z$soc_obs) }))
  cmp <- merge(old, plot_tab[,c("plot_id","year","soc_profile_Mgha")], by=c("plot_id","year"))
  png(file.path(PLOTDIR,"08_new_vs_old_target.png"), 780, 740, res=110)
  ycol <- c(`1985`="#4a9c5d", `2006`="#2C6FB5", `2024`="#c9922b")
  lim <- c(0, max(cmp$old, cmp$soc_profile_Mgha, na.rm=TRUE))
  plot(cmp$old, cmp$soc_profile_Mgha, pch=21, bg=paste0(ycol[as.character(cmp$year)],"88"),
       col="#00000030", xlim=lim, ylim=lim, asp=1,
       xlab="OLD calibration target (inflated, Mg/ha)", ylab="NEW baseline soc_profile (Mg/ha)",
       main="New homogenized target vs old inflated target")
  abline(0,1,lty=2,col="grey40"); abline(0, median(cmp$soc_profile_Mgha/cmp$old,na.rm=TRUE), col="firebrick", lwd=2)
  legend("topleft", c(paste0(names(ycol)," "), sprintf("median new/old = %.2f", median(cmp$soc_profile_Mgha/cmp$old,na.rm=TRUE))),
         col=c(ycol,NA), pch=c(19,19,19,NA), bty="n", cex=.85)
  dev.off()
}

# --- 9. Per-plot SOC change distributions (balanced subset) ------------------
yr_col <- function(y){ d<-pt[pt$year==y, c("plot_id","soc_profile_Mgha")]; names(d)[2]<-paste0("y",y); d }
chg <- merge(merge(yr_col(1985), yr_col(2006), by="plot_id"), yr_col(2024), by="plot_id")
d0685 <- chg$y2006 - chg$y1985
d2406 <- chg$y2024 - chg$y2006
change_hist <- function(x, ttl){ x <- x[is.finite(x)]
  hist(x, breaks=40, col="#9ecae1", border="white", main=ttl, xlab="change in SOC (Mg/ha)")
  abline(v=0, col="grey40", lty=3); abline(v=median(x), col="tomato", lwd=2)
  legend("topright", sprintf("median %+.1f\nmean %+.1f\n%% gaining %.0f",
         median(x), mean(x), 100*mean(x>0)), bty="n", cex=.8) }
png(file.path(PLOTDIR,"09_per_plot_change.png"), 1000, 460, res=110)
par(mfrow=c(1,2))
change_hist(d0685, "1985 -> 2006 (21 yr)")
change_hist(d2406, "2006 -> 2024 (18 yr)")
dev.off()

# --- 10. Maps: plots colored by soc_profile per campaign --------------------
socpal <- function(x, rng){ n<-100; cols<-hcl.colors(n,"YlOrBr",rev=TRUE)
  cols[pmax(1, pmin(n, as.integer((x-rng[1])/(rng[2]-rng[1])*n)+1))] }
rng <- quantile(plot_tab$soc_profile_Mgha, c(.02,.98), na.rm=TRUE)
mk <- tryCatch(if (requireNamespace("sf", quietly=TRUE)) readRDS(file.path(ROOT,"Data/model_inputs/maakunta_sf.rds")) else NULL, error=function(e) NULL)
png(file.path(PLOTDIR,"10_maps.png"), 1150, 560, res=110)
par(mfrow=c(1,3), mar=c(1,1,3,1))
for (y in YRS){ d<-plot_tab[plot_tab$year==y & !is.na(plot_tab$x_ETRS),]
  if (!is.null(mk)) plot(sf::st_geometry(mk), col="grey96", border="grey80", main=paste0(CAMPS[YRS==y]," ",y))
  else plot(d$x_ETRS, d$y_ETRS, type="n", asp=1, axes=FALSE, xlab="", ylab="", main=paste0(CAMPS[YRS==y]," ",y))
  points(d$x_ETRS, d$y_ETRS, pch=21, bg=socpal(d$soc_profile_Mgha,rng), col="#00000040", cex=.8) }
mtext(sprintf("soc_profile (Mg/ha), colour %.0f (pale) -> %.0f (dark)", rng[1], rng[2]), side=1, line=-1, outer=TRUE, cex=.75)
dev.off()

# --- 11. SOC vs drivers (2006 rows; climate joined from site_raw) -----------
s6 <- merge(plot_tab[plot_tab$year==2006,], site[,c("plot_id","mean_temp")], by="plot_id", all.x=TRUE)
png(file.path(PLOTDIR,"11_soc_vs_drivers.png"), 1100, 780, res=110)
par(mfrow=c(2,2), mar=c(4,4,2,1))
sc <- function(x,lab){ ok<-is.finite(x)&is.finite(s6$soc_profile_Mgha)
  plot(x, s6$soc_profile_Mgha, pch=21, bg="#4a9c5d55", col="#00000030", xlab=lab, ylab="soc_profile (Mg/ha)")
  if(sum(ok)>10){ f<-lm(s6$soc_profile_Mgha[ok]~x[ok]); abline(f,col="firebrick",lwd=2)
    mtext(sprintf("slope p=%.3g, r=%.2f", summary(f)$coef[2,4], cor(x[ok],s6$soc_profile_Mgha[ok])), 3, cex=.7)} }
sc(s6$mean_temp,       "mean annual T (degC)")
sc(s6$stand_age_2006,  "stand age 2006 (yr)")
sc(s6$coarse_frag_2006,"coarse fragments (%)")
boxplot(soc_profile_Mgha ~ site_fertility_2024, data=s6, col="#a1d99b",
        xlab="site fertility class (2024)", ylab="soc_profile (Mg/ha)", main="")
dev.off()

# --- 12. Organic-layer share of total SOC, by region x campaign -------------
plot_tab$org_share <- plot_tab$organic_Mgha / plot_tab$soc_0_40_Mgha
png(file.path(PLOTDIR,"12_organic_share.png"), 900, 560, res=110)
plot_tab$grp <- factor(paste(plot_tab$campaign, plot_tab$region_label), levels=as.vector(t(outer(CAMPS,c("South","North"),paste))))
boxplot(org_share ~ grp, data=plot_tab, las=2, col=rep(c("#4a9c5d","#8c4a9c"),3), cex.axis=.7,
        ylab="organic layer / (organic + 0-40 cm)", xlab="", main="Organic-layer share of SOC")
dev.off()

# --- 13. Outliers: extreme stocks and extreme per-plot changes --------------
png(file.path(PLOTDIR,"13_outliers.png"), 1050, 520, res=110)
par(mfrow=c(1,2), mar=c(4,4,3,1))
o <- plot_tab[order(-plot_tab$soc_profile_Mgha),][1:12,]
barplot(rev(o$soc_profile_Mgha), horiz=TRUE, names.arg=rev(paste0(o$plot_id," (",o$year,")")),
        las=1, cex.names=.6, col=ifelse(rev(o$soc_outlier),"#b2182b","#c9922b"),
        main="Highest soc_profile (red = excluded outlier)", xlab="Mg/ha")
if (exists("SOC_OUTLIER_MAX")) abline(v=SOC_OUTLIER_MAX, lty=2, col="grey30")
ch <- data.frame(plot_id=chg$plot_id, d=d2406); ch<-ch[order(-abs(ch$d)),][1:12,]
barplot(rev(ch$d), horiz=TRUE, names.arg=rev(ch$plot_id), las=1, cex.names=.6,
        col=ifelse(rev(ch$d)>0,"#4a9c5d","#b2182b"), main="Largest 2006->2024 change", xlab="Delta SOC (Mg/ha)")
dev.off()

cat("Plots written to", PLOTDIR, "\n")
