# =============================================================================
# soc_vs_input_proportionality.R   (2026-11-05)
#
# THE QUESTION. Is each measured SOC stock proportional to the litter input?
# Steady state gives C = J x MTT, so log SOC = a + b log J with b = 1.
#
# WHAT IS COMPUTED, on the calibrated plots (raw Tupek J, tree litter only; sigma_input
# would only shift the intercept):
#   A  b per campaign x input window (1, 5, 10, 20 yr, all), with and without mean T
#   C  between-plot slope vs within-plot slope (change in SOC vs change in J10)
#   D  b per layer (organic, mineral 0-40, deep tail, profile), 2006 and 2024
#   E  the same regressions on the models' predictions, as the reference
#
# CAVEATS. The input series starts in 1985, so the 1986-95 windows are TRUNCATED.
# A window shorter than the rotation sees the stand stage, not the long-run input,
# so b << 1 rejects plot-level proportionality to THIS product over THIS window,
# not C = J x MTT in the long run. Error in J also attenuates b.
# Run from the repo root.
# =============================================================================
source("manuscript/figures/run_ids.R")
S <- "doublechecks/figures"
d <- read.csv("Data/model_inputs/input_raw_monthly.csv")
cc <- grep("^C_(nwl|fwl|cwl)_", names(d), value=TRUE); d$J <- rowSums(d[,cc])
ann <- aggregate(J ~ plot_id+year, d, sum)
Jm <- with(ann, tapply(J, list(plot_id, year), sum))       # plot x year
Tm <- tapply(d$temp_air, d$plot_id, mean, na.rm=TRUE)
Jw <- function(p, y, w) { yrs <- as.character(max(1985,y-w+1):y)
  mean(Jm[as.character(p), yrs], na.rm=TRUE) }
rd_all <- lapply(FIG_MODELS, function(m){ b <- readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_predictive_%s.rds", m, RID[[m]]))
  x <- b$residuals_df[,c("plot_id","year","soc_obs_tCha","soc_median")]; names(x)[4] <- m; x})
o <- Reduce(function(a,b) merge(a,b[,c("plot_id","year",names(b)[4])],by=c("plot_id","year")), rd_all)
o <- o[is.finite(o$soc_obs_tCha)&o$soc_obs_tCha>0,]
o$camp <- ifelse(o$year<2000,"1986-95",as.character(o$year))
o$T <- Tm[as.character(o$plot_id)]
for (w in c(1,5,10,20)) o[[paste0("J",w)]] <- mapply(Jw,o$plot_id,o$year,MoreArgs=list(w=w))
o$Jall <- mapply(Jw,o$plot_id,o$year,MoreArgs=list(w=100))
o <- o[o$J10>0 & o$J1>0,]
cat("obs",nrow(o),"plots",length(unique(o$plot_id)),"\n"); print(table(o$camp))
cat("J10 quantiles 5/50/95 by campaign\n"); print(sapply(split(o$J10,o$camp),quantile,c(.05,.5,.95)))
# cluster-robust slope
crob <- function(f,dat,term){ m<-lm(f,dat); X<-model.matrix(m); e<-resid(m); cl<-dat[rownames(X),"plot_id"]
  B<-solve(crossprod(X)); U<-rowsum(X*e,cl); V<-B%*%crossprod(U)%*%B; G<-length(unique(cl)); V<-V*G/(G-1)
  b<-coef(m)[term]; se<-sqrt(V[term,term]); c(b=unname(b),lo=unname(b-1.96*se),hi=unname(b+1.96*se),r2=summary(m)$r.squared,n=nrow(X)) }
cat("\n== A. b by campaign x window (no covariate / + mean T) ==\n")
A<-do.call(rbind,lapply(c("1986-95","2006","2024"),function(cp) do.call(rbind,lapply(c("J1","J5","J10","J20","Jall"),function(w){ x<-o[o$camp==cp,]
 a<-crob(as.formula(sprintf("log(soc_obs_tCha)~log(%s)",w)),x,sprintf("log(%s)",w)); t<-crob(as.formula(sprintf("log(soc_obs_tCha)~log(%s)+T",w)),x,sprintf("log(%s)",w))
 data.frame(camp=cp,w=w,b=a["b"],lo=a["lo"],hi=a["hi"],r2=a["r2"],n=a["n"],bT=t["b"],loT=t["lo"],hiT=t["hi"])}))))
rownames(A)<-NULL; print(format(A,digits=3)); write.csv(A,file.path(S,"soc_vs_input_slopes.csv"),row.names=FALSE)
cat("\n== pooled, campaign FE, J10, cluster SE ==\n")
print(round(crob(log(soc_obs_tCha)~log(J10)+camp,o,"log(J10)"),3))
print(round(crob(log(soc_obs_tCha)~log(J10)+camp+T,o,"log(J10)"),3))
cat("level form SOC = a + b J10 (2006,2024): \n"); for(cp in c("2006","2024")) print(round(summary(lm(soc_obs_tCha~J10,o[o$camp==cp,]))$coef[,1:2],2))
cat("\n== C. between vs within ==\n")
bt<-aggregate(cbind(ls=log(soc_obs_tCha),lj=log(J10))~plot_id,o,mean); print(round(c(between=unname(coef(lm(ls~lj,bt))[2]),confint(lm(ls~lj,bt))[2,]),3))
wd<-function(a,b){x<-merge(o[o$camp==a,],o[o$camp==b,],by="plot_id"); x$ds<-log(x$soc_obs_tCha.y/x$soc_obs_tCha.x); x$dj<-log(x$J10.y/x$J10.x); f<-lm(ds~dj,x)
 cat(a,"->",b,"n",nrow(x),"sd dlogJ",round(sd(x$dj),3),"b",round(coef(f)[2],3),"CI",round(confint(f)[2,],3),"r2",round(summary(f)$r.squared,3),"\n")}
wd("2006","2024"); wd("1986-95","2006"); wd("1986-95","2024")
cat("\n== D. by layer, J10 ==\n")
s<-read.csv("Data/SOC_homogeneized/soc_homogenized_plot.csv"); print(unique(s$campaign)); print(unique(s$year))
s<-s[s$plot_id%in%o$plot_id & s$year%in%c(2006,2024),]; s<-merge(s,o[,c("plot_id","year","J10","T","soc_obs_tCha")],by=c("plot_id","year"))
for(v in c("organic_Mgha","mineral_0_40_Mgha","soc_deep_Mgha","soc_profile_Mgha")) for(yy in c(2006,2024)){x<-s[s$year==yy&is.finite(s[[v]])&s[[v]]>0,]
 f<-crob(as.formula(sprintf("log(%s)~log(J10)",v)),x,"log(J10)"); cat(v,yy,"b",round(f["b"],3),"CI",round(f["lo"],3),round(f["hi"],3),"r2",round(f["r2"],3),"n",f["n"],"\n")}
cat("check soc_profile vs obs: ",round(cor(s$soc_profile_Mgha,s$soc_obs_tCha),4),"\n")
cat("\n== E. same regression on model predictions (J10) ==\n")
for(m in FIG_MODELS) cat(m, sapply(c("1986-95","2006","2024"),function(cp){x<-o[o$camp==cp,];round(coef(lm(log(x[[m]])~log(x$J10)))[2],3)}),"\n")
cat("reliability of log J10 needed for true b=1:", round(A$b[A$w=="J10"],3),"\n")
print(sapply(split(o,o$camp),function(x) c(sd_logJ10=sd(log(x$J10)),sd_logSOC=sd(log(x$soc_obs_tCha)),spearman=cor(x$J10,x$soc_obs_tCha,method="s"))))
# === Extra: temperature sign, within-plot slopes of the models ===
cat("\n== T coefficient ==\n")
for(cp in c("2006","2024")){x<-o[o$camp==cp,]; print(round(summary(lm(log(soc_obs_tCha)~log(J10)+T,x))$coef[,c(1,2,4)],4)); cat("cor logJ10,T",round(cor(log(x$J10),x$T),3),"\n")}
x<-merge(o[o$camp=="2006",],o[o$camp=="2024",],by="plot_id"); dj<-log(x$J10.y/x$J10.x)
cat("within 2006-2024 model slopes:\n"); for(m in FIG_MODELS) cat(m,round(coef(lm(log(x[[paste0(m,".y")]]/x[[paste0(m,".x")]])~dj))[2],3),"\n")
# === North / South: LUKE's standard split (~7000 km northing, ~63 N), region.y in
# Boris's tree-litter file. NOT region_label / site_raw$region: those split at ~66 N.
# Region only is read from that file; J is still the current model input.
REG_FILE <- "Data/LitterData/tree_litter_per_site_year_by_component_30.10.26.csv"
reg <- read.csv(pipe(sprintf("awk -F';' 'NR>1{print $13\";\"$NF}' '%s' | sort -u", REG_FILE)),
                sep=";", header=FALSE, col.names=c("plot_id","region_y"), na.strings=c("NA",""))
stopifnot(!anyDuplicated(reg$plot_id[!is.na(reg$plot_id)]) || all(tapply(reg$region_y,reg$plot_id,function(z) length(unique(z)))==1))
o$region <- c("1"="South","2"="North")[as.character(reg$region_y[match(o$plot_id, reg$plot_id)])]
stopifnot(!anyNA(o$region))
cat("\ncalibrated plots by region:\n"); print(table(unique(o[,c("plot_id","region")])$region))
cat("\n== b by campaign x region, J10 (no covariate / + mean T) ==\n")
R <- do.call(rbind, lapply(c("1986-95","2006","2024"), function(cp) do.call(rbind, lapply(sort(unique(na.omit(o$region))), function(g){
  x <- o[o$camp==cp & o$region %in% g,]
  a <- crob(log(soc_obs_tCha)~log(J10),x,"log(J10)"); t <- crob(log(soc_obs_tCha)~log(J10)+T,x,"log(J10)")
  data.frame(camp=cp,region=g,n=a[["n"]],b=a[["b"]],lo=a[["lo"]],hi=a[["hi"]],r2=a[["r2"]],bT=t[["b"]],loT=t[["lo"]],hiT=t[["hi"]],
             medJ=median(x$J10),medSOC=median(x$soc_obs_tCha))}))))
rownames(R) <- NULL; print(format(R,digits=3)); write.csv(R,file.path(S,"soc_vs_input_slopes_region.csv"),row.names=FALSE)

cat("\n== within-plot 2006->2024 slope by region, J10 ==\n")
for(g in c("South","North")){z<-merge(o[o$camp=="2006"&o$region==g,],o[o$camp=="2024"&o$region==g,],by="plot_id"); dj<-log(z$J10.y/z$J10.x); f<-lm(log(z$soc_obs_tCha.y/z$soc_obs_tCha.x)~dj)
 cat(g,"n",nrow(z),"obs b",round(coef(f)[2],3),"CI",round(confint(f)[2,],3),"| models",sapply(FIG_MODELS,function(m) round(coef(lm(log(z[[paste0(m,".y")]]/z[[paste0(m,".x")]])~dj))[2],2)),"\n")}
cat("\n== models by region, J10 ==\n")
for(g in c("South","North")) for(cp in c("1986-95","2006","2024")){z<-o[o$camp==cp&o$region==g,]; cat(g,cp,sapply(FIG_MODELS,function(m) round(coef(lm(log(z[[m]])~log(z$J10)))[2],2)),"\n")}

# === Figure ===
regs <- sort(unique(na.omit(o$region))); cols <- setNames(c("#2166AC","#D9730D"), regs[order(!grepl("north",regs,ignore.case=TRUE))])
png(file.path(S,"soc_vs_input.png"),width=2000,height=1800,res=300,pointsize=8)
par(mfrow=c(2,2),mar=c(3.2,3.2,1.5,.5),mgp=c(2,.5,0),tcl=-.25,las=1,cex=1)
xl<-range(o$J10);yl<-range(o$soc_obs_tCha)
for(cp in c("1986-95","2006","2024")){x<-o[o$camp==cp & !is.na(o$region),];x<-x[sample(nrow(x)),];lx<-log(x$J10);ly<-log(x$soc_obs_tCha);f<-lm(ly~lx)
 plot(x$J10,x$soc_obs_tCha,log="xy",xlim=xl,ylim=yl,pch=16,cex=.55,col=adjustcolor(cols[x$region],.55),bty="l",
  xlab=expression("Mean input, past 10 yr (tC ha"^-1*" yr"^-1*")"),ylab=expression("SOC (tC ha"^-1*")"))
 xs<-seq(min(lx),max(lx),len=50); lines(exp(xs),exp(mean(ly)+xs-mean(lx)),lwd=1.2,lty=2,col="grey25")
 lines(exp(xs),exp(predict(f,data.frame(lx=xs))),lwd=1.5,col="black")
 for(g in regs){i<-x$region==g;fg<-lm(ly[i]~lx[i]);xg<-range(lx[i]);lines(exp(xg),exp(coef(fg)[1]+coef(fg)[2]*xg),lwd=2,col=cols[g])}
 bb<-R[R$camp==cp,]; mtext(sprintf("%s  all %.2f, %s",cp,coef(f)[2],paste(sprintf("%s %.2f",bb$region,bb$b),collapse=", ")),3,.2,adj=0,cex=.9)
 if(cp=="1986-95") legend("topleft",c(names(cols),"all plots","proportional (b = 1)"),pch=c(16,16,NA,NA),lty=c(1,1,1,2),lwd=c(2,2,1.5,1.2),col=c(cols,"black","grey25"),bty="n")}
ww<-c(J1=1,J5=5,J10=10,J20=20); B<-A[A$w%in%names(ww)&A$camp!="1986-95",]; B$x<-log(ww[B$w])+ifelse(B$camp=="2006",-.05,.05)
plot(B$x,B$b,ylim=c(-.05,1.12),xaxt="n",type="n",bty="l",xlab="Input window (years)",ylab="Slope b, all plots"); axis(1,log(ww),ww)
abline(h=c(0,1),lty=c(3,2),col=c("grey60","grey25"))
for(cp in c("2006","2024")){z<-B[B$camp==cp,]; pc<-if(cp=="2006")16 else 1; segments(z$x,z$lo,z$x,z$hi); points(z$x,z$b,pch=pc)}
legend("right",c("2006","2024"),pch=c(16,1),bty="n"); text(log(1),1.07,"proportional",adj=0,col="grey25",cex=.9)
dev.off()
