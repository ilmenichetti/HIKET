source("manuscript/figures/run_ids.R")   # auto-selects current RUN_IDs
setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
# Robustness R4 -- KL(posterior || prior) for the SIMPLE models SP1/TP2/TP3
# (the appendix half of F7, which covers only the operational Yassos). Same faithful
# recompute: saved posterior is physical; prior reproduced via each model's engine
# setup (to_original/best_x/sigma_ppm); engine classify_param/class_cols for colours.
rid <- unlist(RID[c("SP1","TP2","TP3")])

setup_model <- function(MODEL) {
  sp  <- file.path("Calibration_real_data_transient", sprintf("run_%s_transient_calibration.R", MODEL))
  src <- readLines(sp, warn = FALSE); cut <- grep("^t_run <- system.time\\(\\{", src)[1]
  e   <- new.env(parent = globalenv())
  suppressWarnings(suppressMessages(
    source(textConnection(paste(src[seq_len(cut - 1L)], collapse = "\n")), local = e)))
  best_x <- get("best_x", e); sigma_ppm <- get("sigma_ppm", e); to_original <- get("to_original", e)
  set.seed(99)
  pr <- sapply(seq_along(best_x), function(j) rnorm(3000L, best_x[j], sigma_ppm[j]))
  colnames(pr) <- names(best_x); prior <- t(apply(pr, 1, to_original))
  post <- readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_%s.rds", MODEL, rid[[MODEL]]))
  list(FREE=get("FREE_NAMES",e), prior=prior, post=post,
       classify=get("classify_param",e), ccols=get("class_cols",e), clev=get("class_levels",e))
}
kl_from_kde <- function(pv, qv, n=512L){ pv<-pv[is.finite(pv)]; qv<-qv[is.finite(qv)]
  if(length(pv)<10||length(qv)<10) return(NA_real_)
  xl<-c(min(quantile(qv,.005),min(pv)), max(quantile(qv,.995),max(pv)))
  dp<-density(pv,from=xl[1],to=xl[2],n=n); dq<-density(qv,from=xl[1],to=xl[2],n=n)
  sum(pmax(dp$y,1e-10)*log(pmax(dp$y,1e-10)/pmax(dq$y,1e-10)))*diff(dp$x[1:2]) }

cat("sourcing SP1/TP2/TP3 setups...\n")
M <- lapply(names(rid), setup_model); names(M) <- names(rid)

png("manuscript/figures/R4_nonyasso_kl.png", width=10, height=8, units="in", res=200)
par(mfrow=c(3,1), mar=c(6.2,4.8,2.4,0.8), mgp=c(2.9,0.7,0), oma=c(1.5,0,0,0), cex.axis=1.0, cex.lab=1.15)
for (mn in names(rid)) {
  mm <- M[[mn]]; kv <- vapply(mm$FREE, function(nm) kl_from_kde(mm$post[,nm], mm$prior[,nm]), numeric(1))
  cl <- vapply(names(kv), mm$classify, character(1)); bcol <- mm$ccols[cl]
  barplot(kv, col=bcol, border=NA, las=2, cex.names=1.0, ylim=c(0, max(kv,na.rm=TRUE)*1.15),
          ylab="KL (nats)", main=mn, font.main=1, cex.main=1.3)
  abline(h=1, lty=2, col="grey40", lwd=1.3)
  present <- mm$clev[mm$clev %in% cl]
  legend("topright", legend=present, fill=mm$ccols[present], border=NA, bty="n", cex=1.0, title="Parameter class")
}
mtext("KL by parameter, simple models; note y-axes differ",
      side=1, line=0.2, outer=TRUE, cex=0.8, col="grey30")
dev.off(); cat("wrote R4_nonyasso_kl.png\n")
