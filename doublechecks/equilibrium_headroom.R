# =============================================================================
# equilibrium_headroom.R  (2026-09-04)
#
# WHY. Liski et al. (2006) state that Finnish forest soils would keep accumulating
# for "centuries" and settle "at a 38% higher level than in 1922". Our own models
# must be reconciled with that (annotation, storyline v3 p7). This computes the
# same quantity from OUR posteriors: how far below its own equilibrium does each
# model put the soil, at the CURRENT litter flux?
#
# METHOD, using the calibration's own initialiser and no re-implementation:
#   modelled 1985 state   steady_state(mp, lm, xi)                as calibrated
#   equilibrium at J_1985 same call with sigma_init <- 1 and preinit_shape <- 0,
#                         which holds the flux at J_t0*sigma_input = J_1985 for the
#                         whole pre-run, so the routine returns ITS OWN equilibrium.
# Headroom = C_eq/C_1985 - 1, the accumulation still owed at today's input.
# =============================================================================
suppressWarnings(suppressMessages({
  a <- commandArgs(trailingOnly = TRUE)
  N_DRAW <- if (length(a) >= 1) as.integer(a[[1]]) else 5L
  library(BayesianTools)
}))
set.seed(2025)
MODELS <- strsplit(Sys.getenv("HIKET_MODELS","SP1,TP2,TP3,Yasso07,Yasso15,Yasso20"), ",")[[1]]
DIR <- "Calibration_real_data_transient/runs"; NC <- max(1L, parallel::detectCores()-1L)
rid <- function(m){ f <- list.files(DIR, pattern=sprintf("^%s_posterior_[0-9]{8}_[0-9]{6}\\.rds$",m))
  if(!length(f)) NA_character_ else sub(sprintf("^%s_posterior_(.+)\\.rds$",m),"\\1",sort(f,decreasing=TRUE)[1]) }
setup <- function(M){
  src <- readLines(sprintf("Calibration_real_data_transient/run_%s_transient_calibration.R",M),warn=FALSE)
  cut <- grep("^t_run <- system.time",src)[1]; e <- new.env(parent=globalenv())
  invisible(capture.output(suppressMessages(source(textConnection(paste(src[seq_len(cut-1)],collapse="\n")),local=e))))
  st <- grep("ll_fn <- make_likelihood",src)[1]; op <- 0L; en <- NA_integer_
  for(i in seq(st,length(src))){ ch <- strsplit(src[i],"")[[1]]
    op <- op+sum(ch=="(")-sum(ch==")"); if(op==0L){en <- i;break} }
  ml <- as.list(str2lang(sub("^\\s*ll_fn\\s*<-\\s*","",paste(src[seq(st,en)],collapse="\n"))))[-1]
  for(k in c("to_original","assemble_params","compute_xi_mean","steady_state"))
    assign(paste0(".",k), eval(ml[[k]],envir=e), envir=e)
  e
}
cat("\n=== How far below equilibrium do our own models put the soil? ===\n")
cat("Liski et al. 2006, same country and model family: \"centuries later ... a 38% higher level than in 1922\"\n\n")
cat(sprintf("%-8s %10s %10s %9s\n","model","C_1985","C_eq(J_85)","headroom"))
out <- list()
for(M in MODELS){
  r <- rid(M); if(is.na(r)) next
  e <- tryCatch(setup(M), error=function(z) NULL); if(is.null(e)) next
  ch <- file.path(DIR,sprintf("%s_chains_%s.rds",M,r)); if(!file.exists(ch)) next
  s_ <- do.call(rbind, lapply(readRDS(ch), function(z) getSample(z,parametersOnly=FALSE,start=2)))
  free <- names(get("best_x",e)); plots <- get("plots",e)
  cbp <- get("climate_by_plot",e); lms <- get("litter_means",e)
  H <- C0 <- CE <- numeric(0)
  for(k in sample(nrow(s_), min(N_DRAW,nrow(s_)))){
    pf <- tryCatch(e$.to_original(s_[k,free]), error=function(z) NULL); if(is.null(pf)) next
    mp <- e$.assemble_params(pf)
    pf2 <- pf; pf2["sigma_init"] <- 1; mp2 <- e$.assemble_params(pf2)
    v <- do.call(rbind, parallel::mclapply(plots, function(pid){
      xi <- tryCatch(e$.compute_xi_mean(cbp[[pid]],mp), error=function(z) NULL); if(is.null(xi)) return(NULL)
      l <- lms[[pid]]; leq <- l; leq$preinit_shape <- rep(0, length(l$preinit_shape))
      c1 <- tryCatch(sum(e$.steady_state(mp ,l  ,xi)), error=function(z) NA_real_)
      c2 <- tryCatch(sum(e$.steady_state(mp2,leq,xi)), error=function(z) NA_real_)
      if(!is.finite(c1)||!is.finite(c2)) return(NULL); c(c1,c2)
    }, mc.cores=NC))
    if(is.null(v)||!nrow(v)) next
    C0 <- c(C0, mean(v[,1])); CE <- c(CE, mean(v[,2])); H <- c(H, mean(v[,2])/mean(v[,1])-1)
  }
  if(!length(H)) next
  cat(sprintf("%-8s %10.1f %10.1f %8.0f%%\n", M, median(C0), median(CE), 100*median(H)))
  out[[M]] <- data.frame(model=M, run_id=r, n=length(H), C_1985=median(C0),
                         C_eq=median(CE), headroom=median(H))
}
if(length(out)){ df <- do.call(rbind,out); saveRDS(df,"doublechecks/equilibrium_headroom.rds")
  cat("\nsaved doublechecks/equilibrium_headroom.rds\n") }
