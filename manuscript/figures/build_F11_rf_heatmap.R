setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
# F11 (rebuild) -- RF residual-predictor importance heatmap, with NICELY FORMATTED
# predictor names (the pipeline PNG used raw variable names). Relative importance
# (within-model max = 1), rows = union of each model's top-15, sorted by mean.
# Data: <MODEL>_rf_importance_<RUN_ID>.csv + <MODEL>_rf_summary_<RUN_ID>.rds (OOB R2).

rid <- list(SP1="20260710_104903", TP2="20260710_104904", TP3="20260710_104904",
            Yasso07="20260710_104902", Yasso15="20260710_104902", Yasso20="20260710_102431")
MODELS <- names(rid); DG <- "Calibration_real_data_transient/diagnostics"

# --- nice predictor labels (raw variable -> readable); fallback prettifies ---
lab_map <- c(
  basal_area_85="Stand basal area (1985)", dev_class_85="Development class (1985)",
  TotalNitrogen="Soil total N", CN_ratio="Soil C:N ratio", CEC="Cation exchange capacity",
  base_saturation="Base saturation", ExchangeableAl="Exchangeable Al", ExchangeableCa="Exchangeable Ca",
  ExchangeableFe="Exchangeable Fe", pH.H2O.="Soil pH (H2O)", pH.CaCl2.="Soil pH (CaCl2)",
  clay="Clay content", SiltContent="Silt content", SandContent="Sand content",
  CoarseFragments="Coarse fragments", TexturalClass="Textural class",
  EstimatedBulkDensity="Bulk density", profile_depth_cm="Soil profile depth",
  ofh_weight_kgm2="Organic F/H-layer mass", ofh_lower_cm="Organic layer depth", shallow="Shallow soil",
  lat_WGS84="Latitude", lon_WGS84="Longitude", elevation_m="Elevation",
  mean_temp="Mean annual temperature", warmest_month_T="Warmest-month temperature",
  coldest_month_T="Coldest-month temperature", T_seasonality="Temperature seasonality",
  temp_sum_NFI="Temperature sum (NFI)", temp_zone="Temperature zone", GDD5="Growing degree days (>5C)",
  mean_precip_annual="Mean annual precipitation", P_seasonality="Precipitation seasonality",
  aridity_index="Aridity index", koppen_class="Koppen class",
  woody_share="Woody litter share", conifer_share="Conifer share", species_code="Tree species",
  litter_N_frac="Litter N fraction", litter_A_frac="Litter A fraction",
  litter_W_frac="Litter W fraction", litter_E_frac="Litter E fraction",
  kasvup_tyyppi="Site fertility type", kasvyo_ahti="Vegetation zone (Ahti)",
  kasvyo_syke="Vegetation zone (SYKE)", KA="Site type", soil_code="Soil class",
  n_cuts_85_95="No. cuttings (1985-95)", any_cut_85_95="Any cutting (1985-95)",
  any_trt_85_95="Any treatment (1985-95)", soil_prep_pre85="Soil preparation (pre-1985)",
  ojitustilanne="Drainage status", n_soc_obs="No. SOC observations", alaryhma="Site sub-group")
prettify <- function(v) {
  out <- lab_map[v]
  miss <- is.na(out)
  out[miss] <- vapply(v[miss], function(x){
    s <- gsub("_", " ", x); paste0(toupper(substr(s,1,1)), substr(s,2,nchar(s))) }, character(1))
  unname(out)
}

# --- load, relative importance, top-15 union ---
imp <- list(); oob <- setNames(rep(NA_real_, length(MODELS)), MODELS)
for (m in MODELS) {
  d <- read.csv(file.path(DG, m, sprintf("%s_rf_importance_%s.csv", m, rid[[m]])), stringsAsFactors=FALSE)
  d$rel <- d$importance / max(d$importance, na.rm=TRUE); imp[[m]] <- d
  rf <- file.path(DG, m, sprintf("%s_rf_summary_%s.rds", m, rid[[m]]))
  if (file.exists(rf)) oob[m] <- readRDS(rf)$oob_r2
}
top <- unique(unlist(lapply(imp, function(d) d$variable[order(-d$importance)][1:15])))
mat <- sapply(MODELS, function(m){ d<-imp[[m]]; v<-setNames(d$rel, d$variable)[top]; v[is.na(v)]<-0; v })
rownames(mat) <- top
mat <- mat[order(rowMeans(mat)), , drop=FALSE]            # least (bottom) -> most (top)

# --- render ---
np <- nrow(mat); nm <- ncol(mat)
pal <- colorRampPalette(c("#f7fbff","#deebf7","#9ecae1","#4292c6","#08306b"))(100)
collab <- ifelse(is.na(oob[MODELS]), MODELS, sprintf("%s\nR2=%.2f", MODELS, oob[MODELS]))

png("manuscript/figures/F11_rf_heatmap.png", width=9.5, height=max(6, np*0.34+2), units="in", res=200)
layout(matrix(c(1,2), nrow=1), widths=c(1, 0.16))
par(mar=c(3.4, 14.5, 4.4, 0.6), mgp=c(2.4,0.6,0))
image(1:nm, 1:np, t(mat), col=pal, zlim=c(0,1), axes=FALSE, xlab="", ylab="")
axis(2, at=1:np, labels=prettify(rownames(mat)), las=1, tick=FALSE, cex.axis=0.86)
# all 6 column labels via text() (axis() auto-drops overlapping ones)
text(x=1:nm, y=np+0.5, labels=collab, xpd=NA, cex=0.82, font=1, adj=c(0.5,0))
abline(h=seq(0.5, np+0.5, 1), v=seq(0.5, nm+0.5, 1), col="white", lwd=1.2)
mtext("RF residual-predictor importance (relative, within-model max = 1)", side=1, line=1.6, cex=0.95, font=2)
# colour legend
par(mar=c(3.4, 0.5, 3.6, 3.2))
 leg <- seq(0,1,length.out=100)
image(1, leg, matrix(leg, nrow=1), col=pal, axes=FALSE, xlab="", ylab="")
axis(4, at=seq(0,1,0.25), las=1, cex.axis=0.8); mtext("relative importance", side=4, line=2.0, cex=0.8)
dev.off()
cat("wrote manuscript/figures/F11_rf_heatmap.png  (", np, "predictors x", nm, "models )\n")
cat("top predictor overall:", prettify(rownames(mat)[np]), "\n")
