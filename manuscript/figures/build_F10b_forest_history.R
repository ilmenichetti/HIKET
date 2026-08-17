setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
# F10b — Finnish forest history ([History] thread): the physical justification for a
# below-equilibrium 1917 start (sigma_init < 1) and rising, non-stationary litter inputs.
# Total growing stock from the National Forest Inventory (VMI1..VMI13), management eras
# shaded, the two SOC campaigns we calibrate against marked.
#
# SOURCE: Korhonen, Raty et al. (2024) "Forests of Finland 2019-2023 and their
# development 1921-2023", Silva Fennica 58(5) art. 24045 (open access). Total growing
# stock on productive + poorly productive forest land, whole country. Series now read
# from the single sourced data file (endpoints exact from the paper; NFI2-NFI10
# digitized from Fig 10a) so this figure and the C3 pre-run share ONE source of truth.
# Provenance per point: Data/forest_history/nfi_growing_stock.csv + its README.

# --- NFI total growing stock, mill. m3 (inventory midpoint year) -------------
gs  <- read.csv("Data/forest_history/nfi_growing_stock.csv")
vmi <- data.frame(year  = gs$midpoint_year,
                  stock = gs$total_stock_Mm3,
                  label = paste0("VMI", gs$nfi))

# --- management eras (shading) ----------------------------------------------
eras <- list(
  list(x0=1945, x1=1965, col="#c62828", name="Post-war exploitation\n& war reparations"),
  list(x0=1965, x1=1990, col="#f9a825", name="Forest-improvement era\n(drainage, fertilization)"),
  list(x0=1990, x1=2025, col="#2e7d32", name="Modern managed growth")
)

png("manuscript/figures/F10b_forest_history.png", width=10, height=5.2, units="in", res=200)
par(mar=c(4.4, 4.8, 3.2, 1.2), mgp=c(2.7, 0.7, 0))
xl <- c(1915, 2025); yl <- c(1300, 2650)
plot(NA, xlim=xl, ylim=yl, axes=FALSE, xlab="", ylab="")

# era bands
for (e in eras) {
  rect(e$x0, yl[1], e$x1, yl[2], col=adjustcolor(e$col, 0.10), border=NA)
  text(mean(c(e$x0, e$x1)), yl[1] + 60, e$name, cex=0.66, col=e$col, font=2, adj=0.5)
}

axis(1, at=seq(1920, 2020, 20)); axis(2, las=1); box()
mtext("Year", side=1, line=2.6, cex=0.98)
mtext("Total growing stock (mill. m3)", side=2, line=3.2, cex=0.98)

# growing-stock trajectory
lines(vmi$year, vmi$stock, col="#37474f", lwd=2.2)
points(vmi$year, vmi$stock, pch=21, bg="#37474f", col="white", cex=1.3, lwd=1.2)
# label the endpoints and the mid-century low
key <- vmi$label %in% c("VMI1","VMI5","VMI8","VMI13")
text(vmi$year[key], vmi$stock[key], vmi$label[key], pos=c(4,1,2,2), cex=0.62,
     col="#37474f", offset=0.5)

# --- 1917 below-equilibrium anchor (the calibration pre-init year) ----------
abline(v=1917, col="#5e35b1", lwd=1.4, lty=3)
text(1917, yl[2]-40, "1917\npre-init start", cex=0.62, col="#5e35b1", font=2, pos=4, offset=0.3)
arrows(1917, 1560, 1917, 1360, length=0.08, col="#5e35b1", lwd=1.6)

# --- SOC campaigns we calibrate against -------------------------------------
soc <- data.frame(year=c(1985, 2006), name=c("VMI8 SOC\n(1985-86)", "Biosoil SOC\n(2006)"))
for (i in seq_len(nrow(soc))) {
  abline(v=soc$year[i], col="#00695c", lwd=1.2, lty=2)
  text(soc$year[i], 2600, soc$name[i], cex=0.6, col="#00695c", font=2, pos=2, offset=0.25)
}

title(main="Finnish forest history: growing stock",
      cex.main=1.0, font.main=1)
mtext("Source: Korhonen, Raty et al. (2024) Silva Fennica 58(4) art. 24045 (NFI1-NFI13)",
      side=1, line=3.4, cex=0.6, col="grey45", adj=1)
dev.off()
cat("wrote manuscript/figures/F10b_forest_history.png\n")
