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
  list(x0=1945, x1=1965, col="#c62828", name="Post-war\nexploitation"),
  list(x0=1965, x1=1990, col="#f9a825", name="Drainage &\nfertilisation"),
  list(x0=1990, x1=2025, col="#2e7d32", name="Managed\ngrowth")
)

# Print size (0.85 x 6.53 in text width, 9 pt), so nothing prints below ~7 pt.
png("manuscript/figures/F10b_forest_history.png", width=5.55, height=3.23, units="in", res=300, pointsize=9)
par(mar=c(3.2, 4.2, 0.6, 0.6), mgp=c(2.0, 0.5, 0), tcl=-0.3)
xl <- c(1915, 2025); yl <- c(1150, 2650)
plot(NA, xlim=xl, ylim=yl, axes=FALSE, xlab="", ylab="")

# era bands
for (e in eras) {
  rect(e$x0, yl[1], e$x1, yl[2], col=adjustcolor(e$col, 0.10), border=NA)
  text(mean(c(e$x0, e$x1)), yl[1] + 90, e$name, cex=0.85, col=e$col, font=2, adj=0.5)
}

axis(1, at=seq(1920, 2020, 20)); axis(2, las=1); box()
mtext("Year", side=1, line=2.0)
mtext(expression("Growing stock (million m"^3*")"), side=2, line=3.0)

# growing-stock trajectory
lines(vmi$year, vmi$stock, col="#37474f", lwd=2.2)
points(vmi$year, vmi$stock, pch=21, bg="#37474f", col="white", cex=1.3, lwd=1.2)
# label the endpoints and the mid-century low
key <- vmi$label %in% c("VMI1","VMI5","VMI8","VMI13")
text(vmi$year[key], vmi$stock[key], vmi$label[key], pos=c(3,3,2,2), cex=0.85,
     col="#37474f", offset=0.5)

# --- 1917 below-equilibrium anchor (the calibration pre-init year) ----------
abline(v=1917, col="#5e35b1", lwd=1.4, lty=3)
text(1917, yl[2]-40, "1917\npre-run start", cex=0.85, col="#5e35b1", font=2, pos=4, offset=0.3)
arrows(1917, 1560, 1917, 1400, length=0.08, col="#5e35b1", lwd=1.6)

# --- SOC campaigns we calibrate against -------------------------------------
soc <- data.frame(year=c(1985, 2006), name=c("VMI8 SOC\n(1985-86)", "Biosoil SOC\n(2006)"))
for (i in seq_len(nrow(soc))) {
  abline(v=soc$year[i], col="#00695c", lwd=1.2, lty=2)
  text(soc$year[i], 2600, soc$name[i], cex=0.85, col="#00695c", font=2, pos=2, offset=0.25)
}

dev.off()
cat("wrote manuscript/figures/F10b_forest_history.png\n")
