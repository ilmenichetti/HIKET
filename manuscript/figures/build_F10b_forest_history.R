setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
# F10b — Finnish forest history ([History] thread): the physical justification for a
# below-equilibrium 1917 start (sigma_init < 1) and rising, non-stationary litter inputs.
# Total growing stock from the National Forest Inventory (VMI1..VMI13), management eras
# shaded, the two SOC campaigns we calibrate against marked.
#
# SOURCE: Korhonen, Raty et al. (2024) "Forests of Finland 2019-2023 and their
# development 1921-2023", Silva Fennica 58(4) art. 24045 (open access).
# Total growing stock on productive + poorly productive forest land, whole country.
# NUMERIC ANCHORS quoted in that paper: NFI1 (1921-24) ~1400 (recalculated with modern
# volume functions; paper: "1.4 G m3", 84% increase to 2552), NFI11 (2009-13) 2356,
# NFI12 (2014-18) 2475, NFI13 (2019-23) 2552 (SE 13, 0.53%). Intermediate NFI2-NFI10
# read from the paper's development curve (Fig. 10) / the standard published NFI series
# (Tomppo/Henttonen), approximate to ~+-30 mill. m3 — narrative accuracy is sufficient
# here; the anchors and the mid-century-low-then-rise shape are the load-bearing content.

# --- NFI total growing stock, mill. m3 (inventory midpoint year) -------------
vmi <- data.frame(
  year  = c(1922, 1937, 1952, 1962, 1967, 1974, 1980, 1990, 2000, 2006, 2011, 2016, 2021),
  stock = c(1400, 1518, 1508, 1524, 1479, 1519, 1660, 1883, 1937, 2189, 2356, 2475, 2552),
  label = c("VMI1","VMI2","VMI3","VMI4","VMI5","VMI6","VMI7","VMI8",
            "VMI9","VMI10","VMI11","VMI12","VMI13")
)

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
text(1917, 1345, "below-equilibrium\n=> sigma_init < 1", cex=0.6, col="#5e35b1", pos=4, offset=0.3)

# --- SOC campaigns we calibrate against -------------------------------------
soc <- data.frame(year=c(1985, 2006), name=c("VMI8 SOC\n(1985-86)", "Biosoil SOC\n(2006)"))
for (i in seq_len(nrow(soc))) {
  abline(v=soc$year[i], col="#00695c", lwd=1.2, lty=2)
  text(soc$year[i], 2600, soc$name[i], cex=0.6, col="#00695c", font=2, pos=2, offset=0.25)
}

title(main="F10b - Finnish forest history: a rising, non-stationary carbon base",
      cex.main=1.0, font.main=1)
mtext("A depleted, near-stationary base into the 1970s, then a sustained ~70% rise -> the non-equilibrium premise & rising litter inputs",
      side=3, line=0.1, cex=0.72, col="grey35")
mtext("Source: Korhonen, Raty et al. (2024) Silva Fennica 58(4) art. 24045 (NFI1-NFI13)",
      side=1, line=3.4, cex=0.6, col="grey45", adj=1)
dev.off()
cat("wrote manuscript/figures/F10b_forest_history.png\n")
