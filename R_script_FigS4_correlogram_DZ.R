# =============================================================================
# Fig. S4 - Pearson correlation matrix of environmental variables, US-Brw
# Base R only (no extra packages).
#
# Input : TOT_MET_FLX_v2026_11_13_f.csv  (annual values, one row per year)
# Output: FigS4_correlogram.tiff (600 dpi, LZW), FigS4_correlogram.pdf,
#         FigS4_correlogram_stats.csv (r, p and n for every pair)
#
# Correlations use pairwise complete observations (years with missing values
# are dropped only for the pairs involving that variable).
# =============================================================================

infile <- "TOT_MET_FLX_v2026_11_13_f.csv"
d <- read.csv(infile, check.names = FALSE)

# variables (column name = label), in plotting order
vars <- c(
  "SoilT_May"               = "Soil T May",
  "SoilT_June"              = "Soil T June",
  "SoilT_July"              = "Soil T July",
  "SoilT_Aug"               = "Soil T Aug",
  "SoilT_Sept"              = "Soil T Sept",
  "melt_date_NOAA"          = "Melt date",
  "Snow_In_NOAA"            = "Snow-in date",
  "Snow_free_NOAA"          = "Snow-free date",
  "Snow Season Length_NOAA" = "Snow-season length",
  "SWC_May"                 = "SWC May",
  "SWC_June"                = "SWC June",
  "SWC_July"                = "SWC July",
  "SWC_Aug"                 = "SWC Aug",
  "SWC_Sept"                = "SWC Sept",
  "SWC_Oct"                 = "SWC Oct",
  "GDD_June"                = "GDD June",
  "GDD_July"                = "GDD July",
  "GDD_Aug"                 = "GDD Aug",
  "GDD_Sept"                = "GDD Sept",
  "U1_Barrow_mean"          = "Active layer depth"
)
X <- d[, names(vars)]
k <- length(vars)
lab <- unname(vars)

# ---- correlations ------------------------------------------------------------
R <- P <- matrix(NA_real_, k, k)
N <- matrix(0L, k, k)
for (i in 1:k) for (j in i:k) {
  ok <- complete.cases(X[, i], X[, j])
  N[i, j] <- sum(ok)
  if (i == j) { R[i, j] <- 1; P[i, j] <- 0; next }
  ct <- cor.test(X[ok, i], X[ok, j], method = "pearson")
  R[i, j] <- unname(ct$estimate); P[i, j] <- ct$p.value
}

stars <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))

# ---- stats table ---------------------------------------------------------------
idx <- which(upper.tri(R), arr.ind = TRUE)
idx <- idx[order(idx[, 1], idx[, 2]), ]
stats_out <- data.frame(var1 = lab[idx[, 1]], var2 = lab[idx[, 2]],
                        r = round(R[idx], 3), p = signif(P[idx], 3), n = N[idx])
write.csv(stats_out, "FigS4_correlogram_stats.csv", row.names = FALSE)

# ---- colour scale (blue-white-red, -1 to 1) --------------------------------------
pal <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
                          "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B",
                          "#67001F"))(201)
col_of <- function(r) pal[round((r + 1) * 100) + 1]

draw_fig <- function() {
  layout(matrix(1:2, 1), widths = c(1, 0.12))
  # main panel
  par(mar = c(9, 10, 3.5, 0.5), xpd = NA)
  plot(NA, xlim = c(0.5, k + 0.5), ylim = c(k + 0.5, 0.5), axes = FALSE,
       xlab = "", ylab = "", asp = 1, xaxs = "i", yaxs = "i")
  for (i in 1:k) for (j in i:k) {
    rect(j - 0.5, i - 0.5, j + 0.5, i + 0.5, col = col_of(R[i, j]), border = "white", lwd = 0.6)
    text(j, i, paste0(sprintf("%.2f", R[i, j]), stars(P[i, j])), cex = 0.62,
         col = ifelse(abs(R[i, j]) > 0.6, "white", "black"))
  }
  rect(0.5, 0.5, k + 0.5, k + 0.5, border = "black", lwd = 1)
  text(0.3, 1:k, lab, adj = 1, cex = 1)                         # row labels
  text(1:k, k + 0.75, lab, srt = 55, adj = c(1, 1), cex = 1)    # column labels
  title(main = "Pearson correlations", cex.main = 1.6, line = 1.5)
  mtext("* p < 0.05; ** p < 0.01; *** p < 0.001", side = 1, line = 7.5, cex = 1)
  # colour bar
  par(mar = c(9, 0.6, 3.5, 3.8), xpd = NA)
  z <- seq(-1, 1, length.out = 201)
  image(1, z, matrix(z, nrow = 1), col = pal, axes = FALSE, xlab = "", ylab = "")
  box(); axis(4, at = seq(-1, 1, 0.25), las = 1, cex.axis = 0.9)
  mtext("Correlation coefficient", side = 4, line = 2.6, cex = 1)
}

tiff("FigS4_correlogram.tiff", width = 11, height = 10, units = "in", res = 600,
     compression = "lzw", type = "cairo", pointsize = 12)
draw_fig(); dev.off()
pdf("FigS4_correlogram.pdf", width = 11, height = 10, pointsize = 12)
draw_fig(); dev.off()

cat("Done. Wrote FigS4_correlogram.tiff/.pdf and FigS4_correlogram_stats.csv\n")
