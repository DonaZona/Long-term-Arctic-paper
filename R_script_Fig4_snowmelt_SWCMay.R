# =============================================================================
# Fig. 4: Snowmelt date vs May soil water content (SWC), US-Brw (Utqiagvik), 2001-2024
#
# (a) observed values; (b) anomalies after removing the linear trend from each
# variable over its own full record. Pearson correlation with 95% CI band;
# sensitivity tests excluding 2016 and 2015-2016, and Spearman rank correlation.
#
# Outputs (written to the working directory):
#   Fig_4_snowmelt_SWCMay.tiff (600 dpi, LZW), Fig_4_snowmelt_SWCMay.pdf,
#   Fig_4_snowmelt_SWCMay_stats.csv
# =============================================================================

dat <- read.csv("TOT_MET_FLX_v2026_11_13_f.csv", check.names = FALSE)

# ---- detrend each variable over its own full record (residuals of y ~ Year)
detrend <- function(year, y) {
  ok <- !is.na(y)
  out <- rep(NA_real_, length(y))
  out[ok] <- resid(lm(y[ok] ~ year[ok]))
  out
}
dat$melt_detr <- detrend(dat$Year, dat$melt_date_NOAA)
dat$swc_detr  <- detrend(dat$Year, dat$SWC_May)

raw <- na.omit(dat[, c("Year", "melt_date_NOAA", "SWC_May")])
det <- na.omit(dat[, c("Year", "melt_detr", "swc_detr")])

# ---- correlations: all years, and sensitivity to the 2015-2016 extremes
cor_row <- function(df, x, y, label) {
  pr <- cor.test(df[[x]], df[[y]], method = "pearson")
  sp <- suppressWarnings(cor.test(df[[x]], df[[y]], method = "spearman"))
  data.frame(Analysis = label, n = nrow(df),
             Pearson_r = unname(pr$estimate), Pearson_p = pr$p.value,
             R2 = unname(pr$estimate)^2,
             Spearman_rho = unname(sp$estimate), Spearman_p = sp$p.value)
}
stats <- rbind(
  cor_row(raw, "melt_date_NOAA", "SWC_May", "Observed, all years"),
  cor_row(det, "melt_detr", "swc_detr", "Detrended, all years"),
  cor_row(det[det$Year != 2016, ], "melt_detr", "swc_detr", "Detrended, excluding 2016"),
  cor_row(det[!det$Year %in% c(2015, 2016), ], "melt_detr", "swc_detr", "Detrended, excluding 2015-2016")
)
write.csv(stats, "Fig_4_snowmelt_SWCMay_stats.csv", row.names = FALSE)
print(format(stats, digits = 3), row.names = FALSE)

# ---- figure
col_pt <- "#1f77b4"; col_hi <- "#d62728"
fmt_p <- function(p) if (p < 0.001) "p < 0.001" else sprintf("p = %.2f", p)

draw_panel <- function(x, y, years, title, xlab, ylab, zero_lines = FALSE) {
  fit <- lm(y ~ x)
  ct  <- cor.test(x, y)
  xx  <- seq(min(x), max(x), length.out = 100)
  pr  <- predict(fit, newdata = data.frame(x = xx), interval = "confidence")
  plot(x, y, type = "n", bty = "l", xlab = xlab, ylab = ylab,
       ylim = range(c(y, pr[, 2:3])))
  grid(col = "grey90", lty = 1)
  if (zero_lines) abline(h = 0, v = 0, col = "grey60")
  polygon(c(xx, rev(xx)), c(pr[, "lwr"], rev(pr[, "upr"])), col = "grey85", border = NA)
  lines(xx, pr[, "fit"], lwd = 3, lty = if (ct$p.value < 0.05) 1 else 2, col = "grey20")
  hi <- years %in% c(2015, 2016)
  points(x[!hi], y[!hi], pch = 21, bg = col_pt, col = "white", cex = 1.8, lwd = 1.2)
  points(x[hi],  y[hi],  pch = 21, bg = col_hi, col = "white", cex = 2.1, lwd = 1.2)
  text(x[hi], y[hi], years[hi], pos = 4, col = col_hi, font = 2, cex = 0.95)
  legend("topright", bty = "o", bg = "white", box.col = "grey80", cex = 0.95,
         legend = c(sprintf("r = %.2f, %s", ct$estimate, fmt_p(ct$p.value)),
                    as.expression(bquote(R^2 == .(sprintf("%.2f", ct$estimate^2)) * ", n = " * .(length(x))))))
  mtext(title, side = 3, adj = 0, line = 0.6, font = 2, cex = 1.15)
}

draw_fig <- function() {
  par(mfrow = c(1, 2), mar = c(4.8, 5.2, 3, 1.2), mgp = c(3.2, 0.9, 0), las = 1,
      cex.axis = 1.1, cex.lab = 1.2, lwd = 1.5)
  draw_panel(raw$melt_date_NOAA, raw$SWC_May, raw$Year, "(a) Observed",
             "Snowmelt date (day of year)", "May soil water content (%)")
  draw_panel(det$melt_detr, det$swc_detr, det$Year, "(b) Detrended",
             "Snowmelt date anomaly (days)", "May soil water content anomaly (%)",
             zero_lines = TRUE)
}

tiff("Fig_4_snowmelt_SWCMay.tiff", width = 17, height = 7.8, units = "in", res = 600,
     compression = "lzw", type = "cairo", pointsize = 20)
draw_fig(); dev.off()
pdf("Fig_4_snowmelt_SWCMay.pdf", width = 17, height = 7.8, pointsize = 20)
draw_fig(); dev.off()
