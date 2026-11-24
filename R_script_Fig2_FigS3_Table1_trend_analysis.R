#script modified by DZ on 2026-11-24 for the Long term paper
#temporal analysis

# =============================================================================
# Temporal trends in monthly NEE and environmental drivers, US-Brw (Utqiagvik)
# Reproduces the statistics in the Fig. 2 table, Table 1 (all NEE months and
# environmental variables) and the residual ACF figure (Fig. S5)
#
# Methods: OLS regression, Kendall rank correlation, Yue-Pilon trend-free
# pre-whitened Mann-Kendall test (zyp), GLS with continuous-time AR(1) errors
# (nlme), GAM with a smooth term for year (mgcv), lag-1 ACF of OLS residuals,
# and Pearson correlation between NEE and its driver.
#
# Missing years are kept as gaps (NA on a continuous annual grid), so that
# autocorrelation is only estimated between consecutive calendar years.
# =============================================================================

library(zyp)    # Yue-Pilon pre-whitening (Bronaugh and Werner, 2019)
library(nlme)   # GLS with corCAR1 (Pinheiro et al.)
library(mgcv)   # GAM (Wood, 2017)

#below is the path to load in Github
# Set the working directory to the folder containing this script and the data, e.g.
# setwd("path/to/this/folder")
dat <- read.csv("TOT_MET_FLX_v2026_11_13_f.csv", check.names = FALSE)

#below is the path to run on R studio
# data_dir <- "C:/Users/local-admin/Desktop/Long_term_paper/Github_submission/SUBMITTED_GITHUB"
# setwd(data_dir) 
# dat <- read.csv("C:/Users/local-admin/Desktop/Long_term_paper/Github_submission/SUBMITTED_GITHUB/TOT_MET_FLX_v2026_11_13_f.csv", check.names = FALSE)

# Active-layer depth (U1) from the previous year, indexed to the NEE year
dat$U1_prev <- c(NA, head(dat$U1_Barrow_mean, -1))

# Month -> driver used in each Fig. 2 panel.
# drv_plot  = series paired with NEE (lagged for ALD)
# drv_trend = series tested for a temporal trend (full record)
panels <- data.frame(
  month     = c("May", "June", "July", "Aug", "Sept", "Oct"),
  driver    = c("Active-layer depth (previous year)", "Snowmelt date", "Snowmelt date",
                "May SWC", "May SWC", "May SWC"),
  drv_plot  = c("U1_prev", "melt_date_NOAA", "melt_date_NOAA", "SWC_May", "SWC_May", "SWC_May"),
  drv_trend = c("U1_Barrow_mean", "melt_date_NOAA", "melt_date_NOAA", "SWC_May", "SWC_May", "SWC_May"),
  stringsAsFactors = FALSE
)

alpha <- 0.05

# -----------------------------------------------------------------------------
# Trend tests for one series. Years without data are kept as NA on a
# continuous annual grid so lag-1 autocorrelation respects real year spacing.
# -----------------------------------------------------------------------------
trend_tests <- function(year, y) {
  ok   <- !is.na(y)
  grid <- seq(min(year[ok]), max(year[ok]))
  yg   <- y[match(grid, year)]
  df   <- data.frame(Year = grid, y = yg)[!is.na(yg), ]

  # 1. OLS regression
  ols <- lm(y ~ Year, data = df)

  # 2. Kendall rank correlation between value and year (equivalent to the
  #    Mann-Kendall test). Exact p-value when there are no ties; normal
  #    approximation (with a warning, suppressed here) when there are ties.
  ken <- suppressWarnings(cor.test(df$Year, df$y, method = "kendall"))

  # 3. Yue-Pilon trend-free pre-whitening + Mann-Kendall (zyp defaults)
  yp <- zyp.yuepilon(yg, grid)

  # 4. GLS with continuous-time AR(1) errors (handles missing years)
  gls_fit <- gls(y ~ Year, data = df, correlation = corCAR1(form = ~ Year), method = "ML")

  # 5. Lag-1 autocorrelation of OLS residuals on the annual grid
  r <- rep(NA, length(grid)); r[match(df$Year, grid)] <- resid(ols)
  acf1 <- acf(r, lag.max = 1, plot = FALSE, na.action = na.pass)$acf[2]

  list(df = df, grid = grid, resid = r,
       n          = nrow(df),
       ols_slope  = unname(coef(ols)[2]),
       ols_p      = summary(ols)$coefficients["Year", "Pr(>|t|)"],
       kendall_tau = unname(ken$estimate),
       kendall_p  = ken$p.value,
       yuepilon_tau = unname(yp["tau"]),
       sen_slope  = unname(yp["trend"]),
       yuepilon_p = unname(yp["sig"]),
       gls_p      = summary(gls_fit)$tTable["Year", "p-value"],
       acf1       = acf1,
       acf_bound  = qnorm(1 - alpha / 2) / sqrt(nrow(df)))
}

# -----------------------------------------------------------------------------
# Non-linear trend: GAM with a smooth term for year, compared with a linear fit
# -----------------------------------------------------------------------------
gam_test <- function(df) {
  lin <- gam(y ~ Year, data = df, method = "REML")
  sm  <- gam(y ~ s(Year, k = 6), data = df, method = "REML")
  list(fit = sm,
       show  = NA,
       edf   = sum(sm$edf) - 1,                         # edf of the smooth (1 = linear)
       p     = summary(sm)$s.table[1, "p-value"],
       dAIC  = AIC(lin) - AIC(sm))                      # > 2: GAM fits clearly better
}

# -----------------------------------------------------------------------------
# Run for each panel
# -----------------------------------------------------------------------------
results <- list(); acf_store <- list()
for (i in seq_len(nrow(panels))) {
  p   <- panels[i, ]
  nee <- paste0("NEE_", p$month)

  tn <- trend_tests(dat$Year, dat[[nee]])
  gm <- gam_test(tn$df)
  td <- trend_tests(dat$Year, dat[[p$drv_trend]])

  pair <- na.omit(dat[, c(nee, p$drv_plot)])
  pr   <- cor.test(pair[[1]], pair[[2]], method = "pearson")

  show_gam <- gm$edf > 1.5 && gm$dAIC > 2 && gm$p < alpha

  results[[i]] <- data.frame(
    Month            = p$month,
    Driver           = p$driver,
    n                = tn$n,
    OLS_slope        = tn$ols_slope,
    OLS_p            = tn$ols_p,
    Kendall_p        = tn$kendall_p,
    Sen_slope        = tn$sen_slope,
    YuePilon_p       = tn$yuepilon_p,
    GLS_AR1_p        = tn$gls_p,
    GAM_p            = if (show_gam) gm$p else NA,
    GAM_edf          = if (show_gam) gm$edf else NA,
    Resid_lag1_r     = tn$acf1,
    Resid_exceeds_CI = abs(tn$acf1) > tn$acf_bound,
    Drv_OLS_p        = td$ols_p,
    Drv_Kendall_p    = td$kendall_p,
    Drv_YuePilon_p   = td$yuepilon_p,
    Drv_GLS_AR1_p    = td$gls_p,
    Pearson_r        = unname(pr$estimate),
    Pearson_p        = pr$p.value,
    R2               = unname(pr$estimate)^2,
    n_pair           = nrow(pair)
  )
  acf_store[[p$month]] <- tn
}
tab <- do.call(rbind, results)
write.csv(tab, "Table_NEE_trends.csv", row.names = FALSE)
print(format(tab, digits = 3), row.names = FALSE)

# -----------------------------------------------------------------------------
# Fig. SX: ACF of OLS residuals for each month (gaps kept)
# -----------------------------------------------------------------------------
pdf("FigSX_ACF_OLS_residuals.pdf", width = 9, height = 6)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
for (m in names(acf_store)) {
  a <- acf_store[[m]]
  acf(a$resid, lag.max = 8, na.action = na.pass, main = paste("NEE", m),
      xlab = "Lag (years)", ylab = "ACF of OLS residuals")
}
dev.off()

# =============================================================================
# Table 1: temporal trends in monthly NEE and all environmental variables
# Each variable is tested over its own full record (years without data kept
# as gaps). A row is flagged "robust" when OLS, Kendall, Yue-Pilon and GLS
# all give p < 0.05; these rows are shown in bold in Table 1.
# =============================================================================
table1_vars <- data.frame(
  column = c("NEE_May", "NEE_June", "NEE_July", "NEE_Aug", "NEE_Sept", "NEE_Oct",
             "SoilT_June", "SoilT_July", "SoilT_Aug", "SoilT_Sept",
             "melt_date_NOAA", "Snow_In_NOAA", "Snow_free_NOAA",
             "SWC_May", "SWC_June", "SWC_July", "SWC_Aug", "SWC_Sept", "SWC_Oct",
             "GDD_June", "GDD_July", "GDD_Aug", "GDD_Sept", "U1_Barrow_mean"),
  label  = c("May NEE", "June NEE", "July NEE", "Aug NEE", "Sept NEE", "Oct NEE",
             "SoilT_June", "SoilT_July", "SoilT_Aug", "SoilT_Sept",
             "Snow melt", "Snow in", "Snow free",
             "SWC_May", "SWC_June", "SWC_July", "SWC_Aug", "SWC_Sept", "SWC_Oct",
             "GDD_June", "GDD_July", "GDD_Aug", "GDD_Sept", "Thaw depth"),
  units  = c(rep("g C-CO2 m-2 month-1", 6), rep("deg C", 4), rep("day of year", 3),
             rep("%", 6), rep("deg C days", 4), "cm below surface"),
  stringsAsFactors = FALSE
)

t1 <- lapply(seq_len(nrow(table1_vars)), function(i) {
  v  <- table1_vars$column[i]
  tt <- trend_tests(dat$Year, dat[[v]])
  gm <- gam_test(tt$df)
  show_gam <- gm$edf > 1.5 && gm$dAIC > 2 && gm$p < alpha
  data.frame(
    Variable          = table1_vars$label[i],
    Units             = table1_vars$units[i],
    N                 = tt$n,
    Start             = min(tt$df$Year),
    End               = max(tt$df$Year),
    OLS_slope_yr      = tt$ols_slope,
    OLS_slope_decade  = tt$ols_slope * 10,
    OLS_p             = tt$ols_p,
    Kendall_tau       = tt$kendall_tau,
    Kendall_p         = tt$kendall_p,
    YuePilon_tau      = tt$yuepilon_tau,
    YuePilon_p        = tt$yuepilon_p,
    GLS_AR1_p         = tt$gls_p,
    GAM_p             = if (show_gam) gm$p else NA,
    GAM_edf           = if (show_gam) gm$edf else NA,
    Resid_lag1_r      = tt$acf1,
    Resid_exceeds_CI  = abs(tt$acf1) > tt$acf_bound,
    Robust            = all(c(tt$ols_p, tt$kendall_p, tt$yuepilon_p, tt$gls_p) < alpha)
  )
})
table1 <- do.call(rbind, t1)
write.csv(table1, "Table_1_trends.csv", row.names = FALSE)
print(format(table1, digits = 3), row.names = FALSE)

# =============================================================================
# Fig. 2: monthly NEE and its strongest temporally valid driver (z-scores)
#   - NEE z-scored over its own record; driver z-scored over its full record
#     and plotted at NEE years
#   - dashed lines: OLS trends (NEE: red; driver over full record: green)
#   - purple line + band: GAM fit +/- 95% CI, drawn only where the GAM shows
#     significant non-linear change (edf > 1.5, dAIC > 2, p < 0.05)
# Output: Fig_2.tiff (600 dpi, LZW) and Fig_2.pdf
# =============================================================================
zscore <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

col_nee <- "#1f77b4"; col_drv <- "#ff7f0e"; col_nt <- "#d62728"
col_dt  <- "#2ca02c"; col_gam <- "#7b3fa0"
fig2_panels <- data.frame(
  month = c("May", "June", "July", "Aug", "Sept", "Oct"),
  title = c("(a) May", "(b) June", "(c) July", "(d) August", "(e) September", "(f) October"),
  driver_label = c("Active-layer depth (previous year)", "Snowmelt date", "Snowmelt date",
                   "May soil water content", "May soil water content", "May soil water content"),
  driver = c("U1_prev", "melt_date_NOAA", "melt_date_NOAA", "SWC_May", "SWC_May", "SWC_May"),
  stringsAsFactors = FALSE
)

draw_fig2 <- function() {
  layout(rbind(c(7, 7), c(1, 2), c(3, 4), c(5, 6)), heights = c(0.7, 3, 3, 3))
  par(mar = c(4.2, 5.0, 2.6, 1.2), mgp = c(3.1, 0.8, 0), cex.axis = 1.15, cex.lab = 1.2,
      las = 1, tcl = -0.4, lwd = 1.5)
  for (i in seq_len(nrow(fig2_panels))) {
    p   <- fig2_panels[i, ]
    nee <- paste0("NEE_", p$month)

    nd <- na.omit(data.frame(Year = dat$Year, y = dat[[nee]]))
    nd$z <- zscore(nd$y)
    dd <- na.omit(data.frame(Year = dat$Year, x = dat[[p$driver]]))
    dd$z <- zscore(dd$x)
    dsub <- dd[dd$Year %in% nd$Year, ]

    plot(NA, xlim = c(2000, 2025), ylim = c(-3.4, 2.9), xlab = if (i >= 5) "Year" else "",
         ylab = if (i %% 2 == 1) "Standardized value (z-score)" else "", bty = "l", xaxt = "n")
    axis(1, at = seq(2000, 2025, 5))
    abline(v = seq(2000, 2025, 5), h = -3:2, col = "grey90", lwd = 1)
    abline(h = 0, col = "grey60", lwd = 1.5)

    # GAM band and fit (z-scored NEE) where non-linear change is significant
    gm <- gam_test(data.frame(Year = nd$Year, y = nd$z))
    show_gam <- gm$edf > 1.5 && gm$dAIC > 2 && gm$p < alpha
    if (show_gam) {
      yr <- seq(min(nd$Year), max(nd$Year), by = 0.25)
      pr <- predict(gm$fit, newdata = data.frame(Year = yr), se.fit = TRUE)
      polygon(c(yr, rev(yr)), c(pr$fit - 1.96 * pr$se.fit, rev(pr$fit + 1.96 * pr$se.fit)),
              col = adjustcolor(col_gam, 0.15), border = NA)
    }

    # driver
    lines(dsub$Year, dsub$z, col = adjustcolor(col_drv, 0.6), lwd = 2)
    points(dsub$Year, dsub$z, pch = 16, col = col_drv, cex = 1.3)
    dfit <- lm(z ~ Year, dd)
    segments(min(dd$Year), predict(dfit, data.frame(Year = min(dd$Year))),
             max(dd$Year), predict(dfit, data.frame(Year = max(dd$Year))), col = col_dt, lwd = 3.4, lty = 2)

    # NEE
    lines(nd$Year, nd$z, col = adjustcolor(col_nee, 0.6), lwd = 2)
    points(nd$Year, nd$z, pch = 16, col = col_nee, cex = 1.3)
    nfit <- lm(z ~ Year, nd)
    segments(min(nd$Year), predict(nfit, data.frame(Year = min(nd$Year))),
             max(nd$Year), predict(nfit, data.frame(Year = max(nd$Year))), col = col_nt, lwd = 3.4, lty = 2)
    if (show_gam) lines(yr, pr$fit, col = col_gam, lwd = 3.6)

    mtext(p$title, side = 3, adj = 0, line = 0.6, font = 2, cex = 1.15)
    mtext(paste("Driver:", p$driver_label), side = 3, adj = 1, line = 0.6, font = 2, cex = 0.9, col = col_drv)
  }
  # legend panel
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", ncol = 3, bty = "n", cex = 1.25, seg.len = 2.5, x.intersp = 0.6, text.width = NA,
         legend = expression("NEE", "Driver", "NEE trend (OLS)", "Driver trend (OLS, full record)",
                             "NEE GAM" %+-% "95% CI (significant non-linear change)"),
         col = c(col_nee, col_drv, col_nt, col_dt, col_gam), pch = c(16, 16, NA, NA, NA),
         lty = c(1, 1, 2, 2, 1), lwd = c(2, 2, 3.4, 3.4, 3.6))
}

tiff("Fig_2.tiff", width = 17, height = 16, units = "in", res = 600, compression = "lzw", type = "cairo", pointsize = 20)
draw_fig2(); dev.off()
pdf("Fig_2.pdf", width = 17, height = 16, pointsize = 20)
draw_fig2(); dev.off()
