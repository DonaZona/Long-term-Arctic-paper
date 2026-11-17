# Script modified by DZ on 2026-11-16 for the NEE–environmental-driver analysis
# to generate Fig. 2

# Yue–Pilon implementation:
#   Theil–Sen detrending -> lag-1 autocorrelation of detrended residuals
#   -> pre-whitening -> trend restoration -> Kendall test.

library(ggplot2)

# Input and output paths. NEE_DATA_FILE can override the default CSV path.
# Keep the file path as text here; read.csv() is called below after validation.
default_data_file <- "C:/Users/local-admin/Desktop/Long_term_paper/Github_submission/SUBMITTED_GITHUB/TOT_MET_FLX_v2026_11_13_f.csv"
data_file <- Sys.getenv("NEE_DATA_FILE", unset = default_data_file)
data_file <- trimws(data_file)
data_file <- gsub('^"|"$', "", data_file)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1]), mustWork = FALSE))
} else {
  getwd()
}
out_dir <- Sys.getenv("NEE_OUTPUT_DIR", unset = file.path(script_dir, "Temporal_analysis_Fig.2"))

# If the configured file is unavailable, look for the same dataset beside the
# script or in the current working directory (including versioned filenames).
if (!file.exists(data_file)) {
  search_dirs <- unique(c(script_dir, getwd()))
  candidates <- unlist(lapply(search_dirs, function(d) {
    if (!dir.exists(d)) return(character(0))
    c(file.path(d, "TOT_MET_FLX_v2026_11_13_f(5).csv"),
      list.files(d, pattern = "^TOT_MET_FLX_v2026_11_13_f.*\\.csv$",
                 full.names = TRUE, ignore.case = TRUE))
  }), use.names = FALSE)
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates)) {
    data_file <- candidates[1]
    warning("Configured CSV path was not found; using: ", data_file)
  } else if (interactive()) {
    data_file <- file.choose()
  } else {
    stop(
      "Input CSV not found at: ", data_file,
      "\nCheck that the file exists, or set NEE_DATA_FILE to its full path."
    )
  }
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dat <- read.csv(data_file, check.names = FALSE, stringsAsFactors = FALSE,
                na.strings = c("NA", "", "NaN"))
if (!"Year" %in% names(dat)) stop("Required column 'Year' is missing.")
dat$Year <- suppressWarnings(as.numeric(dat$Year))
dat <- dat[is.finite(dat$Year), ]
if (anyDuplicated(dat$Year)) stop("Year must contain one row per year.")
dat <- dat[order(dat$Year), ]

fmt_p <- function(p) {
  if (!is.finite(p)) return("NA")
  if (p < 0.001) return("<0.001")
  sprintf("%.4f", p)
}

# Match Fig. 2: bold statistically significant p-value lines (p < 0.05).
p_line <- function(label, p) {
  list(text = ifelse(is.finite(p), paste0("  ", label, " p = ", fmt_p(p)),
                     paste0("  ", label, " p = NA")),
       font = ifelse(is.finite(p) && p < 0.05, 2, 1))
}

# Theil–Sen slope using the actual time axis.
theil_sen_slope <- function(x, y) {
  z <- is.finite(x) & is.finite(y)
  x <- x[z]; y <- y[z]
  if (length(y) < 2) return(NA_real_)
  slopes <- unlist(lapply(seq_len(length(y) - 1), function(i) {
    (y[(i + 1):length(y)] - y[i]) / (x[(i + 1):length(y)] - x[i])
  }))
  median(slopes, na.rm = TRUE)
}

# Yue–Pilon trend-free pre-whitened Mann–Kendall test.
yue_pilon <- function(x, y) {
  z <- data.frame(x = x, y = y)
  z <- z[is.finite(z$x) & is.finite(z$y), ]
  x <- z$x; y <- z$y
  n <- length(y)
  if (n < 5 || length(unique(y)) < 2) {
    return(list(tau = NA_real_, p = NA_real_, rho = NA_real_, sen_slope = NA_real_))
  }

  sen <- theil_sen_slope(x, y)
  detrended <- y - sen * x
  rho <- suppressWarnings(cor(detrended[-n], detrended[-1], use = "complete.obs"))
  if (!is.finite(rho)) rho <- 0
  rho <- max(0, rho)

  prewhitened_residual <- detrended[-1] - rho * detrended[-n]
  blended <- prewhitened_residual + sen * x[-1]
  kt <- suppressWarnings(cor.test(x[-1], blended, method = "kendall", exact = FALSE))

  list(tau = unname(kt$estimate), p = kt$p.value,
       rho = rho, sen_slope = sen)
}

trend_stats <- function(year, value) {
  z <- data.frame(year = year, value = value)
  z <- z[is.finite(z$year) & is.finite(z$value), ]
  x <- z$year; y <- z$value
  n <- length(y)
  if (n < 3) {
    return(data.frame(n = n, first_year = NA, last_year = NA,
                      slope = NA, change_decade = NA, r2 = NA, ols_p = NA,
                      kendall_tau = NA, kendall_p = NA, yue_tau = NA,
                      yue_p = NA, yue_lag1_rho = NA, mean = NA, min = NA, max = NA))
  }

  ols <- lm(y ~ x)
  ols_coef <- summary(ols)$coefficients
  kt <- suppressWarnings(cor.test(x, y, method = "kendall", exact = FALSE))
  yp <- yue_pilon(x, y)

  data.frame(
    n = n,
    first_year = min(x),
    last_year = max(x),
    slope = unname(coef(ols)[2]),
    change_decade = unname(coef(ols)[2]) * 10,
    r2 = summary(ols)$r.squared,
    ols_p = ols_coef[2, 4],
    kendall_tau = unname(kt$estimate),
    kendall_p = kt$p.value,
    yue_tau = yp$tau,
    yue_p = yp$p,
    yue_lag1_rho = yp$rho,
    mean = mean(y),
    min = min(y),
    max = max(y)
  )
}

units_map <- c(
  NEE_May = "g C-CO2 m^-2 month^-1", NEE_June = "g C-CO2 m^-2 month^-1",
  NEE_July = "g C-CO2 m^-2 month^-1", NEE_Aug = "g C-CO2 m^-2 month^-1",
  NEE_Sept = "g C-CO2 m^-2 month^-1", NEE_Oct = "g C-CO2 m^-2 month^-1",
  SoilT_May = "deg C", SoilT_June = "deg C", SoilT_July = "deg C",
  SoilT_Aug = "deg C", SoilT_Sept = "deg C", SoilT_Oct = "deg C",
  melt_date_NOAA = "day of year", Snow_In_NOAA = "day of year",
  Snow_free_NOAA = "day of year", `Snow Season Length_NOAA` = "days",
  SWC_May = "%", SWC_June = "%", SWC_July = "%", SWC_Aug = "%",
  SWC_Sept = "%", SWC_Oct = "%", SWC_JJA_median = "%",
  GDD_JJA = "deg C days", GDD_June = "deg C days", GDD_July = "deg C days",
  GDD_Aug = "deg C days", GDD_Sept = "deg C days",
  U1_Barrow_mean = "cm below surface"
)

nee_cols <- c("NEE_May", "NEE_June", "NEE_July", "NEE_Aug", "NEE_Sept", "NEE_Oct")
nee_labels <- c(NEE_May = "May NEE", NEE_June = "June NEE", NEE_July = "July NEE",
                NEE_Aug = "August NEE", NEE_Sept = "September NEE", NEE_Oct = "October NEE")

env_cols <- c("SoilT_May", "SoilT_June", "SoilT_July", "SoilT_Aug",
              "SoilT_Sept", "SoilT_Oct",
              "melt_date_NOAA", "Snow_In_NOAA", "Snow_free_NOAA",
              "Snow Season Length_NOAA", "SWC_May", "SWC_June", "SWC_July",
              "SWC_Aug", "SWC_Sept", "SWC_Oct", "SWC_JJA_median",
              "GDD_JJA", "GDD_June", "GDD_July", "GDD_Aug", "GDD_Sept",
              "U1_Barrow_mean")

# Confirm that every requested column exists in the input dataset.
required_cols <- unique(c("Year", nee_cols, env_cols))
missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0) {
  stop(
    "The following required columns are missing from the dataset: ",
    paste(missing_cols, collapse = ", ")
  )
}

make_row <- function(type, label, col, scope) {
  if (!col %in% names(dat)) {
    stop("Column not found in dataset: ", col)
  }

  s <- trend_stats(dat$Year, dat[[col]])
  cbind(data.frame(series_type = type, variable = label, source_column = col,
                   temporal_scope = scope, units = unname(units_map[col])), s)
}

results <- do.call(rbind, c(
  lapply(nee_cols, function(x) make_row("NEE", nee_labels[x], x, "available NEE record")),
  lapply(env_cols, function(x) make_row("environmental variable", x, x,
                                         "entire available environmental record"))
))
rownames(results) <- NULL
results$ols_p_significant_0.05 <- results$ols_p < 0.05
results$kendall_p_significant_0.05 <- results$kendall_p < 0.05
results$yue_p_significant_0.05 <- results$yue_p < 0.05
results$ols_p_display <- vapply(results$ols_p, fmt_p, character(1))
results$kendall_p_display <- vapply(results$kendall_p, fmt_p, character(1))
results$yue_p_display <- vapply(results$yue_p, fmt_p, character(1))

write.csv(results,
          file.path(out_dir, "NEE_and_environmental_temporal_trends_consistent_method.csv"),
          row.names = FALSE, na = "")

# Panel definitions. May active-layer depth is aligned to the following spring.
panels <- data.frame(
  month = c("May", "June", "July", "August", "September", "October"),
  nee = c("NEE_May", "NEE_June", "NEE_July", "NEE_Aug", "NEE_Sept", "NEE_Oct"),
  driver = c("U1_Barrow_mean", "melt_date_NOAA", "melt_date_NOAA", "SWC_May", "SWC_May", "SWC_May"),
  driver_label = c("Active-layer depth (previous year)", "Snowmelt date", "Snowmelt date",
                   "May SWC", "May SWC", "May SWC"),
  lagged = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE), stringsAsFactors = FALSE
)

association_rows <- list()
panel_rows <- list()

png(file.path(out_dir, "Fig.2_FINAL_reproduced.png"),
    width = 5200, height = 8000, res = 300)
layout(matrix(seq_len(nrow(panels) * 2), ncol = 2, byrow = TRUE),
       widths = c(5.5, 1.25))
par(oma = c(0, 0, 0, 0))

for (i in seq_len(nrow(panels))) {
  p <- panels[i, ]
  nee <- dat[, c("Year", p$nee)]
  names(nee) <- c("Year", "NEE")
  nee <- nee[is.finite(nee$NEE), ]

  driver <- dat[, c("Year", p$driver)]
  names(driver) <- c("Year", "driver")
  if (p$lagged) driver$Year <- driver$Year + 1

  pair <- merge(nee, driver, by = "Year")
  pair <- pair[is.finite(pair$NEE) & is.finite(pair$driver), ]

  # Guard correlation and scaling against too few pairs or zero variance.
  assoc <- if (nrow(pair) >= 3 && sd(pair$NEE) > 0 && sd(pair$driver) > 0) {
    suppressWarnings(cor.test(pair$NEE, pair$driver, method = "pearson"))
  } else {
    list(estimate = NA_real_, p.value = NA_real_)
  }
  association_rows[[i]] <- data.frame(
    month = p$month, NEE = p$nee, driver = p$driver_label,
    lagged_driver = p$lagged, n_paired = nrow(pair),
    pearson_r = unname(assoc$estimate), pearson_p = assoc$p.value,
    r_squared = unname(assoc$estimate)^2
  )

  nee_trend <- trend_stats(pair$Year, pair$NEE)
  full_driver <- trend_stats(dat$Year, dat[[p$driver]])

  panel_rows[[i]] <- data.frame(
    month = p$month,
    NEE_variable = p$nee,
    driver = p$driver_label,
    driver_source_column = p$driver,
    driver_trend_scope = "entire available environmental record",
    NEE_n = nee_trend$n,
    NEE_OLS_p = nee_trend$ols_p,
    NEE_Kendall_p = nee_trend$kendall_p,
    NEE_Yue_Pilon_p = nee_trend$yue_p,
    driver_n_full_record = full_driver$n,
    driver_OLS_p_full_record = full_driver$ols_p,
    driver_Kendall_p_full_record = full_driver$kendall_p,
    driver_Yue_Pilon_p_full_record = full_driver$yue_p,
    paired_n = nrow(pair),
    Pearson_r = unname(assoc$estimate),
    Pearson_p = assoc$p.value,
    R_squared = unname(assoc$estimate)^2
  )

  safe_z <- function(x) {
    s <- sd(x, na.rm = TRUE)
    if (!length(x) || !is.finite(s) || s == 0) rep(0, length(x))
    else as.numeric((x - mean(x, na.rm = TRUE)) / s)
  }
  nee_z <- safe_z(pair$NEE)
  drv_z <- safe_z(pair$driver)
  yr <- pair$Year
  yr_range <- range(dat$Year, na.rm = TRUE)

  par(mar = c(3.5, 5.5, 2.5, 1))
  if (length(yr) == 0) {
    plot.new()
    title(main = paste0(p$month, ": NEE and ", p$driver_label))
    text(0.5, 0.5, "No paired observations")
  } else {
  plot(yr, nee_z, type = "o", pch = 16, lwd = 2, col = "#1f77b4",
       xlim = yr_range + c(-0.5, 0.5), ylim = range(c(nee_z, drv_z), na.rm = TRUE) + c(-0.3, 0.3),
       xlab = if (i == nrow(panels)) "Year" else "", ylab = "Standardized value\n(z-score)",
       main = paste0(p$month, ": NEE and ", p$driver_label), cex.main = 1.25,
       cex.lab = 1.1, cex.axis = 1.0)
  lines(yr, drv_z, type = "o", pch = 16, lwd = 2, col = "#ff7f0e")
  abline(h = 0, col = "grey70")
  abline(lm(nee_z ~ yr), col = "#d62728", lty = 2, lwd = 1.5)

  # Project the full-record driver trend onto the paired-year display scale.
  # Keep prediction and scaling vectors aligned to the paired years.
  driver_ok <- is.finite(dat$Year) & is.finite(dat[[p$driver]])
  if (sum(driver_ok) >= 2 && length(yr) >= 2 && sd(pair$driver) > 0) {
    driver_full_df <- data.frame(Year = dat$Year[driver_ok],
                                 driver = dat[[p$driver]][driver_ok])
    full_fit <- lm(driver ~ Year, data = driver_full_df)
    driver_years <- dat$Year[driver_ok]
    pred_raw <- as.numeric(predict(full_fit, newdata = data.frame(Year = driver_years)))
    pred_z <- (pred_raw - mean(pair$driver)) / sd(pair$driver)
    # For May, active-layer depth is aligned to the following year, matching
    # the orange driver series and the paired NEE observations.
    plot_driver_years <- driver_years + if (isTRUE(p$lagged)) 1 else 0
    lines(plot_driver_years, pred_z, col = "#2ca02c", lty = 2, lwd = 1.5)
  }
  legend("topleft", legend = c("NEE", p$driver_label),
         col = c("#1f77b4", "#ff7f0e"), lwd = 2, pch = 16,
         bty = "n", cex = 0.85)

  }

  # Right-hand statistics panel.
  par(mar = c(3.5, 0.2, 2.5, 0.2))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  stat_lines <- list(
    list("NEE trend:", 1),
    p_line("  OLS", nee_trend$ols_p),
    p_line("  Kendall", nee_trend$kendall_p),
    p_line("  Yue–Pilon", nee_trend$yue_p),
    list("Driver trend (full record):", 1),
    p_line("  OLS", full_driver$ols_p),
    p_line("  Kendall", full_driver$kendall_p),
    p_line("  Yue–Pilon", full_driver$yue_p),
    list("NEE–driver association (paired years):", 1),
    list(sprintf("  Pearson r = %.2f; p = %s", assoc$estimate, fmt_p(assoc$p.value)),
         ifelse(is.finite(assoc$p.value) && assoc$p.value < 0.05, 2, 1)),
    list(sprintf("  R² = %.2f; n = %d", assoc$estimate^2, nrow(pair)), 1)
  )
  yy <- seq(0.95, 0.10, length.out = length(stat_lines))
  for (j in seq_along(stat_lines)) {
    text(0.01, yy[j], stat_lines[[j]]$text, adj = c(0, 0.5),
         cex = 0.84, font = stat_lines[[j]]$font)
  }
  box()
}
dev.off()

associations <- do.call(rbind, association_rows)
write.csv(associations,
          file.path(out_dir, "NEE_driver_paired_associations_Fig.2.csv"),
          row.names = FALSE, na = "")

panel_statistics <- do.call(rbind, panel_rows)
write.csv(panel_statistics,
          file.path(out_dir, "NEE_figure_panel_statistics_Fig.2.csv"),
          row.names = FALSE, na = "")

message("Analysis complete. Outputs written to: ", normalizePath(out_dir))
