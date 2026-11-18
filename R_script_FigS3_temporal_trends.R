# Script modified by DZ on 2026-11-16 for the NEE–environmental-driver analysis
# with the original units
# to generate Fig. S3

#below is to load into R studio----
data_file <- "C:/Users/local-admin/Desktop/Long_term_paper/Github_submission/SUBMITTED_GITHUB/TOT_MET_FLX_v2026_11_13_f.csv"
out_dir <- "C:/Users/local-admin/Desktop/Long_term_paper/Github_submission/SUBMITTED_GITHUB"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dat <- read.csv(data_file, check.names = FALSE, stringsAsFactors = FALSE,
                na.strings = c("", "NA", "NaN"))
dat$Year <- as.numeric(dat$Year)
dat <- dat[order(dat$Year), ]

format_p <- function(p) {
  if (!is.finite(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  sprintf("%.4f", p)
}

#below is to load the file on Github----
data_file <- Sys.getenv(
  "NEE_DATA_FILE",
  unset = "TOT_MET_FLX_v2026_11_13_f.csv"
)

out_dir <- Sys.getenv(
  "NEE_OUTPUT_DIR",
  unset = "Figure_S3_outputs"
)

if (!file.exists(data_file)) {
  stop(
    "Input CSV not found: ", data_file,
    "\nEnsure the CSV is in the repository or set NEE_DATA_FILE."
  )
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dat <- read.csv(
  data_file,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA", "NaN")
)

dat$Year <- as.numeric(dat$Year)
dat <- dat[order(dat$Year), ]
#----


# Theil–Sen slope using the actual year axis.
theil_sen_slope <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(y) < 2) return(NA_real_)

  slopes <- unlist(lapply(seq_len(length(y) - 1), function(i) {
    (y[(i + 1):length(y)] - y[i]) /
      (x[(i + 1):length(x)] - x[i])
  }))
  median(slopes, na.rm = TRUE)
}

# Yue–Pilon trend-free pre-whitened Mann–Kendall test.
yue_pilon <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  n <- length(y)

  if (n < 5 || length(unique(y)) < 2) {
    return(list(tau = NA_real_, p = NA_real_, rho = NA_real_))
  }

  sen <- theil_sen_slope(x, y)
  detrended <- y - sen * x
  rho <- suppressWarnings(cor(detrended[-n], detrended[-1],
                              use = "complete.obs"))
  if (!is.finite(rho)) rho <- 0
  rho <- max(0, rho)

  prewhitened <- detrended[-1] - rho * detrended[-n]
  restored <- prewhitened + sen * x[-1]
  test <- suppressWarnings(
    cor.test(x[-1], restored, method = "kendall", exact = FALSE)
  )

  list(tau = unname(test$estimate), p = test$p.value, rho = rho)
}

trend_stats <- function(year, value) {
  keep <- is.finite(year) & is.finite(value)
  x <- year[keep]
  y <- value[keep]
  n <- length(y)

  if (n < 3) {
    return(data.frame(
      n = n, start = NA, end = NA, slope = NA,
      decade_change = NA, r2 = NA, ols_p = NA,
      kendall_p = NA, yue_p = NA
    ))
  }

  fit <- lm(y ~ x)
  fit_summary <- summary(fit)
  kendall <- suppressWarnings(
    cor.test(x, y, method = "kendall", exact = FALSE)
  )
  yue <- yue_pilon(x, y)

  data.frame(
    n = n,
    start = min(x),
    end = max(x),
    slope = unname(coef(fit)[2]),
    decade_change = unname(coef(fit)[2]) * 10,
    r2 = fit_summary$r.squared,
    ols_p = fit_summary$coefficients[2, 4],
    kendall_p = kendall$p.value,
    yue_p = yue$p
  )
}

series <- data.frame(
  label = c("May soil water content", "GDD June", "GDD July",
            "Active layer depth", "NEE May", "NEE August", "NEE September"),
  column = c("SWC_May", "GDD_June", "GDD_July", "U1_Barrow_mean",
             "NEE_May", "NEE_Aug", "NEE_Sept"),
  units = c("%", "°C days", "°C days", "cm below surface",
            "g C-CO₂ m⁻² month⁻¹", "g C-CO₂ m⁻² month⁻¹",
            "g C-CO₂ m⁻² month⁻¹"),
  stringsAsFactors = FALSE
)

missing <- setdiff(c("Year", series$column), names(dat))
if (length(missing) > 0) {
  stop("Missing columns: ", paste(missing, collapse = ", "))
}

# Plot the right-hand statistics box.
draw_stats_box <- function(s) {
  par(mar = c(2.8, 0.1, 1.2, 0.1))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))

  text_lines <- c(
    "OLS trend:",
    sprintf("slope = %.2f unit yr⁻¹", s$slope),
    sprintf("decade = %.2f", s$decade_change),
    paste0("p ", format_p(s$ols_p), "; R² = ", sprintf("%.2f", s$r2)),
    paste0("Kendall: p = ", format_p(s$kendall_p)),
    paste0("Yue–Pilon: p = ", format_p(s$yue_p)),
    sprintf("n = %d (%d–%d)", s$n, s$start, s$end)
  )

  y <- seq(0.91, 0.18, length.out = length(text_lines))
  rect(0.01, 0.08, 0.99, 0.97,
       col = adjustcolor("white", alpha.f = 0.88), border = "grey60")
  for (i in seq_along(text_lines)) {
    text(0.04, y[i], text_lines[i], adj = c(0, 0.5), cex = 0.78,
         font = if (i >= 4 && i <= 6 &&
                   ((i == 4 && s$ols_p < 0.05) ||
                    (i == 5 && s$kendall_p < 0.05) ||
                    (i == 6 && s$yue_p < 0.05))) 2 else 1)
  }
  box()
}

draw_series <- function(i) {
  s <- series[i, ]
  y <- dat[[s$column]]
  keep <- is.finite(dat$Year) & is.finite(y)
  x <- dat$Year[keep]
  y <- y[keep]
  stats <- trend_stats(x, y)

  par(mar = c(2.8, 5.4, 1.6, 0.5))
  ylim <- range(y, na.rm = TRUE)
  padding <- 0.08 * diff(ylim)
  if (!is.finite(padding) || padding == 0) padding <- 1

  # Match the grey late-period shading in the May SWC panel.
  plot(x, y, type = "n", xlim = c(1999.5, 2025.5),
       ylim = ylim + c(-padding, padding),
       xaxt = "n", xlab = "", ylab = s$units,
       cex.lab = 0.95, cex.axis = 0.85,
       panel.first = {
         if (i == 1) {
           rect(2019, -1e6, 2024.5, 1e6,
                col = "grey85", border = NA)
         }
         grid(col = "grey90", lty = 1)
       })

  lines(x, y, type = "o", pch = 16, lwd = 1.6, cex = 0.55,
        col = "#1f77b4")
  abline(lm(y ~ x), col = "red", lty = 2, lwd = 1.4)
  axis(1, at = seq(2000, 2025, by = 5), cex.axis = 0.85)
  box()
  title(main = s$label, cex.main = 1.05, font.main = 2)

  draw_stats_box(stats)
  invisible(stats)
}

png(file.path(out_dir, "Fig.S3_orig_unit_reproduced.png"),
    width = 3200, height = 7200, res = 300)
layout(matrix(seq_len(nrow(series) * 2), ncol = 2, byrow = TRUE),
       widths = c(5.9, 1.25), heights = rep(1, nrow(series)))
par(oma = c(0, 0, 0, 0))

all_stats <- lapply(seq_len(nrow(series)), draw_series)
dev.off()

stats_table <- do.call(rbind, lapply(seq_len(nrow(series)), function(i) {
  s <- series[i, ]
  y <- dat[[s$column]]
  keep <- is.finite(dat$Year) & is.finite(y)
  z <- trend_stats(dat$Year[keep], y[keep])
  data.frame(variable = s$label, column = s$column, units = s$units, z)
}))

write.csv(stats_table,
          file.path(out_dir, "Fig.S3_orig_unit_temporal_statistics.csv"),
          row.names = FALSE)

message("Figure and statistics saved in: ", normalizePath(out_dir))
