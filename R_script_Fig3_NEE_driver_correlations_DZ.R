# Script modified by DZ on 2026-11-16 for the NEE–environmental-driver correlations
# to generate Fig. 3
# Locally, the script assumes that the CSV is in the current working directory.

#to load the file to run on R studio
data_file <- "C:/Users/local-admin/Desktop/Long_term_paper/Github_submission/SUBMITTED_GITHUB/TOT_MET_FLX_v2026_11_13_f.csv"
out_dir <- "C:/Users/local-admin/Desktop/Long_term_paper/Github_submission/SUBMITTED_GITHUB"

#to load the file to run on Github
data_file <- Sys.getenv("NEE_DATA_FILE",
                        unset = "TOT_MET_FLX_v2026_11_13_f.csv")
out_dir <- Sys.getenv("NEE_OUTPUT_DIR", unset = getwd())

if (!file.exists(data_file)) {
  stop("Input CSV not found: ", data_file,
       "\nPlace the CSV beside this script or set NEE_DATA_FILE.")
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dat <- read.csv(data_file, check.names = FALSE, stringsAsFactors = FALSE,
                na.strings = c("", "NA", "NaN"))
dat$Year <- as.numeric(dat$Year)
dat <- dat[order(dat$Year), ]

# Display p-values as in the reference figure.
format_p <- function(p) {
  if (!is.finite(p)) return("NA")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

# Each row is one panel. lag = 1 aligns the previous year's active-layer depth
# with the following year's May NEE.
panels <- data.frame(
  month = c(
    "May", "June", "June", "July", "July", "August", "August",
    "September", "September", "September", "October", "October",
    "October", "October", "October"
  ),
  nee = c(
    "NEE_May", "NEE_June", "NEE_June", "NEE_July", "NEE_July",
    "NEE_Aug", "NEE_Aug", "NEE_Sept", "NEE_Sept", "NEE_Sept",
    "NEE_Oct", "NEE_Oct", "NEE_Oct", "NEE_Oct", "NEE_Oct"
  ),
  driver = c(
    "U1_Barrow_mean", "melt_date_NOAA", "SWC_June", "melt_date_NOAA",
    "SWC_June", "SWC_May", "SWC_June", "SWC_May", "SWC_June",
    "U1_Barrow_mean", "SWC_May", "GDD_Aug", "SWC_July", "SWC_Aug",
    "SWC_June"
  ),
  xlab = c(
    "Active layer depth (cm below surface)\n(previous year)",
    "NOAA melt date\n(day of year)", "June SWC (%)",
    "NOAA melt date\n(day of year)", "June SWC (%)", "May SWC (%)",
    "June SWC (%)", "May SWC (%)", "June SWC (%)",
    "Active layer depth (cm below surface)", "May SWC (%)", "GDD August",
    "July SWC (%)", "August SWC (%)", "June SWC (%)"
  ),
  lag = c(1, rep(0, 14)),
  stringsAsFactors = FALSE
)

# Check all requested columns before plotting.
needed <- unique(c("Year", panels$nee, panels$driver))
missing <- setdiff(needed, names(dat))
if (length(missing) > 0) {
  stop("Missing columns: ", paste(missing, collapse = ", "))
}

make_pairs <- function(nee_col, driver_col, lag_years = 0) {
  nee <- data.frame(Year = dat$Year, NEE = dat[[nee_col]])
  driver <- data.frame(Year = dat$Year + lag_years,
                       driver = dat[[driver_col]])
  z <- merge(nee, driver, by = "Year")
  z <- z[is.finite(z$NEE) & is.finite(z$driver), ]
  z[order(z$Year), ]
}

plot_panel <- function(pair, panel, ylab = "") {
  x <- pair$driver
  y <- pair$NEE

  if (length(x) < 3 || sd(x) == 0 || sd(y) == 0) {
    plot.new()
    box()
    return(invisible(NULL))
  }

  fit <- lm(y ~ x)
  ct <- cor.test(x, y, method = "pearson")
  r <- unname(ct$estimate)
  p <- ct$p.value
  r2 <- r^2

  plot(x, y,
       pch = 16, col = "#377eb8", cex = 0.85,
       xlab = panel$xlab, ylab = ylab,
       axes = TRUE, las = 1,
       cex.lab = 1.0, cex.axis = 0.9,
       panel.first = grid(col = "grey90", lty = 1))
  abline(fit, col = "#377eb8", lty = 2, lwd = 1.2)
  box()

  # Statistics box in the upper-left corner.
  usr <- par("usr")
  dx <- diff(usr[1:2]); dy <- diff(usr[3:4])
  x0 <- usr[1] + 0.035 * dx
  y0 <- usr[4] - 0.055 * dy
  line_gap <- 0.085 * dy
  labels <- c(
    sprintf("r = %.2f", r),
    paste0("p = ", format_p(p)),
    sprintf("R² = %.2f", r2),
    sprintf("n = %d", nrow(pair))
  )
  box_width <- 0.38 * dx
  box_height <- 0.34 * dy
  rect(x0 - 0.02 * dx, y0 - 3.35 * line_gap,
       x0 + box_width, y0 + 0.35 * line_gap,
       col = adjustcolor("white", alpha.f = 0.80), border = "grey60")
  for (j in seq_along(labels)) {
    text(x0, y0 - (j - 1) * line_gap, labels[j],
         adj = c(0, 0.5), cex = 0.86,
         font = if (j == 2 && p < 0.05) 2 else 1)
  }
}

# Five columns are needed for the October row. Blank cells are left empty to
# create the same lower-triangular arrangement as the reference figure.
layout_matrix <- rbind(
  c(1, 0, 0, 0, 0),
  c(2, 3, 0, 0, 0),
  c(4, 5, 0, 0, 0),
  c(6, 7, 0, 0, 0),
  c(8, 9, 10, 0, 0),
  c(11, 12, 13, 14, 15)
)

png(file.path(out_dir, "Fig.3_NEE_driver_correlations.png"),
    width = 5200, height = 6800, res = 300)
layout(layout_matrix, widths = rep(1, 5), heights = rep(1, 6))
par(oma = c(0, 0, 0, 0))

for (i in seq_len(nrow(panels))) {
  pair <- make_pairs(panels$nee[i], panels$driver[i], panels$lag[i])
  # Only the leftmost panel in each row receives a y-axis label.
  row_first <- i %in% c(1, 2, 4, 6, 8, 11)
  ylab <- if (row_first) "NEE (g C-CO₂ m⁻² month⁻¹)" else ""
  plot_panel(pair, panels[i, ], ylab = ylab)
}
dev.off()

# Save the statistics used in the figure.
stats <- do.call(rbind, lapply(seq_len(nrow(panels)), function(i) {
  pair <- make_pairs(panels$nee[i], panels$driver[i], panels$lag[i])
  ct <- cor.test(pair$driver, pair$NEE, method = "pearson")
  data.frame(
    month = panels$month[i],
    NEE = panels$nee[i],
    driver = panels$driver[i],
    lag_years = panels$lag[i],
    n = nrow(pair),
    r = unname(ct$estimate),
    p = ct$p.value,
    R_squared = unname(ct$estimate)^2
  )
}))
write.csv(stats,
          file.path(out_dir, "Fig.3_NEE_driver_correlations_statistics.csv"),
          row.names = FALSE)

message("Figure and statistics saved in: ", normalizePath(out_dir))

