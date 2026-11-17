# Script modified by DZ on 2026-11-16 for the correlogram of
# Fig. S4

data_file <- "C:/Users/local-admin/Desktop/Long_term_paper/Github_submission/SUBMITTED_GITHUB/TOT_MET_FLX_v2026_11_13_f.csv"
out_dir <- "C:/Users/local-admin/Desktop/Long_term_paper/Github_submission/SUBMITTED_GITHUB"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dat <- read.csv(data_file, check.names = FALSE, stringsAsFactors = FALSE,
                na.strings = c("", "NA", "NaN"))

variables <- c(
  "SoilT_June", "SoilT_July", "SoilT_Aug", "SoilT_Sept",
  "melt_date_NOAA", "Snow_In_NOAA", "Snow_free_NOAA",
  "Snow Season Length_NOAA", "SWC_May", "SWC_June", "SWC_July",
  "SWC_Aug", "SWC_Sept", "SWC_Oct", "GDD_June", "GDD_July",
  "GDD_Aug", "GDD_Sept", "U1_Barrow_mean"
)

labels <- c(
  "Soil T June", "Soil T July", "Soil T Aug", "Soil T Sept",
  "Melt date", "Snow-in date", "Snow-free date", "Snow-season length",
  "SWC May", "SWC June", "SWC July", "SWC Aug", "SWC Sept", "SWC Oct",
  "GDD June", "GDD July", "GDD Aug", "GDD Sept", "Active layer depth"
)

missing <- setdiff(c("Year", variables), names(dat))
if (length(missing) > 0) {
  stop("Missing annual-summary columns: ", paste(missing, collapse = ", "),
       "\nThe raw half-hourly AmeriFlux file cannot generate this figure alone.")
}
dat[variables] <- lapply(dat[variables], function(x) as.numeric(as.character(x)))

nvar <- length(variables)
r_mat <- p_mat <- matrix(NA_real_, nvar, nvar,
                          dimnames = list(labels, labels))
n_mat <- matrix(NA_integer_, nvar, nvar,
                dimnames = list(labels, labels))

for (i in seq_len(nvar)) {
  for (j in i:nvar) {
    keep <- is.finite(dat[[variables[i]]]) & is.finite(dat[[variables[j]]])
    n <- sum(keep)
    n_mat[i, j] <- n_mat[j, i] <- n
    if (i == j) {
      r_mat[i, j] <- 1
      p_mat[i, j] <- 0
    } else if (n >= 3 && sd(dat[[variables[i]]][keep]) > 0 &&
               sd(dat[[variables[j]]][keep]) > 0) {
      test <- suppressWarnings(cor.test(
        dat[[variables[i]]][keep], dat[[variables[j]]][keep],
        method = "pearson"
      ))
      r_mat[i, j] <- r_mat[j, i] <- unname(test$estimate)
      p_mat[i, j] <- p_mat[j, i] <- test$p.value
    }
  }
}

stars <- function(p) {
  if (!is.finite(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  ""
}

cell_text <- function(r, p) {
  if (!is.finite(r)) return("")
  paste0(sprintf("%.2f", r), stars(p))
}

palette <- colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(401)
value_to_col <- function(x) {
  index <- round((x + 1) / 2 * (length(palette) - 1)) + 1
  palette[pmax(1, pmin(length(palette), index))]
}

draw_correlogram <- function() {
  par(mar = c(10, 12, 4, 2))
  plot.new()
  plot.window(xlim = c(0.5, nvar + 0.5), ylim = c(0.5, nvar + 0.5), asp = 1)

  for (i in seq_len(nvar)) {
    for (j in i:nvar) {
      r <- r_mat[i, j]
      rect(j - 0.5, nvar - i + 0.5, j + 0.5, nvar - i + 1.5,
           col = value_to_col(r), border = "white", lwd = 0.8)
      text(j, nvar - i + 1, cell_text(r, p_mat[i, j]), cex = 0.78,
           col = if (abs(r) > 0.65) "white" else "black")
    }
  }

  axis(1, at = seq_len(nvar), labels = labels, las = 2,
       cex.axis = 0.78, tick = FALSE)
  axis(2, at = nvar:1, labels = labels, las = 1,
       cex.axis = 0.88, tick = FALSE)
  box(lwd = 1.2)
  title(main = "Pearson correlations", cex.main = 1.45, font.main = 2)
}

draw_colorbar <- function() {
  par(mar = c(10, 1, 4, 5))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(-1, 1))
  y <- seq(-1, 1, length.out = length(palette) + 1)
  for (k in seq_along(palette)) {
    rect(0.15, y[k], 0.55, y[k + 1], col = palette[k], border = NA)
  }
  axis(4, at = seq(-1, 1, by = 0.25), las = 1, cex.axis = 0.82)
  mtext("Correlation coefficient", side = 4, line = 3.2, cex = 1.0)
  box()
}

png(file.path(out_dir, "Fig.S4_correlogram.png"),
    width = 6000, height = 6000, res = 300)
layout(matrix(c(1, 2), nrow = 1), widths = c(10, 1.2))
draw_correlogram()
draw_colorbar()
dev.off()

correlations <- do.call(rbind, lapply(seq_len(nvar), function(i) {
  do.call(rbind, lapply(i:nvar, function(j) {
    data.frame(variable_1 = labels[i], variable_2 = labels[j],
               r = r_mat[i, j], p = p_mat[i, j], n = n_mat[i, j])
  }))
}))

write.csv(correlations,
          file.path(out_dir, "Fig.S4_correlogram_statistics.csv"),
          row.names = FALSE)
message("Correlogram and statistics saved in: ", normalizePath(out_dir))
