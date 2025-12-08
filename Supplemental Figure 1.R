# ============================================================
# Script: HLA_histogram_KDE.R
# Purpose: Generate histogram + KDE plot for HLA methylation data
# This script reproduces the analysis used for Figure X in the manuscript.
# ============================================================

# ========== 1. Load libraries ==========
library(readxl)
library(ggplot2)
library(dplyr)

# ========== 2. File paths (relative paths for reproducibility) ==========
data_path <- "data/HLA_BISseq.xlsx"
save_path <- "figures"

# Create output directory if not present
if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}

# ========== 3. Load data ==========
hla_data <- read_excel(data_path, sheet = 1)

met_ave_data <- hla_data$average
met_ave_data <- met_ave_data[!is.na(met_ave_data)]   # remove NA values

# ========== 4. Basic statistics ==========
cat("Number of samples:", length(met_ave_data), "\n")
cat("Mean:", round(mean(met_ave_data), 2), "\n")
cat("Median:", round(median(met_ave_data), 2), "\n")
cat("SD:", round(sd(met_ave_data), 2), "\n")
cat("Min:", round(min(met_ave_data), 2), "\n")
cat("Max:", round(max(met_ave_data), 2), "\n")

# ========== 5. KDE calculation parameters ==========
bandwidth <- 0.05
grid_step <- 0.01
x_grid <- seq(0, 1, by = grid_step)

kde_result <- density(
  met_ave_data,
  bw = bandwidth,
  from = 0,
  to = 1,
  n = length(x_grid)
)
kde_values <- approx(kde_result$x, kde_result$y, xout = x_grid)$y

# ========== 6. Histogram + KDE plot ==========
p1 <- ggplot(data.frame(met_ave = met_ave_data), aes(x = met_ave)) +
  geom_histogram(aes(y = ..density..), bins = 20,
                 alpha = 0.7, fill = "#69b3a2",
                 color = "white", size = 0.5) +
  geom_line(data = data.frame(x = x_grid, density = kde_values),
            aes(x = x, y = density), color = "#ff6b6b", size = 1.2) +
  labs(title = "Histogram + KDE",
       x = "Average methylation",
       y = "Density") +
  theme_minimal() +
  scale_x_continuous(limits = c(0,1))

# Print plot to screen
print(p1)

# ========== 7. Local minima detection ==========
find_local_minima <- function(x, y) {
  idx <- which(diff(sign(diff(y))) == 2) + 1
  x[idx]
}

local_minima <- sort(find_local_minima(x_grid, kde_values))
local_minima <- local_minima[local_minima >= 0 & local_minima <= 1]

cat("\nDetected local minima:\n")
if (length(local_minima) == 0) {
  cat("None detected.\n")
} else {
  for (i in seq_along(local_minima)) {
    cat("Min", i, ":", round(local_minima[i], 2), "\n")
  }
}

# ========== 8. Plot with minima points ==========
if (length(local_minima) > 0) {
  minima_y <- approx(x_grid, kde_values, xout = local_minima)$y
  
  p2 <- p1 +
    geom_point(data = data.frame(x = local_minima, y = minima_y),
               aes(x = x, y = y),
               shape = 21, color = "red", size = 3, fill = "white") +
    geom_text(data = data.frame(x = local_minima, y = minima_y),
              aes(x = x, y = y, label = paste0("Min:", round(x,1))),
              vjust = -1.5, color = "red", fontface = "bold")
  
  print(p2)
}

# ========== 9. Save plots as PDF ==========
ggsave(file.path(save_path, "HLA_histogram_KDE.pdf"),
       plot = p1, width = 10, height = 6, dpi = 300)

if (exists("p2")) {
  ggsave(file.path(save_path, "HLA_histogram_KDE_minima.pdf"),
         plot = p2, width = 10, height = 6, dpi = 300)
}

# ========== 10. Save results as text file ==========
output_text <- c(
  "HLA methylation KDE analysis",
  "============================",
  paste("Number of samples:", length(met_ave_data)),
  paste("Mean:", round(mean(met_ave_data), 2)),
  paste("Median:", round(median(met_ave_data), 2)),
  paste("SD:", round(sd(met_ave_data), 2)),
  "",
  "Detected local minima:"
)

if (length(local_minima) > 0) {
  for (i in seq_along(local_minima)) {
    output_text <- c(output_text,
                     paste("  Min", i, ":", round(local_minima[i], 2)))
  }
} else {
  output_text <- c(output_text, "  None detected.")
}

writeLines(output_text,
           file.path(save_path, "HLA_KDE_minima_report.txt"))

cat("\n=== Completed ===\n")
cat("Output saved to:", save_path, "\n")
