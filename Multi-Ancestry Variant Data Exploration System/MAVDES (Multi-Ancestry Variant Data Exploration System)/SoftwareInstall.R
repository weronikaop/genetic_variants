# ==============================================================
# Step 0: Install required R packages
# ==============================================================

required_packages <- c(
  "optparse", "GenomicRanges", "ggplot2", "viridis", "tidyverse",
  "readxl", "data.table", "reshape2", "VariantAnnotation", "dplyr",
  "tibble", "RColorBrewer", "tidyr", "scales", "purrr", "cowplot",
  "ggnewscale", "stringr", "patchwork", "ggtext", "forcats",
  "ggrepel", "pheatmap", "Biostrings", "MASS", "haplo.stats", "readr"
)

cat("Installing R dependencies...\n\n")

# Set CRAN mirror
installed <- installed.packages()[, "Package"]

# Install missing packages
for (pkg in required_packages) {
  if (!(pkg %in% installed)) {
    cat("Installing package:", pkg, "...\n")
    install.packages(pkg, dependencies = TRUE)
  }
}

cat("\nR environment successfully configured.\n")
cat("Installed package versions:\n")

installed <- installed.packages()

versions <- installed[rownames(installed) %in% required_packages, "Version"]
print(versions)
