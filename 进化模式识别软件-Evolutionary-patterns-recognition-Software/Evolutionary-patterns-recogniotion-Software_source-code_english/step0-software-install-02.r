# ==============================================================
# Step 1.2: R Package Installation
# ==============================================================

required_packages <- c(
  "ggplot2",
  "ggdendro",
  "igraph",
  "gridExtra",
  "stringr",
  "dplyr",
  "RColorBrewer"
)

cat("Installing R dependencies...\n")

# Check CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

installed <- installed.packages()[, "Package"]

# Install missing packages
for (pkg in required_packages) {
  if (!(pkg %in% installed)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

cat("\nR environment successfully configured.\n")
cat("Installed packages:\n")
print(installed.packages()[required_packages, "Version"])
