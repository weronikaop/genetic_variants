# ==============================================================
# Step 1.2: R包安装
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

# 检查CRAN镜像
options(repos = c(CRAN = "https://cloud.r-project.org"))

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

cat("\nR environment successfully configured.\n")
cat("Installed packages:\n")
print(installed.packages()[required_packages, "Version"])
