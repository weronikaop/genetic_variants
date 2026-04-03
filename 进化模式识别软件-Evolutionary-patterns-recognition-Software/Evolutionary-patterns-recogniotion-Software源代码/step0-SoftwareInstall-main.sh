#!/usr/bin/env bash
# ==============================================================
#  Step 1: 环境与软件安装
#  说明:
#     配置软件环境以进行进化模式分析
# ==============================================================

set -e
set -o pipefail

echo "=============================================================="
echo "Step 1: Software and Environment Installation"
echo "=============================================================="
echo ""

# --- Step 1.1: Python环境配置 ---
echo "[1/2] Setting up Python environment..."
bash ./step0-software-install-01.sh
echo "Python environment setup completed."
echo ""

# --- Step 1.2: R包安装 ---
echo "[2/2] Installing R dependencies..."
Rscript ./step0-software-install-02.r
echo "R packages installed successfully."
echo ""

echo ""
echo "=============================================================="
echo "All software and dependencies have been successfully installed."
echo "You may now proceed to Step 2: ----."
echo "=============================================================="