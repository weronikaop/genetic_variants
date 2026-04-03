#!/usr/bin/env bash
# ==============================================================
#  Step 1: Environment and Software Installation
#  Description:
#     Configure software environment for evolutionary pattern analysis
# ==============================================================

set -e
set -o pipefail

echo "=============================================================="
echo "Step 1: Software and Environment Installation"
echo "=============================================================="
echo ""

# --- Step 1.1: Python environment setup ---
echo "[1/2] Setting up Python environment..."
bash ./step0-software-install-01.sh
echo "Python environment setup completed."
echo ""

# --- Step 1.2: R package installation ---
echo "[2/2] Installing R dependencies..."
Rscript ./step0-software-install-02.r
echo "R packages installed successfully."
echo ""

echo ""
echo "=============================================================="
echo "All software and dependencies have been successfully installed."
echo "You may now proceed to Step 2: ----."
echo "=============================================================="