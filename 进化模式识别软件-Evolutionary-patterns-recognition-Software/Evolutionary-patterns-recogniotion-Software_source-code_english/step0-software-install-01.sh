#!/usr/bin/env bash
# ==============================================================
#  Step 1.1: Python Environment Setup
# ==============================================================

set -e
set -o pipefail

echo "Installing Python dependencies for evolutionary pattern analysis..."
echo ""

# --- Check Python version ---
PYVER=$(python3 -V 2>&1)
echo "Detected Python version: $PYVER"
echo ""

# --- Create virtual environment (recommended) ---
if [ ! -d ".venv" ]; then
    echo "Creating local virtual environment (.venv)..."
    python3 -m venv .venv
fi

source .venv/bin/activate

# --- Upgrade pip ---
pip install --upgrade pip setuptools wheel

# --- Install required Python libraries ---
pip install numpy six ete3

# --- Fix for Python >=3.13 missing cgi module ---
pip install legacy-cgi

echo ""
echo "Python packages installed successfully:"
pip list | grep -E "numpy|six|ete3|legacy-cgi"
echo ""

echo "Python environment setup completed."
deactivate
