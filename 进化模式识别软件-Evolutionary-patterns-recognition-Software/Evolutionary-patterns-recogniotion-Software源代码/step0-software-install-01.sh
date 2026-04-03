#!/usr/bin/env bash
# ==============================================================
#  Step 1.1: Python 环境配置
# ==============================================================

set -e
set -o pipefail

echo "Installing Python dependencies for evolutionary pattern analysis..."
echo ""

# --- 检查Python版本 ---
PYVER=$(python3 -V 2>&1)
echo "Detected Python version: $PYVER"
echo ""

# --- 创建虚拟环境（推荐）---
if [ ! -d ".venv" ]; then
    echo "Creating local virtual environment (.venv)..."
    python3 -m venv .venv
fi

source .venv/bin/activate

# --- 升级pip ---
pip install --upgrade pip setuptools wheel

# ---安装必须Python包 ---
pip install numpy six ete3

# --- 针对Python >=3.13 中缺失cgi模块的修复方案 ---
pip install legacy-cgi

echo ""
echo "Python packages installed successfully:"
pip list | grep -E "numpy|six|ete3|legacy-cgi"
echo ""

echo "Python environment setup completed."
deactivate
