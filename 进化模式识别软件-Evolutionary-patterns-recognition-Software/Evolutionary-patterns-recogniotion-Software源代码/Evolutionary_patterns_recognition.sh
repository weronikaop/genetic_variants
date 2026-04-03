#!/bin/bash
# ================================================
#  进化模式分析流程
# ================================================

set -euo pipefail

# ----------- 环境配置-----------
if [ -d ".venv" ]; then
    echo "[INFO] Activating local Python environment (.venv)"
    source .venv/bin/activate
else
    echo "[WARN] No .venv found — using system Python"
fi

# ----------- 用法-----------
usage() {
    echo ""
    echo "============================================================"
    echo " Evolutionary Pattern Recognition Pipeline"
    echo "============================================================"
    echo ""
    echo "Usage:"
    echo "  bash Evolutionary_patterns_recognition.sh -i <tree_files> -m <motif_file.xml> [options]"
    echo ""
    echo "Required arguments:"
    echo "  -i, --input               One or more phylogenetic tree files (.nwk)"
    echo "  -m, --motif               Motif definition file (.xml, e.g., MEME output)"
    echo ""
    echo "Optional arguments:"
    echo "  -n, --names               List of gene names corresponding to the tree files"
    echo "  -o, --output              Output directory (default: ./output)"
    echo "  -r, --reference           Reference sequence or annotation file (optional)"
    echo "  --format                  Output format: multi | tab | detailed (default: multi)"
    echo "  --validate                Validate motif structure and content before processing"
    echo "  --verbose                 Print detailed progress information"
    echo ""
    echo "Advanced parameters:"
    echo "  -zh, --z-high <float>     Upper Z-score threshold (default: 3.0)"
    echo "  -zl, --z-low <float>      Lower Z-score threshold (default: -2.0)"
    echo "  -ht, --hgt-threshold <f>  Horizontal gene transfer threshold (default: 5.0)"
    echo "  -cs, --converge-similarity <f>  Convergent similarity threshold (default: 0.7)"
    echo "  -ct, --converge-topology <f>   Convergent topology score threshold (default: 3.0)"
    echo ""
    echo "Example (single gene):"
    echo "  bash Evolutionary_patterns_recognition.sh -i treeA.nwk -m motifs.xml -n A"
    echo ""
    echo "Example (multiple genes):"
    echo "  bash Evolutionary_patterns_recognition.sh -i treeA.nwk treeB.nwk -n A B -m motifs.xml -o results/ --validate --verbose"
    echo ""
    echo "============================================================"
    exit 1
}

# ----------- 参数解析 -----------
TREE_FILES=()
GENE_NAMES=()
MOTIF_FILE=""
OUTPUT_DIR="output"
VALIDATE=false
VERBOSE=false
REFERENCE_FILE=""
OUTPUT_FORMAT="multi"
Z_HIGH=3.0
Z_LOW=-2.0
HGT_THRESHOLD=5.0
CONV_SIM=0.7
CONV_TOPO=3.0

while [[ $# -gt 0 ]]; do
  case $1 in
    -i|--input)
      shift
      while [[ $# -gt 0 && ! $1 =~ ^- ]]; do
        TREE_FILES+=("$1")
        shift
      done
      ;;
    -n|--names)
      shift
      while [[ $# -gt 0 && ! $1 =~ ^- ]]; do
        GENE_NAMES+=("$1")
        shift
      done
      ;;
    -m|--motif)
      MOTIF_FILE="$2"
      shift 2
      ;;
    -o|--output)
      OUTPUT_DIR="$2"
      shift 2
      ;;
    --validate)
      VALIDATE=true
      shift
      ;;
    --verbose)
      VERBOSE=true
      shift
      ;;
    -r|--reference)
      REFERENCE_FILE="$2"
      shift 2
      ;;
    --format)
      OUTPUT_FORMAT="$2"
      shift 2
      ;;
        -zh|--z-high)
      Z_HIGH="$2"
      shift 2
      ;;
    -zl|--z-low)
      Z_LOW="$2"
      shift 2
      ;;
    -ht|--hgt-threshold)
      HGT_THRESHOLD="$2"
      shift 2
      ;;
    -cs|--converge-similarity)
      CONV_SIM="$2"
      shift 2
      ;;
    -ct|--converge-topology)
      CONV_TOPO="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      usage
      exit 1
      ;;
  esac
done

# ----------- 检查必需参数 -----------
if [[ ${#TREE_FILES[@]} -eq 0 || -z "$MOTIF_FILE" ]]; then
    echo "[ERROR] Missing required arguments."
    usage
fi

mkdir -p "$OUTPUT_DIR"

echo "==============================================="
echo "Evolutionary Pattern Analysis Pipeline"
echo "Input trees: ${TREE_FILES[*]}"
if [[ ${#GENE_NAMES[@]} -gt 0 ]]; then
  echo "Gene names:  ${GENE_NAMES[*]}"
else
  echo "Gene names:  (inferred from filenames)"
fi
echo "Motif file:  $MOTIF_FILE"
echo "Output dir:  $OUTPUT_DIR"
echo "==============================================="

# ----------- Step 1: 拓扑编码 -----------
echo ""
echo "[STEP 1] Generating topology encoding..."
ENCODER_SCRIPT="step1-topology-encoder.py"

CMD=(python "$ENCODER_SCRIPT" -i "${TREE_FILES[@]}" -o "$OUTPUT_DIR")
if [[ ${#GENE_NAMES[@]} -gt 0 ]]; then
    CMD+=(-n "${GENE_NAMES[@]}")
fi

echo "Running: ${CMD[*]}"
"${CMD[@]}"

SUMMARY_FILE="${OUTPUT_DIR}/topology_codes_summary.txt"
if [[ ! -f "$SUMMARY_FILE" ]]; then
  echo "[ERROR] Topology encoding failed — summary not found at $SUMMARY_FILE"
  deactivate 2>/dev/null || true
  exit 1
fi
echo "Topology encoding complete: $SUMMARY_FILE"

# ----------- Step 2: 模序二进制编码 -----------
echo ""
echo "[STEP 2] Generating motif binary encoding..."
MOTIF_OUT="${OUTPUT_DIR}/motif_binary_encoding.txt"
MOTIF_SCRIPT="step2-meme-encoder.py"

CMD2=(python "$MOTIF_SCRIPT" -i "$MOTIF_FILE" -o "$MOTIF_OUT")
if $VALIDATE; then
    CMD2+=(--validate)
fi
if $VERBOSE; then
    CMD2+=(--verbose)
fi

echo "Running: ${CMD2[*]}"
"${CMD2[@]}"

if [[ ! -f "$MOTIF_OUT" ]]; then
  echo "[ERROR] Motif encoding failed — $MOTIF_OUT not found"
  deactivate 2>/dev/null || true
  exit 1
fi

echo "Motif encoding complete: $MOTIF_OUT"

# ----------- 遍历各基因/进化树 -----------
for idx in "${!TREE_FILES[@]}"; do
    TREE="${TREE_FILES[$idx]}"
    NAME="${GENE_NAMES[$idx]:-$(basename "$TREE" .nwk)}"
    GENE_OUT="${OUTPUT_DIR}/${NAME}"
    mkdir -p "$GENE_OUT"

    echo "[INFO] Processing gene: $NAME"

    TOPO_FILE="${OUTPUT_DIR}/${NAME}_topology_codes.txt"
    if [[ ! -f "$TOPO_FILE" ]]; then
        echo "[ERROR] Topology encoding file not found: $TOPO_FILE"
        continue
    fi

    # Step 3: 合并基因（拓扑代码和模序代码）
    python step3-merge-topology-motif.py \
        -t "$TOPO_FILE" \
        -m "${OUTPUT_DIR}/motif_binary_encoding.txt" \
        -o "${GENE_OUT}/merged_result.txt" \
        $([ -n "$REFERENCE_FILE" ] && echo "-r $REFERENCE_FILE") \
        --format "$OUTPUT_FORMAT" \
        $([ "$VERBOSE" = true ] && echo "--verbose")

    # Step 4: 进化模式识别- 总结
    python step4-summary.py -i "${GENE_OUT}/merged_result.txt" -m 2 3 \
        -zh "$Z_HIGH" -zl "$Z_LOW" -ht "$HGT_THRESHOLD" -cs "$CONV_SIM" -ct "$CONV_TOPO" \
        > "${GENE_OUT}/evolutionary_summary.txt"

    # Step 5: 可视化 (每个基因)
    if command -v Rscript >/dev/null 2>&1; then
        Rscript step5-visualization.R "$GENE_OUT"
    fi
done


# ----------- 完成-----------
echo ""
echo "All analyses and visualizations completed successfully."
echo "Results stored in: $OUTPUT_DIR"
echo "==============================================="

deactivate 2>/dev/null || true
