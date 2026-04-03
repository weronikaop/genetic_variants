#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# Extract variant lines from a VCF (or tab-delimited) file for a specified chromosome region.
# Supports: plain VCF, gzipped VCF (.gz), stdin ("-"), preserves VCF header lines by default.
# Usage examples at the end of the file.

show_help() {
  cat <<'USAGE' >&2
Usage: extract.sh [OPTIONS]

Options:
  -c, --chrom CHROM       Chromosome name (e.g. chr14 or 14)
  -s, --start START       Start position (integer)
  -e, --end END           End position (integer)
  -r, --region REGION     Region in one arg, e.g. chr14:105000-106000 (alternative to -c/-s/-e)
  -i, --input INPUT       Input file (use "-" for stdin). Supports .gz
  -o, --output OUTPUT     Output file (won’t be overwritten unless --force)
  -f, --force             Overwrite output if exists
  --no-header             Do NOT copy header lines starting with '#'
  -h, --help              Show this help
Examples:
  extract.sh -r chr14:105586437-106879843 -i input.vcf.gz -o chr14_snv.vcf.gz
  zcat input.vcf.gz | extract.sh -c 14 -s 1000 -e 2000 -i - -o out.vcf
USAGE
  exit 1
}

# Defaults
CHROM=""
START=""
END=""
REGION=""
INPUT=""
OUTPUT=""
FORCE=0
COPY_HEADER=1

# Parse args (simple loop)
while [[ $# -gt 0 ]]; do
  case "$1" in
    -c|--chrom) CHROM="$2"; shift 2 ;;
    -s|--start) START="$2"; shift 2 ;;
    -e|--end) END="$2"; shift 2 ;;
    -r|--region) REGION="$2"; shift 2 ;;
    -i|--input) INPUT="$2"; shift 2 ;;
    -o|--output) OUTPUT="$2"; shift 2 ;;
    -f|--force) FORCE=1; shift ;;
    --no-header) COPY_HEADER=0; shift ;;
    -h|--help) show_help ;;
    *) echo "Unknown parameter: $1" >&2; show_help ;;
  esac
done

# Parse region if provided
if [[ -n "$REGION" ]]; then
  if [[ "$REGION" =~ ^([^:]+):([0-9]+)-([0-9]+)$ ]]; then
    CHROM="${BASH_REMATCH[1]}"
    START="${BASH_REMATCH[2]}"
    END="${BASH_REMATCH[3]}"
  else
    echo "Error: region must be in format CHR:START-END" >&2
    exit 2
  fi
fi

# Validate required params
if [[ -z "${CHROM}" || -z "${START}" || -z "${END}" || -z "${INPUT}" || -z "${OUTPUT}" ]]; then
  echo "Error: missing required parameters" >&2
  show_help
fi

# Validate numeric positions
if ! [[ "$START" =~ ^[0-9]+$ ]] || ! [[ "$END" =~ ^[0-9]+$ ]]; then
  echo "Error: start and end must be positive integers" >&2
  exit 3
fi
if (( START > END )); then
  echo "Error: start > end" >&2
  exit 4
fi

# Check input existence (unless stdin "-")
if [[ "$INPUT" != "-" && ! -e "$INPUT" ]]; then
  echo "Error: input file does not exist: $INPUT" >&2
  exit 5
fi

# Check output exists
if [[ -e "$OUTPUT" && $FORCE -ne 1 ]]; then
  echo "Error: output file exists. Use --force to overwrite: $OUTPUT" >&2
  exit 6
fi

# Helper to open input (supports gz and stdin)
_open_input_cmd() {
  if [[ "$INPUT" == "-" ]]; then
    cat -
  elif [[ "$INPUT" == *.gz ]]; then
    gzip -cd -- "$INPUT"
  else
    cat -- "$INPUT"
  fi
}

# Main extraction
# - preserves header lines by default
# - compares chrom names tolerant to leading "chr"
# - numeric compare on POS
( _open_input_cmd ) | awk -F'\t' -v chrom="$CHROM" -v start="$START" -v end="$END" -v copy_header="$COPY_HEADER" '
  BEGIN {
    # canonical chrom (no "chr" prefix) used for comparison
    c = chrom
    if (c ~ /^chr/) c = substr(c,4)
  }
  # header handling
  /^#/ {
    if (copy_header == 1) print
    next
  }
  {
    # field 1 may contain "chr" prefix
    f = $1
    if (f ~ /^chr/) nf = substr(f,4); else nf = f
    pos = $2 + 0
    if (nf == c && pos >= start && pos <= end) print
  }
' > "$OUTPUT"

# Report results
if [[ -s "$OUTPUT" ]]; then
  count=$(wc -l < "$OUTPUT" | tr -d ' ')
  echo "OK: extracted $count lines to $OUTPUT" >&2
  exit 0
else
  echo "Warning: no matching lines found; created empty file $OUTPUT" >&2
  exit 0
fi
