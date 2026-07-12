# ==============================================================================
# Script: 02_run_lions.sh
# Purpose:
#   1. Pre-process reference files (remove 'chr' prefix for LIONS compatibility).
#   2. Run LIONS chimericReadSearch.py across multiple sequencing depths.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Define Paths (Dynamic Absolute Paths)
# ------------------------------------------------------------------------------
WORK_DIR="$(pwd)"
BASE="${WORK_DIR}/results"

# 假设 BAM 文件统一存放在 sortbam 目录下
BAM_DIR="${BASE}/sortbam" 
OUTBASE="${BASE}/LIONS_official_by_depth"
REF_DIR="${OUTBASE}/ref"

# IMPORTANT: Path to LIONS installation
LIONS_BIN="${WORK_DIR}/software/LIONS-master/scripts/ChimericReadTool/chimericReadSearch.py"

# Input References
EXON_RAW="${REF_DIR}/official_simulated_reference_90pct.LIONS.exons.clean"
EXON_NOCHR="${REF_DIR}/official_simulated_reference_90pct.LIONS.exons.nochr.clean"
REPEAT="${REF_DIR}/mm10_TE.LIONS.7col"

DEPTHS=("5x" "10x" "25x" "50x" "100x") # 补齐了 100x

mkdir -p "${OUTBASE}" "${REF_DIR}"

# ------------------------------------------------------------------------------
# 1. Clean Chromosome Names in Exon Reference
# ------------------------------------------------------------------------------
echo ">>> Step 1: Removing 'chr' prefix from exon reference..."

test -s "${EXON_RAW}" || { echo "Error: EXON_RAW not found: ${EXON_RAW}"; exit 1; }

awk 'BEGIN{OFS="\t"} {$3=sub(/^chr/,"",$3) ? $3 : $3; print}' \
  "${EXON_RAW}" > "${EXON_NOCHR}"

# Sanity check for 9 columns
awk 'NF != 9 {print "Bad line: ", NR, "NF="NF, $0; exit 1}' "${EXON_NOCHR}"

echo "Exon reference cleaned successfully."

# ------------------------------------------------------------------------------
# 2. Run LIONS Pipeline Loop
# ------------------------------------------------------------------------------
echo ">>> Step 2: Running LIONS chimericReadSearch..."

test -s "${LIONS_BIN}" || { echo "Error: LIONS script not found: ${LIONS_BIN}"; exit 1; }
test -s "${EXON_NOCHR}" || { echo "Error: EXON_NOCHR not found."; exit 1; }
test -s "${REPEAT}" || { echo "Error: REPEAT file not found: ${REPEAT}"; exit 1; }

for depth in "${DEPTHS[@]}"; do
  echo "========== Processing Depth: ${depth} =========="

  BAM="${BAM_DIR}/official_simulated_${depth}_noise0.1.bam"
  OUTDIR="${OUTBASE}/${depth}"
  OUTBED="${OUTDIR}/official_simulated_${depth}_noise0.1.LIONS.bed"
  LOG="${OUTDIR}/official_simulated_${depth}_noise0.1.LIONS.log"

  mkdir -p "${OUTDIR}"

  test -s "${BAM}" || { echo "Error: BAM not found for ${depth}: ${BAM}"; continue; }

  # Remove possible old output to avoid appending issues
  rm -f "${OUTBED}" "${LOG}"

  # Execute LIONS (Note: LIONS is typically a Python 2 script, adjust 'python' command if needed)
  python "${LIONS_BIN}" \
    "${EXON_NOCHR}" \
    "${REPEAT}" \
    "${BAM}" \
    "${OUTBED}" \
    > "${LOG}" 2>&1

  echo "[DONE] ${depth}"
  ls -lh "${OUTBED}"

  # Generate Output Summary
  {
    echo -e "tool\tdepth\toutput_bed\tn_lines"
    echo -e "LIONS\t${depth}\t${OUTBED}\t$(wc -l < "${OUTBED}")"
  } > "${OUTDIR}/LIONS_${depth}_output_summary.tsv"

  cat "${OUTDIR}/LIONS_${depth}_output_summary.tsv"
done

echo ">>> All LIONS runs completed successfully!"