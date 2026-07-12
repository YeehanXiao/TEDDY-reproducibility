# ==============================================================================
# Script: 06_run_teprof2.sh
# Purpose: 
#   1. Prepare TEProf2-specific reference indices and dictionaries (Python 2).
#   2. Run StringTie assembly for each sequencing depth.
#   3. Annotate assembled GTFs using TEProf2 to identify chimeric events.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Define Paths & Parameters (Dynamic Absolute Paths)
# ------------------------------------------------------------------------------
WORK_DIR="$(pwd)"
RESULTS_DIR="${WORK_DIR}/results"
BAM_DIR="${RESULTS_DIR}/sortbam"

OUTBASE="${RESULTS_DIR}/TEProf2_benchmark"
REF_DIR="${OUTBASE}/ref"
LOG_DIR="${OUTBASE}/logs"

# IMPORTANT: Path to TEProf2 installation
TEPROF2_BIN_DIR="${WORK_DIR}/software/TEProf2Paper/bin"

# Reference Files
REFGTF="${RESULTS_DIR}/official_simulated_reference_90pct.gtf"
SIMGTF="${REF_DIR}/official_simulated_reference_90pct.teprof2_input.gtf"
RMSK_BED="${REF_DIR}/rmsk_mm10_from_rds.bed"
RMSK_LST="${REF_DIR}/repeatmasker_description_uniq_mm10_from_rds.lst"

# Arguments file to be generated
ARG_FILE="${REF_DIR}/arguments_mm10_simulation.txt"

# Conda environment names
TEPROF2_ENV="teprof2"
TEPROF2_PY2_ENV="teprof2_py2"

DEPTHS=("5x" "10x" "25x" "50x" "100x")

mkdir -p "${REF_DIR}" "${LOG_DIR}"

# ------------------------------------------------------------------------------
# 1. Reference Preparation (Tabix & Dictionaries)
# ------------------------------------------------------------------------------
echo ">>> Step 1: Preparing TEProf2 Reference Indices..."

test -s "${RMSK_BED}" || { echo "Error: Missing ${RMSK_BED}"; exit 1; }
test -s "${RMSK_LST}" || { echo "Error: Missing ${RMSK_LST}"; exit 1; }
test -s "${SIMGTF}" || { echo "Error: Missing ${SIMGTF}"; exit 1; }

# Sort and index RepeatMasker BED
sort -k1,1 -k2,2n "${RMSK_BED}" > "${REF_DIR}/rmsk_mm10_from_rds.sorted.bed"
bgzip -f "${REF_DIR}/rmsk_mm10_from_rds.sorted.bed"
tabix -f -p bed "${REF_DIR}/rmsk_mm10_from_rds.sorted.bed.gz"

# Format GTF for TEProf2 Python 2 parser
awk 'BEGIN{OFS="\t"} $3=="transcript" || $3=="exon" || $3=="start_codon" {print $0}' "${SIMGTF}" \
| awk -F '; ' '{print $0"\t"$2}' \
> "${REF_DIR}/official_simulated_reference_90pct.teprof2_input.sorted.gtf"

# Generate genecode dictionaries
conda run -n "${TEPROF2_PY2_ENV}" python2 "${TEPROF2_BIN_DIR}/genecode_to_dic.py" \
"${REF_DIR}/official_simulated_reference_90pct.teprof2_input.sorted.gtf"

mv -f "${REF_DIR}/genecode_plus.dic"  "${REF_DIR}/official_simulated_reference_90pct.plus.dic"
mv -f "${REF_DIR}/genecode_minus.dic" "${REF_DIR}/official_simulated_reference_90pct.minus.dic"

# Write arguments configuration file
cat > "${ARG_FILE}" <<EOF
rmsk ${REF_DIR}/rmsk_mm10_from_rds.sorted.bed.gz
rmskannotationfile ${RMSK_LST}
gencodeplusdic ${REF_DIR}/official_simulated_reference_90pct.plus.dic
gencodeminusdic ${REF_DIR}/official_simulated_reference_90pct.minus.dic
EOF

# ------------------------------------------------------------------------------
# 2. Environment Sanity Check
# ------------------------------------------------------------------------------
echo ">>> Step 2: Checking Python 2 and tabix environment..."
conda run -n "${TEPROF2_PY2_ENV}" python2 - <<'PY'
import tabix
print("tabix_ok: Python 2 environment is ready.")
PY

# ------------------------------------------------------------------------------
# 3. Execution Loop (StringTie + TEProf2 Annotation)
# ------------------------------------------------------------------------------
echo ">>> Step 3: Running StringTie and TEProf2..."

for depth in "${DEPTHS[@]}"; do
echo "========== Processing Depth: ${depth} =========="

OUTDIR="${OUTBASE}/${depth}_test"
BAM="${BAM_DIR}/official_simulated_${depth}_noise0.1.bam"
ASM_GTF="${OUTDIR}/assembly/official_simulated_${depth}_noise0.1.stringtie.gtf"

mkdir -p "${OUTDIR}/assembly" "${OUTDIR}/logs"

test -s "${BAM}" || { echo "Error: Missing BAM: ${BAM}"; exit 1; }

# A. StringTie Assembly
conda run -n "${TEPROF2_ENV}" stringtie "${BAM}" \
-o "${ASM_GTF}" \
-p 8 \
-G "${REFGTF}" \
-m 100 \
-c 1 \
-f 0.01 \
-j 1 \
2>&1 | tee "${OUTDIR}/logs/stringtie_${depth}.log"

test -s "${ASM_GTF}" || { echo "Error: StringTie failed at ${depth}"; exit 1; }

# Clean old annotations
rm -f "${ASM_GTF}_annotated_test_all" "${ASM_GTF}_annotated_filtered_test_all"

# B. TEProf2 Annotation
conda run -n "${TEPROF2_PY2_ENV}" python2 "${TEPROF2_BIN_DIR}/rmskhg38_annotate_gtf_update_test_tpm.py" \
"${ASM_GTF}" "${ARG_FILE}" \
2>&1 | tee "${OUTDIR}/logs/teprof2_annotate_${depth}.log"

echo "[DONE] ${depth}"
done

echo ">>> All TEProf2 runs completed successfully!"