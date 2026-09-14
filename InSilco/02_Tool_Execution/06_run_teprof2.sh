#!/usr/bin/env bash
set -euo pipefail

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
TEPROF2_ADAPTER_DIR="${WORK_DIR}/InSilco/02_Tool_Execution/adapters/teprof2"
PREPARE_DICTIONARY_INPUT="${TEPROF2_ADAPTER_DIR}/prepare_teprof2_dictionary_input.py"
GENECODE_TO_DIC="${TEPROF2_ADAPTER_DIR}/genecode_to_dic_simulation.py"
NORMALIZE_DICTIONARIES="${TEPROF2_ADAPTER_DIR}/normalize_teprof2_dictionaries.py"

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
test -s "${REFGTF}" || { echo "Error: Missing ${REFGTF}"; exit 1; }
test -s "${PREPARE_DICTIONARY_INPUT}" || { echo "Error: Missing ${PREPARE_DICTIONARY_INPUT}"; exit 1; }
test -s "${GENECODE_TO_DIC}" || { echo "Error: Missing ${GENECODE_TO_DIC}"; exit 1; }
test -s "${NORMALIZE_DICTIONARIES}" || { echo "Error: Missing ${NORMALIZE_DICTIONARIES}"; exit 1; }

# Sort and index RepeatMasker BED
sort -k1,1 -k2,2n "${RMSK_BED}" > "${REF_DIR}/rmsk_mm10_from_rds.sorted.bed"
bgzip -f "${REF_DIR}/rmsk_mm10_from_rds.sorted.bed"
tabix -f -p bed "${REF_DIR}/rmsk_mm10_from_rds.sorted.bed.gz"

# Format and validate the simulation GTF for the TEProf2 Python 2 parser.
python3 "${PREPARE_DICTIONARY_INPUT}" \
"${REFGTF}" \
"${REF_DIR}/official_simulated_reference_90pct.teprof2_input.sorted.gtf"

# Generate genecode dictionaries
(
cd "${REF_DIR}"
conda run -n "${TEPROF2_PY2_ENV}" python2 "${GENECODE_TO_DIC}" \
"official_simulated_reference_90pct.teprof2_input.sorted.gtf"
)

mv -f "${REF_DIR}/genecode_plus.dic"  "${REF_DIR}/official_simulated_reference_90pct.plus.dic"
mv -f "${REF_DIR}/genecode_minus.dic" "${REF_DIR}/official_simulated_reference_90pct.minus.dic"

conda run -n "${TEPROF2_PY2_ENV}" python2 "${NORMALIZE_DICTIONARIES}" \
"${REF_DIR}/official_simulated_reference_90pct.plus.dic" \
"${REF_DIR}/official_simulated_reference_90pct.minus.dic"

# Write arguments configuration file
printf 'rmsk	%s
rmskannotationfile	%s
gencodeplusdic	%s
gencodeminusdic	%s
' \
"${REF_DIR}/rmsk_mm10_from_rds.sorted.bed.gz" \
"${RMSK_LST}" \
"${REF_DIR}/official_simulated_reference_90pct.plus.dic" \
"${REF_DIR}/official_simulated_reference_90pct.minus.dic" \
> "${ARG_FILE}"

awk -F '	' '
NF != 2 || $1 == "" || $2 == "" {
    print "Invalid TEProf2 argument line:", NR > "/dev/stderr"
    bad = 1
}
END { exit bad }
' "${ARG_FILE}"

while IFS=$'	' read -r key path; do
    test -s "${path}" || {
        echo "Error: Missing ${key} reference: ${path}"
        exit 1
    }
done < "${ARG_FILE}"

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

test -f "${ASM_GTF}_annotated_filtered_test_all" || {
    echo "Error: TEProf2 annotation failed at ${depth}"
    exit 1
}

echo "[DONE] ${depth}"
done

echo ">>> All TEProf2 runs completed successfully!"
