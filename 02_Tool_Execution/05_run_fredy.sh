# ==============================================================================
# Script: 05_run_fredy.sh
# Purpose: 
#   1. Create a CDS-proxy GTF required by FREDY chimeric module.
#   2. Run FREDY star, string, and chimeric modules across simulated depths.
# ==============================================================================

export PATH="$HOME/bin:$PATH"

# ------------------------------------------------------------------------------
# 0. Define Paths (Dynamic Absolute Paths for Docker/Udocker compatibility)
# ------------------------------------------------------------------------------
WORK_DIR="$(pwd)"
DATA_DIR="${WORK_DIR}/data"
RESULTS_DIR="${WORK_DIR}/results"

BASE="${RESULTS_DIR}"
OUTBASE="${BASE}/FREDY_official_by_depth"
REF_DIR="${OUTBASE}/ref"

GENOME="${DATA_DIR}/mm10_no_alt_analysis_set_ENCODE.fasta"
REF_GTF="${BASE}/official_simulated_reference_90pct.gtf"
STAR_INDEX="${BASE}/STAR_index_official_reference_90pct"

FREDY_CDS_GTF="${REF_DIR}/official_simulated_reference_90pct.FREDY_CDSproxy.gtf"
TE_BED="${REF_DIR}/mm10_TE.FREDY.bed4" # From 01c_prepare_fredy_ref.R

THREADS=8
DEPTHS=("5x" "10x" "25x" "50x" "100x")

mkdir -p "${REF_DIR}"

# ------------------------------------------------------------------------------
# 1. Create FREDY CDS proxy GTF
# FREDY requires CDS features for chimeric transcript detection. This step 
# forces 'exon' features to act as 'CDS' proxies.
# ------------------------------------------------------------------------------
echo ">>> Step 1: Generating FREDY CDS proxy GTF..."

awk 'BEGIN{OFS="\t"}
     /^#/ {print; next}
     {print}
     $3=="exon" {
       $3="CDS";
       print
     }' \
  "${REF_GTF}" \
  > "${FREDY_CDS_GTF}"

echo "CDS proxy GTF generated successfully."

# ------------------------------------------------------------------------------
# 2. Run FREDY Pipeline Loop
# ------------------------------------------------------------------------------
echo ">>> Step 2: Running FREDY star, string, and chimeric modules..."

for depth in "${DEPTHS[@]}"; do
  echo "========== Processing Depth: ${depth} =========="
  
  out_dir="${OUTBASE}/${depth}"
  mkdir -p "${out_dir}"

  fq1="${BASE}/fastq/official_simulated_${depth}_noise0.1.R1.fastq.gz"
  fq2="${BASE}/fastq/official_simulated_${depth}_noise0.1.R2.fastq.gz"

  # Create input list for FREDY
  printf "%s\n%s\n" "${fq1}" "${fq2}" > "${out_dir}/fastq_files.txt"

  echo "[$(date)] FREDY star: ${depth}"
  fredy_udocker star \
    -o "${out_dir}" \
    -i "${STAR_INDEX}" \
    -f "${out_dir}/fastq_files.txt" \
    -S -p \
    -t "${THREADS}"

  echo "[$(date)] FREDY string: ${depth}"
  fredy_udocker string \
    -o "${out_dir}" \
    -a "${REF_GTF}" \
    -t "${THREADS}"

  echo "[$(date)] FREDY chimeric: ${depth}"
  fredy_udocker chimeric \
    -o "${out_dir}" \
    -a "${FREDY_CDS_GTF}" \
    -g "${GENOME}" \
    -e "${TE_BED}"

  echo "[$(date)] DONE ${depth}"
done

echo ">>> All FREDY pipeline runs completed successfully!"