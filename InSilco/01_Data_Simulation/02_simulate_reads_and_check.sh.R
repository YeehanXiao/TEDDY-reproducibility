set -euo pipefail

# ==============================================================================
# Script: 02_simulate_reads_and_check.sh
# Purpose:
#   1. Generate STAR and RSEM indices.
#   2. Simulate RNA-seq reads at varying depths (5x, 10x, 25x, 50x, 100x).
#   3. Perform sanity check quantification to validate simulated expressions.
#   4. Compress output FASTQ files for downstream WDL pipelines.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Define Paths and Parameters
# ------------------------------------------------------------------------------
# Directories
DATA_DIR="./data"
RESULTS_DIR="./results"
FASTQ_DIR="${RESULTS_DIR}/fastq"
CHECK_DIR="${RESULTS_DIR}/quant_check"
LOG_DIR="${RESULTS_DIR}/logs"

# Ensure output directories exist
mkdir -p "$FASTQ_DIR" "$CHECK_DIR" "$LOG_DIR"

# Input reference files
GENOME="${DATA_DIR}/mm10_no_alt_analysis_set_ENCODE.fasta"
GTF_REF90="${RESULTS_DIR}/official_simulated_reference_90pct.gtf"
GTF_RSEM_CLEAN="${RESULTS_DIR}/official_simulated_truth_1000_transcripts.rsem_clean.gtf"
ISO="${RESULTS_DIR}/official_simulated_1000.isoforms.results"

# Background model for simulation (Ensure this is provided in your public data repository)
MODEL="${DATA_DIR}/forStat.model"

# Index output paths
STAR_INDEX="${RESULTS_DIR}/STAR_index_official_reference_90pct"
RSEM_REF_PREFIX="${RESULTS_DIR}/mm10_official_simulate_1000_truth_reference"

# Tool paths (Assuming tools are in $PATH, adjust if necessary)
RSEM_BIN_DIR="" # e.g., "/path/to/RSEM/" if not globally installed. Leave empty if in PATH.

# Simulation parameters
READ_LEN=150
THETA0=0.1
DEPTHS=(5 10 25 50 100)

# ------------------------------------------------------------------------------
# 2. Build STAR Index
# ------------------------------------------------------------------------------
echo ">>> Step 1: Building STAR Index..."
rm -rf "$STAR_INDEX"
mkdir -p "$STAR_INDEX"

STAR \
--runThreadN 40 \
--runMode genomeGenerate \
    --limitGenomeGenerateRAM "${STAR_INDEX_RAM:-100000000000}" \
--genomeDir "$STAR_INDEX" \
--genomeFastaFiles "$GENOME" \
--sjdbGTFfile "$GTF_REF90" \
--sjdbOverhang 149

# ------------------------------------------------------------------------------
# 3. Build RSEM Reference
# ------------------------------------------------------------------------------
echo ">>> Step 2: Building RSEM Reference..."
${RSEM_BIN_DIR}rsem-prepare-reference \
--bowtie2 \
--gtf "$GTF_RSEM_CLEAN" \
"$GENOME" \
"$RSEM_REF_PREFIX"

# ------------------------------------------------------------------------------
# 4. Simulate Reads at Multiple Depths
# ------------------------------------------------------------------------------
echo ">>> Step 3: Simulating Reads using RSEM..."

# Calculate total transcript length from the simulated isoform file
total_length=$(awk 'NR>1 {sum += $3} END {print int(sum)}' "$ISO")

for depth in "${DEPTHS[@]}"; do
total_reads=$(echo "$depth * $total_length / $READ_LEN" | bc)
echo "  -> Simulating Depth ${depth}x: ${total_reads} reads"

${RSEM_BIN_DIR}rsem-simulate-reads \
"$RSEM_REF_PREFIX" \
"$MODEL" \
"$ISO" \
"$THETA0" \
"$total_reads" \
"${FASTQ_DIR}/official_simulated_${depth}x_noise0.1" \
--seed "${SIM_SEED:-18}" \
> "${LOG_DIR}/simulate_${depth}x.log" 2>&1
done

# ------------------------------------------------------------------------------
# 5. Sanity Check: Quantify Simulated Reads
# ------------------------------------------------------------------------------
echo ">>> Step 4: Running Sanity Check (Quantification)..."

for depth in "${DEPTHS[@]}"; do
echo "  -> Quantifying ${depth}x"

${RSEM_BIN_DIR}rsem-calculate-expression \
--paired-end \
-p 20 \
--bowtie2 \
"${FASTQ_DIR}/official_simulated_${depth}x_noise0.1_1.fq" \
"${FASTQ_DIR}/official_simulated_${depth}x_noise0.1_2.fq" \
"$RSEM_REF_PREFIX" \
"${CHECK_DIR}/official_simulated_${depth}x" \
> "${LOG_DIR}/quant_check_${depth}x.log" 2>&1
done

# ------------------------------------------------------------------------------
# 6. Check Logs and Print Summaries
# ------------------------------------------------------------------------------
echo ">>> Step 5: Summarizing Sanity Check Results..."

# Check for errors in quantification logs
if grep -i "error\|failed" "${LOG_DIR}"/quant_check_*x.log; then
echo "WARNING: Errors found in quantification logs."
fi

for depth in "${DEPTHS[@]}"; do
echo "===== ${depth}x ====="
awk 'NR>1 && $6 > 0 {n++} END {print "expressed_isoforms_TPM_gt_0:", n+0}' \
"${CHECK_DIR}/official_simulated_${depth}x.isoforms.results"
awk 'NR>1 {sum += $6} END {print "TPM_sum:", sum}' \
"${CHECK_DIR}/official_simulated_${depth}x.isoforms.results"
awk 'NR>1 {sum += $5} END {print "expected_count_sum:", sum}' \
"${CHECK_DIR}/official_simulated_${depth}x.isoforms.results"
done

# ------------------------------------------------------------------------------
# 7. Compress FASTQ Files
# ------------------------------------------------------------------------------
echo ">>> Step 6: Compressing FASTQ Files..."

for depth in "${DEPTHS[@]}"; do
sample="${FASTQ_DIR}/official_simulated_${depth}x_noise0.1"
echo "  -> Compressing ${depth}x"
pigz -c -p 8 "${sample}_1.fq" > "${sample}.R1.fastq.gz"
pigz -c -p 8 "${sample}_2.fq" > "${sample}.R2.fastq.gz"

# Optional: Remove the uncompressed .fq files to save disk space
# rm "${sample}_1.fq" "${sample}_2.fq"
done

echo ">>> All simulations and checks completed successfully!"

# ==============================================================================
# Next Step: Execute Downstream WDL Workflow
# Note: Adjust paths to cromwell and your WDL scripts as needed.
# ==============================================================================
# cd ./pipeline/wdl_RNA
# java -jar /path/to/software/cromwell-71.jar run workflow1.wdl -i official_simulate.json