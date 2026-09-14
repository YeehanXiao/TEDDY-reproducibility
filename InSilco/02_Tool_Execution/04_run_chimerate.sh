#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Script: 04_run_chimerate.sh
# Purpose: 
#   - ChimeraTE mode1 benchmark for official simulation by depth
#   - Uses official stranded-truth simulation FASTQ
#   - Uses ChimeraTE-native mode1 workflow
#   - Runs both supported strand modes because simulated reads are unstranded
#   - Includes pandas >= 2.0 compatibility patch
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Define Paths (Dynamic Absolute Paths to avoid 'cd' issues)
# ------------------------------------------------------------------------------
WORK_DIR="$(pwd)"
DATA_DIR="${WORK_DIR}/data"
RESULTS_DIR="${WORK_DIR}/results"

# IMPORTANT: Path to the cloned ChimeraTE repository
CHIMERATE_DIR="${CHIMERATE_DIR:-${WORK_DIR}/ChimeraTE}"
CHIMERATE_SHA256="${WORK_DIR}/InSilco/02_Tool_Execution/adapters/chimerate/SHA256SUMS"

OUTBASE="${RESULTS_DIR}/ChimeraTE_official_by_depth"

GENOME="${DATA_DIR}/mm10_no_alt_analysis_set_ENCODE.fasta"
TE_GTF="${DATA_DIR}/mm10_TE_annotations.gtf"
REF_GTF="${RESULTS_DIR}/official_simulated_reference_90pct.gtf"
GENE_GTF="${OUTBASE}/ref/official_simulated_reference_90pct.ChimeraTE_gene_exon.gtf"

THREADS=8
DEPTHS=("5x" "10x" "25x" "50x" "100x")
STRANDS=("fwd-stranded" "rf-stranded")

test -s "${CHIMERATE_SHA256}" || { echo "Error: Missing ChimeraTE checksum manifest"; exit 1; }
(cd "${CHIMERATE_DIR}" && sha256sum --check --status "${CHIMERATE_SHA256}") || {
  echo "Error: CHIMERATE_DIR does not contain the validated benchmark fork" >&2
  exit 1
}

mkdir -p "${OUTBASE}/ref"
mkdir -p "${OUTBASE}/logs"
mkdir -p "${OUTBASE}/projects"

rm -f "${OUTBASE}/logs/ChimeraTE_input_tsv_check.txt"
rm -f "${OUTBASE}/logs/ChimeraTE_run_environment.txt"
rm -f "${OUTBASE}/logs/ChimeraTE_strandness_inference.txt"
rm -f "${OUTBASE}/logs/ChimeraTE_remaining_append_calls.txt"

# Enter the ChimeraTE directory as required by its scripts
cd "${CHIMERATE_DIR}"
mkdir -p projects

# ------------------------------------------------------------------------------
# 1. Record environment
# ------------------------------------------------------------------------------
{
  echo "Date: $(date)"
  echo "PWD: $(pwd)"
  echo "Python: $(python3 --version 2>&1)"
  python3 -c 'import numpy, pandas; print("pandas:", pandas.__version__); print("numpy:", numpy.__version__)'
  echo "WORK_DIR=${WORK_DIR}"
  echo "OUTBASE=${OUTBASE}"
  echo "GENOME=${GENOME}"
  echo "TE_GTF=${TE_GTF}"
  echo "REF_GTF=${REF_GTF}"
  echo "GENE_GTF=${GENE_GTF}"
  echo "THREADS=${THREADS}"
  echo "DEPTHS=${DEPTHS[*]}"
  echo "STRANDS=${STRANDS[*]}"
} > "${OUTBASE}/logs/ChimeraTE_run_environment.txt"

# ------------------------------------------------------------------------------
# 2. Prepare ChimeraTE-compatible gene annotation
# ------------------------------------------------------------------------------
test -s "${REF_GTF}" || { echo "Error: REF_GTF not found!"; exit 1; }

awk 'BEGIN{OFS="\t"}
     /^#/ {print; next}
     $3=="gene" || $3=="exon" {print}' \
"${REF_GTF}" \
> "${GENE_GTF}"

{
  echo "ChimeraTE gene annotation generated from:"
  echo "${REF_GTF}"
  echo
  awk '$0 !~ /^#/ && $3=="gene"{n++} END{print "genes:", n+0}' "${GENE_GTF}"
  awk '$0 !~ /^#/ && $3=="exon"{n++} END{print "exons:", n+0}' "${GENE_GTF}"
  echo
  echo "First non-comment lines:"
  awk '$0 !~ /^#/ {print; n++} n==5{exit}' "${GENE_GTF}"
} > "${OUTBASE}/logs/ChimeraTE_gene_gtf_check.txt"

# ------------------------------------------------------------------------------
# 3. Prepare ChimeraTE input TSV
# ------------------------------------------------------------------------------
for depth in "${DEPTHS[@]}"; do
project_input_dir="${OUTBASE}/projects/${depth}"
mkdir -p "${project_input_dir}"

fq1="${RESULTS_DIR}/fastq/official_simulated_${depth}_noise0.1.R1.fastq.gz"
fq2="${RESULTS_DIR}/fastq/official_simulated_${depth}_noise0.1.R2.fastq.gz"

test -s "${fq1}" || { echo "Error: fq1 missing for ${depth}"; exit 1; }
test -s "${fq2}" || { echo "Error: fq2 missing for ${depth}"; exit 1; }

input_tsv="${project_input_dir}/input_mode1_${depth}.tsv"

cat > "${input_tsv}" <<EOF
${fq1}	${fq2}	rep1_${depth}
EOF

{
  echo "Depth: ${depth}"
  echo "Input: ${input_tsv}"
  echo "Rows:"
  wc -l "${input_tsv}"
  echo "Columns per row:"
  awk -F '\t' '{print NF}' "${input_tsv}"
  echo "Preview:"
  cat "${input_tsv}"
  echo
} >> "${OUTBASE}/logs/ChimeraTE_input_tsv_check.txt"

awk -F '\t' 'NF != 3 {exit 1}' "${input_tsv}"
done

# ------------------------------------------------------------------------------
# 4. Record strandness inference
# ------------------------------------------------------------------------------
cat > "${OUTBASE}/logs/ChimeraTE_strandness_inference.txt" <<'EOF'
RSeQC infer_experiment.py result for official_simulated_100x_noise0.1.bam:
  
  This is PairEnd Data
Fraction of reads failed to determine: 0.0000
Fraction of reads explained by "1++,1--,2+-,2-+": 0.5013
Fraction of reads explained by "1+-,1-+,2++,2--": 0.4986

Interpretation:
  The simulated paired-end reads are effectively unstranded.
Because ChimeraTE mode1 only supports fwd-stranded and rf-stranded, both modes
are run and their predictions should be merged/deduplicated for evaluation.
EOF

# ------------------------------------------------------------------------------
# 5. Patch ChimeraTE for pandas >= 2.0
# ------------------------------------------------------------------------------
PATCH_FILE="scripts/mode1_te_exonized.py"
test -s "${PATCH_FILE}" || { echo "Error: Patch file not found! Are you in the ChimeraTE directory?"; exit 1; }

python3 - <<'PY'
from pathlib import Path
import pandas as pd

p = Path("scripts/mode1_te_exonized.py")
s = p.read_text()

replacements = {
  "merging_gene_reads = merging_gene_reads.append(reads_gene_col, ignore_index=True).drop_duplicates()":
    "merging_gene_reads = pd.concat([merging_gene_reads, reads_gene_col], ignore_index=True).drop_duplicates()",
  
  "merging_TE_reads = merging_TE_reads.append(reads_TE_col, ignore_index=True).drop_duplicates()":
    "merging_TE_reads = pd.concat([merging_TE_reads, reads_TE_col], ignore_index=True).drop_duplicates()",
}

for old, new in replacements.items():
  if old in s:
    s = s.replace(old, new)

for old, new in replacements.items():
  if old in s or new not in s:
    raise SystemExit("Unsupported ChimeraTE pandas API in " + str(p))

p.write_text(s)
PY

grep -RIn "\.append" scripts chimTE_mode1.py chimTE_mode2.py \
> "${OUTBASE}/logs/ChimeraTE_remaining_append_calls.txt" || true

# ------------------------------------------------------------------------------
# 6. Check required files
# ------------------------------------------------------------------------------
test -s "${GENOME}"
test -s "${TE_GTF}"
test -s "${GENE_GTF}"

for depth in "${DEPTHS[@]}"; do
test -s "${OUTBASE}/projects/${depth}/input_mode1_${depth}.tsv"
done

# ------------------------------------------------------------------------------
# 7. Run ChimeraTE mode1 across depths and supported strand settings
# ------------------------------------------------------------------------------
for strand in "${STRANDS[@]}"; do
for depth in "${DEPTHS[@]}"; do
echo "[$(date)] Running ChimeraTE mode1 depth=${depth}, strand=${strand}"

input_tsv="${OUTBASE}/projects/${depth}/input_mode1_${depth}.tsv"
project="official_${depth}_${strand}"
result_dir="${OUTBASE}/${depth}_${strand}"
log_file="${OUTBASE}/logs/ChimeraTE_mode1_${depth}_${strand}.log"

rm -rf "projects/${project}"
rm -rf "${result_dir}"

python3 chimTE_mode1.py \
--genome "${GENOME}" \
--input "${input_tsv}" \
--project "${project}" \
--te "${TE_GTF}" \
--gene "${GENE_GTF}" \
--strand "${strand}" \
--threads "${THREADS}" \
--overlap 0.10 \
> "${log_file}" 2>&1

cp -r "projects/${project}" "${result_dir}"

find "${result_dir}" -type f | sort \
> "${OUTBASE}/logs/ChimeraTE_files_${depth}_${strand}.txt"

echo "[$(date)] DONE ChimeraTE mode1 depth=${depth}, strand=${strand}"
done
done

# ------------------------------------------------------------------------------
# 8. Final log summary
# ------------------------------------------------------------------------------
grep -RniE "traceback|error|failed|no such|not found|exception|permission" \
"${OUTBASE}/logs"/ChimeraTE_mode1_*_*.log \
> "${OUTBASE}/logs/ChimeraTE_error_scan.txt" || true

grep -RhiE "ChimeraTE has finished|There are no|TE-initiated|TE-terminated|TE-exonized|Running analysis with" \
"${OUTBASE}/logs"/ChimeraTE_mode1_*_*.log \
> "${OUTBASE}/logs/ChimeraTE_run_status_summary.txt" || true

echo ">>> All ChimeraTE runs and log generation completed successfully!"
