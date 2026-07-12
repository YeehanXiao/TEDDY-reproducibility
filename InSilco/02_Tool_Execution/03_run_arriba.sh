# ==============================================================================
# Script: 03_run_arriba.sh
# Purpose: 
#   1. Convert TE BED to a "fake" GTF and merge with the official reference.
#   2. Build STAR index with the merged GTF.
#   3. Run STAR alignment to generate chimeric BAMs for multiple depths.
#   4. Run Arriba to call fusion/TE-chimeric events.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Paths and parameters (Desensitized)
# ------------------------------------------------------------------------------
DATA_DIR="./data"
RESULTS_DIR="./results"
FASTQ_DIR="${RESULTS_DIR}/fastq"
ARRIBA_DIR="${RESULTS_DIR}/arriba_merge"
REF_DIR="${ARRIBA_DIR}/ref"
STAR_INDEX="${ARRIBA_DIR}/STAR_index_official_mergeTE"

mkdir -p "${REF_DIR}" "${STAR_INDEX}"

# Input files
GENOME_FA="${DATA_DIR}/mm10_no_alt_analysis_set_ENCODE.fasta"
REF_GTF="${RESULTS_DIR}/official_simulated_reference_90pct.gtf"
TE_BED="${REF_DIR}/mm10_TE.arriba.bed6" # From 01_Data_Simulation/01b_prepare_arriba_ref.R

# Output files
TE_GTF="${REF_DIR}/official_simulated_TE_fake.gtf"
MERGE_GTF="${REF_DIR}/official_simulated_mergeTE.gtf"

# Tools
ARRIBA_BIN="arriba" # 假设已加入系统环境变量，如果没有，请写绝对路径，例如 /path/to/arriba
STAR_BIN="STAR"

# Variables
DEPTHS=("5x" "10x" "25x" "50x" "100x") # 统一了所有测序深度

# ------------------------------------------------------------------------------
# 1. Generate "fake" TE GTF
# ------------------------------------------------------------------------------
echo ">>> Step 1: Converting TE BED to fake GTF..."

awk 'BEGIN{
  FS=OFS="\t"
}
function clean(x, y){
  y=x
  gsub(/[^A-Za-z0-9_]/, "_", y)
  return y
}
NF >= 6 {
  chr=$1
  start=$2 + 1
  end=$3
  te_name=$4
  score=$5
  strand=$6

  if (te_name=="" || te_name==".") te_name="TE"
  if (strand!="+" && strand!="-") strand="."

  te_class="unknown"
  te_family="unknown"

  safe_name=clean(te_name)
  safe_chr=clean(chr)
  id=safe_name "_" safe_chr "_" start "_" end

  gene_id="TE_" id
  tx_id=gene_id ".1"

  gene_attr = "gene_id \"" gene_id "\"; gene_name \"" te_name "\"; gene_type \"transposable_element\"; gene_biotype \"transposable_element\"; te_class \"" te_class "\"; te_family \"" te_family "\";"
  tx_attr = "gene_id \"" gene_id "\"; gene_name \"" te_name "\"; transcript_id \"" tx_id "\"; transcript_name \"" tx_id "\"; gene_type \"transposable_element\"; gene_biotype \"transposable_element\"; transcript_type \"transposable_element\"; transcript_biotype \"transposable_element\"; te_class \"" te_class "\"; te_family \"" te_family "\";"
  exon_attr = tx_attr " exon_number \"1\";"

  print chr, "fakeTE", "gene",       start, end, ".", strand, ".", gene_attr
  print chr, "fakeTE", "transcript", start, end, ".", strand, ".", tx_attr
  print chr, "fakeTE", "exon",       start, end, ".", strand, ".", exon_attr
}' "${TE_BED}" > "${TE_GTF}"

awk -F'\t' 'NF!=9{print "bad line", NR, "NF="NF; bad=1} END{exit bad}' "${TE_GTF}"
echo "Fake TE GTF generated."

# ------------------------------------------------------------------------------
# 2. Merge Reference GTF and TE GTF
# ------------------------------------------------------------------------------
echo ">>> Step 2: Merging Reference and TE GTF..."

{
  grep '^#' "${REF_GTF}" || true
  grep -v '^#' "${REF_GTF}"
  grep -v '^#' "${TE_GTF}"
} | awk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} NF==9 {print}' \
| sort -k1,1 -k4,4n -k5,5n > "${MERGE_GTF}"

awk -F'\t' '!/^#/ && NF!=9{print "bad line", NR, "NF="NF; bad=1} END{exit bad}' "${MERGE_GTF}"
echo "Merged GTF generated."

# ------------------------------------------------------------------------------
# 3. Generate STAR Index
# ------------------------------------------------------------------------------
echo ">>> Step 3: Generating STAR Index..."

# Extract READLEN from the first 100x fastq file
R1_TEST="${FASTQ_DIR}/official_simulated_100x_noise0.1.R1.fastq.gz"
if [ -f "$R1_TEST" ]; then
READLEN=$(zcat "${R1_TEST}" | awk 'NR==2{print length($0); exit}')
else
  READLEN=150 # 默认值，防报错
echo "Warning: R1_TEST not found. Using default READLEN=150"
fi
SJDB_OVERHANG=$((READLEN - 1))

echo "Using READLEN=${READLEN}, SJDB_OVERHANG=${SJDB_OVERHANG}"

ulimit -n 65536 || true

${STAR_BIN} --runMode genomeGenerate \
--runThreadN 20 \
--genomeDir "${STAR_INDEX}" \
--genomeFastaFiles "${GENOME_FA}" \
--sjdbGTFfile "${MERGE_GTF}" \
--sjdbOverhang "${SJDB_OVERHANG}"

# ------------------------------------------------------------------------------
# 4. Run STAR and Arriba
# ------------------------------------------------------------------------------
echo ">>> Step 4: Running STAR and Arriba for all depths..."

for depth in "${DEPTHS[@]}"; do
echo "========== Processing Depth: ${depth} =========="
R1="${FASTQ_DIR}/official_simulated_${depth}_noise0.1.R1.fastq.gz"
R2="${FASTQ_DIR}/official_simulated_${depth}_noise0.1.R2.fastq.gz"

PREFIX="${ARRIBA_DIR}/merge_${depth}_"

# Clean old files
rm -f "${PREFIX}Aligned.out.bam" "${PREFIX}Chimeric.out.junction" "${PREFIX}Log.final.out" "${PREFIX}Log.out" "${PREFIX}Log.progress.out" "${PREFIX}SJ.out.tab"
rm -rf "${PREFIX}_STARtmp" "${PREFIX}_STARpass1" "${PREFIX}_STARgenome"

# A. STAR Alignment
# (注意：如果你需要记录时间供 99_timesummary.R 读取，可以在 STAR 前面加上 /usr/bin/time -v -o "${ARRIBA_DIR}/${depth}_arriba.time.log")
${STAR_BIN} --runThreadN 20 \
--genomeDir "${STAR_INDEX}" \
--readFilesIn "${R1}" "${R2}" \
--readFilesCommand zcat \
--sjdbGTFfile "${MERGE_GTF}" \
--twopassMode Basic \
--outSAMtype BAM Unsorted \
--outSAMunmapped Within \
--outBAMcompression 0 \
--outFilterMultimapNmax 50 \
--peOverlapNbasesMin 10 \
--alignSplicedMateMapLminOverLmate 0.5 \
--alignSJstitchMismatchNmax 5 -1 5 5 \
--chimSegmentMin 10 \
--chimJunctionOverhangMin 10 \
--chimOutType WithinBAM HardClip Junctions \
--chimOutJunctionFormat 1 \
--chimScoreDropMax 30 \
--chimScoreJunctionNonGTAG 0 \
--chimScoreSeparation 1 \
--chimSegmentReadGapMax 3 \
--chimMultimapNmax 50 \
--outFileNamePrefix "${PREFIX}"

# B. Arriba Fusion Calling
"${ARRIBA_BIN}" \
-x "${PREFIX}Aligned.out.bam" \
-g "${MERGE_GTF}" \
-a "${GENOME_FA}" \
-f blacklist \
-o "${PREFIX}arriba_fusions.tsv" \
-O "${PREFIX}arriba_fusions.discarded.tsv"

done

echo ">>> All STAR and Arriba runs completed successfully!"