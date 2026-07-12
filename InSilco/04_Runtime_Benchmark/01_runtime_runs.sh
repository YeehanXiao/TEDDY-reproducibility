###############################################################################
# Runtime benchmark commands
#
# This script records the commands used to benchmark the runtime of TEDDY and
# comparison tools on the 100x simulated RNA-seq dataset.
#
# Notes:
# 1. Paths are anonymized and should be replaced by users before running.
# 2. Runtime was measured on the 100x simulated dataset.
# 3. Each tool was run three times with isolated output directories.
# 4. Wall-clock time and maximum resident memory were recorded by /usr/bin/time -v.
###############################################################################

PROJECT_DIR="/path/to/simulation_benchmark"
GENOME_FA="/path/to/mm10.fa"
REFERENCE_GTF="${PROJECT_DIR}/official_simulated_reference_90pct.gtf"
TRUTH_GTF="${PROJECT_DIR}/official_simulated_truth_1000_transcripts.gtf"
TE_ANNOTATION="/path/to/mm10_TE.annotation"
THREADS=8
N_REP=3

FASTQ_R1="${PROJECT_DIR}/fastq/official_simulated_100x_noise0.1.R1.fastq.gz"
FASTQ_R2="${PROJECT_DIR}/fastq/official_simulated_100x_noise0.1.R2.fastq.gz"
BAM_100X="${PROJECT_DIR}/sortbam/official_simulated_100x_noise0.1.bam"

RUNTIME_DIR="${PROJECT_DIR}/runtime_benchmark"
LOG_DIR="${RUNTIME_DIR}/logs"

mkdir -p "${LOG_DIR}"

run_with_time() {
  local tool="$1"
  local rep="$2"
  shift 2

  mkdir -p "${LOG_DIR}"

  /usr/bin/time -v "$@" \
    > "${LOG_DIR}/${tool}_100x_rep${rep}.stdout.log" \
    2> "${LOG_DIR}/${tool}_100x_rep${rep}.time.log"
}

###############################################################################
# TEDDY
###############################################################################

for rep in $(seq 1 "${N_REP}"); do
  OUTDIR="${RUNTIME_DIR}/TEDDY/100x_rep${rep}"
  mkdir -p "${OUTDIR}"

  run_with_time TEDDY "${rep}" \
    Rscript run_TEDDY_100x_runtime.R \
      --bam "${BAM_100X}" \
      --reference_gtf "${REFERENCE_GTF}" \
      --te_annotation "${TE_ANNOTATION}" \
      --outdir "${OUTDIR}" \
      --threads "${THREADS}"
done

###############################################################################
# FREDY
###############################################################################

for rep in $(seq 1 "${N_REP}"); do
  OUTDIR="${RUNTIME_DIR}/FREDY/100x_rep${rep}"
  mkdir -p "${OUTDIR}"

  printf "%s\n%s\n" "${FASTQ_R1}" "${FASTQ_R2}" > "${OUTDIR}/fastq_files.txt"

  run_with_time FREDY "${rep}" \
    bash -c "
      fredy_udocker star \
        -o '${OUTDIR}' \
        -i '/path/to/STAR_index_built_from_reference_gtf' \
        -f '${OUTDIR}/fastq_files.txt' \
        -S -p \
        -t '${THREADS}'

      fredy_udocker string \
        -o '${OUTDIR}' \
        -a '${REFERENCE_GTF}' \
        -t '${THREADS}'

      fredy_udocker chimeric \
        -o '${OUTDIR}' \
        -a '/path/to/FREDY_CDS_proxy_reference.gtf' \
        -g '${GENOME_FA}' \
        -e '/path/to/TE_annotation_for_FREDY.bed'
    "
done

###############################################################################
# TEProf2
###############################################################################

for rep in $(seq 1 "${N_REP}"); do
  OUTDIR="${RUNTIME_DIR}/TEProf2/100x_rep${rep}"
  mkdir -p "${OUTDIR}"

  run_with_time TEProf2 "${rep}" \
    bash -c "
      stringtie '${BAM_100X}' \
        -o '${OUTDIR}/official_simulated_100x.stringtie.gtf' \
        -p '${THREADS}' \
        -G '${REFERENCE_GTF}' \
        -m 100 \
        -c 1 \
        -f 0.01 \
        -j 1

      python /path/to/TEProf2_annotation_script.py \
        '${OUTDIR}/official_simulated_100x.stringtie.gtf' \
        '/path/to/TEProf2_argument_file.txt'
    "
done

###############################################################################
# ChimeraTE
###############################################################################

for rep in $(seq 1 "${N_REP}"); do
  OUTDIR="${RUNTIME_DIR}/ChimeraTE/100x_rep${rep}"
  mkdir -p "${OUTDIR}"

  INPUT_TSV="${OUTDIR}/input_100x.tsv"
  printf "%s\t%s\trep1_100x\n" "${FASTQ_R1}" "${FASTQ_R2}" > "${INPUT_TSV}"

  run_with_time ChimeraTE "${rep}" \
    bash -c "
      python /path/to/ChimeraTE/chimTE_mode1.py \
        --genome '${GENOME_FA}' \
        --input '${INPUT_TSV}' \
        --project 'runtime_100x_rep${rep}' \
        --te '/path/to/TE_annotation_for_ChimeraTE.gtf' \
        --gene '/path/to/gene_annotation_for_ChimeraTE.gtf' \
        --strand fwd-stranded \
        --threads '${THREADS}'

      python /path/to/ChimeraTE/chimTE_mode1.py \
        --genome '${GENOME_FA}' \
        --input '${INPUT_TSV}' \
        --project 'runtime_100x_rep${rep}_rf' \
        --te '/path/to/TE_annotation_for_ChimeraTE.gtf' \
        --gene '/path/to/gene_annotation_for_ChimeraTE.gtf' \
        --strand rf-stranded \
        --threads '${THREADS}'
    "
done

###############################################################################
# LIONS
###############################################################################

for rep in $(seq 1 "${N_REP}"); do
  OUTDIR="${RUNTIME_DIR}/LIONS/100x_rep${rep}"
  mkdir -p "${OUTDIR}"

  run_with_time LIONS "${rep}" \
    bash -c "
      zcat '${FASTQ_R1}' > '${OUTDIR}/R1.fastq'
      paste - - - - < '${OUTDIR}/R1.fastq' | cut -f1,2 | sed 's/^@/>/' | tr '\t' '\n' > '${OUTDIR}/R1.fasta'
      grep '^@' '${OUTDIR}/R1.fastq' | sed 's/^@//' > '${OUTDIR}/R1.ids'

      python /path/to/LIONS/ChimericReadTool/chimericReadSearch.py \
        '${OUTDIR}/R1.fastq' \
        '${OUTDIR}/R1.fasta' \
        '${OUTDIR}/R1.ids' \
        '${OUTDIR}/R1.chimeras'

      zcat '${FASTQ_R2}' > '${OUTDIR}/R2.fastq'
      paste - - - - < '${OUTDIR}/R2.fastq' | cut -f1,2 | sed 's/^@/>/' | tr '\t' '\n' > '${OUTDIR}/R2.fasta'
      grep '^@' '${OUTDIR}/R2.fastq' | sed 's/^@//' > '${OUTDIR}/R2.ids'

      python /path/to/LIONS/ChimericReadTool/chimericReadSearch.py \
        '${OUTDIR}/R2.fastq' \
        '${OUTDIR}/R2.fasta' \
        '${OUTDIR}/R2.ids' \
        '${OUTDIR}/R2.chimeras'
    "
done

###############################################################################
# Arriba
###############################################################################

for rep in $(seq 1 "${N_REP}"); do
  OUTDIR="${RUNTIME_DIR}/ARRIBA/100x_rep${rep}"
  mkdir -p "${OUTDIR}"

  run_with_time ARRIBA "${rep}" \
    /path/to/arriba \
      -x "/path/to/STAR_aligned_100x_Aligned.out.bam" \
      -g "/path/to/official_simulated_mergeTE.gtf" \
      -a "${GENOME_FA}" \
      -f blacklist \
      -o "${OUTDIR}/arriba_fusions.tsv" \
      -O "${OUTDIR}/arriba_fusions.discarded.tsv"
done