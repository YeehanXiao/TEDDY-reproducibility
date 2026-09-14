#!/usr/bin/env bash
set -euo pipefail

WORK_DIR="$(pwd)"
CHIMERATE_DIR="${CHIMERATE_DIR:-${WORK_DIR}/ChimeraTE}"
RESULTS_DIR="${WORK_DIR}/results"
OFFICIAL="${RESULTS_DIR}/ChimeraTE_official_by_depth"
OUT="${RESULTS_DIR}/ChimeraTE_sensitivity_25x"
GENOME="${WORK_DIR}/data/mm10_no_alt_analysis_set_ENCODE.fasta"
TE_GTF="${WORK_DIR}/data/mm10_TE_annotations.gtf"
GENE_GTF="${OFFICIAL}/ref/official_simulated_reference_90pct.ChimeraTE_gene_exon.gtf"
INPUT="${OFFICIAL}/projects/25x/input_mode1_25x.tsv"
THREADS="${THREADS:-8}"
CHIMERATE_SHA256="${WORK_DIR}/InSilco/02_Tool_Execution/adapters/chimerate/SHA256SUMS"

mkdir -p "${OUT}/logs" "${OUT}/results"
cd "${CHIMERATE_DIR}"
mkdir -p projects

for file in "${GENOME}" "${TE_GTF}" "${GENE_GTF}" "${INPUT}"; do
    test -s "${file}" || { echo "Missing required file: ${file}" >&2; exit 1; }
done

(sha256sum --check --status "${CHIMERATE_SHA256}") \
    || { echo "CHIMERATE_DIR does not contain the validated sensitivity fork" >&2; exit 1; }

run_one() {
    local tag="$1" strand="$2" mapq="$3" overlap="$4"
    local project="sensitivity_25x_${tag}_${strand}"
    local result="${OUT}/results/${tag}_${strand}"
    local index="${OFFICIAL}/25x_${strand}/index"

    test -s "${INPUT}" && test -s "${index}/SAindex"
    rm -rf "projects/${project}" "${result}"
    python3 chimTE_mode1.py \
        --genome "${GENOME}" --input "${INPUT}" --project "${project}" \
        --te "${TE_GTF}" --gene "${GENE_GTF}" --strand "${strand}" \
        --threads "${THREADS}" --index "${index}" --mapq "${mapq}" \
        --star_multimap_nmax 10 --overlap "${overlap}" \
        > "${OUT}/logs/${tag}_${strand}.log" 2>&1
    mv "projects/${project}" "${result}"
}

for strand in fwd-stranded rf-stranded; do
    run_one default "${strand}" 255 0.50
    run_one multimap "${strand}" 0 0.50
    run_one overlap10 "${strand}" 255 0.10
    run_one both_relaxed "${strand}" 0 0.10
done

summary="${OUT}/ChimeraTE_sensitivity_25x_counts.tsv"
printf 'configuration\tMAPQ_threshold\tread_TE_overlap_threshold\tstrand_mode\tTE_exonized_rows\tTE_initiated_rows\tTE_terminated_rows\n' > "${summary}"
for config in 'default 255 0.50' 'multimap 0 0.50' 'overlap10 255 0.10' 'both_relaxed 0 0.10'; do
    read -r tag mapq overlap <<< "${config}"
    for strand in fwd-stranded rf-stranded; do
        result="${OUT}/results/${tag}_${strand}/rep1_25x"
        exonized="${result}/TE-exonized-rep1_25x.tsv"
        initiated="${result}/TE-initiated-rep1_25x.tsv"
        terminated="${result}/TE-terminated-rep1_25x.tsv"
        test -f "${exonized}" && test -f "${initiated}" && test -f "${terminated}"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "${tag}" "${mapq}" "${overlap}" "${strand}" \
            "$(wc -l < "${exonized}")" "$(wc -l < "${initiated}")" "$(wc -l < "${terminated}")" \
            >> "${summary}"
    done
done
