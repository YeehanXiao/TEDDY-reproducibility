#!/usr/bin/env bash
set -euo pipefail

# Usage:
# CONDA=/path/to/conda \
# PY2_ENV=teprof2_py2 \
# RSCRIPT=/path/to/Rscript \
# R_LIBS_EXTRA=/path/to/R/library \
# SAMTOOLS_DIR=/path/to/samtools/bin \
# BEDTOOLS_DIR=/path/to/bedtools/bin \
# bash filter_teprof2_reads.sh \
#   TEPROF2_BIN AGGREGATE_DIR BAM_DIR [JOBS]

TEPROF2_BIN=${1:?TEPROF2_BIN is required}
OUT=${2:?AGGREGATE_DIR is required}
BAM_DIR=${3:?BAM_DIR is required}
JOBS=${4:-8}

CONDA=${CONDA:-conda}
PY2_ENV=${PY2_ENV:-teprof2_py2}
RSCRIPT=${RSCRIPT:-Rscript}
R_LIBS_EXTRA=${R_LIBS_EXTRA:-}
SAMTOOLS_DIR=${SAMTOOLS_DIR:-}
BEDTOOLS_DIR=${BEDTOOLS_DIR:-}

[[ -n "$OUT" && "$OUT" != "/" && -d "$OUT" ]] || {
    echo "ERROR: unsafe aggregate directory: $OUT" >&2
    exit 1
}

export PATH="${SAMTOOLS_DIR:+$SAMTOOLS_DIR:}${BEDTOOLS_DIR:+$BEDTOOLS_DIR:}$PATH"

run_py2() {
    env -u PYTHONPATH -u PYTHONHOME \
        "$CONDA" run --no-capture-output \
        -n "$PY2_ENV" python "$@"
}

run_r() {
    if [[ -n "$R_LIBS_EXTRA" ]]; then
        R_LIBS="$R_LIBS_EXTRA${R_LIBS:+:$R_LIBS}" \
            "$RSCRIPT" "$@"
    else
        "$RSCRIPT" "$@"
    fi
}

for file in \
    "$TEPROF2_BIN/commandsmax_speed.py" \
    "$TEPROF2_BIN/rmsk_annotate_bedpe_speed.py" \
    "$TEPROF2_BIN/filterReadCandidates.R" \
    "$OUT/filter_combined_candidates.tsv" \
    "$OUT/Step4.RData"
do
    test -s "$file" || {
        echo "ERROR: missing file: $file" >&2
        exit 1
    }
done

command -v samtools >/dev/null
command -v bedtools >/dev/null
command -v xargs >/dev/null

BAM_DIR="${BAM_DIR%/}/"

# Field 45 is the sample label used by commandsmax_speed.py.
while IFS= read -r sample; do
    bam="${BAM_DIR}${sample}.bam"

    test -s "$bam" || {
        echo "ERROR: missing BAM: $bam" >&2
        exit 1
    }

    samtools quickcheck "$bam"

    if [[ ! -s "${bam}.bai" && ! -s "${bam%.bam}.bai" ]]; then
        samtools index -@ 4 "$bam"
    fi
done < <(
    cut -f45 "$OUT/filter_combined_candidates.tsv" |
    sort -u
)

cd "$OUT"

# Use local copies so the upstream TEProf2 source remains unchanged.
rm -rf read_filter_bin filterreadstats
mkdir -p read_filter_bin filterreadstats

cp "$TEPROF2_BIN/commandsmax_speed.py" \
   read_filter_bin/

cp "$TEPROF2_BIN/rmsk_annotate_bedpe_speed.py" \
   read_filter_bin/

# The helper is invoked directly by generated shell commands.
# Use the Python supplied by the activated legacy conda environment.
sed -i '1c#!/usr/bin/env python' \
    read_filter_bin/rmsk_annotate_bedpe_speed.py

chmod +x read_filter_bin/rmsk_annotate_bedpe_speed.py

rm -f filterreadcommands.txt

run_py2 \
    read_filter_bin/commandsmax_speed.py \
    filter_combined_candidates.tsv \
    "$BAM_DIR"

test -s filterreadcommands.txt

n_commands=$(wc -l < filterreadcommands.txt)
echo "Read-stat commands: $n_commands"

# Run all generated commands inside the verified Python 2 environment.
env -u PYTHONPATH -u PYTHONHOME \
    PATH="$PATH" \
    "$CONDA" run --no-capture-output \
    -n "$PY2_ENV" \
    bash -c '
        set -euo pipefail
        tr "\n" "\0" < filterreadcommands.txt |
            xargs -0 -r -n 1 -P "$1" bash -c
    ' _ "$JOBS"

n_stats=$(
    find filterreadstats \
        -maxdepth 1 \
        -type f \
        -name '*.stats' |
    wc -l
)

echo "Read-stat files: $n_stats"

[[ "$n_stats" -eq "$n_commands" ]] || {
    echo "ERROR: expected $n_commands stats files, found $n_stats" >&2
    exit 1
}

find filterreadstats \
    -maxdepth 1 \
    -type f \
    -name '*.stats' \
    -print0 |
sort -z |
xargs -0 -r grep -H 'e' \
    > resultgrep_filterreadstatsdone.txt

sed 's/:/\t/' \
    resultgrep_filterreadstatsdone.txt \
    > filter_read_stats.txt

[[ "$(wc -l < filter_read_stats.txt)" -eq "$n_commands" ]] || {
    echo "ERROR: incomplete filter_read_stats.txt" >&2
    exit 1
}

# Explicit namespace avoids ArgumentParser name collisions.
cp "$TEPROF2_BIN/filterReadCandidates.R" \
   filterReadCandidates.compat.R

sed -i \
    's/parser <- ArgumentParser$new()/parser <- Xmisc::ArgumentParser$new()/' \
    filterReadCandidates.compat.R

rm -f \
    read_filtered_candidates.tsv \
    candidate_transcripts.gff3 \
    Step6.RData

run_r filterReadCandidates.compat.R \
    -r 10 \
    -s 1 \
    -e 0.15 \
    -d 2500

test -s read_filtered_candidates.tsv
test -s candidate_transcripts.gff3
test -s Step6.RData

{
    printf "metric\tvalue\n"
    printf "read_stat_commands\t%s\n" "$n_commands"
    printf "read_stat_files\t%s\n" "$n_stats"
    printf "read_filtered_candidates\t%s\n" \
        "$(awk 'END{print NR > 0 ? NR-1 : 0}' read_filtered_candidates.tsv)"
    printf "candidate_transcripts\t%s\n" \
        "$(awk '$3=="mRNA"{n++} END{print n+0}' candidate_transcripts.gff3)"
} > read_filter_summary.tsv

cat read_filter_summary.tsv
echo "TEPROF2 READ FILTERING COMPLETE"
