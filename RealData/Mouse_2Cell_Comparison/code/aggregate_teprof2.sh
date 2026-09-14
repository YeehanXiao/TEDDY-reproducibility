#!/usr/bin/env bash
set -euo pipefail

INPUT_ROOT=${1:?INPUT_ROOT is required}
TEPROF2_BIN=${2:?TEPROF2_BIN is required}
ARGUMENT_FILE=${3:?ARGUMENT_FILE is required}
OUTPUT_DIR=${4:?OUTPUT_DIR is required}
shift 4

canonical_path() {
    python3 - "$1" <<'PY'
from pathlib import Path
import sys
print(Path(sys.argv[1]).expanduser().resolve())
PY
}

INPUT_ROOT=$(canonical_path "$INPUT_ROOT")
TEPROF2_BIN=$(canonical_path "$TEPROF2_BIN")
ARGUMENT_FILE=$(canonical_path "$ARGUMENT_FILE")
OUTPUT_DIR=$(canonical_path "$OUTPUT_DIR")

[[ $# -gt 0 ]] || {
    echo "At least one sample is required" >&2
    exit 1
}

SAMPLES=("$@")

CONDA=${CONDA:-conda}
TEPROF2_ENV=${TEPROF2_ENV:-teprof2_py2}
RSCRIPT=${RSCRIPT:-Rscript}
R_LIBS_EXTRA=${R_LIBS_EXTRA:-}
EXPERIMENT_PREFIX=${EXPERIMENT_PREFIX:-2cellrep}

run_py2() {
    env -u PYTHONPATH -u PYTHONHOME \
        "$CONDA" run --no-capture-output \
        -n "$TEPROF2_ENV" python "$@"
}

run_r() {
    if [[ -n "$R_LIBS_EXTRA" ]]; then
        R_LIBS="$R_LIBS_EXTRA${R_LIBS:+:$R_LIBS}" \
            "$RSCRIPT" "$@"
    else
        "$RSCRIPT" "$@"
    fi
}

test -s "$TEPROF2_BIN/annotationtpmprocess.py"
test -s "$TEPROF2_BIN/aggregateProcessedAnnotation.R"
test -s "$ARGUMENT_FILE"

awk -F '\t' '
NF != 2 || $1 == "" || $2 == "" {
    printf "Malformed argument line %d: [%s]\n", NR, $0 > "/dev/stderr"
    bad = 1
}
END {exit bad}
' "$ARGUMENT_FILE"

run_py2 -c '
import sys
import numpy
assert sys.version_info[0] == 2
print("Python: %s" % sys.version.split()[0])
print("NumPy: %s" % numpy.__version__)
'

run_r -e '
stopifnot(requireNamespace("Xmisc", quietly = TRUE))
cat("R:", as.character(getRversion()), "\n")
cat("Xmisc:", as.character(packageVersion("Xmisc")), "\n")
'

for protected in "$INPUT_ROOT" "$TEPROF2_BIN" "$ARGUMENT_FILE"; do
    [[ "$protected/" != "$OUTPUT_DIR/"* && "$OUTPUT_DIR/" != "$protected/"* ]] || {
        echo "OUTPUT_DIR must be disjoint from all inputs: $OUTPUT_DIR" >&2
        exit 1
    }
done

[[ ! -e "$OUTPUT_DIR" || -d "$OUTPUT_DIR" ]] || {
    echo "OUTPUT_DIR exists and is not a directory: $OUTPUT_DIR" >&2
    exit 1
}

[[ ! -d "$OUTPUT_DIR" || -z "$(find "$OUTPUT_DIR" -mindepth 1 -maxdepth 1 -print -quit)" ]] || {
    echo "OUTPUT_DIR must be new or empty: $OUTPUT_DIR" >&2
    exit 1
}

mkdir -p "$OUTPUT_DIR"

for sample in "${SAMPLES[@]}"; do
    input="$INPUT_ROOT/TEProf2/$sample/${sample}.stringtie.gtf_annotated_filtered_test_all"
    unfiltered="${input/_annotated_filtered_test_all/_annotated_test_all}"
    processed="${input}_c"

    test -s "$input"
    test -s "$unfiltered"

    rm -f "$processed"

    echo "[$(date)] Processing $sample"

    run_py2 \
        "$TEPROF2_BIN/annotationtpmprocess.py" \
        "$input"

    test -s "$processed"

    ln -s "$processed" \
        "$OUTPUT_DIR/${sample}.gtf_annotated_filtered_test_all_c"
done

compat="$OUTPUT_DIR/aggregateProcessedAnnotation.compat.R"

cp \
    "$TEPROF2_BIN/aggregateProcessedAnnotation.R" \
    "$compat"

python3 - "$compat" <<'PY'
from pathlib import Path
import re
import sys

path = Path(sys.argv[1])
text = path.read_text()

text = re.sub(
    r'(?m)^\s*library\(\s*["\']argparse["\']\s*\)\s*$\n?',
    '',
    text,
)

if 'Xmisc::ArgumentParser$new()' not in text:
    text, n = re.subn(
        r'parser\s*<-\s*ArgumentParser\$new\(\)',
        'parser <- Xmisc::ArgumentParser$new()',
        text,
        count=1,
    )

    if n != 1:
        raise SystemExit(
            "ArgumentParser$new() was not found exactly once"
        )

path.write_text(text)
PY

cd "$OUTPUT_DIR"

echo "[$(date)] Aggregating ${#SAMPLES[@]} samples"

run_r \
    "$compat" \
    -e "$EXPERIMENT_PREFIX" \
    -n 1 \
    -k no \
    -f yes \
    -a "$ARGUMENT_FILE"

test -s filter_combined_candidates.tsv
test -s initial_candidate_list.tsv
test -s Step4.RData

{
    printf "metric\tvalue\n"
    printf "n_samples\t%s\n" "${#SAMPLES[@]}"
    printf "filter_combined_candidates_lines\t%s\n" \
        "$(wc -l < filter_combined_candidates.tsv)"
    printf "initial_candidate_list_rows\t%s\n" \
        "$(awk 'END{print NR > 0 ? NR - 1 : 0}' initial_candidate_list.tsv)"
} > aggregation_summary.tsv

cat aggregation_summary.tsv
echo "TEPROF2 AGGREGATION COMPLETE"
