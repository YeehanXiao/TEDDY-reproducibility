#!/usr/bin/env bash
set -Eeuo pipefail

[[ $# -eq 4 ]] || { echo "Usage: bash run_lions_simuoptimal_final.sh ROOT LIONS_DIR PYTHON RSCRIPT" >&2; exit 2; }
ROOT=$1
LIONS=$2
PY=$3
RS=$4

BUILDER="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/build_lions_lcsv.py"
CHIMSORT="$LIONS/scripts/ChimericReadTool/chimSort.R"

test -x "$BUILDER"
test -s "$CHIMSORT"

STAMP=$(date +%Y%m%d_%H%M%S)

echo "LIONS final post-processing"
echo "chimSort: $CHIMSORT"
echo "simuOptimal parameters: 1 2 10 2 0.1 2 1.5"
echo

for seed in {101..109}; do
    SD="$ROOT/seed_${seed}"
    TARGET="$SD/results/LIONS_simuOptimal_by_depth"

    if [[ -d "$TARGET" ]]; then
        ARCHIVE="$SD/results/LIONS_simuOptimal_pre_XYfix_${STAMP}"
        mv "$TARGET" "$ARCHIVE"
        echo "[ARCHIVE] seed=$seed -> $ARCHIVE"
    fi

    for depth in 5x 10x 25x 50x 100x; do
        tag="official_simulated_${depth}_noise0.1"

        RAW="$SD/results/LIONS_official_by_depth/$depth/${tag}.LIONS.log"
        BAM="$SD/results/sortbam/${tag}.bam"

        OUT="$TARGET/$depth"
        PREFIX="$OUT/$tag"

        LCSV="${PREFIX}.lcsv"
        MAPPED="${PREFIX}.mapped_reads.txt"
        LION="${PREFIX}.lion"
        LOG="${PREFIX}.chimSort.log"

        test -s "$RAW"
        test -s "$BAM"
        test -s "${BAM}.bai"

        mkdir -p "$OUT"

        echo "[BUILD] seed=$seed depth=$depth"

        "$PY" "$BUILDER" \
            "$RAW" \
            "$BAM" \
            "$LCSV" \
            "$MAPPED"

        mapped=$(cat "$MAPPED")

        "$RS" --vanilla "$CHIMSORT" \
            "$LCSV" \
            "$LION" \
            "$mapped" \
            1 2 10 2 0.1 2 1.5 \
            > "$LOG" 2>&1

        test -s "$LCSV"
        test -s "$LION"
        grep -q "Output Entries" "$LOG"

        events=$(( $(wc -l < "$LION") - 1 ))

        echo "[PASS] seed=$seed depth=$depth events=$events"
    done
done

echo
echo "LIONS final post-processing completed."
