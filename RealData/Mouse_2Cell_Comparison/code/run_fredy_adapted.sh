#!/usr/bin/env bash
set -euo pipefail

DEFAULT_DIR=${1:?DEFAULT_DIR required}
OUTPUT_DIR=${2:?OUTPUT_DIR required}
CALLER=${3:?CALLER required}
ANNOTATION=${4:?ANNOTATION required}
GENOME=${5:?GENOME required}
TE_BED=${6:?TE_BED required}

canonical_path() {
    python3 - "$1" <<'PY'
from pathlib import Path
import sys
print(Path(sys.argv[1]).expanduser().resolve())
PY
}

DEFAULT_DIR=$(canonical_path "$DEFAULT_DIR")
OUTPUT_DIR=$(canonical_path "$OUTPUT_DIR")
CALLER=$(canonical_path "$CALLER")
ANNOTATION=$(canonical_path "$ANNOTATION")
GENOME=$(canonical_path "$GENOME")
TE_BED=$(canonical_path "$TE_BED")

for file in "$CALLER" "$ANNOTATION" "$GENOME" "$TE_BED"; do
    test -s "$file" || {
        echo "Missing input: $file" >&2
        exit 1
    }
done

for protected in "$DEFAULT_DIR" "$CALLER" "$ANNOTATION" "$GENOME" "$TE_BED"; do
    [[ "$protected/" != "$OUTPUT_DIR/"* && "$OUTPUT_DIR/" != "$protected/"* ]] || {
        echo "OUTPUT_DIR must be disjoint from all inputs: $OUTPUT_DIR" >&2
        exit 1
    }
done

for n in 1 2 3 4; do
    test -s "$DEFAULT_DIR/2cellrep${n}/string/merge.gtf" || {
        echo "Missing input: $DEFAULT_DIR/2cellrep${n}/string/merge.gtf" >&2
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
printf "replicate\tstatus\tn_transcripts\n" > "$OUTPUT_DIR/status.tsv"

for n in 1 2 3 4; do
    rep="rep${n}"
    sample="2cellrep${n}"
    source_gtf="$DEFAULT_DIR/$sample/string/merge.gtf"
    out="$OUTPUT_DIR/$sample"
    result="$out/chimeric/protein.gtf"

    mkdir -p "$out/string" "$out/tmp"
    ln -s "$source_gtf" "$out/string/merge.gtf"

    echo "[$(date)] START $rep"

    bash "$CALLER" \
        -o "$out" \
        -a "$ANNOTATION" \
        -g "$GENOME" \
        -e "$TE_BED" \
        -T "$out/tmp" \
        -R

    test -s "$result"

    n_tx=$(awk '$3=="transcript"{n++} END{print n+0}' "$result")
    printf "%s\tcomplete\t%s\n" "$rep" "$n_tx" \
        >> "$OUTPUT_DIR/status.tsv"

    echo "[$(date)] DONE $rep: $n_tx transcripts"
done

cat "$OUTPUT_DIR/status.tsv"
