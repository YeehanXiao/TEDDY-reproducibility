#!/usr/bin/env python3

import re
import sys
from collections import Counter
from pathlib import Path

src = Path(sys.argv[1])
dst = Path(sys.argv[2])

feature_order = {
    "transcript": 0,
    "exon": 1,
    "start_codon": 2,
    "stop_codon": 3
}

records = []
transcript_count = Counter()
feature_tx = set()

with src.open() as handle:
    for line_number, line in enumerate(handle, 1):
        if not line.strip() or line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")

        if len(fields) != 9:
            raise SystemExit(
                f"Invalid GTF line {line_number}: {len(fields)} columns"
            )

        feature = fields[2]

        if feature not in feature_order:
            continue

        match = re.search(
            r'(?:^|;\s*)transcript_id "([^"]+)"',
            fields[8]
        )

        if match is None:
            raise SystemExit(
                f"Missing transcript_id at GTF line {line_number}"
            )

        tx_id = match.group(1)
        records.append((tx_id, fields))

        if feature == "transcript":
            transcript_count[tx_id] += 1
        else:
            feature_tx.add(tx_id)

missing_transcript = sorted(feature_tx - set(transcript_count))

if missing_transcript:
    raise SystemExit(
        "Features without transcript records: "
        + ", ".join(missing_transcript[:10])
    )

duplicate_transcript = sorted(
    tx for tx, count in transcript_count.items()
    if count != 1
)

if duplicate_transcript:
    raise SystemExit(
        "Transcript IDs with non-unique transcript records: "
        + ", ".join(duplicate_transcript[:10])
    )

records.sort(
    key=lambda x: (
        x[0],
        feature_order[x[1][2]],
        int(x[1][3]),
        int(x[1][4])
    )
)

with dst.open("w") as handle:
    for tx_id, fields in records:
        # genecode_to_dic.py expects a tenth field identifying the transcript.
        handle.write(
            "\t".join(
                fields + [f'transcript_id "{tx_id}"']
            )
            + "\n"
        )

print(f"records={len(records)}")
print(f"transcripts={len(transcript_count)}")
print(f"output={dst}")
