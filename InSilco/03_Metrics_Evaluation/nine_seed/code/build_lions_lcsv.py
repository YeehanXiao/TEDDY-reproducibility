#!/usr/bin/env python3

import csv
import os
import sys
import pysam

if len(sys.argv) != 5:
    raise SystemExit(
        "usage: build_lions_lcsv.py RAW_LOG BAM OUT_LCSV OUT_MAPPED"
    )

raw_log, bam_path, out_lcsv, out_mapped = sys.argv[1:]

# ------------------------------------------------------------
# 1. Read every genuine 25-column chimericReadSearch record.
#    Chromosome is deliberately kept as a string (X/Y allowed).
# ------------------------------------------------------------
raw_rows = []

with open(raw_log, errors="replace") as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")

        if len(p) != 25:
            continue

        if ":" not in p[0]:
            continue

        raw_rows.append(p)

if not raw_rows:
    raise SystemExit("No 25-column interaction records found: " + raw_log)

# ------------------------------------------------------------
# 2. Coverage from the existing BAM.
#    q10.F772 corresponds to MAPQ >=10 and exclusion of
#    unmapped / secondary / QC-fail reads.
# ------------------------------------------------------------
bam = pysam.AlignmentFile(bam_path, "rb")
refs = set(bam.references)

EXCLUDE = 4 | 256 | 512

def read_ok(read):
    return (
        read.mapping_quality >= 10
        and (read.flag & EXCLUDE) == 0
    )

def resolve_contig(chrom):
    c = str(chrom)

    candidates = [
        c,
        c.removeprefix("chr"),
        "chr" + c.removeprefix("chr")
    ]

    for x in candidates:
        if x in refs:
            return x

    return None

peak_cache = {}

def peak(chrom, start, end):
    start = max(0, int(start))
    end = max(start, int(end))
    key = (str(chrom), start, end)

    if key in peak_cache:
        return peak_cache[key]

    contig = resolve_contig(chrom)

    if contig is None or end <= start:
        peak_cache[key] = 0
        return 0

    cov = bam.count_coverage(
        contig,
        start,
        end,
        quality_threshold=0,
        read_callback=read_ok
    )

    if len(cov[0]) == 0:
        value = 0
    else:
        value = max(
            a + c + g + t
            for a, c, g, t in zip(
                cov[0], cov[1], cov[2], cov[3]
            )
        )

    peak_cache[key] = int(value)
    return int(value)

# ------------------------------------------------------------
# 3. Construct the chimSort input.
#    RPKM fields are left NA rather than fabricated.
#    Peak-coverage fields are populated explicitly.
# ------------------------------------------------------------
header = [
    "transcriptID",
    "exonRankInTranscript",
    "repeatName",
    "coordinates",
    "ER_Interaction",
    "IsExonic",
    "ExonsOverlappingWithRepeat",
    "ER",
    "DR",
    "DE",
    "DD",
    "Total",
    "Chromosome",
    "EStart",
    "EEnd",
    "RStart",
    "REnd",
    "EStrand",
    "RStrand",
    "RepeatRank",
    "UpExonStart",
    "UpExonEnd",
    "UpThread",
    "DownThread",
    "ExonInGene",
    "ExonRPKM",
    "ExonMax",
    "UpExonRPKM",
    "UpExonMax",
    "RepeatRPKM",
    "RepeatMaxCoverage",
    "UpstreamRepeatRPKM",
    "UpstreamRepeatMaxCoverage"
]

os.makedirs(os.path.dirname(out_lcsv), exist_ok=True)

with open(out_lcsv, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t", lineterminator="\n")
    w.writerow(header)

    for p in raw_rows:
        chrom = p[12]

        estart = int(p[13])
        eend   = int(p[14])
        rstart = int(p[15])
        rend   = int(p[16])
        estrand = int(p[17])

        upstart = int(p[20])
        upend   = int(p[21])

        exon_max   = peak(chrom, estart, eend)
        upexon_max = peak(chrom, upstart, upend)
        repeat_max = peak(chrom, rstart, rend)

        if estrand == 1:
            us = max(0, rstart - 50)
            ue = rstart
        else:
            us = rend
            ue = rend + 50

        upstream_repeat_max = peak(chrom, us, ue)

        extra = [
            "NA", str(exon_max),
            "NA", str(upexon_max),
            "NA", str(repeat_max),
            "NA", str(upstream_repeat_max)
        ]

        w.writerow(p + extra)

mapped = sum(x.mapped for x in bam.get_index_statistics())
bam.close()

with open(out_mapped, "w") as fh:
    fh.write(str(mapped) + "\n")

print(
    f"raw_interactions={len(raw_rows)} "
    f"unique_regions={len(peak_cache)} "
    f"mapped_reads={mapped}"
)
