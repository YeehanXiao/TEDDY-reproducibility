#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Directly assess whether simulated exon-exon TE-host junctions are represented
by split reads containing the expected CIGAR N operation.

This analysis is independent of the intra-exonic breakpoint-support analysis
implemented in scripts 01-03.
"""

import argparse
import csv
import glob
import os
import re
from collections import defaultdict

import pysam


DEPTH_ORDER = {
    "5x": 0,
    "10x": 1,
    "25x": 2,
    "50x": 3,
    "100x": 4
}


def extract_depth(path):
    match = re.search(
        r"(100x|50x|25x|10x|5x)",
        os.path.basename(path)
    )

    return (
        match.group(1)
        if match
        else os.path.basename(path)
    )


def load_tpm(path):
    with open(path) as handle:
        return {
            row["transcript_id"]: float(row["TPM"])
            for row in csv.DictReader(
                handle,
                delimiter="\t"
            )
        }


def load_expressed_truth_junctions(
    path,
    tpm
):
    events = {}

    with open(path) as handle:
        for row in csv.DictReader(
            handle,
            delimiter="\t"
        ):
            if row.get("junction_type") not in {
                "TE_to_host",
                "host_to_TE"
            }:
                continue

            transcript_id = row["transcript_id"]
            transcript_tpm = tpm.get(
                transcript_id,
                0.0
            )

            if transcript_tpm <= 1:
                continue

            chrom = (
                row.get("seqnames")
                or row.get("chrom")
            )

            donor = int(
                float(row["donor_site"])
            )
            acceptor = int(
                float(row["acceptor_site"])
            )

            left = min(donor, acceptor)
            right = max(donor, acceptor)

            event_id = (
                transcript_id,
                chrom,
                left,
                right
            )

            events[event_id] = {
                "transcript_id": transcript_id,
                "gene_id": row.get(
                    "gene_id",
                    "NA"
                ),
                "gene_name": row.get(
                    "gene_name",
                    "NA"
                ),
                "chrom": chrom,
                "strand": row.get(
                    "strand",
                    "NA"
                ),
                "junction_type": row.get(
                    "junction_type",
                    "NA"
                ),
                "donor_site": donor,
                "acceptor_site": acceptor,
                "junction_left": left,
                "junction_right": right,
                "TPM": transcript_tpm
            }

    return list(events.values())


def cigar_n_pairs(read):
    """
    Return exon-boundary pairs represented by CIGAR N operations.

    pysam reference_start is 0-based. For each N operation:

      left exon end, 1-based    = current reference position
      right exon start, 1-based = current reference position + N length + 1
    """

    if read.cigartuples is None:
        return []

    reference_position = read.reference_start
    pairs = []

    for operation, length in read.cigartuples:
        if operation == 3:
            pairs.append((
                reference_position,
                reference_position + length + 1
            ))

        if operation in (0, 2, 3, 7, 8):
            reference_position += length

    return pairs


def process_bam(
    bam_path,
    truth_events
):
    locus_to_events = defaultdict(list)

    for index, event in enumerate(
        truth_events
    ):
        key = (
            event["chrom"],
            event["junction_left"],
            event["junction_right"]
        )

        locus_to_events[key].append(index)

    read_names = [
        set()
        for _ in truth_events
    ]

    mapped_alignments = 0
    alignments_with_n = 0

    with pysam.AlignmentFile(
        bam_path,
        "rb"
    ) as bam:

        for read in bam.fetch(
            until_eof=True
        ):
            if (
                read.is_unmapped
                or read.cigartuples is None
            ):
                continue

            mapped_alignments += 1
            n_pairs = cigar_n_pairs(read)

            if not n_pairs:
                continue

            alignments_with_n += 1
            chrom = read.reference_name

            for left, right in n_pairs:
                key = (
                    chrom,
                    left,
                    right
                )

                for event_index in locus_to_events.get(
                    key,
                    []
                ):
                    read_names[event_index].add(
                        read.query_name
                    )

    depth = extract_depth(bam_path)

    detail = []

    for event, reads in zip(
        truth_events,
        read_names
    ):
        detail.append({
            "Depth": depth,
            **event,
            "CIGAR_N_split_reads": len(reads),
            "CIGAR_N_supported": int(
                bool(reads)
            )
        })

    unique_truth_loci = {
        (
            event["chrom"],
            event["junction_left"],
            event["junction_right"]
        )
        for event in truth_events
    }

    supported_unique_loci = {
        (
            event["chrom"],
            event["junction_left"],
            event["junction_right"]
        )
        for event, reads in zip(
            truth_events,
            read_names
        )
        if reads
    }

    supported_events = sum(
        bool(reads)
        for reads in read_names
    )

    summary = {
        "Depth": depth,
        "BAM": os.path.basename(bam_path),
        "mapped_alignments": mapped_alignments,
        "alignments_with_CIGAR_N": (
            alignments_with_n
        ),
        "CIGAR_N_alignment_rate_pct": round(
            alignments_with_n
            / mapped_alignments
            * 100,
            2
        ) if mapped_alignments else 0,

        "expressed_truth_transcript_junctions": (
            len(truth_events)
        ),
        "exact_CIGAR_N_supported_transcript_junctions": (
            supported_events
        ),
        "exact_CIGAR_N_support_rate_pct": round(
            supported_events
            / len(truth_events)
            * 100,
            2
        ) if truth_events else 0,

        "expressed_truth_unique_genomic_junctions": (
            len(unique_truth_loci)
        ),
        "exact_CIGAR_N_supported_unique_genomic_junctions": (
            len(supported_unique_loci)
        ),
        "exact_CIGAR_N_unique_locus_support_rate_pct": round(
            len(supported_unique_loci)
            / len(unique_truth_loci)
            * 100,
            2
        ) if unique_truth_loci else 0
    }

    return detail, summary


def write_tsv(rows, path):
    if not rows:
        raise ValueError(
            f"No rows available for {path}"
        )

    out_dir = os.path.dirname(path)

    if out_dir:
        os.makedirs(
            out_dir,
            exist_ok=True
        )

    with open(
        path,
        "w",
        newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0].keys()),
            delimiter="\t"
        )

        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--junctions",
        required=True,
        help="Simulated truth-junction TSV"
    )

    parser.add_argument(
        "--isoforms",
        required=True,
        help="Simulated isoform-expression file"
    )

    parser.add_argument(
        "--bam_dir",
        required=True,
        help="Directory containing depth-specific BAM files"
    )

    parser.add_argument(
        "--out_detail",
        required=True
    )

    parser.add_argument(
        "--out_summary",
        required=True
    )

    args = parser.parse_args()

    tpm = load_tpm(args.isoforms)

    truth_events = (
        load_expressed_truth_junctions(
            args.junctions,
            tpm
        )
    )

    print(
        "Expressed truth transcript-junction events:",
        len(truth_events)
    )

    bam_files = sorted(
        glob.glob(
            os.path.join(
                args.bam_dir,
                "*.bam"
            )
        ),
        key=lambda path: DEPTH_ORDER.get(
            extract_depth(path),
            999
        )
    )

    if not bam_files:
        raise FileNotFoundError(
            f"No BAM files found in "
            f"{args.bam_dir}"
        )

    all_detail = []
    all_summary = []

    for bam_path in bam_files:
        print(
            "Scanning:",
            os.path.basename(bam_path)
        )

        detail, summary = process_bam(
            bam_path,
            truth_events
        )

        all_detail.extend(detail)
        all_summary.append(summary)

        print(
            " ",
            summary["Depth"],
            f'{summary["exact_CIGAR_N_supported_unique_genomic_junctions"]}/'
            f'{summary["expressed_truth_unique_genomic_junctions"]}',
            f'({summary["exact_CIGAR_N_unique_locus_support_rate_pct"]}%)'
        )

    write_tsv(
        all_detail,
        args.out_detail
    )

    write_tsv(
        all_summary,
        args.out_summary
    )

    print("Saved:")
    print(args.out_detail)
    print(args.out_summary)


if __name__ == "__main__":
    main()