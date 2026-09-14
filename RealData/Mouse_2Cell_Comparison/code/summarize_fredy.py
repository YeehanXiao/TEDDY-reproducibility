#!/usr/bin/env python3

import argparse
import csv
import re
from collections import defaultdict
from itertools import combinations
from pathlib import Path

REPS = ("rep1", "rep2", "rep3", "rep4")
SAMPLES = ("2cellrep1", "2cellrep2", "2cellrep3", "2cellrep4")
TARGETS = (
    "Nelfa", "Zfp352", "Lmx1a", "Cdk2ap1",
    "Fam172a", "Snai1", "Pou5f1",
)
ATTR = re.compile(r'(\S+)\s+"([^"]*)"')

parser = argparse.ArgumentParser()
parser.add_argument("--fredy-dir", type=Path, required=True)
parser.add_argument("--gtf", type=Path, required=True)
parser.add_argument("--output-dir", type=Path, required=True)
args = parser.parse_args()
args.output_dir.mkdir(parents=True, exist_ok=True)

gene_names = {}

with args.gtf.open() as handle:
    for line in handle:
        fields = line.rstrip().split("\t")

        if len(fields) == 9 and fields[2] == "gene":
            attrs = dict(ATTR.findall(fields[8]))
            gene_id = re.sub(
                r"\.\d+$", "", attrs.get("gene_id", "")
            )

            if gene_id:
                gene_names[gene_id] = attrs.get(
                    "gene_name", gene_id
                )

events = []
native_counts = {}

for replicate, sample in zip(REPS, SAMPLES):
    result_dir = args.fredy_dir / sample / "chimeric"

    initial = []

    with (result_dir / "info.tsv").open() as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) == 5 and row[4] == "Novel Initial":
                initial.append((row[0], row[3]))

    native_counts[replicate] = (
        len(initial),
        len({row[0] for row in initial}),
    )

    gene_map = defaultdict(set)

    with (result_dir / "most_shared.tsv").open() as handle:
        for transcript_id, gene_id, _ in csv.reader(
            handle, delimiter="\t"
        ):
            gene_map[transcript_id].add(
                re.sub(r"\.\d+$", "", gene_id)
            )

    te_map = defaultdict(set)

    with (result_dir / "chimeric.gtf").open() as handle:
        for line in handle:
            fields = line.rstrip().split("\t")

            if len(fields) != 9 or fields[2] != "transcript":
                continue

            attrs = dict(ATTR.findall(fields[8]))
            transcript_id = attrs.get("transcript_id")
            exon_number = attrs.get("chimeric_exon_number")

            if transcript_id and exon_number:
                te_map[(transcript_id, exon_number)].update(
                    family.strip()
                    for family in attrs.get(
                        "chimeric_event", ""
                    ).split(",")
                    if family.strip()
                )

    for transcript_id, exon_number in initial:
        for gene_id in gene_map[transcript_id]:
            events.append({
                "replicate": replicate,
                "transcript_id": transcript_id,
                "chimeric_exon": exon_number,
                "gene_id": gene_id,
                "gene_name": gene_names.get(gene_id, gene_id),
                "n_gene_maps": len(gene_map[transcript_id]),
                "TE_families": ";".join(sorted(
                    te_map[(transcript_id, exon_number)]
                )),
            })

if not events:
    raise ValueError("No FREDY Novel Initial events were mapped")

event_columns = list(events[0])

with (
    args.output_dir / "fredy_novel_initial_events.tsv"
).open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle, event_columns, delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(events)

scopes = {
    "inclusive": events,
    "unambiguous": [
        event for event in events
        if event["n_gene_maps"] == 1
    ],
}

summary = []

for scope, rows in scopes.items():
    gene_sets = {
        replicate: {
            row["gene_id"]
            for row in rows
            if row["replicate"] == replicate
        }
        for replicate in REPS
    }

    union = set().union(*gene_sets.values())
    recurrence = {
        gene_id: sum(
            gene_id in genes
            for genes in gene_sets.values()
        )
        for gene_id in union
    }

    record = {"mapping_scope": scope}

    for replicate in REPS:
        record[f"{replicate}_calls"] = \
            native_counts[replicate][0]
        record[f"{replicate}_transcripts"] = \
            native_counts[replicate][1]
        record[f"{replicate}_genes"] = \
            len(gene_sets[replicate])

    record["union_genes"] = len(union)
    record["at_least_2_reps"] = sum(
        value >= 2 for value in recurrence.values()
    )
    record["at_least_3_reps"] = sum(
        value >= 3 for value in recurrence.values()
    )
    record["all_4_reps"] = sum(
        value == 4 for value in recurrence.values()
    )

    for left, right in combinations(REPS, 2):
        intersection = gene_sets[left] & gene_sets[right]
        pair_union = gene_sets[left] | gene_sets[right]

        record[f"jaccard_{left}_{right}"] = (
            len(intersection) / len(pair_union)
            if pair_union else "NA"
        )

    summary.append(record)

with (
    args.output_dir / "fredy_overlap_summary.tsv"
).open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle, list(summary[0]), delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(summary)

locus_rows = []

for scope, rows in scopes.items():
    for gene_name in TARGETS:
        gene_rows = [
            row for row in rows
            if row["gene_name"] == gene_name
        ]

        for definition in (
            "any_TE",
            "MT2_MERVL_family",
            "exact_MT2_Mm",
        ):
            retained = []

            for row in gene_rows:
                families = (
                    row["TE_families"].split(";")
                    if row["TE_families"] else []
                )

                keep = (
                    definition == "any_TE"
                    or (
                        definition == "MT2_MERVL_family"
                        and any(
                            family.startswith(("MT2", "MERVL"))
                            for family in families
                        )
                    )
                    or (
                        definition == "exact_MT2_Mm"
                        and "MT2_Mm" in families
                    )
                )

                if keep:
                    retained.append(row)

            detected_reps = [
                replicate for replicate in REPS
                if any(
                    row["replicate"] == replicate
                    for row in retained
                )
            ]

            locus_rows.append({
                "gene_name": gene_name,
                "mapping_scope": scope,
                "definition": definition,
                "n_replicates": len(detected_reps),
                "replicates": ",".join(detected_reps),
                "n_transcripts": len({
                    row["transcript_id"]
                    for row in retained
                }),
                "TE_families": ";".join(sorted({
                    family
                    for row in retained
                    for family in row["TE_families"].split(";")
                    if family
                })),
            })

with (
    args.output_dir / "fredy_locus_recovery.tsv"
).open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle, list(locus_rows[0]), delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(locus_rows)
