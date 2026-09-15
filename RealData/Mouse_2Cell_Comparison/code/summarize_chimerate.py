#!/usr/bin/env python3

import argparse
import csv
import hashlib
import re
from itertools import combinations
from pathlib import Path

import pandas as pd


REPS = [f"rep{i}" for i in range(1, 5)]
ATTR = re.compile(r'(\S+)\s+"([^"]*)"')
COORD = re.compile(r"^[^:\t]+:\d+-\d+$")
TARGETS = [
    ("Nelfa", "study_locus"),
    ("Zfp352", "established_control"),
    ("Lmx1a", "study_locus"),
    ("Cdk2ap1", "established_control"),
    ("Fam172a", "figure_example"),
    ("Snai1", "study_locus"),
    ("Pou6f2", "study_locus"),
]


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


parser = argparse.ArgumentParser()
parser.add_argument("--raw-dir", type=Path, required=True)
parser.add_argument("--gtf", type=Path, required=True)
parser.add_argument("--output-dir", type=Path, required=True)
args = parser.parse_args()

raw_dir = args.raw_dir.resolve()
out_dir = args.output_dir.resolve()

out_dir.mkdir(parents=True, exist_ok=True)

if not args.gtf.is_file():
    raise FileNotFoundError(f"GTF not found: {args.gtf}")


# Gene ID -> gene symbol
gene_names = {}

with args.gtf.open() as handle:
    for line in handle:
        if line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")

        if len(fields) != 9 or fields[2] != "gene":
            continue

        attributes = dict(ATTR.findall(fields[8]))
        gene_id = re.sub(r"\.\d+$", "", attributes.get("gene_id", ""))

        if gene_id:
            gene_names[gene_id] = attributes.get("gene_name", gene_id)


events = []
manifest = []

event_columns = [
    "replicate",
    "source_line",
    "gene_id_versioned",
    "gene_id",
    "gene_name",
    "gene_strand",
    "gene_region",
    "TE_family",
    "TE_strand",
    "TE_region",
    "supporting_reads",
    "raw_field_count",
]


# ChimeraTE rows may contain either 6 or 7 fields because TE strand is
# occasionally absent. Parse the fields from the right instead of assuming
# a fixed width.
for replicate in REPS:
    path = raw_dir / f"TE-initiated-{replicate}.tsv"

    if not path.is_file():
        raise FileNotFoundError(f"Raw result not found: {path}")

    n_events = 0

    with path.open() as handle:
        reader = csv.reader(handle, delimiter="\t")

        for line_no, row in enumerate(reader, 1):
            if not row:
                continue

            n_events += 1

            if len(row) == 7:
                (
                    gene_id_versioned,
                    gene_strand,
                    gene_region,
                    te_family,
                    te_strand,
                    te_region,
                    supporting_reads,
                ) = row
            elif len(row) == 6:
                (
                    gene_id_versioned,
                    gene_strand,
                    gene_region,
                    te_family,
                    te_region,
                    supporting_reads,
                ) = row
                te_strand = ""
            else:
                raise ValueError(
                    f"{path}:{line_no}: expected 6 or 7 fields, found {len(row)}"
                )

            if (
                gene_strand not in {"+", "-"}
                or not COORD.fullmatch(gene_region)
                or not te_family
                or (te_strand and te_strand not in {"+", "-"})
                or not COORD.fullmatch(te_region)
                or not supporting_reads.isdigit()
            ):
                raise ValueError(f"{path}:{line_no}: malformed ChimeraTE row")

            gene_id = re.sub(r"\.\d+$", "", gene_id_versioned)

            events.append([
                replicate,
                line_no,
                gene_id_versioned,
                gene_id,
                gene_names.get(gene_id, gene_id),
                gene_strand,
                gene_region,
                te_family,
                te_strand,
                te_region,
                int(supporting_reads),
                len(row),
            ])

    manifest.append([
        replicate,
        n_events,
        sha256(path),
        f"raw/{path.name}",
    ])


events_df = pd.DataFrame(events, columns=event_columns)
events_df["supporting_reads"] = pd.array(
    events_df["supporting_reads"],
    dtype="Int64",
)

events_df.to_csv(
    out_dir / "chimerate_events.tsv",
    sep="\t",
    index=False,
)

pd.DataFrame(
    manifest,
    columns=["replicate", "n_events", "sha256", "file"],
).to_csv(
    out_dir / "chimerate_manifest.tsv",
    sep="\t",
    index=False,
)

categories = {
    "all_TE_initiated": events_df,
    "MT2_Mm_initiated": events_df[
        events_df["TE_family"].eq("MT2_Mm")
    ],
}

summary_rows = []
pairwise_rows = []


for category, category_df in categories.items():
    category_dir = out_dir / category
    category_dir.mkdir(exist_ok=True)

    # Unique genes within each replicate
    genes_by_rep = (
        category_df
        .groupby(
            ["replicate", "gene_id", "gene_name"],
            as_index=False,
        )
        .agg(
            n_events=("gene_id", "size"),
            max_supporting_reads=("supporting_reads", "max"),
            TE_families=(
                "TE_family",
                lambda x: ";".join(sorted(set(x))),
            ),
        )
    )

    for replicate in REPS:
        (
            genes_by_rep[
                genes_by_rep["replicate"].eq(replicate)
            ]
            .drop(columns="replicate")
            .sort_values(["gene_name", "gene_id"])
            .to_csv(
                category_dir / f"{replicate}.genes.tsv",
                sep="\t",
                index=False,
            )
        )

    # Gene-level recurrence across replicates
    recurrence = (
        category_df
        .groupby(
            ["gene_id", "gene_name"],
            as_index=False,
        )
        .agg(
            n_replicates=("replicate", "nunique"),
            replicates=(
                "replicate",
                lambda x: ",".join(
                    rep for rep in REPS if rep in set(x)
                ),
            ),
            total_events=("gene_id", "size"),
            max_supporting_reads=("supporting_reads", "max"),
            TE_families=(
                "TE_family",
                lambda x: ";".join(sorted(set(x))),
            ),
        )
    )

    presence = (
        category_df
        .assign(value=1)
        .pivot_table(
            index=["gene_id", "gene_name"],
            columns="replicate",
            values="value",
            aggfunc="max",
            fill_value=0,
        )
        .reindex(columns=REPS, fill_value=0)
        .rename(columns=lambda x: f"{x}_present")
        .reset_index()
    )

    recurrence = (
        recurrence
        .merge(
            presence,
            on=["gene_id", "gene_name"],
        )
        .sort_values(
            ["n_replicates", "gene_name", "gene_id"],
            ascending=[False, True, True],
        )
    )

    recurrence.to_csv(
        category_dir / "gene_recurrence.tsv",
        sep="\t",
        index=False,
    )

    recurrence[
        recurrence["n_replicates"].ge(2)
    ].to_csv(
        category_dir / "at_least_2_reps.genes.tsv",
        sep="\t",
        index=False,
    )

    recurrence[
        recurrence["n_replicates"].ge(3)
    ].to_csv(
        category_dir / "at_least_3_reps.genes.tsv",
        sep="\t",
        index=False,
    )

    recurrence[
        recurrence["n_replicates"].eq(4)
    ].to_csv(
        category_dir / "all_4_reps.genes.tsv",
        sep="\t",
        index=False,
    )

    gene_sets = {
        replicate: set(
            category_df.loc[
                category_df["replicate"].eq(replicate),
                "gene_id",
            ]
        )
        for replicate in REPS
    }

    summary_rows.append({
        "category": category,
        **{
            f"{replicate}_events":
                int(category_df["replicate"].eq(replicate).sum())
            for replicate in REPS
        },
        **{
            f"{replicate}_genes":
                len(gene_sets[replicate])
            for replicate in REPS
        },
        "union_genes": len(
            set().union(*gene_sets.values())
        ),
        "at_least_2_reps": int(
            recurrence["n_replicates"].ge(2).sum()
        ),
        "at_least_3_reps": int(
            recurrence["n_replicates"].ge(3).sum()
        ),
        "all_4_reps": int(
            recurrence["n_replicates"].eq(4).sum()
        ),
    })

    for rep_a, rep_b in combinations(REPS, 2):
        intersection = gene_sets[rep_a] & gene_sets[rep_b]
        union = gene_sets[rep_a] | gene_sets[rep_b]

        pairwise_rows.append([
            category,
            rep_a,
            len(gene_sets[rep_a]),
            rep_b,
            len(gene_sets[rep_b]),
            len(intersection),
            len(union),
            len(intersection) / len(union) if union else pd.NA,
        ])


summary_df = pd.DataFrame(summary_rows)

summary_df.to_csv(
    out_dir / "chimerate_overlap_summary.tsv",
    sep="\t",
    index=False,
)

pd.DataFrame(
    pairwise_rows,
    columns=[
        "category",
        "rep_a",
        "n_genes_a",
        "rep_b",
        "n_genes_b",
        "intersection",
        "union",
        "jaccard",
    ],
).to_csv(
    out_dir / "chimerate_pairwise_overlap.tsv",
    sep="\t",
    index=False,
)

locus_rows = []
for gene_name, panel in TARGETS:
    gene_events = events_df[events_df["gene_name"].eq(gene_name)]
    definitions = {
        "any_TE": gene_events,
        "MT2_MERVL_family": gene_events[
            gene_events["TE_family"].str.startswith(("MT2", "MERVL"))
        ],
        "exact_MT2_Mm": gene_events[gene_events["TE_family"].eq("MT2_Mm")],
    }
    for definition, subset in definitions.items():
        detected = set(subset["replicate"])
        locus_rows.append({
            "gene_name": gene_name,
            "panel": panel,
            "definition": definition,
            "n_replicates": len(detected),
            "replicates": ",".join(rep for rep in REPS if rep in detected),
            "n_events": len(subset),
            "TE_families": ";".join(sorted(set(subset["TE_family"]))),
            "max_supporting_reads": (
                subset["supporting_reads"].max() if not subset.empty else 0
            ),
            **{f"{rep}_present": int(rep in detected) for rep in REPS},
        })

pd.DataFrame(locus_rows).to_csv(
    out_dir / "chimerate_locus_recovery.tsv",
    sep="\t",
    index=False,
)

pd.DataFrame(
    [
        ["tool", "ChimeraTE"],
        ["version", "1.2"],
        ["mode", "mode1"],
        ["analysis", "TE-initiated"],
        ["stage", "mouse_2cell"],
        ["replicates", ",".join(REPS)],
        ["strand", "rf-stranded"],
        ["threads", "8"],
        ["genome", "mm10_no_alt_analysis_set_ENCODE.fasta"],
        ["gene_annotation", "gencode.vM25.annotation.gtf"],
        ["TE_annotation", "mm10_TE_annotations.gtf"],
        ["overlap_unit", "version-stripped Ensembl gene_id"],
        ["MT2_Mm_rule", "exact TE_family == MT2_Mm"],
    ],
    columns=["key", "value"],
).to_csv(
    out_dir / "chimerate_run_metadata.tsv",
    sep="\t",
    index=False,
)

print(summary_df.to_string(index=False))
