#!/usr/bin/env python3

import argparse
import itertools
import os
import re
from pathlib import Path

import pandas as pd


REPS = ("rep1", "rep2", "rep3", "rep4")
TARGETS = (
    "Nelfa",
    "Zfp352",
    "Lmx1a",
    "Cdk2ap1",
    "Fam172a",
    "Snai1",
    "Pou5f1",
)


def pick_column(df, candidates, fallback=None):
    for col in candidates:
        if col in df.columns:
            return col
    if fallback is not None:
        return df.columns[fallback]
    raise ValueError(
        "None of the expected columns were found: "
        + ", ".join(candidates)
    )


parser = argparse.ArgumentParser()
parser.add_argument("--final", type=Path, required=True)
parser.add_argument("--read-stats", type=Path, required=True)
parser.add_argument("--output-dir", type=Path, required=True)
args = parser.parse_args()

args.output_dir.mkdir(parents=True, exist_ok=True)

final = pd.read_csv(args.final, sep="\t", dtype=str)

uid_col = pick_column(
    final,
    ("uniqid", "uniqueid", "candidate_id"),
    fallback=0,
)

te_col = pick_column(
    final,
    ("repName", "repeatName", "repeat_name", "TE_name", "te_name"),
    fallback=1,
)

gene_col = pick_column(
    final,
    ("gene2", "gene_name", "gene", "splice_gene"),
)

final = final.rename(
    columns={
        uid_col: "candidate_id",
        te_col: "TE_family",
        gene_col: "gene_name",
    }
)

final = final[
    final["gene_name"].notna()
    & ~final["gene_name"].isin(("None", "NA", ""))
].copy()

# TEProf2 filter_read_stats.txt:
# filename.stats  read  startread  endread  filetype
stats = pd.read_csv(
    args.read_stats,
    sep="\t",
    header=None,
    names=[
        "stats_file",
        "read",
        "startread",
        "endread",
        "filetype",
    ],
)

for col in ("read", "startread", "endread"):
    stats[col] = pd.to_numeric(
        stats[col],
        errors="coerce",
    ).fillna(0).astype(int)

def parse_stats_name(value):
    name = os.path.basename(str(value))
    name = re.sub(r"\.stats$", "", name)

    match = re.match(
        r"^(.*)--(2cellrep[1-4])$",
        name,
    )

    if match is None:
        raise ValueError(
            f"Cannot parse TEProf2 stats filename: {value}"
        )

    candidate_id, sample = match.groups()

    return pd.Series(
        {
            "candidate_id": candidate_id,
            "replicate": sample.replace("2cell", ""),
        }
    )

parsed = stats["stats_file"].apply(parse_stats_name)
stats = pd.concat([stats, parsed], axis=1)

# Restrict to final Step-6 candidates.
stats = stats[
    stats["candidate_id"].isin(final["candidate_id"])
].copy()

# Conservative replicate-level support:
# the same replicate must independently satisfy both TE-read and
# TE-to-gene start-read requirements.
stats["supported"] = (
    (stats["read"] >= 10)
    & (stats["startread"] >= 1)
)

# Candidate-level annotated event table.
support = (
    stats
    .pivot_table(
        index="candidate_id",
        columns="replicate",
        values="supported",
        aggfunc="max",
        fill_value=False,
    )
    .reindex(columns=REPS, fill_value=False)
    .reset_index()
)

support = support.rename(
    columns={rep: f"{rep}_supported" for rep in REPS}
)

events = (
    final
    .merge(
        support,
        on="candidate_id",
        how="left",
    )
)

for rep in REPS:
    col = f"{rep}_supported"
    events[col] = events[col].fillna(False).astype(bool)

events["n_replicates"] = events[
    [f"{rep}_supported" for rep in REPS]
].sum(axis=1)

events.to_csv(
    args.output_dir / "teprof2_final_events.tsv",
    sep="\t",
    index=False,
)

# Gene-level recurrence.
gene_presence = (
    events[
        ["gene_name"]
        + [f"{rep}_supported" for rep in REPS]
    ]
    .groupby("gene_name", as_index=False)
    .max()
)

gene_te = (
    events
    .groupby("gene_name")["TE_family"]
    .agg(lambda x: ";".join(sorted(set(x))))
    .rename("TE_families")
    .reset_index()
)

gene_recurrence = gene_presence.merge(
    gene_te,
    on="gene_name",
    how="left",
)

gene_recurrence["n_replicates"] = gene_recurrence[
    [f"{rep}_supported" for rep in REPS]
].sum(axis=1)

gene_recurrence = gene_recurrence.sort_values(
    ["n_replicates", "gene_name"],
    ascending=[False, True],
)

gene_recurrence.to_csv(
    args.output_dir / "teprof2_gene_recurrence.tsv",
    sep="\t",
    index=False,
)

# Overall summary.
rep_gene_sets = {
    rep: set(
        gene_recurrence.loc[
            gene_recurrence[f"{rep}_supported"],
            "gene_name",
        ]
    )
    for rep in REPS
}

summary = pd.DataFrame(
    [
        {
            "final_candidates": len(final),
            "final_unique_genes": final["gene_name"].nunique(),
            **{
                f"{rep}_supported_candidates":
                    int(events[f"{rep}_supported"].sum())
                for rep in REPS
            },
            **{
                f"{rep}_supported_genes":
                    len(rep_gene_sets[rep])
                for rep in REPS
            },
            "union_supported_genes":
                len(set().union(*rep_gene_sets.values())),
            "at_least_2_reps":
                int((gene_recurrence["n_replicates"] >= 2).sum()),
            "at_least_3_reps":
                int((gene_recurrence["n_replicates"] >= 3).sum()),
            "all_4_reps":
                int((gene_recurrence["n_replicates"] == 4).sum()),
        }
    ]
)

summary.to_csv(
    args.output_dir / "teprof2_overlap_summary.tsv",
    sep="\t",
    index=False,
)

pairwise = []
for rep_a, rep_b in itertools.combinations(REPS, 2):
    intersection = rep_gene_sets[rep_a] & rep_gene_sets[rep_b]
    union = rep_gene_sets[rep_a] | rep_gene_sets[rep_b]
    pairwise.append({
        "replicate_1": rep_a,
        "replicate_2": rep_b,
        "intersection": len(intersection),
        "union": len(union),
        "jaccard": len(intersection) / len(union) if union else 0,
    })

pd.DataFrame(pairwise).to_csv(
    args.output_dir / "teprof2_pairwise_overlap.tsv",
    sep="\t",
    index=False,
)

# Locus recovery.
rows = []

for gene in TARGETS:
    x = events[events["gene_name"] == gene]

    row = {
        "gene_name": gene,
        "final_candidate": int(not x.empty),
        "n_candidates": len(x),
        "TE_families": (
            ";".join(sorted(set(x["TE_family"])))
            if not x.empty else ""
        ),
    }

    for rep in REPS:
        families = sorted(set(x.loc[x[f"{rep}_supported"], "TE_family"]))
        row[f"{rep}_supported"] = int(bool(families))
        row[f"{rep}_TE_families"] = ";".join(families)

    row["n_replicates"] = sum(
        row[f"{rep}_supported"]
        for rep in REPS
    )

    rows.append(row)

pd.DataFrame(rows).to_csv(
    args.output_dir / "teprof2_locus_recovery.tsv",
    sep="\t",
    index=False,
)

print("===== TEPROF2 SUMMARY =====")
print(summary.to_string(index=False))

print("\n===== LOCUS RECOVERY =====")
print(
    pd.DataFrame(rows).to_string(index=False)
)
