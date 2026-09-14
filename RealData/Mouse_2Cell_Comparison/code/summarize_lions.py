#!/usr/bin/env python3

import argparse
import itertools
from pathlib import Path

import pandas as pd

REPS = ["rep1", "rep2", "rep3", "rep4"]

TARGETS = [
    "Nelfa",
    "Zfp352",
    "Lmx1a",
    "Cdk2ap1",
    "Fam172a",
    "Snai1",
    "Pou5f1",
]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--work-dir", type=Path, required=True)
    p.add_argument("--filter-qc", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    all_calls = []
    for i, rep in enumerate(REPS, 1):
        f = args.work_dir / f"2cellrep{i}.lion"

        x = pd.read_csv(f, sep="\t", dtype=str)

        if "transcriptID" not in x.columns:
            raise RuntimeError(f"{f}: transcriptID column missing")

        if "repeatName" not in x.columns:
            raise RuntimeError(f"{f}: repeatName column missing")

        x["gene_name"] = (
            x["transcriptID"]
            .astype(str)
            .str.split(":")
            .str[0]
        )

        x["TE_name"] = (
            x["repeatName"]
            .astype(str)
            .str.split(":")
            .str[0]
        )

        x["replicate"] = rep

        all_calls.append(x)

    calls = pd.concat(all_calls, ignore_index=True)

    qc = pd.read_csv(args.filter_qc, sep="\t")
    required_qc = {
        "replicate",
        "raw_interactions",
        "mapped_reads",
        "effective_scREADS",
        "calls",
        "genes",
    }
    if set(qc.columns) != required_qc or set(qc["replicate"]) != set(REPS):
        raise RuntimeError(f"{args.filter_qc}: invalid LIONS filter QC")
    observed = calls.groupby("replicate").agg(
        calls=("gene_name", "size"),
        genes=("gene_name", "nunique"),
    )
    expected = qc.set_index("replicate")[["calls", "genes"]].astype(int)
    if not observed.equals(expected.loc[observed.index]):
        raise RuntimeError("LIONS calls do not match filter QC")

    calls.to_csv(
        args.output_dir / "lions_final_calls_by_replicate.tsv",
        sep="\t",
        index=False,
    )

    rep_sets = {
        rep: set(
            calls.loc[
                calls["replicate"] == rep,
                "gene_name"
            ]
        )
        for rep in REPS
    }

    union = set().union(*rep_sets.values())

    rows = []

    for gene in sorted(union):
        y = calls[calls["gene_name"] == gene]

        row = {
            "gene_name": gene,
            **{
                f"{rep}_supported":
                    gene in rep_sets[rep]
                for rep in REPS
            },
            "TE_families":
                ";".join(sorted(set(y["TE_name"]))),
        }

        row["n_replicates"] = sum(
            row[f"{rep}_supported"]
            for rep in REPS
        )

        rows.append(row)

    rec = pd.DataFrame(rows)

    rec = rec.sort_values(
        ["n_replicates", "gene_name"],
        ascending=[False, True],
    )

    rec.to_csv(
        args.output_dir / "lions_gene_recurrence.tsv",
        sep="\t",
        index=False,
    )

    summary = pd.DataFrame([{
        "final_calls": len(calls),
        **{
            f"{rep}_calls":
                int((calls["replicate"] == rep).sum())
            for rep in REPS
        },
        **{
            f"{rep}_genes":
                len(rep_sets[rep])
            for rep in REPS
        },
        "union_genes": len(union),
        "at_least_2_reps":
            int((rec["n_replicates"] >= 2).sum()),
        "at_least_3_reps":
            int((rec["n_replicates"] >= 3).sum()),
        "all_4_reps":
            int((rec["n_replicates"] == 4).sum()),
    }])

    summary.to_csv(
        args.output_dir / "lions_overlap_summary.tsv",
        sep="\t",
        index=False,
    )

    pairwise = []

    for a, b in itertools.combinations(REPS, 2):
        inter = rep_sets[a] & rep_sets[b]
        uni = rep_sets[a] | rep_sets[b]

        pairwise.append({
            "replicate_1": a,
            "replicate_2": b,
            "intersection": len(inter),
            "union": len(uni),
            "jaccard":
                len(inter) / len(uni)
                if len(uni) else 0,
        })

    pd.DataFrame(pairwise).to_csv(
        args.output_dir / "lions_pairwise_overlap.tsv",
        sep="\t",
        index=False,
    )

    locus_rows = []

    for gene in TARGETS:
        y = calls[calls["gene_name"] == gene]

        row = {
            "gene_name": gene,
            "final_candidate": int(not y.empty),
            "n_calls": len(y),
            "TE_families":
                ";".join(sorted(set(y["TE_name"])))
                if not y.empty else "",
        }

        for rep in REPS:
            families = sorted(set(y.loc[y["replicate"].eq(rep), "TE_name"]))
            row[f"{rep}_supported"] = int(bool(families))
            row[f"{rep}_TE_families"] = ";".join(families)

        row["n_replicates"] = sum(
            row[f"{rep}_supported"]
            for rep in REPS
        )

        locus_rows.append(row)

    pd.DataFrame(locus_rows).to_csv(
        args.output_dir /
        "lions_locus_recovery.tsv",
        sep="\t",
        index=False,
    )

    pd.DataFrame(qc).to_csv(
        args.output_dir / "lions_qc.tsv",
        sep="\t",
        index=False,
    )

    print("===== QC =====")
    print(pd.DataFrame(qc).to_string(index=False))

    print("\n===== OVERLAP =====")
    print(summary.to_string(index=False))

    print("\n===== LOCUS RECOVERY =====")
    print(pd.DataFrame(locus_rows).to_string(index=False))


if __name__ == "__main__":
    main()
