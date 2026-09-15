#!/usr/bin/env python3

import argparse
import itertools
import math
import statistics
from pathlib import Path

import pandas as pd


REPS = tuple(f"rep{i}" for i in range(1, 5))
MASTER_LOCI = ("Nelfa", "Zfp352", "Lmx1a", "Cdk2ap1", "Snai1", "Pou6f2")
DISPLAY_LOCI = ("Cdk2ap1", "Zfp352", "Lmx1a", "Snai1", "Nelfa", "Pou6f2")

METHODS = (
    {
        "method": "ChimeraTE",
        "analysis": "Primary",
        "summary": "chimerate_overlap_summary.tsv",
        "summary_filter": {"category": "all_TE_initiated"},
        "gene_columns": tuple(f"rep{i}_genes" for i in range(1, 5)),
        "union_column": "union_genes",
        "pairwise": "chimerate_pairwise_overlap.tsv",
        "pairwise_filter": {"category": "all_TE_initiated"},
        "locus": "chimerate_locus_recovery.tsv",
        "locus_filter": {"definition": "MT2_MERVL_family"},
    },
    {
        "method": "FREDY",
        "analysis": "Primary",
        "summary": "fredy_overlap_summary.tsv",
        "summary_filter": {"mapping_scope": "unambiguous"},
        "gene_columns": tuple(f"rep{i}_genes" for i in range(1, 5)),
        "union_column": "union_genes",
        "pairwise": None,
        "pairwise_filter": {},
        "locus": "fredy_locus_recovery.tsv",
        "locus_filter": {
            "mapping_scope": "unambiguous",
            "definition": "MT2_MERVL_family",
        },
    },
    {
        "method": "FREDY-adapted",
        "analysis": "Sensitivity analysis",
        "summary": "fredy_adapted_overlap_summary.tsv",
        "summary_filter": {},
        "gene_columns": tuple(f"rep{i}_genes" for i in range(1, 5)),
        "union_column": "union_genes",
        "pairwise": "fredy_adapted_pairwise_overlap.tsv",
        "pairwise_filter": {},
        "locus": "fredy_adapted_locus_recovery.tsv",
        "locus_filter": {},
    },
    {
        "method": "TEProf2",
        "analysis": "Primary",
        "summary": "teprof2_overlap_summary.tsv",
        "summary_filter": {},
        "gene_columns": tuple(f"rep{i}_supported_genes" for i in range(1, 5)),
        "union_column": "union_supported_genes",
        "pairwise": "teprof2_pairwise_overlap.tsv",
        "pairwise_filter": {},
        "locus": "teprof2_locus_recovery.tsv",
        "locus_filter": {},
    },
    {
        "method": "LIONS",
        "analysis": "Primary",
        "summary": "lions_overlap_summary.tsv",
        "summary_filter": {},
        "gene_columns": tuple(f"rep{i}_genes" for i in range(1, 5)),
        "union_column": "union_genes",
        "pairwise": "lions_pairwise_overlap.tsv",
        "pairwise_filter": {},
        "locus": "lions_locus_recovery.tsv",
        "locus_filter": {},
    },
)


def read_tsv(path):
    if not path.is_file():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep="\t", dtype=str)


def select(frame, filters, label, one=False):
    result = frame
    for column, value in filters.items():
        if column not in result.columns:
            raise ValueError(f"{label}: missing column {column}")
        result = result[result[column].eq(value)]
    if one and len(result) != 1:
        raise ValueError(f"{label}: expected one row, found {len(result)}")
    return result.iloc[0] if one else result


def number(value):
    if pd.isna(value) or str(value).strip() in {"", "NA", "None"}:
        return 0.0
    return float(value)


def mean_sd(values, digits):
    values = tuple(map(float, values))
    return f"{statistics.mean(values):.{digits}f} ± {statistics.stdev(values):.{digits}f}"


def count_percent(value, denominator):
    value, denominator = int(value), int(denominator)
    percent = 100 * value / denominator if denominator else 0
    return f"{value} ({percent:.1f}%)"


def family_tokens(value):
    if pd.isna(value):
        return set()
    return {
        token.strip()
        for token in str(value).replace(",", ";").split(";")
        if token.strip() not in {"", "NA", "None", "nan"}
    }


def matching_families(value):
    return {
        family
        for family in family_tokens(value)
        if family.startswith("MT2") or "MERVL" in family.upper()
    }


def pairwise_jaccards(results, config, summary):
    if config["pairwise"] is None:
        columns = (
            f"jaccard_{left}_{right}"
            for left, right in itertools.combinations(REPS, 2)
        )
        values = [number(summary[column]) for column in columns]
    else:
        path = results / config["pairwise"]
        frame = select(
            read_tsv(path),
            config["pairwise_filter"],
            path.name,
        )
        if "jaccard" not in frame.columns:
            raise ValueError(f"{path.name}: missing jaccard")
        pair_columns = next(
            (
                columns
                for columns in (("rep_a", "rep_b"), ("replicate_1", "replicate_2"))
                if all(column in frame.columns for column in columns)
            ),
            None,
        )
        expected_pairs = {
            frozenset(pair) for pair in itertools.combinations(REPS, 2)
        }
        observed_pairs = [
            frozenset(pair)
            for pair in frame.loc[:, pair_columns].itertuples(index=False, name=None)
        ] if pair_columns else []
        if len(observed_pairs) != 6 or set(observed_pairs) != expected_pairs:
            raise ValueError(f"{path.name}: incomplete or duplicate replicate pairs")
        values = [number(value) for value in frame["jaccard"]]
    if len(values) != 6 or any(
        not math.isfinite(value) or value < 0 or value > 1 for value in values
    ):
        raise ValueError(f"{config['method']}: invalid pairwise Jaccards")
    return values


def locus_index(results, config):
    path = results / config["locus"]
    frame = select(read_tsv(path), config["locus_filter"], path.name)
    if "gene_name" not in frame.columns or frame["gene_name"].duplicated().any():
        raise ValueError(f"{path.name}: invalid gene_name keys")
    return frame.set_index("gene_name")


def locus_value(row):
    if row is None:
        return 0, set()
    per_rep = tuple(f"{rep}_TE_families" for rep in REPS)
    if all(column in row.index for column in per_rep):
        families = set().union(*(matching_families(row[column]) for column in per_rep))
        count = sum(bool(matching_families(row[column])) for column in per_rep)
    else:
        families = matching_families(row.get("TE_families", ""))
        count = int(number(row.get("n_replicates", 0))) if families else 0
    return count, families


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, required=True)
    args = parser.parse_args()

    rows = []
    for config in METHODS:
        summary_path = args.results / config["summary"]
        summary = select(
            read_tsv(summary_path),
            config["summary_filter"],
            summary_path.name,
            one=True,
        )
        gene_counts = [number(summary[column]) for column in config["gene_columns"]]
        union = int(number(summary[config["union_column"]]))
        recurrent = {
            threshold: int(number(summary[column]))
            for threshold, column in (
                (2, "at_least_2_reps"),
                (3, "at_least_3_reps"),
                (4, "all_4_reps"),
            )
        }
        if not union >= recurrent[2] >= recurrent[3] >= recurrent[4] >= 0:
            raise ValueError(f"{config['method']}: invalid recurrence counts")

        row = {
            "Method": config["method"],
            "Analysis": config["analysis"],
            "Genes/replicate (mean ± SD)": mean_sd(gene_counts, 1),
            "Pairwise Jaccard (mean ± SD)": mean_sd(
                pairwise_jaccards(args.results, config, summary), 3
            ),
            "Union genes": union,
            "Genes in ≥2/4, n (%)": count_percent(recurrent[2], union),
            "Genes in ≥3/4, n (%)": count_percent(recurrent[3], union),
            "Genes in 4/4, n (%)": count_percent(recurrent[4], union),
        }
        loci = locus_index(args.results, config)
        missing_loci = set(MASTER_LOCI) - set(loci.index)
        if missing_loci:
            raise ValueError(
                f"{config['locus']}: missing loci {', '.join(sorted(missing_loci))}"
            )
        for gene in MASTER_LOCI:
            count, families = locus_value(loci.loc[gene])
            suffix = f" [{','.join(sorted(families))}]" if families and count else ""
            row[gene] = f"{count}/4{suffix}"
        rows.append(row)

    master = pd.DataFrame(rows)
    master.to_csv(args.results / "method_comparison.tsv", sep="\t", index=False)

    reproducibility_columns = (
        "Method",
        "Genes/replicate (mean ± SD)",
        "Pairwise Jaccard (mean ± SD)",
        "Union genes",
        "Genes in ≥2/4, n (%)",
        "Genes in 4/4, n (%)",
    )
    master.loc[:, reproducibility_columns].to_csv(
        args.results / "reproducibility.tsv", sep="\t", index=False
    )

    known = master.set_index("Method").loc[:, DISPLAY_LOCI]
    known = known.map(lambda value: value.split(" ", 1)[0]).reset_index()
    known.to_csv(args.results / "known_loci.tsv", sep="\t", index=False)

    print(master.to_string(index=False))


if __name__ == "__main__":
    main()
