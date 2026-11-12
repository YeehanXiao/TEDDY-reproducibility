# Final nine-seed benchmark

Seeds **101-109** were evaluated independently at 5x, 10x, 25x, 50x and
100x, giving 270 seed/depth/tool combinations. Configuration choices were
established using seed18, which is excluded from this nine-seed summary.
The ten-seed TEDDY structural analysis and Fig2A are in
[`../../07_Stratification/`](../../07_Stratification/).

## Replay retained final scoring

From this directory, with R packages `dplyr` and `readr` installed:

```sh
mkdir -p results/metrics
cp inputs/metrics_101_109.rds metrics_101_109.tsv results/metrics/
Rscript code/update_lions_final_metrics.R results inputs/LIONS_final_calls
Rscript code/export_sheet3_final.R results results
Rscript code/check_metrics_final.R results
```

The compact inputs contain the retained metrics/gene sets and 45 final LIONS
`.lion` files. This replay recalculates LIONS host-gene metrics from final
calls, preserves the other tools' stored metrics, and exports the two Sheet3
tables. It does not regenerate simulations or rerun callers. Existing metrics
are archived before LIONS updates; generated files stay in ignored `results/`.

The original base scorer reads pre-chimSort LIONS log candidates. Its LIONS
rows must be replaced with `update_lions_final_metrics.R` before final export.
Final LIONS uses the official `chimSort.R` classifier with built-in
`simuOptimal` parameters `1 2 10 2 0.1 2 1.5`. The other final settings are
TEDDY original; Arriba final+discarded with support >=1; ChimeraTE overlap
0.10 with both stranded modes; FREDY all-transcript TE50 before CDS filters;
and TEProf2 with the simulation dictionary adapter. Truth scopes remain in
`metrics_101_109.tsv`.

## Score the original result tree

`SEED_ROOT` is an existing independent-seed result directory. Each
`seed_101` through `seed_109` must contain:

```text
data/mm10_TE.rds
results/
  official_simulated_transcript_truth_status.tsv
  official_simulated_exon_truth_status.tsv
  official_simulated_1000.isoforms.results
  official_simulated_truth_1000_transcripts.gtf
  TEDDY_official_by_depth/DEPTH/meta/official_chi_GTF_DEPTH.rds
  arriba_merge/merge_DEPTH_arriba_fusions{,.discarded}.tsv
  ChimeraTE_official_by_depth/DEPTH_{fwd,rf}-stranded/rep1_DEPTH/TE-*-rep1_DEPTH.tsv
  LIONS_official_by_depth/DEPTH/official_simulated_DEPTH_noise0.1.LIONS.log
  LIONS_simuOptimal_by_depth/DEPTH/official_simulated_DEPTH_noise0.1.lion
  TEProf2_benchmark/DEPTH_test/assembly/*_annotated_filtered_test_all
  FREDY_relaxed_by_depth/DEPTH/chimeric/protein.gtf
```

With those existing files and R packages `rtracklayer`, `GenomicRanges`,
`dplyr`, `purrr` and `readr`, run:

```sh
Rscript code/score_seed_benchmarks_base.R "$SEED_ROOT"
Rscript code/update_lions_final_metrics.R "$SEED_ROOT"
Rscript code/export_sheet3_final.R "$SEED_ROOT" results
Rscript code/check_metrics_final.R "$SEED_ROOT"
```

This writes `SEED_ROOT/metrics/metrics_101_109.{tsv,rds}`. The base scorer
overwrites those metrics files; use a separate result-tree copy to preserve
a previous scoring run. Caller outputs are read only.

`run_lions_simuoptimal_final.sh ROOT LIONS_DIR PYTHON RSCRIPT` optionally
rebuilds the LCSV files and reruns final LIONS classification from existing
logs and indexed BAMs at `seed_N/results/sortbam/`. It requires Python >=3.9
with `pysam` and the original LIONS installation. It archives an existing
`LIONS_simuOptimal_by_depth` directory before rebuilding. The helper
`audit_lions_final.py ROOT` checks interaction counts, including X/Y.

## Validation and provenance

The final updater was tested after replacing all 45 LIONS metric rows and
prediction sets with incorrect empty results. Updating from retained calls
restored all 270 metric rows and both frozen Sheet3 tables (270 raw rows,
30 summary rows), with all three TSV files byte-identical to the frozen
references. Count identities and final
`simuOptimal` settings also passed. Raw caller execution was not rerun.

Source: the retained `benchmark_9seed_final` bundle dated 2026-09-12.
Baseline comment edits preserve its R AST.
The updater adds an optional compact-call directory; the shell runner now
requires explicit machine paths. Other numerical logic is unchanged.
Original-source SHA256 values:

| Original file | SHA256 |
| --- | --- |
| `score_seed_benchmarks_base.R` | `29d3aabd80ea88d37abe0753a69f74f5e3fae7822322d4e0de8a789ca4eae6a7` |
| `update_lions_final_metrics.R` | `1a1f1c42a765d672248e15a04cc57c4d11c4091a5d10f84e55bb83fa9c33ee86` |
| `export_sheet3_final.R` | `f5ba7822f71cab90b07a5a73c07daf3a16f62903cd85540f880536272c81e312` |
| `check_metrics_final.R` | `62aaec4d98d3929fd33acba475cbe3bfe6ebddd464ffb0554dfea6f4ead1a781` |
| `run_lions_simuoptimal_final.sh` | `e0f6168ea2294a3bbf1bd4849981550bc85956be8385a7d865b0db9bd46a1068` |
| `build_lions_lcsv.py` | `2dbb3366e8ef32603f07a4b8c5f0efcc78cf491e1c8984b7fbde62d6d91c71f8` |
| `audit_lions_final.py` | `adedb50e9b340581e1d0872413b8807601cefb809522127663449393a974f140` |
