# Ten-seed structural recovery and contribution stratification

This module reproduces TEDDY scoring across seeds **18 and 101-109** at
**5x, 10x, 25x, 50x and 100x** sequencing depth. It includes the three-panel
Figure 2A and the final contribution-stratified recall figure.

## Run

Requirements: R >= 4.1, `data.table`, `ggplot2`, `patchwork`, and an R build
with Cairo support for PDF export. No genome annotation packages are needed
to replay these retained annotations.

```sh
Rscript -e 'install.packages(c("data.table", "ggplot2", "patchwork"), repos="https://cloud.r-project.org")'
cd InSilco/07_Stratification
Rscript 01_score_stratification.R
Rscript 04_verify.R
Rscript 02_plot_fig2a.R
Rscript 03_plot_stratification.R
```

The scripts resolve default inputs relative to their own location. Generated
tables and PDF/SVG figures go to `results/`, which is ignored by Git.
To use another output directory, supply these positional arguments:

```sh
Rscript 01_score_stratification.R inputs/structural_annotations.rds results
Rscript 04_verify.R results inputs/reference
Rscript 02_plot_fig2a.R results/TE_structural_recovery_10seeds.tsv results
Rscript 03_plot_stratification.R results/stratified_performance.tsv results
```

## Files and outputs

| File | Purpose |
| --- | --- |
| `inputs/structural_annotations.rds` | Retained truth exon-TE annotations, expressed transcript IDs and predicted annotations for all ten seeds and five depths |
| `inputs/reference/stratified_performance.tsv` | Frozen 3,640-row gene/locus result table |
| `inputs/reference/structural_stratified_performance.tsv` | Frozen 5,850-row structural result table |
| `inputs/reference/TE_structural_recovery_10seeds.rds` | Archived 50-row recovery object with all four historical metrics |
| `inputs/reference/03_plot_TE_structural_recovery_10seeds.R` | Archived four-panel plot source; trailing whitespace removed |
| `01_score_stratification.R` | Recompute the two complete tables and `TE_structural_recovery_10seeds.tsv` |
| `02_plot_fig2a.R` | Plot locus, category and context-aware relative-position recovery |
| `03_plot_stratification.R` | Plot locus recall by overlap, exon fraction and TE incorporation fraction |
| `04_verify.R` | Check every frozen table field, Fig2A denominators/rates and count conservation |

## Definitions

Truth is restricted to the retained expressed-transcript set (TPM > 1).
Annotations use same-strand exon overlap with a single annotated TE interval
of at least 5 bp. Original structural labels retain a 5-bp annotation tolerance.

The matching units are assigned host gene; host gene plus exact TE locus
(`chromosome:TE_start:TE_end:TE_strand`); or that locus key plus the original
category, position or boundary label. Repeated exon/isoform associations
collapse only within a complete key. Each contribution measure independently
takes its maximum within that key.

For exon E and annotated TE interval T, overlap is |E intersect T| in bp,
exon fraction is overlap / |E|, and TE incorporation fraction is overlap / |T|.
Bins retain their original left-inclusive boundaries:

| Axis | Bins |
| --- | --- |
| Overlap | <50, 50-<75, 75-<150, >=150 bp |
| Exon fraction | <10%, 10-<25%, 25-<50%, >=50% |
| TE incorporation fraction | <50%, 50-<95%, >=95% |

Recall is stratified using truth attributes and includes missed truth units.
Precision is separately stratified using all prediction attributes, including
false positives. No stratified F1 is calculated. Tables retain overall F1 for
gene/locus scoring, explicit denominator counts, empty-stratum NA values and
`n_nonempty_seeds`. Summaries are equal-weight seed means and sample SD.
`all10` includes all ten seeds; `formal9` contains 101-109; `seed18` retains
the original benchmark seed separately.

## Figure 2A field mapping

| Manuscript panel | Archived field |
| --- | --- |
| TE-chimeric locus | `TE_locus` |
| TE-chimeric category | `TE_category` |
| Relative TE position | `TE_boundary` |

The manuscript's relative-position panel incorporates exon-chain context,
so it uses `TE_boundary`. The historical `TE_position` field describes plain
exon-relative position and remains in the recovery table and archived object.
The three-panel plotting script recreates the manuscript layout from the
archived four-panel source; it is not the lost original layout-edit script.
Grey points are individual seeds; red curves show mean +/- sample SD.

## Provenance and validation

The portable annotation input was extracted from the retained annotation
cache on 2026-11-06. Frozen locus/gene stratification is from 2026-11-02;
structural stratification is from 2026-11-06. The original recovery object
and four-panel plotting script were archived on 2026-09-09.

The compact scorer preserves the original scoring and aggregation rules.
Validated with R 4.5.1, data.table 1.18.4, ggplot2 4.0.3 and patchwork 1.3.2.
Verification compares all keys, fields and NA positions in all 3,640 and
5,850 frozen rows, plus 50 seed/depth rows containing four denominator/rate
pairs each. Numerical tolerance is 1e-12. TP + FN = N_truth and
pred_TP + FP = N_pred are also checked.

This is a replay from stored annotations. It does not run simulation,
alignment, assembly, annotation or callers. Final independent-seed gene-level
scoring sources and retained LIONS calls are in
[`../03_Metrics_Evaluation/nine_seed/`](../03_Metrics_Evaluation/nine_seed/).
The original simulation-generation runner for seed replicates is not included.
