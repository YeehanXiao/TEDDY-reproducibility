# TEDDY reproducibility materials

This repository contains benchmark scripts used in the TEDDY manuscript.

Corresponding TEDDY version: v1.2.8  
Software DOI: https://doi.org/10.5281/zenodo.21302347

## Contents

- `InSilco/`: simulated-data accuracy, sensitivity, and runtime benchmarks.
- [`InSilco/03_Metrics_Evaluation/nine_seed/`](InSilco/03_Metrics_Evaluation/nine_seed/):
  final independent-seed scoring sources, including the LIONS correction and
  retained inputs for replaying the final metrics.
- [`InSilco/07_Stratification/`](InSilco/07_Stratification/): ten-seed TEDDY
  contribution stratification and the three-panel Figure 2A, with portable
  annotations and verification against archived results.
- `RealData/Mouse_2Cell_Comparison/`: four-replicate mouse 2-cell comparison
  of external TE-initiated-transcript callers.
- `CrossPlatform/`: annotation-version, platform, and cell-line analyses.
- `Drosophila/`: illustrative Drosophila analysis.
- `MTA_Analysis/`: locus-level read-support analysis.
