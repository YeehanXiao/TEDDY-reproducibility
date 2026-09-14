# Sensitivity analysis

`run_chimerate_sensitivity_25x.sh` reproduces the 2 x 2 MAPQ and read-TE
overlap comparison for both supported strand modes. It writes run products and
a deterministic row-count summary under `results/ChimeraTE_sensitivity_25x`.
`CHIMERATE_DIR` must point to the validated ChimeraTE 1.2 sensitivity fork;
the runner verifies all seven modified Mode 1 files by SHA256 before execution.

`ChimeraTE_seed18_25x_chat_counts.tsv` is the retained seed-18 count record used
to select overlap 0.10. Its verification column preserves the original
provenance status instead of claiming independent revalidation.

`FREDY_filter_sensitivity_25x.tsv` is the retained gene-level comparison used
to select the frozen all-transcript TE50 pre-CDS caller. Its SHA256 is
`360fca2c5fcfc03b96db69e4393d190a1fcbb66ec9745da45f73c0dd6c450fb3`.
