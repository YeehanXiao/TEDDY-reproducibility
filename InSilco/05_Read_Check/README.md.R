# Read-support quality checks

This folder contains the scripts used to evaluate the read-level quality of the simulated RNA-seq data reported in Supplementary Table 1.

The analyses assess whether simulated TE-chimeric transcript structures are supported by appropriate genome-coordinate alignment evidence, including:
  
  - continuous alignments spanning intra-exonic TE-host boundaries;

- paired-end read linkage between TE-side and host-side anchor regions;

- exact CIGAR-N split-read support for simulated exon-exon TE-host splice junctions.

Scripts:
  
  1. `01_prepare_simulated_breakpoints.R`  

Defines evaluable TE-host breakpoint regions and corresponding TE-side and host-side anchors from the simulated truth annotation.

2. `02_quantify_simulated_support.py`  

Quantifies continuous-alignment and paired-anchor support across sequencing depths. Soft-clipped reads are recorded only as an auxiliary diagnostic.

3. `03_summarize_simulated_support.R`  

Summarizes read-support rates across sequencing depths and transcript abundance groups for Supplementary Table 1.

4. `04_check_CIGAR_N_split_reads.py`  

Independently checks whether expressed simulated TE-host splice junctions are supported by reads containing an exact matching CIGAR `N` operation.

5. `run_read_support_checks.sh`  

Provides the command sequence used to run the analyses.