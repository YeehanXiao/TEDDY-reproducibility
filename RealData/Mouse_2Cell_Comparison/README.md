# Mouse 2-cell comparator-method analysis

This directory contains the compact reproducibility code and summary tables for
external TE-initiated transcript callers evaluated against the study loci.
TEDDY is not included in the cross-method table.

## Data

- GEO accession: GSE97778
- Samples: four mouse 2-cell embryo biological replicates
- Reference genome: mm10
- Gene annotation: GENCODE vM25

Raw FASTQ, BAM, bigWig, reference, log, and tool-intermediate files are not
included. All input and software paths are supplied at runtime.

## Methods

- ChimeraTE v1.2: mode 1, `rf-stranded`.
- FREDY v1.0.0: `Novel Initial` calls with unambiguous gene mapping.
- FREDY-adapted sensitivity analysis: the frozen all-transcript TE50 pre-CDS
  caller, followed by a same-strand overlap of at least 5 bp between the TE and
  the transcript first exon. The frozen caller is shared with
  `InSilco/02_Tool_Execution/adapters/fredy/`.
- TEProf2: final Step 6 candidates; replicate support requires at least 10 TE
  reads and at least one TE-to-gene start read in the same replicate.
- LIONS: existing forward/reverse bigWigs supply the official real-data
  coverage fields before `chimSort.R` filtering with base `scREADS=3`,
  `scTHREAD=10`, `scDownThread=10`, `scRPKM=1`, `scCONTR=0.10`, `scUPCOV=2`,
  and `scUPEXON=1.5`. Effective `scREADS` is
  `max(3, round(mapped_reads / 20,000,000))` and equals 3 in every replicate.
  The BAM-based LCSV adapter under `InSilco/` is specific to the simulated
  benchmark, which does not use these real-data bigWigs.

FREDY-adapted is a sensitivity analysis, not an independent method.

Within-method recurrence uses the unambiguous gene identifier available from
each workflow: stable GENCODE gene IDs for ChimeraTE and FREDY, and gene symbols
when the native result lacks a stable gene ID. No gene sets are merged across
methods.

## Outputs

- `results/method_comparison.tsv`: master comparator-method audit table.
- `results/reproducibility.tsv`: compact within-method reproducibility table;
  the full table also reports the number of genes recovered in at least three
  replicates.
- `results/known_loci.tsv`: recovery of six established TE-initiated loci.
- Tool-specific summary, pairwise-overlap, locus-recovery, and QC tables are
  retained as the inputs to the final tables.

The six-locus panel contains four loci characterized in this study (`Nelfa`,
`Lmx1a`, `Snai1`, and `Pou6f2`) and two established controls (`Zfp352` and
`Cdk2ap1`). Tool-specific diagnostic tables also retain `Fam172a`, which is not
part of the cross-method panel. Locus recovery requires an MT2/MERVL-compatible
call and is reported as the number of supported replicates out of four.
`Pou6f2` was resolved as GENCODE vM25 gene `ENSMUSG00000009734`; neither
ChimeraTE nor unambiguously mapped FREDY calls recovered an MT2/MERVL-compatible
event for this gene in any replicate.

Rebuild the three final tables after regenerating the tool-specific summaries:

```bash
python3 code/build_comparison_tables.py --results results
```

The scripts require Python 3 with pandas. Tool-specific preprocessing also uses
the dependencies required by the original tools, including R, samtools,
bedtools, pyBigWig, and the legacy TEProf2 Python 2 environment.
