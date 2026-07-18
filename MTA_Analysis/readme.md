1. `01_prepare_strictly_supported_loci.R`

- Identifies genomic coordinates for target TEs and their neighboring host transcript exons.

- Outputs a coordinate support table (`MTA_locus_neighbor_support_table.tsv`) used by downstream validation scripts.



2. `02_quantify_transcript_level_chimeric_rpm.py`

- Quantifies the macroscopic intensity of chimeric junctions (Transcript-level RPM).

- Counts reads where one end anchors within the TE (>= 15 bp) and the mate extends outside the annotated TE boundary.



3. `03_quantify_strictly_supported_loci.py`

- Implements an ultra-strict dual-anchor physical constraint via query-name intersection.

- Verifies strictly supported loci by ensuring reads or their pairs anchor securely on *both* the TE side (>= 15 bp) and the host gene side (>= 15 bp). 