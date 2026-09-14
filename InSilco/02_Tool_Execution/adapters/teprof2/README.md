# TEProf2 simulation adapters

SHA256 comparison across seeds 101 through 109 found one common version of
each recovered script:

| File | Recovered SHA256 |
| --- | --- |
| `prepare_teprof2_dictionary_input.py` | `dc8ae4c4630e8afebc8272b435494584d025a37b177e302f4b109c4208bef7b2` |
| `normalize_teprof2_dictionaries.py` | `e18252b60c66f912400c782e2af8ed4ab2570fd71d27f7bc21c54ec29c90b673` |
| `genecode_to_dic_simulation.py` | `0d9c31bbc264fc8c3c3c79143b71662239bea93583eb87b7c502cb2579e87c05` |
| `run_teprof2_all_depths.sh` | `a5c49b95219bdb3dbc1ff09d1d253595d879eb5af2e576394b986caed4094f45` |

The input-preparation and normalization adapters are byte-identical to the
recovered files. The local dictionary builder fixes the recovered builder's
`stop_codon` key removal typo; its SHA256 is
`024ea67e3346351a7678a127ac40188bf9795f33f225cb28ca6901fd2c5bbd48`.
This does not affect the validated transcript/exon-only simulation input.

The server-specific orchestration from the recovered shell script was
integrated into `../../06_run_teprof2.sh`; the local shell file is a portable
entrypoint with SHA256
`267d47ba74af9b2d23879206aec7fbfc76dd6d397dad2f76a3b50ac34e36b414`.
