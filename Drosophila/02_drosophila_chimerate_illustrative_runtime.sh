PROJECT_DIR="/path/to/benchmark_project"
CHIMERATE_DIR="/path/to/ChimeraTE"

GENOME_FA="${CHIMERATE_DIR}/example_data/mode1/dmel_genome_sample.fa"
INPUT_TSV="${CHIMERATE_DIR}/example_data/mode1/input_mode1.tsv"
TE_GTF="${CHIMERATE_DIR}/example_data/mode1/dmel_TEs_sample.gtf"
GENE_GTF="${CHIMERATE_DIR}/example_data/mode1/dmel_genes_sample.gtf"

OUTDIR="${PROJECT_DIR}/runtime_drosophila_chimerate"
LOGDIR="${OUTDIR}/logs"
STRAND="rf-stranded"

mkdir -p "${OUTDIR}" "${LOGDIR}"

cd "${CHIMERATE_DIR}"

for rep in 1 2 3; do
PROJECT_NAME="example_mode1_runtime_rep${rep}"

rm -rf "projects/${PROJECT_NAME}"

/usr/bin/time -v python3 chimTE_mode1.py \
--genome "${GENOME_FA}" \
--input "${INPUT_TSV}" \
--project "${PROJECT_NAME}" \
--te "${TE_GTF}" \
--gene "${GENE_GTF}" \
--strand "${STRAND}" \
> "${LOGDIR}/ChimeraTE_drosophila_rep${rep}.stdout.log" \
2> "${LOGDIR}/ChimeraTE_drosophila_rep${rep}.time.log"

cp -r "projects/${PROJECT_NAME}" "${OUTDIR}/"
done