suppressPackageStartupMessages({
  library(Teddy)
  library(parallel)
  library(GenomicRanges)
  library(GenomeInfoDb)
})

###############################################################################
# Drosophila illustrative runtime example: TEDDY
#
# Runtime can be measured externally with:
#   /usr/bin/time -v Rscript 03_drosophila_teddy_illustrative_runtime.R
#
# Replace PROJECT_DIR according to the local installation.
###############################################################################

PROJECT_DIR <- "/path/to/benchmark_project"

chim_dir <- file.path(PROJECT_DIR, "ChimeraTE/example_data/mode1")
work_dir <- file.path(PROJECT_DIR, "runtime_teddy_drosophila")

dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(work_dir, "gtf"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(work_dir, "GTF"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(work_dir, "meta"), recursive = TRUE, showWarnings = FALSE)

bamfiles <- list.files(
  chim_dir,
  pattern = "\\.bam$",
  full.names = TRUE
)

reference <- file.path(chim_dir, "test.gtf")
te_file <- file.path(chim_dir, "dmel_TEs_sample.gtf")

stopifnot(length(bamfiles) > 0)
stopifnot(file.exists(reference))
stopifnot(file.exists(te_file))

mclapply(
  bamfiles,
  FUN = function(x, reference, outdir) {
    outfile <- file.path(outdir, sub("\\.bam$", ".gtf", basename(x)))
    Teddy::stringtieAssembly(
      bam = x,
      reference = reference,
      outfile = outfile,
      params = "-p 1"
    )
  },
  reference = reference,
  outdir = file.path(work_dir, "gtf"),
  mc.cores = 1
)

gtffiles <- list.files(
  file.path(work_dir, "gtf"),
  pattern = "\\.gtf$",
  full.names = TRUE
)

stopifnot(length(gtffiles) > 0)

merged_gtf <- file.path(work_dir, "GTF", "fly.gtf")

Teddy::stringtieMerge(
  reference = reference,
  gtfFiles = gtffiles,
  outfile = merged_gtf,
  params = "-p 1"
)

anno_gtf <- file.path(work_dir, "GTF", "fly.annotated.gtf")

Teddy::gffcompareAnno(
  reference = reference,
  gtffile = merged_gtf,
  outfile = anno_gtf,
  overwrite = TRUE
)

te <- read.delim(te_file, header = FALSE, sep = "\t")

te_G <- GRanges(
  seqnames = te[, 1],
  ranges = IRanges(start = te[, 4], end = te[, 5]),
  strand = te[, 7],
  names = te[, 9],
  family = te[, 9],
  class = te[, 9]
)

anno <- Teddy::prepareAnno(
  gtffile = anno_gtf,
  transposon = te_G,
  minoverlap = 0,
  cores = 1
)

se <- Teddy::countAnno(
  annotation = anno,
  bamfile = bamfiles,
  nthreads = 1
)

combineSE <- Teddy::stringtieCombine(
  reference = anno_gtf,
  params = "-p 1",
  gtfFiles = gtffiles,
  bamFiles = bamfiles,
  longRead = FALSE,
  cores = 1
)

CHI_GTF <- Teddy::processGTF(
  te = te_G,
  combineSE = combineSE,
  minoverlap = 0,
  threads = 1
)

saveRDS(anno, file = file.path(work_dir, "meta", "anno.rds"))
saveRDS(se, file = file.path(work_dir, "meta", "se.rds"))
saveRDS(combineSE, file = file.path(work_dir, "meta", "combineSE.rds"))
saveRDS(CHI_GTF, file = file.path(work_dir, "meta", "CHI_GTF.rds"))

summary <- data.frame(
  dataset = "Drosophila illustrative example",
  tool = "TEDDY",
  n_bam = length(bamfiles),
  n_gtf = length(gtffiles),
  n_TE_chimeric_entries = length(CHI_GTF)
)

write.table(
  summary,
  file.path(work_dir, "drosophila_teddy_illustrative_runtime_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

print(summary)