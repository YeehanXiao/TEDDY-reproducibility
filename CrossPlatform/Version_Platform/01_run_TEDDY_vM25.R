suppressPackageStartupMessages({
  library(Teddy)
  library(parallel)
})

## =========================
## 1. Paths
## =========================

base_dir <- "path/to/project"
work_dir <- file.path(base_dir, "vM25_2Clike")

dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(work_dir, "gtf"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(work_dir, "GTF"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(work_dir, "meta"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(work_dir, "logs"), recursive = TRUE, showWarnings = FALSE)

reference <- "path/to/gencode.vM25.annotation.gtf"
te_rds <- "path/to/mm10_TE.rds"

es_bam_dir <- "path/to/ESC_2Clike_BAMs"
embryo_bam_dir <- "path/to/preimplantation_embryo_BAMs"

## =========================
## 2. Sample information
## =========================

sample_tbl <- data.frame(
  sample = c(
    paste0("bioneg_", 1:3),
    paste0("biopos_", 1:3),
    paste0("pos_", 1:3),
    paste0("2cellrep", 1:4),
    paste0("4cellrep", 1:4),
    paste0("8cellrep", 1:3),
    paste0("ICMrep", 1:4),
    paste0("TErep", 1:4)
  ),
  group = c(
    rep("ESCs_ZM", 3),
    rep("intermediate_ZpMneg", 3),
    rep("2Clike_ZpMp", 3),
    rep("2cell", 4),
    rep("4cell", 4),
    rep("8cell", 3),
    rep("ICM", 4),
    rep("TE", 4)
  ),
  bam_dir = c(
    rep(es_bam_dir, 9),
    rep(embryo_bam_dir, 19)
  ),
  stringsAsFactors = FALSE
)

sample_tbl$bam <- file.path(
  sample_tbl$bam_dir,
  paste0(sample_tbl$sample, ".bam")
)

sample_tbl$gtf <- file.path(
  work_dir,
  "gtf",
  paste0(sample_tbl$sample, ".gtf")
)

stopifnot(file.exists(reference))
stopifnot(file.exists(te_rds))

missing_bam <- sample_tbl$bam[!file.exists(sample_tbl$bam)]

if (length(missing_bam) > 0) {
  stop(
    "Missing BAM files:\n",
    paste(missing_bam, collapse = "\n")
  )
}

## =========================
## 3. StringTie assembly
## =========================

message("Step 1: StringTie assembly per BAM")

parallel::mclapply(
  seq_len(nrow(sample_tbl)),
  FUN = function(i) {
    Teddy::stringtieAssembly(
      bam = sample_tbl$bam[i],
      reference = reference,
      outfile = sample_tbl$gtf[i],
      params = "-p 8"
    )
  },
  mc.cores = 12
)

gtf_files <- sample_tbl$gtf

missing_gtf <- gtf_files[!file.exists(gtf_files)]

if (length(missing_gtf) > 0) {
  stop(
    "Missing assembled GTF files:\n",
    paste(missing_gtf, collapse = "\n")
  )
}

## =========================
## 4. Merge and annotate GTF
## =========================

message("Step 2: StringTie merge")

merged_gtf <- file.path(
  work_dir,
  "GTF",
  "vM25_2Clike.merged.gtf"
)

annotated_gtf <- file.path(
  work_dir,
  "GTF",
  "vM25_2Clike.annotated.gtf"
)

Teddy::stringtieMerge(
  reference = reference,
  gtfFiles = gtf_files,
  outfile = merged_gtf,
  params = "-p 16"
)

message("Step 3: gffcompare annotation")

Teddy::gffcompareAnno(
  reference = reference,
  gtffile = merged_gtf,
  outfile = annotated_gtf,
  overwrite = TRUE
)

## =========================
## 5. TEDDY analysis
## =========================

message("Step 4: Import TE annotation")

te_gr <- readRDS(te_rds)

message("Step 5: prepareAnno")

anno_compare_vM25 <- Teddy::prepareAnno(
  gtffile = annotated_gtf,
  transposon = te_gr,
  cores = 40
)

message("Step 6: countAnno")

se_vM25 <- Teddy::countAnno(
  annotation = anno_compare_vM25,
  bamfile = sample_tbl$bam
)

colnames(se_vM25) <- sample_tbl$sample

message("Step 7: stringtieCombine")

combineSE_vM25 <- Teddy::stringtieCombine(
  reference = annotated_gtf,
  params = "-p 70",
  gtfFiles = gtf_files,
  bamFiles = sample_tbl$bam,
  longRead = FALSE,
  cores = 68
)

colnames(combineSE_vM25) <- sample_tbl$sample

message("Step 8: processGTF")

chi_GTF_vM25 <- Teddy::processGTF(
  te = te_gr,
  combineSE = combineSE_vM25,
  threads = 8,
  minoverlap = 5
)

## =========================
## 6. Save outputs
## =========================

saveRDS(
  anno_compare_vM25,
  file.path(work_dir, "meta", "anno_compare.vM25_2Clike.rds")
)

saveRDS(
  se_vM25,
  file.path(work_dir, "meta", "se.vM25_2Clike.rds")
)

saveRDS(
  combineSE_vM25,
  file.path(work_dir, "meta", "combineSE.vM25_2Clike.rds")
)

saveRDS(
  chi_GTF_vM25,
  file.path(work_dir, "meta", "chi_GTF.vM25_2Clike.rds")
)

write.table(
  sample_tbl[, c("sample", "group")],
  file = file.path(work_dir, "meta", "sample_info.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Completed.")
