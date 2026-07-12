#!/usr/bin/env Rscript

# ==============================================================================
# Script: 01_run_teddy.R
# Purpose: Execute TEDDY pipeline (Assembly, Annotation, Quantification, Chimeric)
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Paths and Setup (Desensitized)
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Teddy)
  library(parallel)
  library(GenomeInfoDb)
  library(rtracklayer)
  library(dplyr)
})

data_dir <- "./data"
results_dir <- "./results"
bam_dir <- file.path(results_dir, "sortbam")
work_dir <- file.path(results_dir, "TEDDY_official_by_depth")

reference <- file.path(results_dir, "official_simulated_reference_90pct.gtf")
mm10_TE_path <- file.path(data_dir, "mm10_TE.rds")

if (dir.exists(work_dir)) unlink(work_dir, recursive = TRUE, force = TRUE)
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

depths <- c("5x", "10x", "25x", "50x", "100x")
bamfiles_all <- sort(list.files(bam_dir, pattern = "\\.bam$", full.names = TRUE))
bamfiles_all <- bamfiles_all[!grepl("\\.bai$", bamfiles_all)]

stopifnot(file.exists(reference), length(bamfiles_all) > 0)

mm10_TE <- readRDS(mm10_TE_path)
if (exists("NCBI_check")) mm10_TE <- NCBI_check(mm10_TE, ncbi_style = FALSE)

# ------------------------------------------------------------------------------
# 1. Execution Loop
# ------------------------------------------------------------------------------
for (depth in depths) {
  message("========== Running TEDDY depth: ", depth, " ==========")
  
  depth_dir <- file.path(work_dir, depth)
  gtf_dir   <- file.path(depth_dir, "gtf")
  GTF_dir   <- file.path(depth_dir, "GTF")
  meta_dir  <- file.path(depth_dir, "meta")
  
  dir.create(gtf_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(GTF_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(meta_dir, recursive = TRUE, showWarnings = FALSE)
  
  bamfiles <- grep(paste0("_", depth, "_noise0\\.1\\.bam$"), bamfiles_all, value = TRUE)
  stopifnot(length(bamfiles) == 1)
  
  # Step A: Assembly
  assembled_gtf <- file.path(gtf_dir, sub("\\.bam$", ".gtf", basename(bamfiles)))
  Teddy::stringtieAssembly(bam = bamfiles, reference = reference, outfile = assembled_gtf, params = "-p 8")
  
  # Step B: Annotate
  annotated_gtf <- file.path(GTF_dir, paste0("official_simulated_", depth, ".assembled.annotated.gtf"))
  gffcompareAnno(reference = reference, gtffile = assembled_gtf, outfile = annotated_gtf, overwrite = TRUE)
  
  # Step C: TEDDY core algorithms
  anno_compare <- prepareAnno(gtffile = annotated_gtf, transposon = mm10_TE, cores = 40)
  official_se <- countAnno(annotation = anno_compare, bamfile = bamfiles)
  
  official_combineSE <- stringtieCombine(
    reference = annotated_gtf, params = "-p 40", gtfFiles = assembled_gtf,
    bamFiles = bamfiles, longRead = FALSE, cores = 4
  )
  
  chi_GTF <- processGTF(te = mm10_TE, combineSE = official_combineSE, minoverlap = 5, threads = 8)
  
  # Save clean outputs for evaluation downstream
  saveRDS(anno_compare, file = file.path(meta_dir, paste0("anno_compare_official_", depth, ".rds")))
  saveRDS(official_se, file = file.path(meta_dir, paste0("official_se_", depth, ".rds")))
  saveRDS(official_combineSE, file = file.path(meta_dir, paste0("official_combineSE_", depth, ".rds")))
  saveRDS(chi_GTF, file = file.path(meta_dir, paste0("official_chi_GTF_", depth, ".rds")))
}
message("--- TEDDY Execution Completed ---")