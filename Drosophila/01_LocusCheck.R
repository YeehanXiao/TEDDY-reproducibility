#!/usr/bin/env Rscript

# ==============================================================================
# Script: 01_LocusCheck.R
# Purpose: Benchmark TEDDY against ChimeraTE using Drosophila example data.
#          Includes TEDDY processing, comparative evaluation, and Gviz visualization.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. User Configuration (Modify these to your local paths)
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Rstringtie)
  library(GenomicRanges)
  library(rtracklayer)
  library(GenomicFeatures)
  library(Gviz)
  library(eulerr)
  library(dplyr)
})

base_dir <- "/mnt/datadisk/xiaoyihan/manuals/compare/ChimeraTE/example_data/mode1"
benchmark_dir <- "/mnt/datadisk/xiaoyihan/manuals/compare/ChimeraTE"

# Set output directories
gtf_out_dir <- "../gtf" # Relative to script execution dir
asm_out_dir <- file.path(gtf_out_dir, "assembly")
dir.create(asm_out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. Assembly with Loose Parameters (Fair Comparison)
# ------------------------------------------------------------------------------
start_time <- Sys.time()
bamfiles <- list.files(base_dir, pattern = "*.bam$", full.names = TRUE)
reference_gtf <- file.path(base_dir, "test.gtf")
params <- "-p 8 -f 0.02 -c 0.5 -j 1 -g 20"

# Parallel Assembly (Write outputs to assembly directory, not the raw data dir)
mclapply(bamfiles, FUN = function(x) {
  outfile <- file.path(asm_out_dir, paste0(gsub("\\.bam$", ".loose.gtf", basename(x))))
  Teddy::stringtieAssembly(bam = x, reference = reference_gtf, outfile = outfile, params = params)
}, mc.cores = 8)

# Merge
gtffiles <- list.files(asm_out_dir, pattern = "loose.gtf", full.names = TRUE)
merged_gtf <- file.path(gtf_out_dir, "fly.gtf")
stringtieMerge(reference = reference_gtf, gtfFiles = gtffiles, outfile = merged_gtf, params = "-p 8")

# ------------------------------------------------------------------------------
# 2. TE Processing and Annotation
# ------------------------------------------------------------------------------
te_gtf <- read.delim(file.path(base_dir, "dmel_TEs_sample.gtf"), header = FALSE, sep = "\t")
te_G <- GRanges(seqnames = te_gtf[,1], IRanges(start = te_gtf[,4], end = te_gtf[,5]), 
                strand = te_gtf[,7], names = te_gtf[,9], class = te_gtf[,9])

annotated_gtf <- file.path(gtf_out_dir, "fly.annotated.gtf")
gffcompareAnno(reference = reference_gtf, gtffile = merged_gtf, outfile = annotated_gtf)

anno <- prepareAnno(gtffile = annotated_gtf, transposon = te_G)
se <- countAnno(annotation = anno, bamfile = bamfiles)
combineSE <- stringtieCombine(reference = annotated_gtf, bamFiles = bamfiles, 
                              params = "-p 70", gtfFiles = gtffiles)

# ------------------------------------------------------------------------------
# 3. Benchmark Evaluation
# ------------------------------------------------------------------------------
# Load ChimeraTE baseline results
expected_data <- readRDS(file.path(benchmark_dir, "intermediate/ChimeraTE_results.rds"))
chimera_genes <- unique(c(expected_data$init_gr$gene$gene_id, expected_data$exon_gr$gene$gene_id, expected_data$term_gr$gene$gene_id))

CHI_GTF <- processGTF(te = te_G, combineSE = combineSE, minoverlap = 0)

# ------------------------------------------------------------------------------
# 4. Visualization (Gviz Tracks)
# ------------------------------------------------------------------------------
# FIX: Import the GTF generated in Step 2 so plotting doesn't fail
GTF <- rtracklayer::import(annotated_gtf)

# ---- Target Gene 1: FBgn0031188 ----
target_gene1 <- "FBgn0031188"
pad <- 5000
gene_exons1 <- GTF[mcols(GTF)$gene_name == target_gene1]
gene_span1  <- range(gene_exons1)
plot_from1  <- start(gene_span1) - pad
plot_to1    <- end(gene_span1) + pad
plot_chr1   <- as.character(seqnames(gene_span1))

# Tracks for Gene 1
genomeAxis <- GenomeAxisTrack(name = "axis", col = "black")
bw_rev1 <- import.bw(file.path(benchmark_dir, "Benchmark/fly/CR1_reverse.bw"), which = GRanges(plot_chr1, IRanges(plot_from1, plot_to1)))
data_dt1 <- DataTrack(range = bw_rev1, genome = "", chromosome = plot_chr1, name = "CR_reverse", type = "polygon", col = "black", fill = "black")

Gene_track1 <- GeneRegionTrack(GTF[GTF$transcript_id %in% c("MSTRG.82.1")], transcript = "MSTRG.82.1", fill = "darkblue", col = NA)
Gene_track_2_1 <- GeneRegionTrack(GTF[GTF$transcript_id %in% c("FBtr0335047")], transcript = "FBtr0335047", fill = "darkblue", col = NA)

te_annotationTrack <- AnnotationTrack(te_G[grep("S2", te_G$names)], feature = "TE", name = "TE", fill = "darkblue", col = NA)

plotTracks(list(genomeAxis, te_annotationTrack, Gene_track1, Gene_track_2_1, data_dt1), 
           from = plot_from1, to = plot_to1, chromosome = plot_chr1, scale = 0.2, geneSymbols = TRUE)

# ---- Target Gene 2: FBgn0262731 ----
target_gene2 <- "FBgn0262731"
gene_exons2 <- GTF[mcols(GTF)$gene_name == target_gene2]
gene_span2  <- range(gene_exons2)
plot_from2  <- start(gene_span2) - pad
plot_to2    <- end(gene_span2) + pad
plot_chr2   <- as.character(seqnames(gene_span2))

bw_rev2 <- import.bw(file.path(benchmark_dir, "Benchmark/fly/CR2_reverse.bw"), which = GRanges(plot_chr2, IRanges(plot_from2, plot_to2)))
data_dt2 <- DataTrack(range = bw_rev2, genome = "", chromosome = plot_chr2, name = "CR2_reverse", type = "polygon", col = "black", fill = "black")

Gene_track2 <- GeneRegionTrack(GTF[GTF$transcript_id %in% c("FBtr0089254")], transcript = "FBtr0089254", fill = "darkblue", col = NA)
FB_highlight <- HighlightTrack(trackList = list(data_dt2), range = te_G[grep("FB", te_G$names)][2], col = NA)

plotTracks(list(genomeAxis, te_annotationTrack, Gene_track2, FB_highlight), 
           from = plot_from2, to = plot_to2, chromosome = plot_chr2, scale = 0.2, geneSymbols = TRUE)

# ------------------------------------------------------------------------------
# 5. Venn Diagram & Time Summary
# ------------------------------------------------------------------------------
Teddy_genes <- unique(CHI_GTF$gene_name)
sets <- list(TEDDY = Teddy_genes, ChimeraTE = chimera_genes)
plot(euler(sets), fills = list(fill = c("#EBC9C7", "#918579"), alpha = 0.55), quantities = TRUE)

elapsed_time <- Sys.time() - start_time
print(paste("Total execution time: ", elapsed_time, " seconds"))