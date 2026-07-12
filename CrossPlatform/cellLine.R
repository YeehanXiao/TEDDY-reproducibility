library(Teddy)
library(parallel)
library(GenomeInfoDb)
library(readr)
library(dplyr)
library(stringr)
library(tidyr)

# ------------------------------------------------------------------------------
# 1. Illumina platform: assembly, merge, annotation, counting, combine
# ------------------------------------------------------------------------------

reference <- "input/reference/Homo_sapiens.GRCh38.91.gtf"

bamfiles <- list.files(
  "input/illumina/mapping/sortbam",
  pattern = "*.bam$",
  full.names = TRUE
)

bamfiles_HepG2 <- grep(x = bamfiles, pattern = "Hep", value = TRUE)
bamfiles_K562  <- grep(x = bamfiles, pattern = "K562", value = TRUE)

setwd("work/illumina")

mclapply(
  bamfiles,
  FUN = function(x, reference, outfile) {
    bam <- basename(x)
    outfile <- gsub(".bam$", ".gtf", bam)
    
    Teddy::stringtieAssembly(
      bam = x,
      reference = reference,
      outfile = outfile,
      params = "-p 8"
    )
  },
  reference = reference,
  mc.cores = 20
)

gtffiles_HepG2 <- list.files(".", full.names = TRUE, pattern = "Hep*")
gtffiles_K562  <- list.files(".", full.names = TRUE, pattern = "K5*")

stringtieMerge(
  reference = reference,
  gtfFiles = gtffiles_HepG2,
  outfile = "GTF/HepG2.gtf",
  params = "-p 8"
)

stringtieMerge(
  reference = reference,
  gtfFiles = gtffiles_K562,
  outfile = "GTF/K562.gtf",
  params = "-p 8"
)

gffcompareAnno(
  reference = reference,
  gtffile = "GTF/HepG2.gtf",
  outfile = "GTF/gffcompare_HepG2.gtf"
)

gffcompareAnno(
  reference = reference,
  gtffile = "GTF/K562.gtf",
  outfile = "GTF/gffcompare_K562.gtf"
)

te_G_H <- readRDS("input/meta/te_G_H.rds")

seqlevelsStyle(te_G_H) <- "NCBI"
seqlevels(te_G_H) <- sub("^chr", "", seqlevels(te_G_H))
seqnames(te_G_H) <- sub("^chr", "", seqnames(te_G_H))

anno_compare_HepG2 <- prepareAnno(
  gtffile = "GTF/gffcompare_HepG2.annotated.gtf",
  transposon = te_G_H,
  cores = 80
)

anno_compare_K562 <- prepareAnno(
  gtffile = "GTF/gffcompare_K562.annotated.gtf",
  transposon = te_G_H,
  cores = 80
)

HepG2_se <- countAnno(
  annotation = anno_compare_HepG2,
  bamfile = bamfiles_HepG2
)

K562_se <- countAnno(
  annotation = anno_compare_K562,
  bamfile = bamfiles_K562
)

HepG2_combineSE <- stringtieCombine(
  reference = "GTF/gffcompare_HepG2.annotated.gtf",
  params = "-p 70",
  gtfFiles = gtffiles_HepG2,
  bamFiles = bamfiles_HepG2,
  longRead = FALSE,
  cores = 4
)

K562_combineSE <- stringtieCombine(
  reference = "GTF/gffcompare_K562.annotated.gtf",
  params = "-p 70",
  gtfFiles = gtffiles_K562,
  bamFiles = bamfiles_K562,
  longRead = FALSE,
  cores = 4
)

saveRDS(anno_compare_HepG2, file = "meta/anno_compare_HepG2.rds")
saveRDS(anno_compare_K562, file = "meta/anno_compare_K562.rds")
saveRDS(HepG2_se, file = "meta/HepG2_se.rds")
saveRDS(K562_se, file = "meta/K562_se.rds")
saveRDS(HepG2_combineSE, file = "meta/HepG2_combineSE.rds")
saveRDS(K562_combineSE, file = "meta/K562_combineSE.rds")

# ------------------------------------------------------------------------------
# 2. PacBio platform: assembly, annotation, counting, combine
# ------------------------------------------------------------------------------

setwd("work/pacbio")

bamfiles <- list.files("rmdup_bam", pattern = "*.bam$", full.names = TRUE)

mclapply(
  bamfiles,
  FUN = function(x, reference, outfile) {
    bam <- basename(x)
    outfile <- gsub(".bam$", ".gtf", bam)
    
    Teddy::stringtieAssembly(
      bam = x,
      reference = reference,
      outfile = outfile,
      params = "-p 20",
      longRead = TRUE
    )
  },
  reference = reference,
  mc.cores = 2
)

gffcompareAnno(
  reference = reference,
  gtffile = "GTF/HepG2_PacBio.gtf",
  outfile = "GTF/gff_HepG2_PacBio.gtf"
)

gffcompareAnno(
  reference = reference,
  gtffile = "GTF/K562_PacBio.gtf",
  outfile = "GTF/gff_K562_PacBio.gtf"
)

te_G_H <- readRDS("input/meta/te_G_H.rds")
seqlevelsStyle(te_G_H) <- "NCBI"

anno_compare_HepG2_pacbio <- prepareAnno(
  gtffile = "GTF/gff_HepG2_PacBio.annotated.gtf",
  transposon = te_G_H,
  cores = 40
)

anno_compare_K562_pacbio <- prepareAnno(
  gtffile = "GTF/gff_K562_PacBio.annotated.gtf",
  transposon = te_G_H,
  cores = 40
)

HepG2_se_pacbio <- countAnno(
  annotation = anno_compare_HepG2_pacbio,
  bamfile = bamfiles[1],
  isLongRead = TRUE,
  strandSpecific = 0,
  isPairedEnd = FALSE,
  nthreads = 8,
  annot.inbuilt = NULL
)

K562_se_pacbio <- countAnno(
  annotation = anno_compare_K562_pacbio,
  bamfile = bamfiles[2],
  isLongRead = TRUE,
  strandSpecific = 0,
  isPairedEnd = FALSE,
  nthreads = 8,
  annot.inbuilt = NULL
)

HepG2_combineSE_pacbio <- stringtieCombine(
  reference = "GTF/gff_HepG2_PacBio.annotated.gtf",
  bamFiles = bamfiles[1],
  params = "-p 70",
  gtfFiles = "GTF/HepG2_PacBio.gtf"
)

K562_combineSE_pacbio <- stringtieCombine(
  reference = "GTF/gff_K562_PacBio.annotated.gtf",
  bamFiles = bamfiles[2],
  params = "-p 70",
  gtfFiles = "GTF/K562_PacBio.gtf"
)

saveRDS(anno_compare_HepG2_pacbio, file = "meta/anno_compare_HepG2_pacbio.rds")
saveRDS(anno_compare_K562_pacbio, file = "meta/anno_compare_K562_pacbio.rds")
saveRDS(HepG2_se_pacbio, file = "meta/HepG2_se_pacbio.rds")
saveRDS(K562_se_pacbio, file = "meta/K562_se_pacbio.rds")
saveRDS(HepG2_combineSE_pacbio, file = "meta/HepG2_combineSE_pacbio.rds")
saveRDS(K562_combineSE_pacbio, file = "meta/K562_combineSE_pacbio.rds")

# ------------------------------------------------------------------------------
# 3. Cross-platform pooled annotation and mapping back to pooled reference
# ------------------------------------------------------------------------------

setwd("work/cross_platform")

gtffiles_HepG2_platforms <- c(
  "../pacbio/GTF/gff_HepG2_PacBio.annotated.gtf",
  "../illumina/GTF/gffcompare_HepG2.annotated.gtf",
  "../ont/GTF/gff_HepG2_ontcDNA_all.annotated.gtf"
)

stringtieMerge(
  reference = reference,
  gtfFiles = gtffiles_HepG2_platforms,
  outfile = "GTF/HepG2_pool.gtf",
  params = "-p 18"
)

gtffiles_K562_platforms <- c(
  "../pacbio/GTF/gff_K562_PacBio.annotated.gtf",
  "../illumina/GTF/gffcompare_K562.annotated.gtf",
  "../ont/GTF/gff_K562_ontcDNA_all.annotated.gtf"
)

stringtieMerge(
  reference = reference,
  gtfFiles = gtffiles_K562_platforms,
  outfile = "GTF/K562_pool.gtf",
  params = "-p 18"
)

gffcompareAnno(
  reference = "GTF/HepG2_pool.gtf",
  gtffile = "../pacbio/GTF/gff_HepG2_PacBio.annotated.gtf",
  outfile = "GTF/HepG2_map_PacBio"
)

gffcompareAnno(
  reference = "GTF/HepG2_pool.gtf",
  gtffile = "../illumina/GTF/gffcompare_HepG2.annotated.gtf",
  outfile = "GTF/HepG2_map_Illumina"
)

gffcompareAnno(
  reference = "GTF/HepG2_pool.gtf",
  gtffile = "../ont/GTF/gff_HepG2_ontcDNA_all.annotated.gtf",
  outfile = "GTF/HepG2_map_ONT"
)

gffcompareAnno(
  reference = "GTF/K562_pool.gtf",
  gtffile = "../pacbio/GTF/gff_K562_PacBio.annotated.gtf",
  outfile = "GTF/K562_map_PacBio"
)

gffcompareAnno(
  reference = "GTF/K562_pool.gtf",
  gtffile = "../illumina/GTF/gffcompare_K562.annotated.gtf",
  outfile = "GTF/K562_map_Illumina"
)

gffcompareAnno(
  reference = "GTF/K562_pool.gtf",
  gtffile = "../ont/GTF/gff_K562_ontcDNA_all.annotated.gtf",
  outfile = "GTF/K562_map_ONT"
)

# ------------------------------------------------------------------------------
# 4. Recovery summary based on gffcompare tracking files
# ------------------------------------------------------------------------------

extract_query_tx <- function(x) {
  x <- as.character(x)
  x[x == "-" | is.na(x)] <- NA_character_
  x <- na.omit(x)
  if (length(x) == 0) return(character(0))
  out <- str_extract(x, "(?<=\\|)[^|,;]+$|(?<=:)[^|,;]+")
  out[is.na(out)] <- x[is.na(out)]
  unique(out)
}

# Example: K562 recovery evaluation
k_plat_col_map <- c(
  PacBio   = "gff_K562_PacBio.annotated",
  Illumina = "gffcompare_K562.annotated",
  ONT      = "gff_K562_ontcDNA_all.annotated"
)

k_plat_trk_map <- c(
  PacBio   = "GTF/K562_map_PacBio.tracking",
  Illumina = "GTF/K562_map_Illumina.tracking",
  ONT      = "GTF/K562_map_ONT.tracking"
)

# wide2 should be generated from pooled mapping summary before this step
# and contain one column per platform-specific transcript identifier.

coverage_list_K <- lapply(names(k_plat_trk_map), function(plat) {
  trk_path <- k_plat_trk_map[[plat]]
  trk <- read_tsv(trk_path, col_names = FALSE, comment = "#", show_col_types = FALSE)
  stopifnot(ncol(trk) >= 5)
  
  base_set   <- unique(na.omit(extract_query_tx(trk[[5]])))
  mapped_col <- k_plat_col_map[[plat]]
  stopifnot(mapped_col %in% colnames(wide2))
  mapped_set <- unique(na.omit(wide2[[mapped_col]]))
  matched <- sum(base_set %in% mapped_set)
  
  data.frame(
    cell_line = "K562",
    platform = plat,
    base_total = length(base_set),
    mapped_in_pool = matched,
    coverage_rate = if (length(base_set) > 0) matched / length(base_set) else NA_real_,
    stringsAsFactors = FALSE
  )
})

K562_platform_coverage_by_wideclean <- do.call(rbind, coverage_list_K) %>%
  mutate(coverage_pct = round(100 * coverage_rate, 2))

print(K562_platform_coverage_by_wideclean)

coverage_list_eq_K <- lapply(names(k_plat_trk_map), function(plat) {
  trk_path <- k_plat_trk_map[[plat]]
  trk <- read_tsv(trk_path, col_names = FALSE, comment = "#", show_col_types = FALSE)
  stopifnot(ncol(trk) >= 5)
  
  eq_rows <- trk[[4]] == "="
  base_set_eq <- unique(na.omit(extract_query_tx(trk[[5]][eq_rows])))
  mapped_col  <- k_plat_col_map[[plat]]
  mapped_set  <- unique(na.omit(wide2[[mapped_col]]))
  matched_eq  <- sum(base_set_eq %in% mapped_set)
  
  data.frame(
    cell_line = "K562",
    platform = plat,
    base_total_eq = length(base_set_eq),
    mapped_in_pool_eq = matched_eq,
    coverage_rate_eq = if (length(base_set_eq) > 0) matched_eq / length(base_set_eq) else NA_real_,
    stringsAsFactors = FALSE
  )
})

K562_platform_coverage_EQ_by_wideclean <- do.call(rbind, coverage_list_eq_K) %>%
  mutate(coverage_pct_eq = round(100 * coverage_rate_eq, 2))

