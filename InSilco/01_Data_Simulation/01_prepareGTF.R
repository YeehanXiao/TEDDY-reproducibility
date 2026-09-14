# ==============================================================================
# Script: 01_prepareGTF.R
# Purpose:
#   1. Generate simulated benchmark GTF from official mouse annotation.
#   2. Define strand-aware TE-overlap ground truth.
#   3. Generate RSEM-compatible clean GTF.
#   4. Generate 90% partial reference GTF.
#   5. Generate transcript abundance table for RSEM simulation.
# ==============================================================================

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(IRanges)
  library(dplyr)
})

# ==============================================================================
# User parameters
# ==============================================================================

data_dir <- "./data"
output_dir <- "./results"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

baseline_gtf_path <- file.path(data_dir, "gencode.vM7.annotation.gtf")
te_annotation_rds <- file.path(data_dir, "mm10_TE.rds")
old_iso_path <- file.path(data_dir, "forStat.isoforms.results")

output_gtf_path <- file.path(output_dir, "official_simulated_truth_1000_transcripts.gtf")
output_rsem_gtf_path <- file.path(output_dir, "official_simulated_truth_1000_transcripts.rsem_clean.gtf")
output_ref90_gtf_path <- file.path(output_dir, "official_simulated_reference_90pct.gtf")
output_truth_junction_path <- file.path(output_dir, "official_simulated_truth_junctions.tsv")

n_genes_to_sample <- 1000L
min_te_overlap <- 5L
reference_fraction <- 0.90

min_downstream_if_sparse <- 1L
max_downstream_if_sparse <- 2L
min_downstream_if_rich <- 3L
max_downstream_if_rich <- 5L

filter_core_genes <- TRUE

set.seed(as.integer(Sys.getenv("SIM_SEED", "18")))

# ==============================================================================
# Generate simulated benchmark GTF from official mouse annotation
#
# Main truth:
#   A transcript is TE-chimeric if any exon in its simulated exon-chain
#   overlaps a TE on the same strand by at least min_te_overlap bp.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load GTF and TE
# ------------------------------------------------------------------------------

message("Loading official GTF and TE annotation...")

gtf_all <- import(baseline_gtf_path)
te <- readRDS(te_annotation_rds)

if (!"name" %in% names(mcols(te)) && "names" %in% names(mcols(te))) {
  te$name <- te$names
}
if (!"class" %in% names(mcols(te))) {
  te$class <- NA_character_
}
if (!"family" %in% names(mcols(te))) {
  te$family <- NA_character_
}
if (!"name" %in% names(mcols(te))) {
  stop("TE annotation does not contain 'name' or 'names'.")
}

common_seq <- intersect(seqlevels(gtf_all), seqlevels(te))
if (length(common_seq) == 0) {
  stop("No common seqlevels between GTF and TE.")
}

gtf_all <- keepSeqlevels(gtf_all, common_seq, pruning.mode = "coarse")
te <- keepSeqlevels(te, common_seq, pruning.mode = "coarse")

message("Common seqlevels: ", length(common_seq))

# ------------------------------------------------------------------------------
# 2. Prepare official exon pool
# ------------------------------------------------------------------------------

message("Preparing official exon pool...")

gtf_exon <- gtf_all[gtf_all$type == "exon"]
gtf_exon <- gtf_exon[strand(gtf_exon) != "*"]

if (!"gene_name" %in% names(mcols(gtf_exon))) {
  gtf_exon$gene_name <- gtf_exon$gene_id
}
if (!"gene_type" %in% names(mcols(gtf_exon))) {
  gtf_exon$gene_type <- NA_character_
}

gtf_exon$gene_id <- as.character(gtf_exon$gene_id)
gtf_exon$gene_name <- as.character(gtf_exon$gene_name)
gtf_exon$gene_type <- as.character(gtf_exon$gene_type)
gtf_exon$transcript_id <- as.character(gtf_exon$transcript_id)

gtf_df <- as.data.frame(gtf_exon)
gtf_df$row_id <- seq_len(nrow(gtf_df))

gtf_df <- gtf_df %>%
  group_by(transcript_id) %>%
  arrange(start, end, .by_group = TRUE) %>%
  mutate(
    exon_number = row_number(),
    tx_exon_rank = if (as.character(strand[1]) == "-") {
      rev(seq_len(n()))
    } else {
      seq_len(n())
    }
  ) %>%
  ungroup() %>%
  arrange(row_id)

gtf_exon$exon_number <- gtf_df$exon_number
gtf_exon$tx_exon_rank <- gtf_df$tx_exon_rank

# ------------------------------------------------------------------------------
# 3. Mark official strand-aware TE-overlap exons
# ------------------------------------------------------------------------------

message("Marking strand-aware TE-overlap exons...")

gtf_exon$is_TE_overlap_exon <- FALSE
gtf_exon$TE_name <- "none"
gtf_exon$TE_class <- "none"
gtf_exon$TE_family <- "none"
gtf_exon$TE_overlap_width <- 0L

hits <- findOverlaps(
  gtf_exon,
  te,
  ignore.strand = FALSE,
  minoverlap = min_te_overlap
)

message("Strand-aware TE-overlap exon hits: ", length(hits))

if (length(hits) > 0) {
  ov_width <- width(pintersect(
    gtf_exon[queryHits(hits)],
    te[subjectHits(hits)],
    ignore.strand = FALSE
  ))
  
  hit_tbl <- data.frame(
    exon_idx = queryHits(hits),
    TE_name = as.character(te$name[subjectHits(hits)]),
    TE_class = as.character(te$class[subjectHits(hits)]),
    TE_family = as.character(te$family[subjectHits(hits)]),
    TE_overlap_width = ov_width,
    stringsAsFactors = FALSE
  )
  
  hit_tbl2 <- hit_tbl %>%
    group_by(exon_idx) %>%
    summarise(
      TE_name = paste(unique(TE_name), collapse = ","),
      TE_class = paste(unique(TE_class), collapse = ","),
      TE_family = paste(unique(TE_family), collapse = ","),
      TE_overlap_width = max(TE_overlap_width, na.rm = TRUE),
      .groups = "drop"
    )
  
  gtf_exon$is_TE_overlap_exon[hit_tbl2$exon_idx] <- TRUE
  gtf_exon$TE_name[hit_tbl2$exon_idx] <- hit_tbl2$TE_name
  gtf_exon$TE_class[hit_tbl2$exon_idx] <- hit_tbl2$TE_class
  gtf_exon$TE_family[hit_tbl2$exon_idx] <- hit_tbl2$TE_family
  gtf_exon$TE_overlap_width[hit_tbl2$exon_idx] <- hit_tbl2$TE_overlap_width
}

# ------------------------------------------------------------------------------
# 4. Define clean gene pool
# ------------------------------------------------------------------------------

first_exons <- gtf_exon[gtf_exon$tx_exon_rank == 1]
all_exons <- gtf_exon

candidate_gene_tbl <- as.data.frame(first_exons) %>%
  transmute(
    gene_id = as.character(gene_id),
    gene_name = as.character(gene_name),
    gene_type = as.character(gene_type)
  ) %>%
  distinct() %>%
  filter(!is.na(gene_id), gene_id != "")

if (filter_core_genes) {
  candidate_gene_tbl <- candidate_gene_tbl %>%
    filter(
      gene_type %in% c("protein_coding", "lncRNA", "lincRNA"),
      !grepl("^Gm[0-9]+", gene_name),
      !grepl("^RP[0-9]", gene_name),
      !grepl("Rik$", gene_name),
      !grepl("-ps$", gene_name)
    )
}

candidate_gene_ids <- unique(candidate_gene_tbl$gene_id)

message("Candidate genes after filtering: ", length(candidate_gene_ids))

if (length(candidate_gene_ids) < n_genes_to_sample) {
  warning("Only ", length(candidate_gene_ids), " candidate genes available; using all.")
  n_genes_to_sample <- length(candidate_gene_ids)
}

# ------------------------------------------------------------------------------
# 5. Simulate exon-chain
# ------------------------------------------------------------------------------

simulate_one_isoform <- function(gene_id, iso_idx, first_exons, all_exons, transcript_id_prefix) {
  gene_first_exons <- first_exons[first_exons$gene_id == gene_id]
  if (length(gene_first_exons) == 0) return(NULL)
  
  initiated_exon <- gene_first_exons[sample(seq_along(gene_first_exons), 1)]
  gene_strand <- as.character(strand(initiated_exon))
  
  gene_exons <- all_exons[all_exons$gene_id == gene_id]
  if (length(gene_exons) == 0) return(NULL)
  
  if (gene_strand == "+") {
    downstream_exons <- gene_exons[start(gene_exons) > end(initiated_exon)]
  } else {
    downstream_exons <- gene_exons[end(gene_exons) < start(initiated_exon)]
  }
  
  if (length(downstream_exons) > 0) {
    downstream_exons <- downstream_exons[order(start(downstream_exons), end(downstream_exons))]
    
    n_downstream <- if (length(downstream_exons) < 5) {
      sample(min_downstream_if_sparse:max_downstream_if_sparse, 1)
    } else {
      sample(min_downstream_if_rich:max_downstream_if_rich, 1)
    }
    
    n_downstream <- min(n_downstream, length(downstream_exons))
    
    selected_exons <- c(
      initiated_exon,
      downstream_exons[sample(seq_along(downstream_exons), n_downstream)]
    )
  } else {
    selected_exons <- initiated_exon
  }
  
  selected_exons <- selected_exons[order(start(selected_exons), end(selected_exons))]
  
  selected_df <- as.data.frame(selected_exons)
  keep_idx <- !duplicated(
    paste(selected_df$seqnames, selected_df$start, selected_df$end, selected_df$strand, sep = ":")
  )
  selected_exons <- selected_exons[keep_idx]
  
  if (length(selected_exons) == 0) return(NULL)
  
  new_tx_id <- paste0(transcript_id_prefix, gene_id, "_iso", iso_idx)
  
  exon_out <- selected_exons
  exon_out$type <- "exon"
  exon_out$gene_id <- gene_id
  exon_out$gene_name <- selected_exons$gene_name[1]
  exon_out$gene_type <- selected_exons$gene_type[1]
  exon_out$transcript_id <- new_tx_id
  exon_out$transcript_name <- new_tx_id
  
  exon_df <- as.data.frame(exon_out)
  exon_df$row_id <- seq_len(nrow(exon_df))
  
  exon_df <- exon_df %>%
    arrange(start, end) %>%
    mutate(
      exon_number = row_number(),
      tx_exon_rank = if (gene_strand == "-") {
        rev(seq_len(n()))
      } else {
        seq_len(n())
      }
    ) %>%
    arrange(row_id)
  
  exon_out$exon_number <- exon_df$exon_number
  exon_out$tx_exon_rank <- exon_df$tx_exon_rank
  
  tx_range <- range(exon_out)
  
  tx_out <- GRanges(
    seqnames = seqnames(tx_range),
    ranges = ranges(tx_range),
    strand = strand(tx_range),
    type = "transcript",
    gene_id = gene_id,
    gene_name = exon_out$gene_name[1],
    gene_type = exon_out$gene_type[1],
    transcript_id = new_tx_id,
    transcript_name = new_tx_id,
    is_TE_chimeric_tx = any(exon_out$is_TE_overlap_exon),
    n_TE_overlap_exons = sum(exon_out$is_TE_overlap_exon),
    max_TE_overlap_width = max(exon_out$TE_overlap_width, na.rm = TRUE)
  )
  
  c(tx_out, exon_out)
}

simulate_one_gene <- function(gene_id, first_exons, all_exons, transcript_id_prefix) {
  gene_exons <- all_exons[all_exons$gene_id == gene_id]
  if (length(gene_exons) == 0) return(NULL)
  
  n_unique_exons <- length(unique(
    paste(seqnames(gene_exons), start(gene_exons), end(gene_exons), strand(gene_exons), sep = ":")
  ))
  
  n_isoforms <- if (n_unique_exons < 5) {
    sample(1:2, 1)
  } else {
    sample(2:3, 1)
  }
  
  iso_list <- lapply(
    seq_len(n_isoforms),
    simulate_one_isoform,
    gene_id = gene_id,
    first_exons = first_exons,
    all_exons = all_exons,
    transcript_id_prefix = transcript_id_prefix
  )
  
  iso_list <- Filter(Negate(is.null), iso_list)
  if (length(iso_list) == 0) return(NULL)
  
  do.call(c, iso_list)
}

message("Sampling genes and simulating transcript models...")

sampled_gene_ids <- sample(candidate_gene_ids, n_genes_to_sample, replace = FALSE)

simulated_list <- lapply(
  sampled_gene_ids,
  simulate_one_gene,
  first_exons = first_exons,
  all_exons = all_exons,
  transcript_id_prefix = "OFFICIAL_SIMTX_"
)

simulated_list <- Filter(Negate(is.null), simulated_list)

if (length(simulated_list) == 0) {
  stop("No simulated transcripts were generated.")
}

simulated_raw <- sort(do.call(c, simulated_list))
simulated_exon_raw <- simulated_raw[simulated_raw$type == "exon"]

# ------------------------------------------------------------------------------
# 6. Standard GTF builder
# ------------------------------------------------------------------------------

make_standard_gtf <- function(exon_gr) {
  stopifnot(all(exon_gr$type == "exon"))
  
  exon_df <- as.data.frame(exon_gr)
  
  needed <- c(
    "gene_id", "gene_name", "gene_type", "transcript_id", "transcript_name",
    "is_TE_overlap_exon", "TE_name", "TE_class", "TE_family", "TE_overlap_width"
  )
  for (cc in needed) {
    if (!cc %in% colnames(exon_df)) exon_df[[cc]] <- NA
  }
  
  exon_df <- exon_df %>%
    mutate(
      seqnames = as.character(seqnames),
      strand = as.character(strand),
      gene_id = as.character(gene_id),
      gene_name = as.character(gene_name),
      gene_type = as.character(gene_type),
      transcript_id = as.character(transcript_id),
      transcript_name = as.character(transcript_id),
      is_TE_overlap_exon = as.logical(is_TE_overlap_exon),
      is_TE_overlap_exon = ifelse(is.na(is_TE_overlap_exon), FALSE, is_TE_overlap_exon),
      TE_name = ifelse(is.na(TE_name), "none", as.character(TE_name)),
      TE_class = ifelse(is.na(TE_class), "none", as.character(TE_class)),
      TE_family = ifelse(is.na(TE_family), "none", as.character(TE_family)),
      TE_overlap_width = ifelse(is.na(TE_overlap_width), 0L, as.integer(TE_overlap_width))
    ) %>%
    arrange(seqnames, gene_id, transcript_id, start, end)
  
  exon_df <- exon_df %>%
    group_by(transcript_id) %>%
    arrange(start, end, .by_group = TRUE) %>%
    mutate(
      exon_number = row_number(),
      tx_exon_rank = if (as.character(strand[1]) == "-") {
        rev(seq_len(n()))
      } else {
        seq_len(n())
      }
    ) %>%
    ungroup()
  
  tx_df <- exon_df %>%
    group_by(seqnames, strand, gene_id, gene_name, gene_type, transcript_id) %>%
    summarise(
      start = min(start),
      end = max(end),
      is_TE_chimeric_tx = any(is_TE_overlap_exon),
      n_TE_overlap_exons = sum(is_TE_overlap_exon),
      max_TE_overlap_width = max(TE_overlap_width, na.rm = TRUE),
      .groups = "drop"
    )
  
  gene_df <- tx_df %>%
    group_by(seqnames, strand, gene_id, gene_name, gene_type) %>%
    summarise(
      start = min(start),
      end = max(end),
      is_TE_chimeric_gene = any(is_TE_chimeric_tx),
      n_TE_chimeric_tx = sum(is_TE_chimeric_tx),
      n_TE_overlap_exons = sum(n_TE_overlap_exons),
      max_TE_overlap_width = max(max_TE_overlap_width, na.rm = TRUE),
      .groups = "drop"
    )
  
  gene_gr <- GRanges(
    seqnames = gene_df$seqnames,
    ranges = IRanges(gene_df$start, gene_df$end),
    strand = gene_df$strand,
    type = "gene",
    gene_id = gene_df$gene_id,
    gene_name = gene_df$gene_name,
    gene_type = gene_df$gene_type,
    is_TE_chimeric_gene = gene_df$is_TE_chimeric_gene,
    n_TE_chimeric_tx = gene_df$n_TE_chimeric_tx,
    n_TE_overlap_exons = gene_df$n_TE_overlap_exons,
    max_TE_overlap_width = gene_df$max_TE_overlap_width
  )
  
  tx_gr <- GRanges(
    seqnames = tx_df$seqnames,
    ranges = IRanges(tx_df$start, tx_df$end),
    strand = tx_df$strand,
    type = "transcript",
    gene_id = tx_df$gene_id,
    gene_name = tx_df$gene_name,
    gene_type = tx_df$gene_type,
    transcript_id = tx_df$transcript_id,
    transcript_name = tx_df$transcript_id,
    is_TE_chimeric_tx = tx_df$is_TE_chimeric_tx,
    n_TE_overlap_exons = tx_df$n_TE_overlap_exons,
    max_TE_overlap_width = tx_df$max_TE_overlap_width
  )
  
  exon_gr2 <- GRanges(
    seqnames = exon_df$seqnames,
    ranges = IRanges(exon_df$start, exon_df$end),
    strand = exon_df$strand,
    type = "exon",
    gene_id = exon_df$gene_id,
    gene_name = exon_df$gene_name,
    gene_type = exon_df$gene_type,
    transcript_id = exon_df$transcript_id,
    transcript_name = exon_df$transcript_id,
    exon_number = exon_df$exon_number,
    tx_exon_rank = exon_df$tx_exon_rank,
    is_TE_overlap_exon = exon_df$is_TE_overlap_exon,
    TE_name = exon_df$TE_name,
    TE_class = exon_df$TE_class,
    TE_family = exon_df$TE_family,
    TE_overlap_width = exon_df$TE_overlap_width
  )
  
  out <- c(gene_gr, tx_gr, exon_gr2)
  out$feature_order <- match(out$type, c("gene", "transcript", "exon"))
  
  out <- out[order(
    as.character(seqnames(out)),
    start(out),
    out$feature_order,
    as.character(out$gene_id),
    as.character(out$transcript_id),
    start(out),
    end(out)
  )]
  
  out$feature_order <- NULL
  out
}

make_rsem_clean_gtf <- function(gtf_gr) {
  keep_cols <- c(
    "type", "gene_id", "gene_name", "gene_type",
    "transcript_id", "transcript_name", "exon_number"
  )
  keep_cols <- intersect(keep_cols, names(mcols(gtf_gr)))
  mcols(gtf_gr) <- mcols(gtf_gr)[, keep_cols, drop = FALSE]
  gtf_gr
}

simulated_gtf <- make_standard_gtf(simulated_exon_raw)

export(simulated_gtf, output_gtf_path, format = "gtf")
export(make_rsem_clean_gtf(simulated_gtf), output_rsem_gtf_path, format = "gtf")

message("Done truth GTF: ", output_gtf_path)
message("Done RSEM-clean GTF: ", output_rsem_gtf_path)

# ------------------------------------------------------------------------------
# 7. Build 90% reference GTF
# ------------------------------------------------------------------------------

message("Generating 90% reference GTF...")

tx_tbl <- as.data.frame(simulated_gtf[simulated_gtf$type == "transcript"]) %>%
  transmute(
    transcript_id = as.character(transcript_id),
    gene_id = as.character(gene_id)
  ) %>%
  distinct()

target_n_ref <- floor(nrow(tx_tbl) * reference_fraction)

anchor_tx <- tx_tbl %>%
  group_by(gene_id) %>%
  slice_sample(n = 1) %>%
  ungroup()

remaining_tx <- tx_tbl %>%
  filter(!transcript_id %in% anchor_tx$transcript_id)

n_extra <- target_n_ref - nrow(anchor_tx)

extra_tx <- if (n_extra > 0) {
  remaining_tx %>% slice_sample(n = min(n_extra, nrow(remaining_tx)))
} else {
  remaining_tx[0, ]
}

ref_tx_ids <- bind_rows(anchor_tx, extra_tx) %>%
  distinct(transcript_id) %>%
  pull(transcript_id)

ref_exon <- simulated_gtf[
  simulated_gtf$type == "exon" &
    simulated_gtf$transcript_id %in% ref_tx_ids
]

reference_90pct_gtf <- make_standard_gtf(ref_exon)

export(reference_90pct_gtf, output_ref90_gtf_path, format = "gtf")

message("Done reference 90% GTF: ", output_ref90_gtf_path)

# ------------------------------------------------------------------------------
# 8. Truth tables
# ------------------------------------------------------------------------------

message("Exporting truth tables...")

sim_exon <- simulated_gtf[simulated_gtf$type == "exon"]

exon_truth_tbl <- as.data.frame(sim_exon) %>%
  transmute(
    transcript_id = as.character(transcript_id),
    gene_id = as.character(gene_id),
    gene_name = as.character(gene_name),
    gene_type = as.character(gene_type),
    seqnames = as.character(seqnames),
    start,
    end,
    strand = as.character(strand),
    exon_number = as.integer(exon_number),
    tx_exon_rank = as.integer(tx_exon_rank),
    is_TE_overlap_exon = as.logical(is_TE_overlap_exon),
    TE_name = as.character(TE_name),
    TE_class = as.character(TE_class),
    TE_family = as.character(TE_family),
    TE_overlap_width = as.integer(TE_overlap_width)
  )

tx_truth_tbl <- exon_truth_tbl %>%
  group_by(transcript_id, gene_id, gene_name, gene_type) %>%
  summarise(
    n_exons = n(),
    n_TE_overlap_exons = sum(is_TE_overlap_exon),
    is_TE_chimeric_tx = any(is_TE_overlap_exon),
    max_TE_overlap_width = max(TE_overlap_width, na.rm = TRUE),
    TE_names = paste(unique(TE_name[TE_name != "none"]), collapse = ","),
    TE_classes = paste(unique(TE_class[TE_class != "none"]), collapse = ","),
    TE_families = paste(unique(TE_family[TE_family != "none"]), collapse = ","),
    .groups = "drop"
  ) %>%
  mutate(
    TE_names = ifelse(TE_names == "", "none", TE_names),
    TE_classes = ifelse(TE_classes == "", "none", TE_classes),
    TE_families = ifelse(TE_families == "", "none", TE_families)
  )

gene_truth_tbl <- tx_truth_tbl %>%
  group_by(gene_id, gene_name, gene_type) %>%
  summarise(
    n_simulated_tx = n_distinct(transcript_id),
    n_TE_chimeric_tx = sum(is_TE_chimeric_tx),
    is_TE_chimeric_gene = any(is_TE_chimeric_tx),
    n_TE_overlap_exons = sum(n_TE_overlap_exons),
    max_TE_overlap_width = max(max_TE_overlap_width, na.rm = TRUE),
    .groups = "drop"
  )

write.table(
  exon_truth_tbl,
  file.path(output_dir, "official_simulated_exon_truth_status.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  tx_truth_tbl,
  file.path(output_dir, "official_simulated_transcript_truth_status.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  gene_truth_tbl,
  file.path(output_dir, "official_simulated_gene_truth_status.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 9. Secondary TE-host junction table
# ------------------------------------------------------------------------------

sim_df <- exon_truth_tbl %>%
  arrange(transcript_id, tx_exon_rank)

truth_junction_all <- sim_df %>%
  group_by(transcript_id) %>%
  arrange(tx_exon_rank, .by_group = TRUE) %>%
  mutate(
    next_start = lead(start),
    next_end = lead(end),
    next_is_TE_overlap_exon = lead(is_TE_overlap_exon),
    next_TE_name = lead(TE_name),
    next_TE_class = lead(TE_class),
    next_TE_family = lead(TE_family),
    next_TE_overlap_width = lead(TE_overlap_width)
  ) %>%
  ungroup() %>%
  filter(!is.na(next_start)) %>%
  mutate(
    left_is_TE = is_TE_overlap_exon,
    right_is_TE = next_is_TE_overlap_exon,
    junction_type = case_when(
      left_is_TE & !right_is_TE ~ "TE_to_host",
      !left_is_TE & right_is_TE ~ "host_to_TE",
      left_is_TE & right_is_TE ~ "TE_to_TE",
      TRUE ~ "host_to_host"
    ),
    donor_site = ifelse(strand == "+", end, start),
    acceptor_site = ifelse(strand == "+", next_start, next_end)
  )

truth_junction_detail <- truth_junction_all %>%
  filter(junction_type %in% c("TE_to_host", "host_to_TE")) %>%
  transmute(
    transcript_id,
    gene_id,
    gene_name,
    gene_type,
    seqnames,
    strand,
    left_start = start,
    left_end = end,
    left_is_TE,
    left_TE_name = TE_name,
    left_TE_class = TE_class,
    left_TE_family = TE_family,
    left_TE_overlap_width = TE_overlap_width,
    right_start = next_start,
    right_end = next_end,
    right_is_TE,
    right_TE_name = next_TE_name,
    right_TE_class = next_TE_class,
    right_TE_family = next_TE_family,
    right_TE_overlap_width = next_TE_overlap_width,
    junction_type,
    donor_site,
    acceptor_site
  )

write.table(
  truth_junction_detail,
  output_truth_junction_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  truth_junction_all,
  file.path(output_dir, "official_simulated_all_junctions.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 10. Summary
# ------------------------------------------------------------------------------

gene_event_summary <- data.frame(
  simulated_genes = n_distinct(gene_truth_tbl$gene_id),
  simulated_transcripts = n_distinct(tx_truth_tbl$transcript_id),
  simulated_exons = nrow(exon_truth_tbl),
  TE_chimeric_genes = sum(gene_truth_tbl$is_TE_chimeric_gene),
  non_TE_chimeric_genes = sum(!gene_truth_tbl$is_TE_chimeric_gene),
  TE_chimeric_transcripts = sum(tx_truth_tbl$is_TE_chimeric_tx),
  non_TE_chimeric_transcripts = sum(!tx_truth_tbl$is_TE_chimeric_tx),
  TE_overlap_exons = sum(exon_truth_tbl$is_TE_overlap_exon),
  TE_host_truth_junctions = nrow(truth_junction_detail),
  all_exon_exon_junctions = nrow(truth_junction_all),
  min_te_overlap = min_te_overlap,
  reference_fraction = reference_fraction,
  filter_core_genes = filter_core_genes,
  sampled_Gm_gene_fraction = mean(grepl("^Gm[0-9]+", gene_truth_tbl$gene_name)),
  sampled_RP_gene_fraction = mean(grepl("^RP[0-9]", gene_truth_tbl$gene_name)),
  sampled_Rik_gene_fraction = mean(grepl("Rik$", gene_truth_tbl$gene_name)),
  sampled_ps_gene_fraction = mean(grepl("-ps$", gene_truth_tbl$gene_name)),
  sampled_protein_coding_gene_fraction = mean(gene_truth_tbl$gene_type == "protein_coding", na.rm = TRUE)
)

write.table(
  gene_event_summary,
  file.path(output_dir, "official_simulated_gene_event_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nTruth GTF feature types:\n")
print(table(simulated_gtf$type))

cat("\nReference 90% GTF feature types:\n")
print(table(reference_90pct_gtf$type))

cat("\nJunction type table:\n")
print(table(truth_junction_all$junction_type))

cat("\nGene event summary:\n")
print(gene_event_summary)

# ==============================================================================
# Generate RSEM isoforms.results from simulated truth GTF
# ==============================================================================

truth_gtf <- rtracklayer::import(output_rsem_gtf_path)

tx_map <- as.data.frame(truth_gtf[truth_gtf$type == "exon"]) %>%
  transmute(
    transcript_id = as.character(transcript_id),
    gene_id = as.character(gene_id),
    exon_length = as.integer(end - start + 1)
  ) %>%
  group_by(transcript_id, gene_id) %>%
  summarise(
    length = sum(exon_length),
    .groups = "drop"
  )

old_iso <- read.delim(
  old_iso_path,
  check.names = FALSE
)

expr_pool <- old_iso %>%
  filter(!is.na(TPM), TPM >= 0) %>%
  pull(TPM)

set.seed(as.integer(Sys.getenv("SIM_SEED", "18")))

sim_tpm <- sample(expr_pool, nrow(tx_map), replace = TRUE)

if (sum(sim_tpm) == 0) {
  sim_tpm <- sample(old_iso$TPM[old_iso$TPM > 0], nrow(tx_map), replace = TRUE)
}

sim_tpm <- sim_tpm / sum(sim_tpm) * 1e6

out <- tx_map %>%
  mutate(
    effective_length = pmax(length - 150 + 1, 1),
    TPM = sim_tpm,
    FPKM = TPM,
    expected_count = TPM * effective_length
  )

out$expected_count <- out$expected_count / sum(out$expected_count) * 1e6
out$IsoPct <- 100

out <- out %>%
  select(
    transcript_id,
    gene_id,
    length,
    effective_length,
    expected_count,
    TPM,
    FPKM,
    IsoPct
  )

write.table(
  out,
  file.path(output_dir, "official_simulated_1000.isoforms.results"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nGenerated isoforms.results:\n")
cat("n_tx:", nrow(out), "\n")
cat("TPM_sum:", sum(out$TPM), "\n")
cat("n_TPM_gt_0:", sum(out$TPM > 0), "\n")
cat("n_TPM_gt_1:", sum(out$TPM > 1), "\n")
cat("overlap_with_truth:", length(intersect(out$transcript_id, tx_map$transcript_id)), "\n")
print(summary(out$length))
print(summary(out$TPM))

message("All done.")