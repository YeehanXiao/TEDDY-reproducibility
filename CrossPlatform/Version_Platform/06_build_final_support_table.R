#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(openxlsx)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4L) {
  stop(
    paste0(
      "Usage:\n",
      "Rscript 04_build_final_support_table.R ",
      "<twoC_candidates.rds> ",
      "<nanopore_tx_support.tsv> ",
      "<annotation_best.rds> ",
      "<output_prefix>"
    )
  )
}

candidate_file <- args[1]
support_file <- args[2]
annotation_file <- args[3]
output_prefix <- args[4]

dir.create(
  dirname(output_prefix),
  recursive = TRUE,
  showWarnings = FALSE
)

stopifnot(
  file.exists(candidate_file),
  file.exists(support_file),
  file.exists(annotation_file)
)

## ============================================================
## 1. Load inputs
## ============================================================

candidate_tbl <- readRDS(candidate_file)

support_tbl <- readr::read_tsv(
  support_file,
  show_col_types = FALSE
)

annotation_best <- readRDS(annotation_file)

candidate_tbl <- as.data.frame(candidate_tbl) |>
  tibble::as_tibble()

support_tbl <- as.data.frame(support_tbl) |>
  tibble::as_tibble()

annotation_best <- as.data.frame(annotation_best) |>
  tibble::as_tibble()

## ============================================================
## 2. Validate required columns
## ============================================================

required_candidate_cols <- c(
  "transcript_id",
  "gene_name",
  "TE_name",
  "TE_class"
)

required_support_cols <- c(
  "transcript_id",
  "eval_type",
  "max_TE_cov_frac",
  "TE_host_junction_support_10bp",
  "TE_host_anchor75_support"
)

required_annotation_cols <- c(
  "TEDDY_vM7_tx_id",
  "TEDDY_vM25_traceable",
  "GENCODE_vM25_matched_tx_id"
)

missing_candidate_cols <- setdiff(
  required_candidate_cols,
  colnames(candidate_tbl)
)

missing_support_cols <- setdiff(
  required_support_cols,
  colnames(support_tbl)
)

missing_annotation_cols <- setdiff(
  required_annotation_cols,
  colnames(annotation_best)
)

if (length(missing_candidate_cols) > 0L) {
  stop(
    "Missing columns in candidate table: ",
    paste(missing_candidate_cols, collapse = ", ")
  )
}

if (length(missing_support_cols) > 0L) {
  stop(
    "Missing columns in Nanopore support table: ",
    paste(missing_support_cols, collapse = ", ")
  )
}

if (length(missing_annotation_cols) > 0L) {
  stop(
    "Missing columns in annotation-version table: ",
    paste(missing_annotation_cols, collapse = ", ")
  )
}

## Ensure one row per candidate transcript.
candidate_tbl <- candidate_tbl |>
  dplyr::mutate(
    transcript_id = as.character(transcript_id)
  )

if (anyDuplicated(candidate_tbl$transcript_id) > 0L) {
  stop("Candidate table contains duplicated transcript IDs.")
}

## ============================================================
## 3. Annotation-version status
## ============================================================

annotation_status <- annotation_best |>
  dplyr::transmute(
    transcript_id = as.character(TEDDY_vM7_tx_id),
    
    Teddy_vM25_reconstructed = dplyr::if_else(
      tidyr::replace_na(
        as.logical(TEDDY_vM25_traceable),
        FALSE
      ),
      "Yes",
      "-"
    ),
    
    Catalog_vM25 = dplyr::if_else(
      !is.na(GENCODE_vM25_matched_tx_id) &
        as.character(GENCODE_vM25_matched_tx_id) != "" &
        as.character(GENCODE_vM25_matched_tx_id) != "-",
      "Yes",
      "-"
    )
  ) |>
  dplyr::group_by(transcript_id) |>
  dplyr::summarise(
    Teddy_vM25_reconstructed = dplyr::if_else(
      any(Teddy_vM25_reconstructed == "Yes"),
      "Yes",
      "-"
    ),
    
    Catalog_vM25 = dplyr::if_else(
      any(Catalog_vM25 == "Yes"),
      "Yes",
      "-"
    ),
    
    .groups = "drop"
  )

## ============================================================
## 4. Nanopore support status
## ============================================================

support_status <- support_tbl |>
  dplyr::transmute(
    transcript_id = as.character(transcript_id),
    eval_type = as.character(eval_type),
    
    Nanopore_JunctionSupported = dplyr::if_else(
      tidyr::replace_na(
        TE_host_junction_support_10bp == "Yes",
        FALSE
      ),
      "Yes",
      "-"
    ),
    
    Nanopore_ExonchainSupported = dplyr::case_when(
      eval_type == "multi_exon" &
        tidyr::replace_na(
          TE_host_anchor75_support == "Yes",
          FALSE
        ) ~ "Yes",
      
      eval_type == "single_exon" &
        tidyr::replace_na(
          TE_host_junction_support_10bp == "Yes",
          FALSE
        ) &
        tidyr::replace_na(
          as.numeric(max_TE_cov_frac) >= 0.75,
          FALSE
        ) ~ "Yes",
      
      TRUE ~ "-"
    )
  ) |>
  dplyr::group_by(transcript_id) |>
  dplyr::summarise(
    Nanopore_JunctionSupported = dplyr::if_else(
      any(Nanopore_JunctionSupported == "Yes"),
      "Yes",
      "-"
    ),
    
    Nanopore_ExonchainSupported = dplyr::if_else(
      any(Nanopore_ExonchainSupported == "Yes"),
      "Yes",
      "-"
    ),
    
    .groups = "drop"
  )

## ============================================================
## 5. Short single-exon candidates
## ============================================================

## Five short single-exon candidates lacked a definable >=15-bp host
## flank and were confirmed by inspection of the Nanopore-derived
## TE-chimeric transcript annotation at the corresponding loci.

short_single_exon_junction_supported <- c(
  "ENSMUST00000196294.1",
  "ENSMUST00000116700.1",
  "ENSMUST00000116704.1",
  "ENSMUST00000116730.1",
  "ENSMUST00000116769.1"
)

stopifnot(
  all(
    short_single_exon_junction_supported %in%
      candidate_tbl$transcript_id
  )
)

## These candidates receive junction support only.
## Their exon-chain support status is not altered.

## ============================================================
## 6. Build Supplementary Table 4
## ============================================================

final_tbl <- candidate_tbl |>
  dplyr::left_join(
    annotation_status,
    by = "transcript_id"
  ) |>
  dplyr::left_join(
    support_status,
    by = "transcript_id"
  ) |>
  dplyr::mutate(
    Teddy_vM25_reconstructed = dplyr::if_else(
      Teddy_vM25_reconstructed == "Yes",
      "Yes",
      "-"
    ),
    
    Catalog_vM25 = dplyr::if_else(
      Catalog_vM25 == "Yes",
      "Yes",
      "-"
    ),
    
    Nanopore_JunctionSupported = dplyr::case_when(
      Nanopore_JunctionSupported == "Yes" ~ "Yes",
      
      transcript_id %in%
        short_single_exon_junction_supported ~ "Yes",
      
      TRUE ~ "-"
    ),
    
    Nanopore_ExonchainSupported = dplyr::if_else(
      Nanopore_ExonchainSupported == "Yes",
      "Yes",
      "-"
    )
  )

## ============================================================
## 7. Arrange final columns
## ============================================================

preferred_columns <- c(
  "seqnames",
  "start",
  "end",
  "strand",
  "transcript_id",
  "gene_name",
  "TE_name",
  "TE_class",
  "Teddy_vM25_reconstructed",
  "Catalog_vM25",
  "Nanopore_JunctionSupported",
  "Nanopore_ExonchainSupported"
)

existing_preferred_columns <- intersect(
  preferred_columns,
  colnames(final_tbl)
)

remaining_columns <- setdiff(
  colnames(final_tbl),
  existing_preferred_columns
)

final_tbl <- final_tbl |>
  dplyr::select(
    dplyr::all_of(existing_preferred_columns),
    dplyr::all_of(remaining_columns)
  )

## Rename the genomic chromosome column for the final table.
if ("seqnames" %in% colnames(final_tbl)) {
  final_tbl <- final_tbl |>
    dplyr::rename(
      Chromosome = seqnames
    )
}

final_tbl <- final_tbl |>
  dplyr::rename(
    Transcript_id = transcript_id,
    Gene_name = gene_name
  )

## ============================================================
## 8. Sanity checks
## ============================================================

stopifnot(
  nrow(final_tbl) == nrow(candidate_tbl),
  
  dplyr::n_distinct(final_tbl$Transcript_id) ==
    nrow(final_tbl),
  
  !anyNA(final_tbl$Teddy_vM25_reconstructed),
  
  !anyNA(final_tbl$Catalog_vM25),
  
  !anyNA(final_tbl$Nanopore_JunctionSupported),
  
  !anyNA(final_tbl$Nanopore_ExonchainSupported),
  
  all(
    final_tbl$Teddy_vM25_reconstructed %in% c("Yes", "-")
  ),
  
  all(
    final_tbl$Catalog_vM25 %in% c("Yes", "-")
  ),
  
  all(
    final_tbl$Nanopore_JunctionSupported %in% c("Yes", "-")
  ),
  
  all(
    final_tbl$Nanopore_ExonchainSupported %in% c("Yes", "-")
  )
)

## ============================================================
## 9. Summary
## ============================================================

summary_tbl <- final_tbl |>
  dplyr::summarise(
    n_transcripts = dplyr::n_distinct(Transcript_id),
    
    n_vM25_reconstructed = sum(
      Teddy_vM25_reconstructed == "Yes"
    ),
    
    pct_vM25_reconstructed = round(
      100 * n_vM25_reconstructed / n_transcripts,
      2
    ),
    
    n_catalog_vM25 = sum(
      Catalog_vM25 == "Yes"
    ),
    
    pct_catalog_vM25 = round(
      100 * n_catalog_vM25 / n_transcripts,
      2
    ),
    
    n_nanopore_junction_supported = sum(
      Nanopore_JunctionSupported == "Yes"
    ),
    
    pct_nanopore_junction_supported = round(
      100 * n_nanopore_junction_supported /
        n_transcripts,
      2
    ),
    
    n_nanopore_exonchain_supported = sum(
      Nanopore_ExonchainSupported == "Yes"
    ),
    
    pct_nanopore_exonchain_supported = round(
      100 * n_nanopore_exonchain_supported /
        n_transcripts,
      2
    )
  )

## Expected values for the current 633-transcript candidate set.
expected_summary <- tibble::tibble(
  metric = c(
    "n_transcripts",
    "n_nanopore_junction_supported",
    "n_nanopore_exonchain_supported"
  ),
  expected = c(
    633L,
    615L,
    515L
  ),
  observed = c(
    summary_tbl$n_transcripts,
    summary_tbl$n_nanopore_junction_supported,
    summary_tbl$n_nanopore_exonchain_supported
  )
)

expected_summary <- expected_summary |>
  dplyr::mutate(
    matched = observed == expected
  )

if (!all(expected_summary$matched)) {
  warning(
    paste0(
      "One or more output counts differ from the expected ",
      "633-transcript analysis. Check the summary table."
    )
  )
}

## ============================================================
## 10. Write outputs
## ============================================================

readr::write_tsv(
  final_tbl,
  paste0(output_prefix, ".tsv")
)

readr::write_tsv(
  summary_tbl,
  paste0(output_prefix, ".summary.tsv")
)

readr::write_tsv(
  expected_summary,
  paste0(output_prefix, ".expected_counts_check.tsv")
)

saveRDS(
  final_tbl,
  paste0(output_prefix, ".rds")
)

workbook <- openxlsx::createWorkbook()

openxlsx::addWorksheet(
  workbook,
  "Supplementary Table 4"
)

openxlsx::writeData(
  workbook,
  sheet = "Supplementary Table 4",
  x = final_tbl,
  withFilter = TRUE
)

openxlsx::freezePane(
  workbook,
  sheet = "Supplementary Table 4",
  firstRow = TRUE
)

openxlsx::setColWidths(
  workbook,
  sheet = "Supplementary Table 4",
  cols = seq_len(ncol(final_tbl)),
  widths = "auto"
)

header_style <- openxlsx::createStyle(
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  wrapText = TRUE
)

openxlsx::addStyle(
  workbook,
  sheet = "Supplementary Table 4",
  style = header_style,
  rows = 1,
  cols = seq_len(ncol(final_tbl)),
  gridExpand = TRUE,
  stack = TRUE
)

openxlsx::addWorksheet(
  workbook,
  "Summary"
)

openxlsx::writeData(
  workbook,
  sheet = "Summary",
  x = summary_tbl
)

openxlsx::addWorksheet(
  workbook,
  "Expected counts check"
)

openxlsx::writeData(
  workbook,
  sheet = "Expected counts check",
  x = expected_summary
)

openxlsx::saveWorkbook(
  workbook,
  paste0(output_prefix, ".xlsx"),
  overwrite = TRUE
)

## ============================================================
## 11. Print summary
## ============================================================

final_tbl
summary_tbl
expected_summary
