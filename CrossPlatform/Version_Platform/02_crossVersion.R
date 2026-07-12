## =========================
## 1. Paths
## =========================
suppressPackageStartupMessages({
  library(Teddy)
  library(dplyr)
  library(tidyr)
})

project_dir <- "path/to/project"
work_dir <- file.path(project_dir, "annotation_version_comparison")

dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

gtf_vM7_teddy  <- "path/to/TEDDY_vM7.annotated.gtf"
gtf_vM25_teddy <- "path/to/TEDDY_vM25.annotated.gtf"
gtf_vM7_gene   <- "path/to/gencode.vM7.annotation.gtf"
gtf_vM25_gene  <- "path/to/gencode.vM25.annotation.gtf"

required_files <- c(
  gtf_vM7_teddy,
  gtf_vM25_teddy,
  gtf_vM7_gene,
  gtf_vM25_gene
)

stopifnot(all(file.exists(required_files)))

annotation_gtfs <- c(
  TEDDY_vM7    = gtf_vM7_teddy,
  TEDDY_vM25   = gtf_vM25_teddy,
  GENCODE_vM7  = gtf_vM7_gene,
  GENCODE_vM25 = gtf_vM25_gene
)

## =========================
## 2. Build a pooled annotation space
## =========================

pooled_gtf <- file.path(
  work_dir,
  "annotation_version_pool.gtf"
)

Teddy::stringtieMerge(
  reference = gtf_vM7_gene,
  gtfFiles = annotation_gtfs,
  outfile = pooled_gtf,
  params = "-p 16"
)

stopifnot(file.exists(pooled_gtf))

## Compare each annotation against the same pooled reference.
## buildIsoformSupport also generates the corresponding tracking files.

annotation_version_support <- Teddy::buildIsoformSupport(
  reference = pooled_gtf,
  gtffiles = annotation_gtfs,
  conditions = names(annotation_gtfs),
  out_dir = work_dir,
  out_prefix = "annotation_version_pool",
  cores = 16
)

## =========================
## 3. Parse gffcompare tracking files
## =========================

parse_tracking <- function(tracking_file, prefix) {
  x <- read.delim(
    tracking_file,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE,
    quote = ""
  )
  
  ref_raw <- as.character(x$V3)
  query_raw <- as.character(x$V5)
  
  pooled_gene_id <- sub("\\|.*$", "", ref_raw)
  pooled_tx_id <- sub("^.*\\|", "", ref_raw)
  
  matched <- !is.na(query_raw) &
    query_raw != "" &
    query_raw != "-"
  
  query_clean <- sub("^q[0-9]+:", "", query_raw)
  query_parts <- strsplit(query_clean, "\\|")
  
  get_part <- function(x, i) {
    if (length(x) >= i) x[i] else NA_character_
  }
  
  matched_gene_id <- vapply(
    query_parts,
    get_part,
    character(1),
    i = 1
  )
  
  matched_tx_id <- vapply(
    query_parts,
    get_part,
    character(1),
    i = 2
  )
  
  matched_gene_id[!matched] <- NA_character_
  matched_tx_id[!matched] <- NA_character_
  
  out <- data.frame(
    pooled_gene_id = pooled_gene_id,
    pooled_tx_id = pooled_tx_id,
    class_code = as.character(x$V4),
    matched_gene_id = matched_gene_id,
    matched_tx_id = matched_tx_id,
    stringsAsFactors = FALSE
  )
  
  out$exact <- !is.na(out$class_code) &
    out$class_code == "="
  
  out$traceable <- !is.na(out$class_code) &
    !out$class_code %in% c("u", "-")
  
  colnames(out)[3:7] <- paste0(
    prefix,
    "_",
    colnames(out)[3:7]
  )
  
  out
}

tracking_tables <- lapply(
  names(annotation_gtfs),
  function(prefix) {
    parse_tracking(
      file.path(
        work_dir,
        paste0(
          "annotation_version_pool_",
          prefix,
          ".tracking"
        )
      ),
      prefix
    )
  }
)

names(tracking_tables) <- names(annotation_gtfs)

## =========================
## 4. Collapse duplicate pooled mappings
## =========================

collapse_tracking <- function(df, prefix) {
  class_col <- paste0(prefix, "_class_code")
  gene_col <- paste0(prefix, "_matched_gene_id")
  tx_col <- paste0(prefix, "_matched_tx_id")
  exact_col <- paste0(prefix, "_exact")
  trace_col <- paste0(prefix, "_traceable")
  
  collapse_character <- function(x) {
    x <- unique(x[!is.na(x) & x != "" & x != "-"])
    
    if (length(x) == 0L) {
      NA_character_
    } else {
      paste(x, collapse = ";")
    }
  }
  
  df |>
    dplyr::group_by(pooled_gene_id, pooled_tx_id) |>
    dplyr::group_modify(function(.x, .y) {
      exact_rows <- tidyr::replace_na(.x[[exact_col]], FALSE)
      traceable_rows <- tidyr::replace_na(.x[[trace_col]], FALSE)
      
      keep_rows <- if (any(exact_rows)) {
        exact_rows
      } else if (any(traceable_rows)) {
        traceable_rows
      } else {
        rep(TRUE, nrow(.x))
      }
      
      data.frame(
        exact = any(exact_rows),
        traceable = any(traceable_rows),
        class_code = collapse_character(
          .x[[class_col]][keep_rows]
        ),
        matched_gene_id = collapse_character(
          .x[[gene_col]][keep_rows]
        ),
        matched_tx_id = collapse_character(
          .x[[tx_col]][keep_rows]
        ),
        stringsAsFactors = FALSE
      )
    }) |>
    dplyr::ungroup() |>
    dplyr::rename_with(
      ~ paste0(prefix, "_", .x),
      c(
        exact,
        traceable,
        class_code,
        matched_gene_id,
        matched_tx_id
      )
    )
}

tracking_collapsed <- lapply(
  names(tracking_tables),
  function(prefix) {
    collapse_tracking(
      tracking_tables[[prefix]],
      prefix
    )
  }
)

names(tracking_collapsed) <- names(tracking_tables)

## =========================
## 5. Build pooled-level comparison table
## =========================

annotation_version_final <- annotation_version_support |>
  dplyr::rename(
    pooled_gene_id = gene_id,
    pooled_tx_id = tx_id,
    TEDDY_vM7_exact_raw = TEDDY_vM7,
    TEDDY_vM25_exact_raw = TEDDY_vM25,
    GENCODE_vM7_exact_raw = GENCODE_vM7,
    GENCODE_vM25_exact_raw = GENCODE_vM25
  )

for (prefix in names(tracking_collapsed)) {
  annotation_version_final <- annotation_version_final |>
    dplyr::left_join(
      tracking_collapsed[[prefix]],
      by = c("pooled_gene_id", "pooled_tx_id")
    )
}

annotation_version_final <- annotation_version_final |>
  dplyr::mutate(
    dplyr::across(
      dplyr::ends_with("_exact"),
      ~ tidyr::replace_na(.x, FALSE)
    ),
    dplyr::across(
      dplyr::ends_with("_traceable"),
      ~ tidyr::replace_na(.x, FALSE)
    ),
    
    TEDDY_vM7_exact =
      TEDDY_vM7_exact_raw | TEDDY_vM7_exact,
    
    TEDDY_vM25_exact =
      TEDDY_vM25_exact_raw | TEDDY_vM25_exact,
    
    GENCODE_vM7_exact =
      GENCODE_vM7_exact_raw | GENCODE_vM7_exact,
    
    GENCODE_vM25_exact =
      GENCODE_vM25_exact_raw | GENCODE_vM25_exact,
    
    TEDDY_vM7_traceable =
      TEDDY_vM7_traceable | TEDDY_vM7_exact,
    
    TEDDY_vM25_traceable =
      TEDDY_vM25_traceable | TEDDY_vM25_exact,
    
    GENCODE_vM7_traceable =
      GENCODE_vM7_traceable | GENCODE_vM7_exact,
    
    GENCODE_vM25_traceable =
      GENCODE_vM25_traceable | GENCODE_vM25_exact
  )

stopifnot(
  nrow(annotation_version_final) ==
    nrow(annotation_version_support)
)

## =========================
## 6. Map each TEDDY vM7 transcript to pooled space
## =========================

vM7_tracking <- tracking_tables[["TEDDY_vM7"]]

vM7_TEDDY_annotation_version_expanded <- vM7_tracking |>
  dplyr::filter(
    !is.na(TEDDY_vM7_matched_tx_id),
    TEDDY_vM7_matched_tx_id != "",
    TEDDY_vM7_matched_tx_id != "-"
  ) |>
  dplyr::transmute(
    TEDDY_vM7_gene_id = TEDDY_vM7_matched_gene_id,
    TEDDY_vM7_tx_id = TEDDY_vM7_matched_tx_id,
    pooled_gene_id,
    pooled_tx_id,
    TEDDY_vM7_class_code_to_pool =
      TEDDY_vM7_class_code,
    TEDDY_vM7_exact_to_pool =
      TEDDY_vM7_class_code == "=",
    TEDDY_vM7_traceable_to_pool =
      TEDDY_vM7_traceable
  ) |>
  dplyr::left_join(
    annotation_version_final |>
      dplyr::select(
        pooled_gene_id,
        pooled_tx_id,
        
        GENCODE_vM7_exact,
        GENCODE_vM7_traceable,
        GENCODE_vM7_class_code,
        GENCODE_vM7_matched_tx_id,
        
        TEDDY_vM25_exact,
        TEDDY_vM25_traceable,
        TEDDY_vM25_class_code,
        TEDDY_vM25_matched_tx_id,
        
        GENCODE_vM25_exact,
        GENCODE_vM25_traceable,
        GENCODE_vM25_class_code,
        GENCODE_vM25_matched_tx_id
      ),
    by = c("pooled_gene_id", "pooled_tx_id")
  ) |>
  dplyr::mutate(
    mapping_priority = dplyr::case_when(
      TEDDY_vM7_exact_to_pool &
        GENCODE_vM25_exact ~ 1L,
      
      TEDDY_vM7_exact_to_pool &
        GENCODE_vM25_traceable ~ 2L,
      
      TEDDY_vM7_exact_to_pool &
        TEDDY_vM25_exact ~ 3L,
      
      TEDDY_vM7_exact_to_pool &
        TEDDY_vM25_traceable ~ 4L,
      
      TEDDY_vM7_exact_to_pool ~ 5L,
      TEDDY_vM7_traceable_to_pool ~ 6L,
      TRUE ~ 7L
    )
  ) |>
  dplyr::group_by(TEDDY_vM7_tx_id) |>
  dplyr::arrange(
    mapping_priority,
    pooled_gene_id,
    pooled_tx_id,
    .by_group = TRUE
  ) |>
  dplyr::mutate(
    n_pooled_mapping = dplyr::n(),
    mapping_rank = dplyr::row_number(),
    is_best_mapping = mapping_rank == 1L
  ) |>
  dplyr::ungroup()

vM7_to_vM25_annotation_support <-
  vM7_TEDDY_annotation_version_expanded |>
  dplyr::filter(is_best_mapping)

stopifnot(
  dplyr::n_distinct(
    vM7_to_vM25_annotation_support$TEDDY_vM7_tx_id
  ) ==
    nrow(vM7_to_vM25_annotation_support)
)

## =========================
## 7. Save outputs
## =========================

saveRDS(
  annotation_version_final,
  file.path(
    work_dir,
    "annotation_version_pooled_comparison.rds"
  )
)

saveRDS(
  vM7_TEDDY_annotation_version_expanded,
  file.path(
    work_dir,
    "vM7_TEDDY_annotation_version_expanded.rds"
  )
)

saveRDS(
  vM7_to_vM25_annotation_support,
  file.path(
    work_dir,
    "vM7_to_vM25_annotation_support.rds"
  )
)

write.table(
  annotation_version_final,
  file = file.path(
    work_dir,
    "annotation_version_pooled_comparison.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Annotation-version comparison completed.")
