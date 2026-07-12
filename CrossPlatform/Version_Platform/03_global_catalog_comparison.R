#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(SummarizedExperiment)
})

## =========================
## 1. Paths
## =========================

project_dir <- "path/to/project"

annotation_best_file <- file.path(
  project_dir,
  "annotation_version_comparison",
  "vM7_to_vM25_annotation_support.rds"
)

combineSE_vM7_file <- "path/to/combineSE_vM7.rds"
chi_GTF_vM7_file <- "path/to/chi_GTF_vM7.rds"

output_dir <- file.path(
  project_dir,
  "global_catalog_comparison"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

stopifnot(
  file.exists(annotation_best_file),
  file.exists(combineSE_vM7_file),
  file.exists(chi_GTF_vM7_file)
)

## =========================
## 2. Load inputs
## =========================

annotation_best <- readRDS(annotation_best_file)
combineSE_vM7 <- readRDS(combineSE_vM7_file)
chi_GTF_vM7 <- readRDS(chi_GTF_vM7_file)

fpkm_vM7 <- SummarizedExperiment::assay(
  combineSE_vM7,
  "FPKM"
)

## =========================
## 3. Collapse redundant transcript structures
## =========================

CollapseRedundantTx <- function(
    chi_GTF,
    expr_mat = NULL,
    tx_keep = NULL,
    min_fpkm = 0.3,
    min_fraction = NULL,
    mode = c("tss", "structure"),
    scope = c("all", "first_exon_LINE1", "first_exon_TE"),
    terminal_tolerance = 15,
    internal_tolerance = 0,
    include_first_exon_TE = TRUE,
    representative_by = c("expression", "annotated", "longest", "tx_id")
) {
  mode <- match.arg(mode)
  scope <- match.arg(scope)
  representative_by <- match.arg(representative_by)
  
  make_te_signature <- function(TE_name, is_TE) {
    x <- TE_name[is_TE]
    x <- x[!is.na(x) & x != "" & x != "none"]
    
    if (length(x) == 0) {
      "no_TE"
    } else {
      paste(sort(unique(x)), collapse = ";")
    }
  }
  
  gtf_df <- as.data.frame(chi_GTF) |>
    dplyr::mutate(
      seqnames = as.character(seqnames),
      strand = as.character(strand),
      transcript_id = as.character(transcript_id),
      gene_id = as.character(gene_id),
      gene_name = as.character(gene_name),
      TE_name = ifelse(
        is.na(TE_name),
        "",
        as.character(TE_name)
      ),
      TE_class = ifelse(
        is.na(TE_class),
        "",
        as.character(TE_class)
      ),
      TE_family = ifelse(
        is.na(TE_family),
        "",
        as.character(TE_family)
      ),
      tx_exon_rank = as.integer(tx_exon_rank),
      is_TE = TE_name != "" & TE_name != "none",
      is_LINE1 =
        stringr::str_detect(TE_family, "L1") |
        stringr::str_detect(TE_name, "(^|[,;])L1") |
        stringr::str_detect(TE_class, "LINE"),
      is_novel_tx = stringr::str_detect(
        transcript_id,
        "^MSTRG"
      )
    )
  
  all_tx <- unique(gtf_df$transcript_id)
  all_tx <- all_tx[
    !is.na(all_tx) &
      all_tx != ""
  ]
  
  expr_score <- tibble::tibble(
    transcript_id = all_tx,
    expr_score = NA_real_
  )
  
  if (!is.null(expr_mat)) {
    tx_in_expr <- intersect(
      all_tx,
      rownames(expr_mat)
    )
    
    expr_score <- tibble::tibble(
      transcript_id = tx_in_expr,
      expr_score = rowMeans(
        log2(
          expr_mat[
            tx_in_expr,
            ,
            drop = FALSE
          ] + 1
        ),
        na.rm = TRUE
      )
    )
    
    if (is.null(tx_keep)) {
      min_samples <- if (is.null(min_fraction)) {
        1L
      } else {
        ceiling(
          min_fraction * ncol(expr_mat)
        )
      }
      
      tx_keep <- tx_in_expr[
        rowSums(
          expr_mat[
            tx_in_expr,
            ,
            drop = FALSE
          ] > min_fpkm,
          na.rm = TRUE
        ) >= min_samples
      ]
    }
  }
  
  if (is.null(tx_keep)) {
    tx_keep <- all_tx
  }
  
  gtf_df <- gtf_df |>
    dplyr::filter(
      transcript_id %in% tx_keep
    )
  
  exon_unique <- gtf_df |>
    dplyr::distinct(
      transcript_id,
      gene_id,
      gene_name,
      seqnames,
      strand,
      tx_exon_rank,
      start,
      end,
      .keep_all = TRUE
    )
  
  tx_info <- exon_unique |>
    dplyr::group_by(
      transcript_id,
      gene_id,
      gene_name,
      seqnames,
      strand
    ) |>
    dplyr::summarise(
      n_exons = dplyr::n_distinct(
        tx_exon_rank
      ),
      total_exon_width = sum(
        end - start + 1,
        na.rm = TRUE
      ),
      tss = ifelse(
        dplyr::first(strand) == "+",
        min(
          start[
            tx_exon_rank ==
              min(tx_exon_rank, na.rm = TRUE)
          ],
          na.rm = TRUE
        ),
        max(
          end[
            tx_exon_rank ==
              min(tx_exon_rank, na.rm = TRUE)
          ],
          na.rm = TRUE
        )
      ),
      tes = ifelse(
        dplyr::first(strand) == "+",
        max(
          end[
            tx_exon_rank ==
              max(tx_exon_rank, na.rm = TRUE)
          ],
          na.rm = TRUE
        ),
        min(
          start[
            tx_exon_rank ==
              max(tx_exon_rank, na.rm = TRUE)
          ],
          na.rm = TRUE
        )
      ),
      .groups = "drop"
    ) |>
    dplyr::left_join(
      expr_score,
      by = "transcript_id"
    ) |>
    dplyr::mutate(
      is_annotated = !stringr::str_detect(
        transcript_id,
        "^MSTRG"
      ),
      expr_score = ifelse(
        is.na(expr_score),
        -Inf,
        expr_score
      )
    )
  
  if (mode == "tss") {
    collapse_input <- gtf_df |>
      dplyr::filter(
        tx_exon_rank == 1
      ) |>
      dplyr::group_by(
        transcript_id,
        gene_id,
        gene_name,
        seqnames,
        strand
      ) |>
      dplyr::summarise(
        tss = ifelse(
          dplyr::first(strand) == "+",
          min(start, na.rm = TRUE),
          max(end, na.rm = TRUE)
        ),
        first_exon_TE_signature = if (
          include_first_exon_TE
        ) {
          make_te_signature(
            TE_name,
            is_TE
          )
        } else {
          "all"
        },
        is_LINE1 = any(is_LINE1),
        is_TE = any(is_TE),
        .groups = "drop"
      )
    
    if (scope == "first_exon_LINE1") {
      collapse_input <- collapse_input |>
        dplyr::filter(is_LINE1)
    }
    
    if (scope == "first_exon_TE") {
      collapse_input <- collapse_input |>
        dplyr::filter(is_TE)
    }
    
    collapse_input <- collapse_input |>
      dplyr::mutate(
        feature_key =
          first_exon_TE_signature
      ) |>
      dplyr::group_by(
        gene_id,
        gene_name,
        seqnames,
        strand,
        feature_key
      ) |>
      dplyr::arrange(
        tss,
        transcript_id,
        .by_group = TRUE
      ) |>
      dplyr::mutate(
        cluster_index = cumsum(
          dplyr::if_else(
            dplyr::row_number() == 1 |
              (tss - dplyr::lag(tss)) >
              terminal_tolerance,
            1L,
            0L
          )
        ),
        cluster_id = paste(
          "tss",
          gene_id,
          seqnames,
          strand,
          feature_key,
          cluster_index,
          sep = "|"
        )
      ) |>
      dplyr::ungroup() |>
      dplyr::select(
        transcript_id,
        gene_id,
        gene_name,
        seqnames,
        strand,
        tss,
        first_exon_TE_signature,
        feature_key,
        cluster_id
      )
  }
  
  if (mode == "structure") {
    intron_key <- exon_unique |>
      dplyr::arrange(
        transcript_id,
        tx_exon_rank
      ) |>
      dplyr::group_by(transcript_id) |>
      dplyr::mutate(
        next_start = dplyr::lead(start),
        next_end = dplyr::lead(end),
        intron_left =
          pmin(
            end,
            next_end,
            na.rm = TRUE
          ) + 1,
        intron_right =
          pmax(
            start,
            next_start,
            na.rm = TRUE
          ) - 1
      ) |>
      dplyr::filter(
        !is.na(next_start),
        !is.na(next_end)
      ) |>
      dplyr::mutate(
        intron_left_key = ifelse(
          internal_tolerance > 0,
          round(
            intron_left /
              internal_tolerance
          ) * internal_tolerance,
          intron_left
        ),
        intron_right_key = ifelse(
          internal_tolerance > 0,
          round(
            intron_right /
              internal_tolerance
          ) * internal_tolerance,
          intron_right
        ),
        intron_piece = paste(
          intron_left_key,
          intron_right_key,
          sep = "-"
        )
      ) |>
      dplyr::summarise(
        intron_chain_key = paste(
          intron_piece,
          collapse = "|"
        ),
        .groups = "drop"
      )
    
    first_exon_key <- gtf_df |>
      dplyr::filter(
        tx_exon_rank == 1
      ) |>
      dplyr::group_by(transcript_id) |>
      dplyr::summarise(
        first_exon_TE_signature = if (
          include_first_exon_TE
        ) {
          make_te_signature(
            TE_name,
            is_TE
          )
        } else {
          "all"
        },
        has_first_exon_LINE1 =
          any(is_LINE1),
        has_first_exon_TE =
          any(is_TE),
        .groups = "drop"
      )
    
    collapse_input <- tx_info |>
      dplyr::left_join(
        intron_key,
        by = "transcript_id"
      ) |>
      dplyr::left_join(
        first_exon_key,
        by = "transcript_id"
      ) |>
      dplyr::mutate(
        intron_chain_key = ifelse(
          is.na(intron_chain_key),
          "single_exon",
          intron_chain_key
        ),
        first_exon_TE_signature = ifelse(
          is.na(first_exon_TE_signature),
          "no_TE",
          first_exon_TE_signature
        ),
        has_first_exon_LINE1 = ifelse(
          is.na(has_first_exon_LINE1),
          FALSE,
          has_first_exon_LINE1
        ),
        has_first_exon_TE = ifelse(
          is.na(has_first_exon_TE),
          FALSE,
          has_first_exon_TE
        )
      )
    
    if (scope == "first_exon_LINE1") {
      collapse_input <- collapse_input |>
        dplyr::filter(
          has_first_exon_LINE1
        )
    }
    
    if (scope == "first_exon_TE") {
      collapse_input <- collapse_input |>
        dplyr::filter(
          has_first_exon_TE
        )
    }
    
    collapse_input <- collapse_input |>
      dplyr::mutate(
        feature_key =
          first_exon_TE_signature
      ) |>
      dplyr::group_by(
        gene_id,
        gene_name,
        seqnames,
        strand,
        n_exons,
        intron_chain_key,
        feature_key
      ) |>
      dplyr::arrange(
        tss,
        tes,
        transcript_id,
        .by_group = TRUE
      ) |>
      dplyr::mutate(
        tss_cluster = cumsum(
          dplyr::if_else(
            dplyr::row_number() == 1 |
              abs(
                tss -
                  dplyr::lag(tss)
              ) > terminal_tolerance,
            1L,
            0L
          )
        )
      ) |>
      dplyr::group_by(
        gene_id,
        gene_name,
        seqnames,
        strand,
        n_exons,
        intron_chain_key,
        feature_key,
        tss_cluster
      ) |>
      dplyr::arrange(
        tes,
        transcript_id,
        .by_group = TRUE
      ) |>
      dplyr::mutate(
        tes_cluster = cumsum(
          dplyr::if_else(
            dplyr::row_number() == 1 |
              abs(
                tes -
                  dplyr::lag(tes)
              ) > terminal_tolerance,
            1L,
            0L
          )
        ),
        cluster_id = paste(
          "structure",
          gene_id,
          seqnames,
          strand,
          n_exons,
          intron_chain_key,
          feature_key,
          tss_cluster,
          tes_cluster,
          sep = "|"
        )
      ) |>
      dplyr::ungroup() |>
      dplyr::select(
        transcript_id,
        gene_id,
        gene_name,
        seqnames,
        strand,
        tss,
        first_exon_TE_signature,
        feature_key,
        cluster_id
      )
  }
  
  if (nrow(collapse_input) == 0) {
    return(
      list(
        collapsed_gtf = gtf_df[0, ],
        tx_map = tibble::tibble(),
        cluster_summary = tibble::tibble(),
        parameters = list(
          mode = mode,
          scope = scope,
          terminal_tolerance =
            terminal_tolerance,
          internal_tolerance =
            internal_tolerance,
          min_fpkm = min_fpkm,
          min_fraction = min_fraction,
          representative_by =
            representative_by
        )
      )
    )
  }
  
  rep_candidates <- collapse_input |>
    dplyr::left_join(
      tx_info |>
        dplyr::select(
          transcript_id,
          expr_score,
          is_annotated,
          total_exon_width,
          n_exons
        ),
      by = "transcript_id"
    ) |>
    dplyr::mutate(
      representative_score =
        dplyr::case_when(
          representative_by ==
            "expression" ~ expr_score,
          representative_by ==
            "annotated" ~
            as.numeric(is_annotated),
          representative_by ==
            "longest" ~ total_exon_width,
          TRUE ~ 0
        )
    )
  
  rep_tbl <- rep_candidates |>
    dplyr::group_by(cluster_id) |>
    dplyr::arrange(
      dplyr::desc(
        representative_score
      ),
      dplyr::desc(expr_score),
      dplyr::desc(is_annotated),
      dplyr::desc(total_exon_width),
      dplyr::desc(n_exons),
      transcript_id,
      .by_group = TRUE
    ) |>
    dplyr::slice(1) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      cluster_id,
      representative_tx =
        transcript_id
    )
  
  tx_map <- rep_candidates |>
    dplyr::left_join(
      rep_tbl,
      by = "cluster_id"
    ) |>
    dplyr::mutate(
      is_representative =
        transcript_id ==
        representative_tx
    )
  
  cluster_summary <- tx_map |>
    dplyr::group_by(
      cluster_id,
      representative_tx,
      gene_id,
      gene_name,
      seqnames,
      strand,
      first_exon_TE_signature,
      feature_key
    ) |>
    dplyr::summarise(
      n_members =
        dplyr::n_distinct(
          transcript_id
        ),
      n_novel_MSTRG =
        dplyr::n_distinct(
          transcript_id[
            stringr::str_detect(
              transcript_id,
              "^MSTRG"
            )
          ]
        ),
      n_annotated =
        dplyr::n_distinct(
          transcript_id[
            !stringr::str_detect(
              transcript_id,
              "^MSTRG"
            )
          ]
        ),
      has_novel_MSTRG =
        n_novel_MSTRG > 0,
      has_annotated =
        n_annotated > 0,
      cluster_class =
        dplyr::case_when(
          has_novel_MSTRG &
            has_annotated ~
            "mixed_MSTRG_ENST",
          has_novel_MSTRG &
            !has_annotated ~
            "novel_MSTRG_only",
          !has_novel_MSTRG &
            has_annotated ~
            "annotated_ENST_only",
          TRUE ~ "unknown"
        ),
      .groups = "drop"
    )
  
  collapsed_gtf <- gtf_df |>
    dplyr::semi_join(
      rep_tbl |>
        dplyr::transmute(
          transcript_id =
            representative_tx
        ),
      by = "transcript_id"
    )
  
  list(
    collapsed_gtf = collapsed_gtf,
    tx_map = tx_map,
    cluster_summary = cluster_summary,
    parameters = list(
      mode = mode,
      scope = scope,
      terminal_tolerance =
        terminal_tolerance,
      internal_tolerance =
        internal_tolerance,
      min_fpkm = min_fpkm,
      min_fraction = min_fraction,
      representative_by =
        representative_by
    )
  )
}

## =========================
## 4. Expression summary
## =========================

expr_vM7_tbl <- tibble::tibble(
  TEDDY_vM7_tx_id = rownames(fpkm_vM7),
  FPKM_max = apply(
    fpkm_vM7,
    1,
    max,
    na.rm = TRUE
  ),
  FPKM_mean = rowMeans(
    fpkm_vM7,
    na.rm = TRUE
  ),
  n_sample_FPKM_gt1 = rowSums(
    fpkm_vM7 > 1,
    na.rm = TRUE
  )
)

## =========================
## 5. Select expressed TEDDY-novel
##    TE-chimeric MSTRG transcripts
## =========================

chi_tx_tbl <- as.data.frame(
  chi_GTF_vM7
) |>
  dplyr::transmute(
    TEDDY_vM7_tx_id =
      as.character(transcript_id),
    gene_id =
      as.character(gene_id),
    gene_name =
      as.character(gene_name)
  ) |>
  dplyr::filter(
    !is.na(TEDDY_vM7_tx_id),
    TEDDY_vM7_tx_id != ""
  ) |>
  dplyr::distinct()

novel_expr_tbl <- annotation_best |>
  dplyr::left_join(
    expr_vM7_tbl,
    by = "TEDDY_vM7_tx_id"
  ) |>
  dplyr::inner_join(
    chi_tx_tbl,
    by = "TEDDY_vM7_tx_id"
  ) |>
  dplyr::filter(
    grepl(
      "^MSTRG",
      TEDDY_vM7_tx_id
    ),
    !tidyr::replace_na(
      GENCODE_vM7_exact,
      FALSE
    ),
    FPKM_max > 1
  ) |>
  dplyr::distinct(
    TEDDY_vM7_tx_id,
    .keep_all = TRUE
  )

novel_expr_tx <- unique(
  novel_expr_tbl$TEDDY_vM7_tx_id
)

## =========================
## 6. Collapse redundant structures
## =========================

vM7_novel_expr_collapse <- CollapseRedundantTx(
  chi_GTF = chi_GTF_vM7,
  expr_mat = fpkm_vM7,
  tx_keep = novel_expr_tx,
  mode = "structure",
  scope = "all",
  terminal_tolerance = 15,
  internal_tolerance = 0,
  representative_by = "expression"
)

representative_tx <- vM7_novel_expr_collapse$tx_map |>
  dplyr::filter(
    is_representative
  ) |>
  dplyr::pull(
    representative_tx
  ) |>
  unique()

## =========================
## 7. High-confidence representative set
## =========================

global_catalog_table <- novel_expr_tbl |>
  dplyr::filter(
    TEDDY_vM7_tx_id %in%
      representative_tx,
    FPKM_max > 5
  ) |>
  dplyr::mutate(
    GENCODE_vM25_matched =
      dplyr::if_else(
        !is.na(
          GENCODE_vM25_matched_tx_id
        ) &
          GENCODE_vM25_matched_tx_id != "" &
          GENCODE_vM25_matched_tx_id != "-",
        "Yes",
        "-"
      ),
    
    GENCODE_vM25_exact =
      dplyr::if_else(
        tidyr::replace_na(
          GENCODE_vM25_class_code == "=",
          FALSE
        ),
        "Yes",
        "-"
      ),
    
    Recovered_by_TEDDY_vM25 =
      dplyr::if_else(
        tidyr::replace_na(
          TEDDY_vM25_traceable,
          FALSE
        ),
        "Yes",
        "-"
      )
  ) |>
  dplyr::select(
    TEDDY_vM7_tx_id,
    gene_id,
    gene_name,
    FPKM_max,
    FPKM_mean,
    n_sample_FPKM_gt1,
    
    GENCODE_vM25_matched_tx_id,
    GENCODE_vM25_class_code,
    GENCODE_vM25_matched,
    GENCODE_vM25_exact,
    
    TEDDY_vM25_matched_tx_id,
    TEDDY_vM25_class_code,
    Recovered_by_TEDDY_vM25,
    
    pooled_gene_id,
    pooled_tx_id
  ) |>
  dplyr::arrange(
    gene_name,
    TEDDY_vM7_tx_id
  )

## =========================
## 8. Summary
## =========================

global_catalog_summary <- global_catalog_table |>
  dplyr::summarise(
    category = paste(
      "High-confidence representative",
      "vM7 TEDDY-novel TE-chimeric",
      "MSTRG structures with FPKM_max > 5"
    ),
    
    n_total = dplyr::n(),
    
    n_GENCODE_vM25_matched =
      sum(
        GENCODE_vM25_matched == "Yes"
      ),
    
    pct_GENCODE_vM25_matched =
      round(
        100 *
          n_GENCODE_vM25_matched /
          n_total,
        2
      ),
    
    n_GENCODE_vM25_exact =
      sum(
        GENCODE_vM25_exact == "Yes"
      ),
    
    pct_GENCODE_vM25_exact =
      round(
        100 *
          n_GENCODE_vM25_exact /
          n_total,
        2
      ),
    
    n_TEDDY_vM25_recovered =
      sum(
        Recovered_by_TEDDY_vM25 ==
          "Yes"
      ),
    
    pct_TEDDY_vM25_recovered =
      round(
        100 *
          n_TEDDY_vM25_recovered /
          n_total,
        2
      )
  )

global_catalog_summary

## Values reported in the response letter
stopifnot(
  global_catalog_summary$n_total == 4521,
  global_catalog_summary$n_GENCODE_vM25_matched == 220,
  global_catalog_summary$n_GENCODE_vM25_exact == 86,
  global_catalog_summary$n_TEDDY_vM25_recovered == 4395
)

## =========================
## 9. Save outputs
## =========================

readr::write_tsv(
  global_catalog_table,
  file.path(
    output_dir,
    "global_representative_annotation_version_table.tsv"
  )
)

readr::write_tsv(
  global_catalog_summary,
  file.path(
    output_dir,
    "global_representative_annotation_version_summary.tsv"
  )
)

readr::write_tsv(
  vM7_novel_expr_collapse$tx_map,
  file.path(
    output_dir,
    "redundant_transcript_collapse_map.tsv"
  )
)

readr::write_tsv(
  vM7_novel_expr_collapse$cluster_summary,
  file.path(
    output_dir,
    "redundant_transcript_cluster_summary.tsv"
  )
)

saveRDS(
  global_catalog_table,
  file.path(
    output_dir,
    "global_representative_annotation_version_table.rds"
  )
)

saveRDS(
  vM7_novel_expr_collapse,
  file.path(
    output_dir,
    "vM7_novel_TE_chimeric_transcript_collapse.rds"
  )
)

message("Global annotation-version comparison completed.")
