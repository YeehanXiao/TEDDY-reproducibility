#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(dplyr)
  library(purrr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript score_seed_benchmarks.R ROOT")
root <- normalizePath(args[1], mustWork = TRUE)
depths <- c("5x", "10x", "25x", "50x", "100x")

scope <- c(
  TEDDY = "expressed_TE_chimeric",
  Arriba = "expressed_TE_chimeric",
  ChimeraTE = "expressed_TE_chimeric",
  FREDY = "expressed_TE50_unstranded",
  LIONS = "expressed_5prime_exon_TE",
  TEProf2 = "expressed_5prime_exon_TE"
)
setting <- c(
  TEDDY = "original",
  Arriba = "final+discarded;support>=1",
  ChimeraTE = "overlap=0.10;fwd+rf_union;per_sample",
  FREDY = "all_transcript_TE50_preCDS",
  LIONS = "original",
  TEProf2 = "simulation_dictionary_adapter"
)

need <- function(f, nonempty = TRUE) {
  if (!file.exists(f) || isTRUE(file.info(f)$isdir) ||
      (nonempty && file.info(f)$size == 0)) {
    stop("Missing/empty input: ", f)
  }
  f
}

read_tab <- function(f, header = TRUE) {
  need(f, FALSE)
  if (file.info(f)$size == 0) return(tibble())
  x <- read_tsv(
    f, col_names = header, col_types = cols(.default = col_character()),
    quote = "", progress = FALSE, show_col_types = FALSE
  )
  if (nrow(problems(x))) stop("Malformed TSV: ", f)
  x
}

clean <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "None", ".")] <- NA_character_
  x
}
genes <- function(x) sort(unique(na.omit(clean(x))))

require_cols <- function(x, cols, label) {
  if (!all(cols %in% names(x))) {
    stop(label, ": missing columns: ",
         paste(setdiff(cols, names(x)), collapse = ", "))
  }
}

score_seed <- function(seed) {
  message("[READ] seed", seed)
  d <- file.path(root, paste0("seed_", seed), "results")

  tx <- read_tab(file.path(d, "official_simulated_transcript_truth_status.tsv"))
  ex <- read_tab(file.path(d, "official_simulated_exon_truth_status.tsv"))
  iso <- read_tab(file.path(d, "official_simulated_1000.isoforms.results"))

  require_cols(
    tx, c("transcript_id", "gene_id", "gene_name", "is_TE_chimeric_tx"),
    "transcript truth"
  )
  require_cols(
    ex, c("transcript_id", "gene_name", "seqnames", "start", "end",
          "strand", "is_TE_overlap_exon"),
    "exon truth"
  )
  require_cols(iso, c("transcript_id", "TPM"), "isoforms")

  stopifnot(
    !anyDuplicated(tx$transcript_id),
    !anyDuplicated(iso$transcript_id),
    !anyNA(as.numeric(iso$TPM)),
    setequal(tx$transcript_id, iso$transcript_id)
  )

  expressed <- iso$transcript_id[as.numeric(iso$TPM) > 1]
  full <- genes(tx$gene_name[
    tx$is_TE_chimeric_tx == "TRUE" & tx$transcript_id %in% expressed
  ])
  universe <- genes(tx$gene_name)

  ex <- ex |> mutate(start = as.integer(start), end = as.integer(end))
  stopifnot(!anyNA(ex$start), !anyNA(ex$end), all(ex$start <= ex$end))

  # Find the first exon using all truth exons, then test TE overlap.
  first <- ex |>
    group_by(transcript_id) |>
    filter(
      (strand == "+" & start == min(start)) |
      (strand == "-" & end == max(end))
    ) |>
    ungroup() |>
    filter(is_TE_overlap_exon == "TRUE", transcript_id %in% expressed)
  first <- genes(first$gene_name)

  gtf <- import(need(file.path(
    d, "official_simulated_truth_1000_transcripts.gtf"
  )))
  gx <- gtf[gtf$type == "exon"]
  stopifnot(setequal(as.character(gx$transcript_id), tx$transcript_id))

  te <- readRDS(need(file.path(
    root, paste0("seed_", seed), "data/mm10_TE.rds"
  )))
  h <- findOverlaps(gx, te, ignore.strand = TRUE)
  frac <- width(pintersect(
    gx[queryHits(h)], te[subjectHits(h)], ignore.strand = TRUE
  )) / width(te)[subjectHits(h)]
  hit <- queryHits(h)[frac >= 0.5]
  fredy_truth <- genes(
    gx$gene_name[hit][gx$transcript_id[hit] %in% expressed]
  )

  truths <- list(
    expressed_TE_chimeric = full,
    expressed_5prime_exon_TE = first,
    expressed_TE50_unstranded = fredy_truth
  )
  stopifnot(length(full) > 0, length(first) > 0, length(fredy_truth) > 0)

  idmap <- tx |> distinct(gene_id, gene_name)
  stopifnot(!anyDuplicated(idmap$gene_id))

  # Retain the validated FREDY exon-overlap to host-gene mapping.
  # Report unmapped transcripts and ties for the largest overlap.
  map_fredy <- function(f) {
    p <- import(need(f))
    p <- p[p$type == "exon"]
    all_tx <- unique(as.character(p$transcript_id))
    h <- findOverlaps(p, gx, ignore.strand = FALSE)

    m <- tibble(
      pred_tx = as.character(p$transcript_id[queryHits(h)]),
      gene_name = clean(gx$gene_name[subjectHits(h)]),
      bp = width(pintersect(
        p[queryHits(h)], gx[subjectHits(h)], ignore.strand = FALSE
      ))
    ) |>
      filter(!is.na(gene_name)) |>
      group_by(pred_tx, gene_name) |>
      summarise(bp = sum(bp), .groups = "drop") |>
      group_by(pred_tx) |>
      mutate(tied_best = sum(bp == max(bp)) > 1L) |>
      arrange(desc(bp), gene_name, .by_group = TRUE) |>
      slice_head(n = 1) |>
      ungroup()

    tibble(pred_tx = all_tx) |> left_join(m, by = "pred_tx")
  }

  map_dfr(depths, function(depth) {
    tag <- paste0("official_simulated_", depth, "_noise0.1")

    teddy <- as.data.frame(readRDS(need(file.path(
      d, "TEDDY_official_by_depth", depth, "meta",
      paste0("official_chi_GTF_", depth, ".rds")
    ))))
    require_cols(teddy, "gene_name", "TEDDY")

    # Retain the original relaxed-main Arriba candidate definition.
    arriba <- map(c("", ".discarded"), function(suffix) {
      f <- file.path(
        d, "arriba_merge",
        paste0("merge_", depth, "_arriba_fusions", suffix, ".tsv")
      )
      x <- read_tab(need(f))
      names(x) <- sub("^#", "", names(x))
      require_cols(
        x, c("gene1", "gene2", "split_reads1",
             "split_reads2", "discordant_mates"), f
      )

      support <- as.numeric(x$split_reads1) +
        as.numeric(x$split_reads2) + as.numeric(x$discordant_mates)
      if (anyNA(support)) stop("Invalid Arriba support: ", f)

      x <- x[support >= 1, ]
      g <- unlist(strsplit(c(x$gene1, x$gene2), ",", fixed = TRUE))
      intersect(genes(sub("\\(.*\\)$", "", trimws(g))), universe)
    }) |> unlist(use.names = FALSE) |> genes()

    # Read only the six formal ChimeraTE sample TSVs per depth; empty files are valid.
    chim_ids <- map(c("fwd-stranded", "rf-stranded"), function(strand) {
      map(c("initiated", "terminated", "exonized"), function(kind) {
        f <- file.path(
          d, "ChimeraTE_official_by_depth", paste0(depth, "_", strand),
          paste0("rep1_", depth),
          paste0("TE-", kind, "-rep1_", depth, ".tsv")
        )
        x <- read_tab(f, FALSE)
        if (nrow(x) && ncol(x) < 6L) {
          stop("Unexpected ChimeraTE columns: ", f)
        }
        if (nrow(x)) as.character(x[[1]]) else character()
      }) |> unlist(use.names = FALSE)
    }) |> unlist(use.names = FALSE) |> genes()

    chim <- clean(idmap$gene_name[match(chim_ids, idmap$gene_id)])
    chim[is.na(chim)] <- paste0(
      "UNMAPPED_GENE_ID:", chim_ids[is.na(chim)]
    )

    lf <- need(file.path(
      d, "LIONS_official_by_depth", depth, paste0(tag, ".LIONS.log")
    ))
    lines <- readLines(lf, warn = FALSE)
    if (any(grepl("Traceback|Error:|ERROR:|FATAL|Killed|Aborted|Disk quota exceeded|No space left", lines)) ||
        !any(grepl("Chimeric processing complete", lines, fixed = TRUE))) {
      stop("LIONS run failed/incomplete: ", lf)
    }
    lines <- sub("\t.*$", "", lines)
    lions <- genes(sub(
      ":.*$", "", grep(":OFFICIAL_SIMTX", lines, value = TRUE)
    ))

    pf <- file.path(
      d, "TEProf2_benchmark", paste0(depth, "_test"), "assembly",
      paste0(tag, ".stringtie.gtf_annotated_filtered_test_all")
    )
    x <- read_tab(pf, FALSE)
    if (nrow(x)) {
      if (ncol(x) < 16) stop("Unexpected TEProf2 columns: ", pf)
      tp <- coalesce(clean(x[[16]]), clean(x[[3]]))
    } else {
      tp <- character()
    }

    fm <- map_fredy(file.path(
      d, "FREDY_relaxed_by_depth", depth, "chimeric/protein.gtf"
    ))

    pred <- list(
      TEDDY = genes(teddy$gene_name),
      Arriba = arriba,
      ChimeraTE = genes(chim),
      FREDY = genes(fm$gene_name),
      LIONS = lions,
      TEProf2 = genes(tp)
    )

    imap_dfr(pred, function(g, tool) {
      tibble(
        seed = seed, depth = depth, tool = tool,
        setting = unname(setting[tool]),
        truth_scope = unname(scope[tool]),
        scoring_unit = "assigned_host_gene_name",
        truth = list(truths[[scope[[tool]]]]),
        pred = list(g),
        universe = list(universe),
        unmapped_tx = if (tool == "FREDY") {
          sum(is.na(fm$gene_name))
        } else NA_integer_,
        tied_best_tx = if (tool == "FREDY") {
          sum(fm$tied_best %in% TRUE)
        } else NA_integer_,
        mapping = list(if (tool == "FREDY") fm else NULL)
      )
    })
  })
}

sets <- map_dfr(101:109, score_seed)

metrics <- sets |>
  mutate(
    truth_gene = lengths(truth),
    pred_gene = lengths(pred),
    TP = map2_int(pred, truth, function(p, t) length(intersect(p, t))),
    FP = pred_gene - TP,
    FN = truth_gene - TP,
    precision = if_else(pred_gene > 0, TP / pred_gene, NA_real_),
    recall = if_else(truth_gene > 0, TP / truth_gene, NA_real_),
    F1 = if_else(
      2 * TP + FP + FN > 0,
      2 * TP / (2 * TP + FP + FN),
      NA_real_
    )
  ) |>
  dplyr::select(
    seed, depth, tool, setting, truth_scope, scoring_unit,
    truth_gene, pred_gene, TP, FP, FN, precision, recall, F1,
    unmapped_tx, tied_best_tx
  ) |>
  arrange(seed, match(depth, depths), match(tool, names(scope)))

stopifnot(
  nrow(metrics) == 270L,
  !anyDuplicated(metrics[c("seed", "depth", "tool")]),
  all(metrics$TP + metrics$FN == metrics$truth_gene),
  all(metrics$TP + metrics$FP == metrics$pred_gene)
)

out <- file.path(root, "metrics")
dir.create(out, showWarnings = FALSE)

# Write summaries after all scoring completes; preserve caller outputs.
write_tsv(metrics, file.path(out, "metrics_101_109.tsv"))
saveRDS(
  list(metrics = metrics, gene_sets = sets),
  file.path(out, "metrics_101_109.rds")
)

print(
  metrics |>
    filter(seed == 101) |>
    dplyr::select(
      depth, tool, truth_gene, pred_gene, TP, FP, FN,
      precision, recall, F1, unmapped_tx
    ),
  n = 30, width = Inf
)
cat("\nCOMPLETE: 270 rows. Output: ", out, "\n", sep = "")
