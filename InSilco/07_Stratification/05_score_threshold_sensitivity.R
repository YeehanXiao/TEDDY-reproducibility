#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

event_rows <- function(annotation, expressed_tx = NULL) {
  x <- as.data.table(annotation)
  required <- c("seqnames", "start", "end", "gene_name", "transcript_id", "TE_start",
                "TE_end", "TE_strand", "overlap_bp")
  stopifnot(all(required %in% names(x)))
  if (!is.null(expressed_tx)) x <- x[transcript_id %chin% expressed_tx]
  x[, `:=`(gene_name = as.character(gene_name),
           TE_key = paste(seqnames, TE_start, TE_end, TE_strand, sep = ":"))]
  stopifnot(!anyNA(x[, .(start, end, gene_name, TE_key, overlap_bp)]),
            all(nzchar(x$gene_name)), all(x$overlap_bp >= 5))
  x
}

key_score <- function(truth, pred, matched = NULL) {
  truth_keys <- unique(truth[, .(gene_name, TE_key)])
  pred_keys <- unique(pred[, .(gene_name, TE_key)])
  tp <- if (is.null(matched)) nrow(merge(truth_keys, pred_keys, by = c("gene_name", "TE_key"))) else nrow(matched)
  c(N_truth = nrow(truth_keys), N_pred = nrow(pred_keys), TP = tp,
    precision = if (nrow(pred_keys)) tp / nrow(pred_keys) else NA_real_,
    recall = if (nrow(truth_keys)) tp / nrow(truth_keys) else NA_real_)
}

boundary_scores <- function(truth, pred, seed, depth) {
  keys <- c("gene_name", "TE_key")
  baseline <- key_score(truth, pred)
  pairs <- merge(
    unique(truth[, .(gene_name, TE_key, truth_start = start, truth_end = end)]),
    unique(pred[, .(gene_name, TE_key, pred_start = start, pred_end = end)]),
    by = keys, allow.cartesian = TRUE
  )
  settings <- c("original_locus", paste0("tolerance_", c(0, 5, 10, 25), "bp"))
  thresholds <- c(NA_real_, 0, 5, 10, 25)
  rbindlist(lapply(seq_along(settings), function(i) {
    score <- if (is.na(thresholds[i])) baseline else key_score(
      truth, pred,
      unique(pairs[abs(pred_start - truth_start) <= thresholds[i] &
                   abs(pred_end - truth_end) <= thresholds[i], ..keys])
    )
    out <- data.table(analysis = "boundary_tolerance", setting = settings[i],
                      threshold_bp = thresholds[i], seed = seed, depth = depth,
                      boundary_rule = if (is.na(thresholds[i])) "gene_and_exact_TE_locus" else "both_exon_boundaries_within_tolerance",
                      annotation_boundary_tolerance_bp = 5)
    out[, names(score) := as.list(score)]
    out
  }))
}

overlap_scores <- function(truth, pred, seed, depth) {
  thresholds <- c(5, 10, 20, 50)
  rbindlist(lapply(thresholds, function(threshold) {
    score <- key_score(truth, pred[overlap_bp >= threshold])
    out <- data.table(analysis = "minimum_overlap", setting = paste0("minimum_", threshold, "bp"),
                      threshold_bp = threshold, seed = seed, depth = depth,
                      boundary_rule = "not_applicable", annotation_boundary_tolerance_bp = 5)
    out[, names(score) := as.list(score)]
    out
  }))
}

args <- commandArgs(trailingOnly = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
if (length(args) > 2L) stop("Usage: Rscript 05_score_threshold_sensitivity.R INPUT_RDS OUTPUT_DIR")
input_rds <- if (length(args)) args[1L] else file.path(script_dir, "inputs/structural_annotations.rds")
output_dir <- if (length(args) > 1L) args[2L] else file.path(script_dir, "results")
performance_tsv <- file.path(output_dir, "stratified_performance.tsv")
stopifnot(file.exists(input_rds), file.exists(performance_tsv))
depths <- c("5x", "10x", "25x", "50x", "100x")
cache <- readRDS(input_rds)$seed_annotations
stopifnot(setequal(as.integer(names(cache)), c(18L, 101:109)))
per_seed <- rbindlist(lapply(cache, function(seed_data) {
  stopifnot(setequal(names(seed_data$predictions), depths))
  truth <- event_rows(seed_data$truth_annotation, seed_data$expressed_tx)
  rbindlist(lapply(depths, function(depth) {
    pred <- event_rows(seed_data$predictions[[depth]]$annotation)
    rbindlist(list(boundary_scores(truth, pred, seed_data$seed, depth),
                   overlap_scores(truth, pred, seed_data$seed, depth)))
  }))
}))

baseline_check <- per_seed[analysis == "boundary_tolerance" & setting == "original_locus"]
stored <- fread(performance_tsv)[
  row_type == "seed" & cohort == "per_seed" &
  evaluation_unit == "gene_plus_exact_annotated_TE_locus" &
  stratification == "overall" & performance_metric %chin% c("recall", "precision"),
  .(seed, depth, performance_metric, stored_estimate = estimate)
]
check <- merge(melt(baseline_check, id.vars = c("seed", "depth"),
                    measure.vars = c("precision", "recall"),
                    variable.name = "performance_metric", value.name = "estimate"),
               stored, by = c("seed", "depth", "performance_metric"))
stopifnot(nrow(check) == 100L, max(abs(check$estimate - check$stored_estimate)) < 1e-12)

summary_table <- per_seed[, .(
  n_seeds = uniqueN(seed), N_truth_mean = mean(N_truth), N_pred_mean = mean(N_pred),
  TP_mean = mean(TP), precision = mean(precision, na.rm = TRUE),
  precision_sd = sd(precision, na.rm = TRUE), recall = mean(recall), recall_sd = sd(recall)
), by = .(analysis, setting, threshold_bp, depth, boundary_rule,
          annotation_boundary_tolerance_bp)]
summary_table[, `:=`(
  baseline_setting = fifelse(analysis == "boundary_tolerance", "original_locus", "minimum_5bp"),
  precision_delta_pp = 100 * (precision - precision[setting == fifelse(
    analysis[1] == "boundary_tolerance", "original_locus", "minimum_5bp")]),
  recall_delta_pp = 100 * (recall - recall[setting == fifelse(
    analysis[1] == "boundary_tolerance", "original_locus", "minimum_5bp")])
), by = .(analysis, depth)]
summary_table[, depth := factor(depth, depths)]
setorder(summary_table, analysis, depth, threshold_bp, na.last = FALSE)
summary_table[, depth := as.character(depth)]
stopifnot(nrow(per_seed) == 450L, nrow(summary_table) == 45L, all(summary_table$n_seeds == 10L))
fwrite(summary_table, file.path(output_dir, "TEDDY_bp_threshold_sensitivity.tsv"),
       sep = "\t", na = "NA", quote = FALSE)
message("Created 45 threshold sensitivity rows in ", output_dir)
