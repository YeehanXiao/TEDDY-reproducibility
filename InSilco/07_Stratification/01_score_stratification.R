#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

depths <- c("5x", "10x", "25x", "50x", "100x")
seeds <- c(18L, 101:109)
counts <- c("N_truth", "TP", "FN", "N_pred", "pred_TP", "FP")
axes <- list(
  O = list("overlap_bp", c(-Inf, 50, 75, 150, Inf), c("<50 bp", "50–<75 bp", "75–<150 bp", "≥150 bp")),
  F_E = list("overlap_fraction_exon", c(-Inf, .10, .25, .50, Inf), c("<10%", "10–<25%", "25–<50%", "≥50%")),
  F_inc = list("overlap_fraction_TE", c(-Inf, .50, .95, Inf), c("<50%", "50–<95%", "≥95%"))
)
labels <- c(chimeric_category = "TE_category", relative_TE_exon_position = "TE_position",
            context_aware_TE_boundary = "TE_boundary")
mean_na <- function(x) if (any(is.finite(x))) mean(x, na.rm = TRUE) else NA_real_
sd_na <- function(x) if (sum(is.finite(x)) > 1L) sd(x, na.rm = TRUE) else NA_real_

event_frame <- function(annotation, expressed = NULL) {
  x <- as.data.table(annotation)
  required <- c("gene_name", "transcript_id", "seqnames", "TE_start", "TE_end", "TE_strand",
                "overlap_bp", "overlap_fraction_exon", "overlap_fraction_TE", unname(labels))
  stopifnot(all(required %in% names(x)))
  if (!is.null(expressed)) x <- x[transcript_id %chin% expressed]
  stopifnot(!anyNA(x[, ..required]), all(nzchar(x$gene_name)), all(x$overlap_bp >= 5))
  x[, `:=`(gene_name = as.character(gene_name),
           TE_key = paste(seqnames, TE_start, TE_end, TE_strand, sep = ":"))]
  x
}

units <- function(x, key) {
  if (any(!nzchar(as.matrix(x[, ..key])))) stop("Empty scoring key")
  x[, lapply(.SD, max), by = key,
    .SDcols = vapply(axes, `[[`, "", 1L)]
}

score <- function(truth, pred, key) {
  tp <- nrow(merge(truth[, ..key], pred[, ..key], by = key))
  c(N_truth = nrow(truth), TP = tp, FN = nrow(truth) - tp,
    N_pred = nrow(pred), pred_TP = tp, FP = nrow(pred) - tp)
}

stratified_rows <- function(truth, pred, key, include_unassigned) {
  rbindlist(lapply(names(axes), function(axis) {
    spec <- axes[[axis]]
    truth_bin <- as.character(cut(truth[[spec[[1]]]], spec[[2]], spec[[3]], right = FALSE))
    pred_bin <- as.character(cut(pred[[spec[[1]]]], spec[[2]], spec[[3]], right = FALSE))
    truth_bin[is.na(truth_bin)] <- "unassigned"
    pred_bin[is.na(pred_bin)] <- "unassigned"
    rbindlist(lapply(c(spec[[3]], if (include_unassigned) "unassigned"), function(bin) {
      recall <- score(truth[truth_bin == bin], pred, key)
      precision <- score(truth, pred[pred_bin == bin], key)
      data.table(stratification = axis, stratum = bin,
                 grouping_basis = c("truth", "prediction"), performance_metric = c("recall", "precision"),
                 N_truth = c(recall[["N_truth"]], NA_real_), TP = c(recall[["TP"]], NA_real_),
                 FN = c(recall[["FN"]], NA_real_), N_pred = c(NA_real_, precision[["N_pred"]]),
                 pred_TP = c(NA_real_, precision[["TP"]]), FP = c(NA_real_, precision[["FP"]]),
                 estimate = c(if (recall[["N_truth"]]) recall[["TP"]] / recall[["N_truth"]] else NA_real_,
                              if (precision[["N_pred"]]) precision[["TP"]] / precision[["N_pred"]] else NA_real_))
    }))
  }))
}

add_metadata <- function(rows, seed, depth) {
  rows[, `:=`(row_type = "seed", cohort = "per_seed", seed = as.integer(seed), depth = depth,
              tool = "TEDDY", configuration = "original", truth_scope = "expressed_TE_chimeric",
              n_seeds = 1L, n_nonempty_seeds = as.integer(is.finite(estimate)), sd = NA_real_)]
  rows
}

baseline_rows <- function(truth, pred, key, seed, depth) {
  gene <- identical(key, "gene_name")
  rows <- stratified_rows(truth, pred, key, gene)
  total <- score(truth, pred, key)
  overall <- as.data.table(as.list(total))[rep(1L, 3L)]
  overall[, `:=`(stratification = "overall", stratum = "all", grouping_basis = "overall",
                 performance_metric = c("recall", "precision", "F1"),
                 estimate = c(TP[1] / N_truth[1], TP[1] / N_pred[1],
                              2 * TP[1] / (2 * TP[1] + FP[1] + FN[1])))]
  rows <- rbindlist(list(rows, overall), use.names = TRUE)
  rows[, `:=`(evaluation_unit = if (gene) "assigned_host_gene_name" else "gene_plus_exact_annotated_TE_locus",
              aggregation = if (gene) "independent_gene_max" else "independent_locus_max")]
  add_metadata(rows, seed, depth)
}

structural_rows <- function(truth, pred, label, level, seed, depth) {
  key <- c("gene_name", "TE_key", label)
  rows <- stratified_rows(truth, pred, key, TRUE)
  total <- score(truth, pred, key)
  overall <- data.table(stratification = "overall", stratum = "all", grouping_basis = c("truth", "prediction"),
                        performance_metric = c("recall", "precision"),
                        N_truth = c(total[["N_truth"]], NA_real_), TP = c(total[["TP"]], NA_real_),
                        FN = c(total[["FN"]], NA_real_), N_pred = c(NA_real_, total[["N_pred"]]),
                        pred_TP = c(NA_real_, total[["TP"]]), FP = c(NA_real_, total[["FP"]]),
                        estimate = c(total[["TP"]] / total[["N_truth"]], total[["TP"]] / total[["N_pred"]]))
  rows <- rbindlist(list(rows, overall), use.names = TRUE)
  rows[, `:=`(scoring_level = level, matching_label = label,
              evaluation_unit = paste0("gene_plus_exact_annotated_TE_locus_plus_", label),
              aggregation = "independent_full_structural_key_max")]
  add_metadata(rows, seed, depth)
}

summaries <- function(rows, structural = FALSE) {
  grouping <- c("depth", "tool", "configuration", "truth_scope",
                if (structural) c("scoring_level", "matching_label"), "evaluation_unit", "aggregation",
                "stratification", "stratum", "grouping_basis", "performance_metric")
  cohorts <- list(all10 = seeds, formal9 = 101:109, seed18 = 18L)
  rbindlist(lapply(names(cohorts), function(cohort_name) {
    rows[seed %in% cohorts[[cohort_name]], c(
      lapply(.SD, mean_na), setNames(lapply(.SD, sd_na), paste0(counts, "_sd")),
      list(estimate = mean_na(estimate), sd = sd_na(estimate), n_seeds = uniqueN(seed),
           n_nonempty_seeds = sum(is.finite(estimate)), row_type = "summary", cohort = cohort_name, seed = NA_integer_)
    ), by = grouping, .SDcols = counts]
  }))
}

run_analysis <- function(input_rds, output_dir) {
  annotations <- readRDS(input_rds)$seed_annotations
  stopifnot(setequal(as.integer(names(annotations)), seeds))
  scored <- lapply(annotations, function(s) {
    stopifnot(setequal(names(s$predictions), depths))
    truth <- event_frame(s$truth_annotation, s$expressed_tx)
    lapply(depths, function(depth) {
      pred <- event_frame(s$predictions[[depth]]$annotation)
      baseline <- rbindlist(lapply(list("gene_name", c("gene_name", "TE_key")), function(key)
        baseline_rows(units(truth, key), units(pred, key), key, s$seed, depth)), use.names = TRUE)
      structural <- rbindlist(lapply(names(labels), function(level) {
        key <- c("gene_name", "TE_key", labels[[level]])
        structural_rows(units(truth, key), units(pred, key), labels[[level]], level, s$seed, depth)
      }), use.names = TRUE)
      recovery <- lapply(c(locus = "", category = "TE_category", position = "TE_position", boundary = "TE_boundary"), function(label) {
        key <- c("gene_name", "TE_key", if (nzchar(label)) label)
        score(units(truth, key), units(pred, key), key)
      })
      fig2a <- data.table(seed = s$seed, depth = depth)
      fig2a[, paste0("n_", names(recovery)) := lapply(recovery, `[[`, "N_truth")]
      fig2a[, paste0("TE_", names(recovery)) := lapply(recovery, function(x) x[["TP"]] / x[["N_truth"]])]
      list(baseline = baseline, structural = structural, fig2a = fig2a)
    })
  }) |> unlist(recursive = FALSE)
  baseline <- rbindlist(lapply(scored, `[[`, "baseline"))
  structural <- rbindlist(lapply(scored, `[[`, "structural"))
  baseline <- rbindlist(list(baseline, summaries(baseline)), use.names = TRUE, fill = TRUE)
  structural <- rbindlist(list(structural, summaries(structural, TRUE)), use.names = TRUE, fill = TRUE)
  first <- c("row_type", "cohort", "seed", "depth", "tool", "configuration", "truth_scope")
  rest <- c("evaluation_unit", "aggregation", "stratification", "stratum", "grouping_basis", "performance_metric",
            counts, "estimate", "sd", "n_seeds", "n_nonempty_seeds", paste0(counts, "_sd"))
  setcolorder(baseline, c(first, rest))
  setcolorder(structural, c(first, "scoring_level", "matching_label", rest))
  baseline[, `:=`(depth_order = match(depth, depths), axis_order = match(stratification, c("overall", names(axes))))]
  setorderv(baseline, c("row_type", "cohort", "seed", "depth_order", "evaluation_unit", "axis_order", "grouping_basis", "stratum"))
  baseline[, c("depth_order", "axis_order") := NULL]
  setorderv(structural, c("row_type", "cohort", "seed", "depth", "scoring_level", "stratification", "stratum", "grouping_basis"))
  stopifnot(nrow(baseline) == 3640L, nrow(structural) == 5850L)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(baseline, file.path(output_dir, "stratified_performance.tsv"), sep = "\t", na = "NA", quote = FALSE)
  fwrite(structural, file.path(output_dir, "structural_stratified_performance.tsv"), sep = "\t", na = "NA", quote = FALSE)
  fwrite(rbindlist(lapply(scored, `[[`, "fig2a")), file.path(output_dir, "TE_structural_recovery_10seeds.tsv"), sep = "\t", quote = FALSE)
  message("Created 3,640 baseline rows, 5,850 structural rows, and 50 Fig2A seed/depth rows in ", output_dir)
}

cli <- commandArgs(trailingOnly = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
if (length(cli) > 2L) stop("Usage: Rscript 01_score_stratification.R INPUT_RDS OUTPUT_DIR")
run_analysis(if (length(cli)) cli[1L] else file.path(script_dir, "inputs/structural_annotations.rds"),
             if (length(cli) > 1L) cli[2L] else file.path(script_dir, "results"))
