#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

compare_table <- function(actual, expected, keys, label) {
  stopifnot(setequal(names(actual), names(expected)), nrow(actual) == nrow(expected),
            !anyDuplicated(actual[, ..keys]), !anyDuplicated(expected[, ..keys]))
  setcolorder(actual, names(expected))
  setorderv(actual, keys)
  setorderv(expected, keys)
  mismatches <- vapply(names(expected), function(column) {
    a <- actual[[column]]
    b <- expected[[column]]
    if (is.numeric(b)) {
      sum(xor(is.na(a), is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) > 1e-12))
    } else sum(xor(is.na(a), is.na(b)) | (!is.na(a) & !is.na(b) & as.character(a) != as.character(b)))
  }, integer(1))
  if (any(mismatches)) stop(label, " differs in: ", paste(names(mismatches)[mismatches > 0], collapse = ", "))
  message(label, ": ", nrow(actual), " rows; all keys, values and NA cells match")
  invisible(actual)
}

check_counts <- function(x) {
  seed_rows <- x[row_type == "seed"]
  stopifnot(all(seed_rows[!is.na(N_truth), N_truth == TP + FN]),
            all(seed_rows[!is.na(N_pred), N_pred == pred_TP + FP]),
            all(seed_rows[!is.na(N_truth), N_truth >= 0 & TP >= 0 & FN >= 0]),
            all(seed_rows[!is.na(N_pred), N_pred >= 0 & pred_TP >= 0 & FP >= 0]),
            all(seed_rows[performance_metric == "recall" & N_truth > 0, abs(estimate - TP / N_truth) < 1e-12]),
            all(seed_rows[performance_metric == "precision" & N_pred > 0, abs(estimate - pred_TP / N_pred) < 1e-12]),
            !any(x[stratification != "overall", performance_metric == "F1"]))
}

cli <- commandArgs(trailingOnly = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
if (length(cli) > 2L) stop("Usage: Rscript 04_verify.R RESULTS_DIR REFERENCE_DIR")
results <- if (length(cli)) cli[1L] else file.path(script_dir, "results")
reference <- if (length(cli) > 1L) cli[2L] else file.path(script_dir, "inputs/reference")
keys <- c("row_type", "cohort", "seed", "depth", "evaluation_unit", "stratification", "stratum", "performance_metric")
invisible(lapply(c("stratified_performance.tsv", "structural_stratified_performance.tsv"), function(filename) {
  actual <- fread(file.path(results, filename))
  compare_table(actual, fread(file.path(reference, filename)), keys, filename)
  check_counts(actual)
}))

recovery <- fread(file.path(results, "TE_structural_recovery_10seeds.tsv"))
archived <- as.data.table(readRDS(file.path(reference, "TE_structural_recovery_10seeds.rds"))$recovery)
archived[, depth := as.character(depth)]
compare_table(recovery, archived, c("seed", "depth"), "Fig2A recovery")
stopifnot(nrow(recovery) == 50L, setequal(recovery$seed, c(18L, 101:109)),
          setequal(recovery$depth, c("5x", "10x", "25x", "50x", "100x")))
sensitivity <- fread(file.path(results, "TEDDY_bp_threshold_sensitivity.tsv"))
compare_table(sensitivity, fread(file.path(reference, "TEDDY_bp_threshold_sensitivity.tsv")),
              c("analysis", "setting", "depth"), "Threshold sensitivity")
stopifnot(nrow(sensitivity) == 45L, all(sensitivity$n_seeds == 10L),
          all(sensitivity$N_truth_mean > 0), all(sensitivity$N_pred_mean >= 0),
          all(sensitivity$TP_mean >= 0),
          all(sensitivity$TP_mean <= sensitivity$N_truth_mean),
          all(sensitivity$TP_mean <= sensitivity$N_pred_mean),
          all(sensitivity$precision >= 0 & sensitivity$precision <= 1),
          all(sensitivity$recall >= 0 & sensitivity$recall <= 1),
          all(sensitivity[, uniqueN(N_truth_mean), by = depth]$V1 == 1L))
message("Verified both stratification tables, 200 Fig2A denominator/rate pairs, 45 threshold rows, and count conservation.")
