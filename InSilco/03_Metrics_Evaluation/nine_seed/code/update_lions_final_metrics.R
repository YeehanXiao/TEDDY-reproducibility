suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || length(args) > 2L) {
    stop("Usage: Rscript update_lions_final_metrics.R ROOT [LIONS_CALLS_ROOT]")
}
root <- args[[1]]

metric_dir <- file.path(root, "metrics")
rds_file <- file.path(metric_dir, "metrics_101_109.rds")
tsv_file <- file.path(metric_dir, "metrics_101_109.tsv")

obj <- readRDS(rds_file)
metrics <- obj$metrics
gene_sets <- obj$gene_sets

stopifnot(
    nrow(metrics) == 270L,
    nrow(gene_sets) == 270L,
    n_distinct(
        metrics$seed,
        metrics$depth,
        metrics$tool
    ) == 270L
)

archive_dir <- file.path(metric_dir, "archive")
dir.create(
    archive_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

file.copy(
    rds_file,
    file.path(
        archive_dir,
        paste0(
            "metrics_101_109_before_final_LIONS_",
            stamp,
            ".rds"
        )
    ),
    overwrite = FALSE
)

file.copy(
    tsv_file,
    file.path(
        archive_dir,
        paste0(
            "metrics_101_109_before_final_LIONS_",
            stamp,
            ".tsv"
        )
    ),
    overwrite = FALSE
)

idx <- which(gene_sets$tool == "LIONS")
stopifnot(length(idx) == 45L)

for (i in idx) {

    seed_i <- gene_sets$seed[[i]]
    depth_i <- gene_sets$depth[[i]]

    lion_file <- if (length(args) == 2L) {
        file.path(
            args[[2]], paste0("seed_", seed_i), depth_i,
            paste0("official_simulated_", depth_i, "_noise0.1.lion")
        )
    } else file.path(
        root,
        paste0("seed_", seed_i),
        "results",
        "LIONS_simuOptimal_by_depth",
        depth_i,
        paste0(
            "official_simulated_",
            depth_i,
            "_noise0.1.lion"
        )
    )

    if (!file.exists(lion_file)) {
        stop(
            "Missing final LIONS output: ",
            lion_file
        )
    }

    x <- read_tsv(
        lion_file,
        show_col_types = FALSE,
        col_types = cols(.default = col_character())
    )

    if (!"transcriptID" %in% names(x)) {
        stop(
            "transcriptID column missing: ",
            lion_file
        )
    }

    pred <- if (nrow(x) == 0L) {
        character()
    } else {
        sort(
            unique(
                sub(
                    ":.*$",
                    "",
                    x$transcriptID
                )
            )
        )
    }

    truth <- sort(
        unique(
            gene_sets$truth[[i]]
        )
    )

    TP <- length(
        intersect(pred, truth)
    )

    FP <- length(
        setdiff(pred, truth)
    )

    FN <- length(
        setdiff(truth, pred)
    )

    n_pred <- length(pred)
    n_truth <- length(truth)

    precision <- if (n_pred == 0L) {
        NA_real_
    } else {
        TP / n_pred
    }

    recall <- TP / n_truth

    F1 <- if ((n_pred + n_truth) == 0L) {
        NA_real_
    } else {
        2 * TP / (n_pred + n_truth)
    }

    gene_sets$pred[[i]] <- pred

    if ("setting" %in% names(gene_sets)) {
        gene_sets$setting[[i]] <- "simuOptimal"
    }

    j <- which(
        metrics$tool == "LIONS" &
        metrics$seed == seed_i &
        metrics$depth == depth_i
    )

    stopifnot(length(j) == 1L)

    metrics$pred_gene[[j]] <- n_pred
    metrics$truth_gene[[j]] <- n_truth
    metrics$TP[[j]] <- TP
    metrics$FP[[j]] <- FP
    metrics$FN[[j]] <- FN
    metrics$precision[[j]] <- precision
    metrics$recall[[j]] <- recall
    metrics$F1[[j]] <- F1

    if ("setting" %in% names(metrics)) {
        metrics$setting[[j]] <- "simuOptimal"
    }
}

stopifnot(
    nrow(metrics) == 270L,
    n_distinct(
        metrics$seed,
        metrics$depth,
        metrics$tool
    ) == 270L,
    all(
        metrics$TP + metrics$FP ==
        metrics$pred_gene
    ),
    all(
        metrics$TP + metrics$FN ==
        metrics$truth_gene
    )
)

obj$metrics <- metrics
obj$gene_sets <- gene_sets

saveRDS(obj, rds_file)
write_tsv(metrics, tsv_file)

summary <- metrics |>
    filter(tool == "LIONS") |>
    group_by(depth) |>
    summarise(
        Predicted_Genes = mean(pred_gene),
        Truth_Genes = mean(truth_gene),
        TP = mean(TP),
        FP = mean(FP),
        FN = mean(FN),
        Precision_mean = mean(precision),
        Precision_SD = sd(precision),
        Recall_mean = mean(recall),
        Recall_SD = sd(recall),
        F1_mean = mean(F1),
        F1_SD = sd(F1),
        .groups = "drop"
    ) |>
    arrange(
        match(
            depth,
            c(
                "5x",
                "10x",
                "25x",
                "50x",
                "100x"
            )
        )
    )

print(
    summary,
    n = Inf,
    width = Inf
)

cat(
    "\nPASS: canonical metrics updated ",
    "from final LIONS .lion outputs.\n",
    sep = ""
)
