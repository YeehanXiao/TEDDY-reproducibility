suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
root <- args[[1]]
outdir <- args[[2]]

dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = FALSE
)

m <- read_tsv(
    file.path(root, "metrics", "metrics_101_109.tsv"),
    show_col_types = FALSE
)

tools <- c(
    "TEDDY",
    "TEProf2",
    "FREDY",
    "LIONS",
    "Arriba",
    "ChimeraTE"
)

depths <- c(
    "5x",
    "10x",
    "25x",
    "50x",
    "100x"
)

stopifnot(
    nrow(m) == 270L,
    dplyr::n_distinct(m$seed, m$depth, m$tool) == 270L,
    setequal(m$seed, 101:109)
)

summary_tbl <- m |>
    group_by(
        tool,
        depth,
        setting,
        truth_scope
    ) |>
    summarise(
        n_seed = n_distinct(seed),
        pred_gene = mean(pred_gene),
        truth_gene = mean(truth_gene),
        TP = mean(TP),
        FP = mean(FP),
        FN = mean(FN),
        precision_mean = mean(precision),
        precision_SD = sd(precision),
        recall_mean = mean(recall),
        recall_SD = sd(recall),
        F1_mean = mean(F1),
        F1_SD = sd(F1),
        .groups = "drop"
    )

stopifnot(
    nrow(summary_tbl) == 30L,
    all(summary_tbl$n_seed == 9L)
)

summary_out <- summary_tbl |>
    arrange(
        match(tool, tools),
        match(depth, depths)
    ) |>
    transmute(
        Tool = tool,
        Depth = depth,
        `Predicted Genes` = sprintf("%.1f", pred_gene),
        `Truth Genes` = sprintf("%.1f", truth_gene),
        TP = sprintf("%.1f", TP),
        FP = sprintf("%.1f", FP),
        FN = sprintf("%.1f", FN),
        Precision = sprintf(
            "%.3f ± %.3f",
            precision_mean,
            precision_SD
        ),
        Recall = sprintf(
            "%.3f ± %.3f",
            recall_mean,
            recall_SD
        ),
        F1 = sprintf(
            "%.3f ± %.3f",
            F1_mean,
            F1_SD
        )
    )

raw_out <- m |>
    arrange(
        seed,
        match(tool, tools),
        match(depth, depths)
    ) |>
    transmute(
        Seed = seed,
        Tool = tool,
        Depth = depth,
        `Predicted Genes` = pred_gene,
        `Truth Genes` = truth_gene,
        TP = TP,
        FP = FP,
        FN = FN,
        Precision = precision,
        Recall = recall,
        F1 = F1
    )

write_tsv(
    summary_out,
    file.path(
        outdir,
        "Sheet3_summary_9seeds.tsv"
    )
)

write_tsv(
    raw_out,
    file.path(
        outdir,
        "Sheet3_raw_9seeds.tsv"
    )
)

cat("Sheet 3 summary: 30 rows.\n")
cat("Per-seed metrics: 270 rows.\n")
cat("PASS: Sheet 3 files exported successfully.\n")
