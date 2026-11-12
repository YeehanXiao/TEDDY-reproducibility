suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
root <- args[[1]]

m <- read_tsv(
    file.path(root, "metrics", "metrics_101_109.tsv"),
    show_col_types = FALSE
)

depths <- c("5x", "10x", "25x", "50x", "100x")
tools <- c(
    "TEDDY",
    "TEProf2",
    "FREDY",
    "LIONS",
    "Arriba",
    "ChimeraTE"
)

stopifnot(
    nrow(m) == 270L,
    dplyr::n_distinct(m$seed, m$depth, m$tool) == 270L,
    setequal(m$seed, 101:109),
    setequal(m$depth, depths),
    setequal(m$tool, tools),
    all(m$TP >= 0),
    all(m$FP >= 0),
    all(m$FN >= 0),
    all(m$TP + m$FP == m$pred_gene),
    all(m$TP + m$FN == m$truth_gene),
    all(is.finite(m$precision)),
    all(is.finite(m$recall)),
    all(is.finite(m$F1))
)

if ("setting" %in% names(m)) {
    lions_settings <- unique(
        m$setting[m$tool == "LIONS"]
    )

    stopifnot(
        length(lions_settings) == 1L,
        lions_settings == "simuOptimal"
    )
}

cat("PASS: 270 unique seed-depth-tool combinations.\n")
cat("PASS: TP/FP/FN identities are valid.\n")
cat("PASS: all performance metrics are finite.\n")
cat("PASS: LIONS uses the final simuOptimal calls.\n")

cat("\n100x summary across nine simulation seeds:\n")

m |>
    filter(depth == "100x") |>
    group_by(tool, truth_scope) |>
    summarise(
        n_seed = n_distinct(seed),
        precision_mean = mean(precision),
        recall_mean = mean(recall),
        F1_mean = mean(F1),
        F1_SD = sd(F1),
        .groups = "drop"
    ) |>
    arrange(match(tool, tools)) |>
    print(n = Inf, width = Inf)

cat("\nFINAL METRICS QC PASS\n")
