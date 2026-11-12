#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
if (length(args) > 2L) stop("Usage: Rscript 03_plot_stratification.R INPUT_TSV OUTPUT_DIR")
input <- if (length(args)) args[1L] else file.path(script_dir, "results/stratified_performance.tsv")
output <- if (length(args) > 1L) args[2L] else file.path(script_dir, "results")
depths <- c("5x", "10x", "25x", "50x", "100x")
palette <- c("#6C6C9D", "#E5C2A4", "#8C3837", "#383867")

plot_data <- fread(input)[
  row_type == "summary" & cohort == "all10" &
  evaluation_unit == "gene_plus_exact_annotated_TE_locus" &
  stratification %chin% c("O", "F_E", "F_inc") & stratum != "unassigned" &
  performance_metric == "recall"
]
stopifnot(nrow(plot_data) == 55L, all(plot_data$n_seeds == 10L),
          all(is.finite(plot_data$estimate)), all(is.finite(plot_data$sd)))
plot_data[, depth := factor(depth, depths)]

panel_plot <- function(metric, title, levels, colours) {
  x <- copy(plot_data[stratification == metric])
  stopifnot(setequal(x$stratum, levels), nrow(x) == 5L * length(levels))
  x[, stratum := factor(stratum, levels)]
  ggplot(x, aes(depth, estimate, colour = stratum, group = stratum)) +
    geom_hline(yintercept = seq(0, 1, 0.25), colour = "#E5E5E5", linewidth = 0.35) +
    geom_errorbar(aes(ymin = estimate - sd, ymax = estimate + sd),
                  width = 0.10, linewidth = 0.45, position = position_dodge(0.24)) +
    geom_line(linewidth = 0.75, position = position_dodge(0.24)) +
    geom_point(size = 2.2, position = position_dodge(0.24)) +
    scale_colour_manual(values = colours, breaks = levels, drop = FALSE, name = NULL) +
    guides(colour = guide_legend(nrow = if (length(levels) > 3) 2 else 1, byrow = TRUE)) +
    scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.25),
                       labels = function(x) paste0(round(x * 100), "%"), expand = c(0, 0)) +
    scale_x_discrete(expand = expansion(add = 0.25)) +
    labs(title = title, x = "Sequencing depth", y = "Recall") +
    theme_classic(base_size = 10, base_family = "sans") +
    theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5, margin = margin(b = 5)),
          axis.title = element_text(size = 9.5), axis.text = element_text(size = 9, colour = "black"),
          axis.line = element_line(linewidth = 0.4, colour = "black"),
          axis.ticks = element_line(linewidth = 0.4, colour = "black"),
          legend.position = "bottom", legend.text = element_text(size = 8.5),
          legend.key.width = grid::unit(0.65, "cm"), legend.spacing.x = grid::unit(0.05, "cm"),
          plot.margin = margin(4, 7, 2, 7))
}

figure <- panel_plot("O", "Overlap (bp)",
                     c("<50 bp", "50–<75 bp", "75–<150 bp", "≥150 bp"), palette) +
  panel_plot("F_E", "Exon fraction", c("<10%", "10–<25%", "25–<50%", "≥50%"), palette) +
  panel_plot("F_inc", "TE incorporation fraction", c("<50%", "50–<95%", "≥95%"), palette[1:3]) +
  plot_layout(nrow = 1)

dir.create(output, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(output, "TEDDY_contribution_stratified_recall.pdf"), figure,
       device = cairo_pdf, width = 11.4, height = 4.1, bg = "white")
ggsave(file.path(output, "TEDDY_contribution_stratified_recall.svg"), figure,
       device = grDevices::svg, width = 11.4, height = 4.1, bg = "white")
message("Created contribution-stratified recall PDF and SVG in ", normalizePath(output))
