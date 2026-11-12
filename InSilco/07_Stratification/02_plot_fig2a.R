#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
if (length(args) > 2L) stop("Usage: Rscript 02_plot_fig2a.R INPUT_TSV OUTPUT_DIR")
input <- if (length(args)) args[1L] else file.path(script_dir, "results/TE_structural_recovery_10seeds.tsv")
output <- if (length(args) > 1L) args[2L] else file.path(script_dir, "results")
depths <- c("5x", "10x", "25x", "50x", "100x")

# Recreate the manuscript's three-panel layout from the archived four-metric
# source. Its relative-position panel uses TE_boundary, including exon context.
metrics <- c(TE_locus = "TE-chimeric locus", TE_category = "Chimeric category",
             TE_boundary = "Relative TE position")
recovery <- fread(input)
stopifnot(nrow(recovery) == 50L, !anyDuplicated(recovery[, .(seed, depth)]),
          setequal(recovery$seed, c(18L, 101:109)), setequal(recovery$depth, depths),
          all(vapply(recovery[, names(metrics), with = FALSE], function(x) all(is.finite(x) & x >= 0 & x <= 1), logical(1))))
x <- melt(recovery, id.vars = c("seed", "depth"), measure.vars = names(metrics),
          variable.name = "metric", value.name = "recovery")
x[, `:=`(depth = factor(depth, depths), metric = factor(metric, names(metrics), metrics))]
s <- x[, .(mean = mean(recovery), sd = sd(recovery)), by = .(metric, depth)]

figure <- ggplot(s, aes(depth, mean, group = 1)) +
  geom_point(data = x, aes(y = recovery),
             position = position_jitter(width = 0.05, height = 0, seed = 1),
             colour = "grey75", size = 1.2, alpha = 0.65) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                width = 0.12, linewidth = 0.6, colour = "#8B3530") +
  geom_line(linewidth = 0.85, colour = "#8B3530") +
  geom_point(size = 2.6, colour = "#8B3530") +
  facet_wrap(~metric, nrow = 1) +
  scale_x_discrete(labels = c("5×", "10×", "25×", "50×", "100×"), expand = expansion(add = 0.35)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                     labels = function(x) paste0(round(100 * x), "%"),
                     expand = expansion(mult = c(0, 0.03))) +
  labs(x = "Sequencing depth", y = "Recovery rate") +
  theme_classic(base_size = 12, base_family = "sans") +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", hjust = 0, size = 11),
        axis.title = element_text(size = 11), axis.text = element_text(size = 9),
        panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.3),
        panel.spacing.x = grid::unit(0.8, "lines"),
        plot.margin = margin(3, 4, 3, 3), legend.position = "none")

dir.create(output, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(output, "Fig2A_TE_structural_recovery.pdf"), figure,
       device = cairo_pdf, width = 19, height = 5.5, units = "cm", bg = "white")
ggsave(file.path(output, "Fig2A_TE_structural_recovery.svg"), figure,
       device = grDevices::svg, width = 19, height = 5.5, units = "cm", bg = "white")
message("Created Fig2A PDF and SVG in ", normalizePath(output))
