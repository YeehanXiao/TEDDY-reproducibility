plotTERecovery <-
function(recovery) {

    x <- recovery |>
        tidyr::pivot_longer(
            cols = c(TE_locus, TE_category, TE_position, TE_boundary),
            names_to = "metric", values_to = "recovery"
        ) |>
        dplyr::mutate(
            depth = factor(depth, c("5x", "10x", "25x", "50x", "100x")),
            metric = factor(
                metric,
                c("TE_locus", "TE_category", "TE_position", "TE_boundary"),
                c("TE locus", "TE category", "TE position", "TE boundary")
            )
        )

    s <- x |>
        dplyr::group_by(metric, depth) |>
        dplyr::summarise(
            mean = mean(recovery),
            sd = stats::sd(recovery),
            .groups = "drop"
        )

    ggplot2::ggplot(s, ggplot2::aes(depth, mean, group = 1)) +
        ggplot2::geom_point(
            data = x, ggplot2::aes(y = recovery),
            position = ggplot2::position_jitter(
                width = 0.05, height = 0, seed = 1
            ),
            colour = "grey75", size = 1.2, alpha = 0.65
        ) +
        ggplot2::geom_errorbar(
            ggplot2::aes(ymin = mean - sd, ymax = mean + sd),
            width = 0.12, linewidth = 0.6, colour = "#8B3530"
        ) +
        ggplot2::geom_line(linewidth = 0.85, colour = "#8B3530") +
        ggplot2::geom_point(size = 2.6, colour = "#8B3530") +
        ggplot2::facet_wrap(~metric, nrow = 1) +
        ggplot2::scale_x_discrete(
            labels = c("5×", "10×", "25×", "50×", "100×"),
            expand = ggplot2::expansion(add = 0.35)
        ) +
        ggplot2::scale_y_continuous(
            limits = c(0, 1), breaks = seq(0, 1, 0.25),
            labels = scales::label_percent(accuracy = 1),
            expand = ggplot2::expansion(mult = c(0, 0.03))
        ) +
        ggplot2::labs(x = "Sequencing depth", y = "Recovery rate") +
        ggplot2::theme_classic(base_size = 12) +
        ggplot2::theme(
            aspect.ratio = NULL,
            strip.background = ggplot2::element_blank(),
            strip.text = ggplot2::element_text(
                face = "bold", hjust = 0, size = 11
            ),
            axis.title = ggplot2::element_text(size = 11),
            axis.text = ggplot2::element_text(size = 9),
            panel.grid.major.y = ggplot2::element_line(
                colour = "grey93", linewidth = 0.3
            ),
            panel.spacing.x = grid::unit(0.8, "lines"),
            plot.margin = ggplot2::margin(3, 4, 3, 3),
            legend.position = "none"
        )
}

# Run from this directory, or use source(..., chdir = TRUE).
# Red: equal-weight mean +/- SD across 10 seeds; grey: individual seeds.
p <- plotTERecovery(
    readRDS("TE_structural_recovery_10seeds.rds")$recovery
)
print(p)
ggplot2::ggsave(
    "TE_structural_recovery_10seeds.svg",
    plot = p, device = svglite::svglite,
    width = 19, height = 5.5, units = "cm", bg = "white"
)
