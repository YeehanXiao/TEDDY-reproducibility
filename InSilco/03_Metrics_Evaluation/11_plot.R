## ============================================================
## Benchmark accuracy plot
## Precision / Recall / F1 across sequencing depth
## ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

## ---- 1. Paths ----

project_dir <- "/path/to/simulation_benchmark"

accuracy_summary_file <- file.path(
  project_dir,
  "03_Metrics_Evaluation",
  "accuracy_summary.rds"
)

out_dir <- file.path(
  project_dir,
  "03_Metrics_Evaluation",
  "figures"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---- 2. Load accuracy summary ----

accuracy_summary <- readRDS(accuracy_summary_file)

depths <- c("5x", "10x", "25x", "50x", "100x")

tool_levels <- c(
  "TEDDY",
  "TEProf2",
  "FREDY",
  "LIONS",
  "Arriba",
  "ChimeraTE"
)

## ---- 3. Prepare plotting table ----

plot_df <- accuracy_summary |>
  select(tool, depth, precision, recall, F1) |>
  pivot_longer(
    cols = c(precision, recall, F1),
    names_to = "Metric",
    values_to = "Value"
  ) |>
  mutate(
    Metric = factor(
      Metric,
      levels = c("precision", "recall", "F1"),
      labels = c("Precision", "Recall", "F1")
    ),
    depth = factor(as.character(depth), levels = depths),
    tool = factor(as.character(tool), levels = tool_levels),
    TE_category = if_else(
      tool %in% c("TEProf2", "LIONS"),
      "TE-initiated only",
      "All TE-chimeric categories"
    ),
    TE_category = factor(
      TE_category,
      levels = c("All TE-chimeric categories", "TE-initiated only")
    )
  )

## ---- 4. Colors and line types ----

tool_cols <- c(
  "TEDDY"      = "#8c3837",
  "TEProf2"   = "#383867",
  "FREDY"     = "#705296",
  "LIONS"     = "#F7D070",
  "Arriba"    = "#EBC9C7",
  "ChimeraTE" = "#666666"
)

category_linetypes <- c(
  "All TE-chimeric categories" = "solid",
  "TE-initiated only" = "dashed"
)

## ---- 5. Single-panel plotting function ----

make_metric_plot <- function(metric_name, show_y = FALSE) {
  pdat <- plot_df |>
    filter(Metric == metric_name)
  
  ggplot(
    pdat,
    aes(
      x = depth,
      y = Value,
      group = tool,
      color = tool,
      linetype = TE_category
    )
  ) +
    geom_line(linewidth = 0.9, alpha = 0.95) +
    geom_point(size = 2.7, shape = 16, alpha = 1) +
    
    geom_line(
      data = pdat |> filter(tool == "TEDDY"),
      linewidth = 1.15,
      alpha = 1
    ) +
    geom_point(
      data = pdat |> filter(tool == "TEDDY"),
      size = 3.0,
      shape = 16,
      alpha = 1
    ) +
    
    scale_color_manual(values = tool_cols, drop = FALSE) +
    scale_linetype_manual(values = category_linetypes, drop = FALSE) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.25),
      expand = expansion(mult = c(0.02, 0.04))
    ) +
    labs(
      title = metric_name,
      x = "Sequencing depth",
      y = if (show_y) "Score" else NULL,
      color = "Tool",
      linetype = "TE category"
    ) +
    guides(
      color = guide_legend(
        order = 1,
        nrow = 1,
        override.aes = list(
          linetype = "solid",
          linewidth = 1.0,
          size = 2.8,
          alpha = 1
        )
      ),
      linetype = guide_legend(
        order = 2,
        nrow = 1,
        override.aes = list(
          color = "black",
          linewidth = 0.9
        )
      )
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      axis.title.x = element_text(size = 11.5),
      axis.title.y = element_text(size = 11.5),
      axis.text = element_text(size = 10.5, color = "black"),
      panel.grid.major = element_line(linewidth = 0.32, color = "grey88"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.45),
      plot.margin = margin(5, 8, 5, 8),
      legend.position = "top",
      legend.title = element_text(size = 10.5, face = "bold"),
      legend.text = element_text(size = 10),
      legend.key.width = unit(0.75, "cm")
    )
}

## ---- 6. Generate figure ----

p_precision <- make_metric_plot("Precision", show_y = TRUE)
p_recall <- make_metric_plot("Recall", show_y = FALSE)
p_f1 <- make_metric_plot("F1", show_y = FALSE)

p_accuracy <- (
  p_precision + p_recall + p_f1 +
    plot_layout(nrow = 1, guides = "collect")
) &
  theme(
    legend.position = "top",
    legend.justification = "center",
    legend.box = "vertical"
  )

print(p_accuracy)

