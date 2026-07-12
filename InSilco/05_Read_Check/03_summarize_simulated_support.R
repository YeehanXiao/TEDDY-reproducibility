# ==============================================================================
# Summarize read-level support for simulated TE-host breakpoint structures.
#
# Main output:
#   output/ReadSupport_for_Supplementary_Table.tsv
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})


script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)

script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_arg)))
} else {
  normalizePath(getwd())
}

project_dir <- normalizePath(file.path(script_dir, ".."))
input_dir <- file.path(project_dir, "01_Data_Simulation", "input")
output_dir <- file.path(script_dir, "output")

support_file <- file.path(
  output_dir,
  "truth_TE_host_breakpoint_read_support_by_depth.tsv"
)

iso_file <- file.path(
  input_dir,
  "official_simulated_1000.isoforms.results"
)

stopifnot(file.exists(support_file))
stopifnot(file.exists(iso_file))

support <- read_tsv(
  support_file,
  show_col_types = FALSE
)

iso <- read.delim(
  iso_file,
  check.names = FALSE
) |>
  as_tibble() |>
  select(
    transcript_id,
    TPM,
    expected_count
  )

support_anno <- support |>
  left_join(
    iso,
    by = "transcript_id"
  ) |>
  mutate(
    TPM = coalesce(TPM, 0),
    expected_count = coalesce(
      expected_count,
      0
    ),
    Depth = str_extract(
      depth,
      "(5x|10x|25x|50x|100x)"
    ),
    Depth = factor(
      Depth,
      levels = c(
        "5x",
        "10x",
        "25x",
        "50x",
        "100x"
      )
    ),
    Evidence_type = case_when(
      evidence_type == "short_full_span" ~
        "Short TE full-span",

      evidence_type %in% c(
        "long_boundary_left",
        "long_boundary_right"
      ) ~
        "Long TE boundary",

      TRUE ~ NA_character_
    ),
    Combined_supported = as.integer(
      Anchor_supported == 1 |
        CIGAR_supported == 1
    )
  ) |>
  filter(!is.na(Evidence_type))

summary_by_depth <- support_anno |>
  filter(TPM > 1) |>
  group_by(
    Depth,
    Evidence_type
  ) |>
  summarise(
    Truth_loci = n(),
    Anchor_pair_supported_loci = sum(
      Anchor_supported
    ),
    CIGAR_supported_loci = sum(
      CIGAR_supported
    ),
    Supported_loci = sum(
      Combined_supported
    ),
    Support_rate = (
      Supported_loci / Truth_loci * 100
    ),
    .groups = "drop"
  ) |>
  arrange(
    Depth,
    Evidence_type
  )

summary_by_TPM_bin <- support_anno |>
  mutate(
    Expression_group = cut(
      TPM,
      breaks = c(
        -Inf,
        1,
        10,
        50,
        Inf
      ),
      labels = c(
        "Very low TPM (<=1)",
        "Low TPM (1-10)",
        "Medium TPM (10-50)",
        "High TPM (>50)"
      )
    )
  ) |>
  group_by(
    Depth,
    Evidence_type,
    Expression_group
  ) |>
  summarise(
    Truth_loci = n(),
    Anchor_pair_supported_loci = sum(
      Anchor_supported
    ),
    CIGAR_supported_loci = sum(
      CIGAR_supported
    ),
    Supported_loci = sum(
      Combined_supported
    ),
    Support_rate = (
      Supported_loci / Truth_loci * 100
    ),
    .groups = "drop"
  ) |>
  arrange(
    Depth,
    Evidence_type,
    Expression_group
  )

read_support_depth <- summary_by_depth |>
  transmute(
    Loci_subset = "Expressed loci across depths",
    Depth = as.character(Depth),
    Evidence_type,
    Expression_group = "TPM > 1",
    Truth_loci,
    Anchor_pair_supported_loci,
    CIGAR_supported_loci,
    Supported_loci,
    Support_rate
  )

read_support_100x_TPM <- summary_by_TPM_bin |>
  filter(Depth == "100x") |>
  transmute(
    Loci_subset = "100x loci by TPM group",
    Depth = as.character(Depth),
    Evidence_type,
    Expression_group = as.character(
      Expression_group
    ),
    Truth_loci,
    Anchor_pair_supported_loci,
    CIGAR_supported_loci,
    Supported_loci,
    Support_rate
  )

read_support_table <- bind_rows(
  read_support_depth,
  read_support_100x_TPM
) |>
  mutate(
    Loci_subset = factor(
      Loci_subset,
      levels = c(
        "Expressed loci across depths",
        "100x loci by TPM group"
      )
    ),
    Depth = factor(
      Depth,
      levels = c(
        "5x",
        "10x",
        "25x",
        "50x",
        "100x"
      )
    ),
    Evidence_type = factor(
      Evidence_type,
      levels = c(
        "Long TE boundary",
        "Short TE full-span"
      )
    ),
    Expression_group = factor(
      Expression_group,
      levels = c(
        "TPM > 1",
        "Very low TPM (<=1)",
        "Low TPM (1-10)",
        "Medium TPM (10-50)",
        "High TPM (>50)"
      )
    ),
    Support_rate = round(
      Support_rate,
      2
    )
  ) |>
  arrange(
    Loci_subset,
    Depth,
    Evidence_type,
    Expression_group
  ) |>
  mutate(
    Loci_subset = as.character(
      Loci_subset
    ),
    Depth = as.character(Depth),
    Evidence_type = as.character(
      Evidence_type
    ),
    Expression_group = as.character(
      Expression_group
    )
  )

write_tsv(
  read_support_table,
  file.path(
    output_dir,
    "ReadSupport_for_Supplementary_Table.tsv"
  )
)

saveRDS(
  list(
    support_anno = support_anno,
    summary_by_depth = summary_by_depth,
    summary_by_TPM_bin = summary_by_TPM_bin,
    read_support_table = read_support_table
  ),
  file.path(
    output_dir,
    "truth_breakpoint_support_sanity_results.rds"
  )
)

message("Saved read-support summaries to: ", output_dir)

read_support_table