# ==============================================================================
# Script: 02_timesummary.R
# Purpose: Parse and summarize wall-clock time and memory (RSS) usage 
#          from benchmarking log files.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. User parameters (Desensitized Paths)
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

log_dir <- "./results/logs"
output_dir <- "./results/runtime_benchmark"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 2. Parse Log Files
# ------------------------------------------------------------------------------
files <- list.files(
  log_dir,
  pattern = "\\.time\\.log$",
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No .time.log files found in the specified log directory.")
}

parse_elapsed_to_sec <- function(x) {
  x <- trimws(x)
  parts <- strsplit(x, ":", fixed = TRUE)[[1]]
  nums <- as.numeric(parts)
  
  if (length(nums) == 3) {
    nums[1] * 3600 + nums[2] * 60 + nums[3]
  } else if (length(nums) == 2) {
    nums[1] * 60 + nums[2]
  } else {
    as.numeric(x)
  }
}

parse_one_log <- function(f) {
  b <- basename(f)
  b <- sub("\\.time\\.log$", "", b)
  
  tool <- sub("_.*$", "", b)
  depth_rep <- sub("^[^_]+_", "", b)
  depth <- sub("_rep[0-9]+$", "", depth_rep)
  replicate <- sub("^.*_rep", "rep", depth_rep)
  
  lines <- readLines(f, warn = FALSE)
  
  elapsed_line <- grep("Elapsed \\(wall clock\\) time", lines, value = TRUE)
  rss_line <- grep("Maximum resident set size", lines, value = TRUE)
  
  elapsed <- if (length(elapsed_line) > 0) {
    sub(".*: ", "", elapsed_line[length(elapsed_line)])
  } else {
    NA_character_
  }
  
  max_rss_kb <- if (length(rss_line) > 0) {
    as.numeric(sub(".*: ", "", rss_line[length(rss_line)]))
  } else {
    NA_real_
  }
  
  tibble(
    tool = tool,
    depth = depth,
    replicate = replicate,
    wall_time = elapsed,
    wall_seconds = parse_elapsed_to_sec(elapsed),
    wall_minutes = wall_seconds / 60,
    max_rss_kb = max_rss_kb,
    max_rss_gb = max_rss_kb / 1024 / 1024,
    log_file = f
  )
}

# ------------------------------------------------------------------------------
# 3. Summarize and Export
# ------------------------------------------------------------------------------
runtime_raw <- bind_rows(lapply(files, parse_one_log)) |>
  filter(depth == "100x", grepl("^rep[0-9]+$", replicate)) |>
  arrange(tool, replicate)

runtime_summary <- runtime_raw |>
  group_by(tool, depth) |>
  summarise(
    n_runs = n(),
    mean_wall_seconds = mean(wall_seconds, na.rm = TRUE),
    sd_wall_seconds = sd(wall_seconds, na.rm = TRUE),
    mean_wall_minutes = mean_wall_seconds / 60,
    sd_wall_minutes = sd_wall_seconds / 60,
    mean_max_rss_gb = mean(max_rss_gb, na.rm = TRUE),
    sd_max_rss_gb = sd(max_rss_gb, na.rm = TRUE),
    max_rss_gb = max(max_rss_gb, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(mean_wall_seconds)

write_tsv(runtime_raw, file.path(output_dir, "runtime_raw.tsv"))
write_tsv(runtime_summary, file.path(output_dir, "runtime_summary.tsv"))

message("Runtime metrics parsed and saved to: ", output_dir)