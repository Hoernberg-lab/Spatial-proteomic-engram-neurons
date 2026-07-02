# ===========================================================
# GCT Comparison Splitter: Per-Comparison Forward and Reverse Output
# Robust version for Metric.Comparison formatted GCT files
# Handles duplicate column names by keeping first occurrence only
# Author: Tobias Pohl
# ===========================================================

# -------------------------------
# Library Setup
# -------------------------------

required_pkgs <- c("dplyr", "readr", "stringr", "purrr", "fs", "tibble")

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
library(pacman)
pacman::p_load(char = required_pkgs)
source(file.path("R", "analysis_labels.R"))

# -------------------------------
# Parameters
# -------------------------------

labels <- source_analysis_labels()

# -------------------------------
# Input
# -------------------------------

setwd("S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/")

input_file <- "pg.matrix_Two-sample_mod_T_2025-12-15-transformed-p-val_n120x5349"

gct_path <- file.path(
  "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Datasets/gct/data",
  paste0(input_file, ".gct")
)

# -------------------------------
# Helper Functions
# -------------------------------

safe_name <- function(x) {
  x |>
    stringr::str_replace_all("[^A-Za-z0-9._-]", "_") |>
    stringr::str_replace_all("_+", "_")
}

format_side <- function(unit, sample_class, condition_code) {
  condition <- labels$condition_code_map[[condition_code]]
  if (is.null(condition) || is.na(condition)) condition <- condition_code
  paste(unit, sample_class, condition, sep = "_")
}

format_comparison <- function(case_unit, case_class, case_code, ref_unit, ref_class, ref_code) {
  paste0(
    format_side(case_unit, case_class, case_code),
    "_vs_",
    format_side(ref_unit, ref_class, ref_code)
  )
}

swap_comparison <- function(comp_key) {
  parts <- stringr::str_split(comp_key, "\\.over\\.", simplify = TRUE)

  if (ncol(parts) != 2) return(NA_character_)

  paste0(parts[2], ".over.", parts[1])
}

split_col <- function(col) {
  m <- stringr::str_match(
    col,
    "^([A-Za-z0-9\\.]+)\\.([A-Za-z0-9_]+\\.over\\.[A-Za-z0-9_]+)$"
  )

  if (is.na(m[1, 1])) {
    return(list(metric = NA_character_, comparison = NA_character_))
  }

  list(
    metric = m[1, 2],
    comparison = m[1, 3]
  )
}

parse_compkey <- function(key) {
  m <- stringr::str_match(
    key,
    "^([A-Za-z0-9]+)_([a-z]+)_([1234])\\.over\\.([A-Za-z0-9]+)_([a-z]+)_([1234])$"
  )

  if (!is.na(m[1, 1])) {
    r1 <- m[1, 2]
    g1 <- m[1, 3]
    l1 <- m[1, 4]

    r2 <- m[1, 5]
    g2 <- m[1, 6]
    l2 <- m[1, 7]

    return(format_comparison(r1, g1, l1, r2, g2, l2))
  }

  key2 <- stringr::str_replace_all(key, "\\.over\\.", "_")
  for (code in names(labels$condition_code_map)) {
    key2 <- stringr::str_replace_all(key2, paste0("_", code, "(?=($|_))"), paste0("_", labels$condition_code_map[[code]]))
  }
  key2 <- stringr::str_replace_all(key2, "[^A-Za-z0-9_]", "")

  key2
}

# -------------------------------
# Read GCT File
# -------------------------------

raw <- utils::read.delim(
  gct_path,
  header = FALSE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  comment.char = ""
)

raw_clean <- raw[-c(1:2), ]

header_row <- which(tolower(trimws(raw_clean[[1]])) == "id")[1]

if (is.na(header_row)) {
  stop("Could not find GCT header row. Expected first column to contain 'id'.")
}

col_names <- as.character(unlist(raw_clean[header_row, ]))

data <- raw_clean[(header_row + 1):nrow(raw_clean), , drop = FALSE]
colnames(data) <- col_names

# -------------------------------
# Remove Duplicate Columns
# Keep first occurrence only
# -------------------------------

dup_cols <- duplicated(names(data))

if (any(dup_cols)) {
  message(
    "Removing duplicate columns, keeping first occurrence: ",
    paste(unique(names(data)[dup_cols]), collapse = ", ")
  )

  data <- data[, !dup_cols, drop = FALSE]
}

if (!"id" %in% names(data)) {
  stop("No 'id' column found after assigning GCT header.")
}

if (anyDuplicated(names(data))) {
  stop(
    "Duplicate columns remain after cleanup: ",
    paste(unique(names(data)[duplicated(names(data))]), collapse = ", ")
  )
}

# -------------------------------
# Remove Annotation Rows
# -------------------------------

annotation_rows <- c(
  "AnimalID", "ReplicateGroup", "sample_class",
  "condition_code", "condition", "sample_class_condition", "plate",
  "sampleNumber", "shortname"
)

data <- data |>
  dplyr::filter(!id %in% annotation_rows, id != "na")

# -------------------------------
# Parse Feature Columns
# -------------------------------

feature_cols <- setdiff(names(data), "id")
split_info <- lapply(feature_cols, split_col)

comparison_keys <- setNames(
  vapply(split_info, `[[`, character(1), "comparison"),
  feature_cols
)

metric_keys <- setNames(
  vapply(split_info, `[[`, character(1), "metric"),
  feature_cols
)

valid_cols <- names(comparison_keys)[!is.na(comparison_keys)]

if (length(valid_cols) == 0) {
  stop("No valid Metric.Comparison columns detected.")
}

by_comparison <- split(valid_cols, comparison_keys[valid_cols])

# -------------------------------
# Convert Numeric Columns
# -------------------------------

data[valid_cols] <- lapply(data[valid_cols], readr::parse_number)

# -------------------------------
# Output Folder Structure
# -------------------------------

outdir_base <- "raw"

fname <- basename(gct_path)
subfolder <- basename(fs::path_dir(gct_path))

if (is.na(subfolder) || subfolder == "") {
  subfolder <- "unknown-comparison"
}

outdir <- file.path(outdir_base, subfolder)
outdir_fwd <- file.path(outdir, "forward")
outdir_rev <- file.path(outdir, "reverse")

fs::dir_create(outdir_fwd)
fs::dir_create(outdir_rev)

# -------------------------------
# Metric Rename Map
# -------------------------------

recode_map <- c(
  "adj.P.Val"   = "padj",
  "P.Value"     = "pval",
  "logFC"       = "log2fc",
  "RawlogFC"    = "rawlog2fc",
  "Log.P.Value" = "logpval",
  "AveExpr"     = "aveExpr",
  "t"           = "t"
)

# -------------------------------
# Main Loop: Write Forward & Reverse CSVs
# -------------------------------

written_index <- purrr::imap_dfr(by_comparison, function(cols, comp_key) {

  df_out <- data |>
    dplyr::select(id, dplyr::all_of(cols)) |>
    dplyr::mutate(id = as.character(id))

  names(df_out)[1] <- "gene_symbol"

  metrics <- metric_keys[cols]

  new_names <- vapply(metrics, function(m) {
    if (m %in% names(recode_map)) {
      recode_map[[m]]
    } else {
      m
    }
  }, character(1))

  new_names <- make.unique(new_names)
  names(df_out)[-1] <- new_names

  comp2 <- parse_compkey(comp_key)

  fwd_file <- file.path(
    outdir_fwd,
    paste0(safe_name(comp2), ".csv")
  )

  utils::write.csv(
    df_out,
    fwd_file,
    row.names = FALSE,
    quote = TRUE
  )

  message("Wrote: ", fwd_file)

  rev_file <- NA_character_
  rev_comp <- NA_character_

  if (stringr::str_detect(comp_key, "\\.over\\.")) {

    df_rev <- df_out

    log_cols <- names(df_rev)[
      stringr::str_detect(
        names(df_rev),
        stringr::regex("log.*fc", ignore_case = TRUE)
      )
    ]

    for (col in log_cols) {
      df_rev[[col]] <- suppressWarnings(as.numeric(df_rev[[col]]) * -1)
    }

    m <- stringr::str_match(
      comp_key,
      "^([A-Za-z0-9]+)_([a-z]+)_([1234])\\.over\\.([A-Za-z0-9]+)_([a-z]+)_([1234])$"
    )

    if (!is.na(m[1, 1])) {
      r1 <- m[1, 2]
      g1 <- m[1, 3]
      l1 <- m[1, 4]

      r2 <- m[1, 5]
      g2 <- m[1, 6]
      l2 <- m[1, 7]

      rev_comp <- format_comparison(r2, g2, l2, r1, g1, l1)

    } else {
      rev_comp <- swap_comparison(comp_key)
      for (code in names(labels$condition_code_map)) {
        rev_comp <- stringr::str_replace_all(rev_comp, paste0("_", code, "(?=($|_|\\.))"), paste0("_", labels$condition_code_map[[code]]))
      }
      rev_comp <- stringr::str_replace_all(rev_comp, "[^A-Za-z0-9_]", "")
    }

    rev_file <- file.path(
      outdir_rev,
      paste0(safe_name(rev_comp), ".csv")
    )

    utils::write.csv(
      df_rev,
      rev_file,
      row.names = FALSE,
      quote = TRUE
    )

    message("Wrote reversed: ", rev_file)
  }

  tibble::tibble(
    comparison = comp_key,
    parsed_forward_comparison = comp2,
    parsed_reverse_comparison = rev_comp,
    n_columns = length(cols),
    columns_used = paste(cols, collapse = ";"),
    forward_file = fwd_file,
    reverse_file = rev_file
  )
})

# -------------------------------
# Index File
# -------------------------------

readr::write_csv(
  written_index,
  file.path(outdir, "indexComparisons.csv")
)

message("Finished successfully.")
