# ================================================================
# Normalize raw tracing/cFos input files before downstream analysis
# ================================================================
# Purpose:
#   Fix inconsistent upstream region metadata in the raw input files used by
#   tracing_cfos_correlation.r, especially duplicated/truncated Class/Level
#   metadata for the same Allen-style region abbreviation.
#
# Input examples:
#   CellCountNoOutliersFin.csv
#   IntensityNoOutliersFin.xlsx
#
# Output:
#   raw_input_normalized/CellCountNoOutliersFin_normalized.csv
#   raw_input_normalized/IntensityNoOutliersFin_normalized.xlsx
#   raw_input_normalized/IntensityNoOutliersFin_normalized.csv
#   raw_input_normalized/QC_region_metadata_corrections.csv
#   raw_input_normalized/QC_abbreviation_canonical_reference.csv
#   raw_input_normalized/QC_abbreviation_ambiguous_candidates.csv
#
# Recommended downstream use:
#   In tracing_cfos_correlation.r, point the raw input paths to the normalized
#   outputs instead of the original raw files.
# ================================================================

packages <- c("tidyverse", "readr", "openxlsx", "janitor", "stringr", "fs")

install_if_missing <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) install.packages(missing)
}
install_if_missing(packages)
invisible(lapply(packages, library, character.only = TRUE))

# -----------------------------
# 1. User settings
# -----------------------------
input_dir <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/Results Final tracing"

cell_count_raw_file <- file.path(input_dir, "CellCountNoOutliersFin.csv")
intensity_raw_file  <- file.path(input_dir, "IntensityNoOutliersFin.xlsx")

out_dir <- file.path(input_dir, "raw_input_normalized")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Safety default: never overwrite originals unless explicitly changed.
overwrite_originals <- FALSE

# -----------------------------
# 2. Helpers
# -----------------------------
clean_chr <- function(x) {
  x <- as.character(x)
  x <- stringr::str_squish(x)
  dplyr::na_if(x, "")
}

clean_abbrev <- function(x) {
  clean_chr(x) %>%
    stringr::str_replace_all("[^A-Za-z0-9]", "")
}

is_bad_class <- function(x) {
  x_low <- stringr::str_to_lower(clean_chr(x))
  is.na(x_low) |
    x_low %in% c("unknown", "na", "nan", "null", "cortic", "cortical") |
    stringr::str_length(x_low) < 5
}

metadata_quality_score <- function(annotation, class, level) {
  annotation <- clean_chr(annotation)
  class <- clean_chr(class)

  dplyr::case_when(
    TRUE ~
      100 * !is.na(level) +
      50 * !is_bad_class(class) +
      0.10 * stringr::str_length(dplyr::coalesce(class, "")) +
      0.01 * stringr::str_length(dplyr::coalesce(annotation, ""))
  )
}

first_nonmissing <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) NA else x[[1]]
}

read_tracing_raw_file <- function(path, file_role = c("cell_count", "intensity")) {
  file_role <- match.arg(file_role)
  if (!fs::file_exists(path)) stop("Input file not found: ", path)

  ext <- tolower(fs::path_ext(path))
  message("Reading ", file_role, " file: ", path)

  x <- switch(
    ext,
    csv  = readr::read_csv(path, show_col_types = FALSE, guess_max = 100000),
    txt  = readr::read_csv(path, show_col_types = FALSE, guess_max = 100000),
    xlsx = openxlsx::read.xlsx(path),
    stop("Unsupported file extension: ", ext, ". Expected csv, txt, or xlsx.")
  )

  x <- x %>% janitor::clean_names()

  # Handle minimal two-column CSVs such as Annotation,Cell_Count if needed.
  # For full files, this leaves all existing columns intact.
  names(x) <- stringr::str_replace_all(names(x), "^cellcount$", "cell_count")

  if (!"annotation" %in% names(x)) {
    stop("Missing required column 'annotation' after clean_names() in ", fs::path_file(path))
  }

  if (!"abbreviation" %in% names(x)) {
    # Fall back to Annotation when the raw table only has region abbreviations.
    x <- x %>% mutate(abbreviation = annotation)
  }
  if (!"class" %in% names(x)) x <- x %>% mutate(class = NA_character_)
  if (!"level" %in% names(x)) x <- x %>% mutate(level = NA_integer_)

  x %>%
    mutate(
      .source_file = fs::path_file(path),
      .file_role = file_role,
      annotation = clean_chr(annotation),
      abbreviation = clean_chr(abbreviation),
      class = clean_chr(class),
      level = suppressWarnings(as.integer(level)),
      .abbreviation_clean = clean_abbrev(abbreviation),
      .annotation_clean = stringr::str_to_lower(clean_chr(annotation)),
      .metadata_score_before = metadata_quality_score(annotation, class, level)
    )
}

make_canonical_reference <- function(combined_raw) {
  candidates <- combined_raw %>%
    filter(!is.na(.abbreviation_clean), .abbreviation_clean != "") %>%
    group_by(.abbreviation_clean, annotation, abbreviation, class, level) %>%
    summarise(
      n_rows = n(),
      n_files = n_distinct(.source_file),
      score = first(metadata_quality_score(annotation, class, level)),
      .groups = "drop"
    ) %>%
    mutate(
      total_score = score + log1p(n_rows) + 0.5 * n_files
    )

  canonical <- candidates %>%
    arrange(.abbreviation_clean, desc(total_score), desc(n_rows), desc(score)) %>%
    group_by(.abbreviation_clean) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(
      .abbreviation_clean,
      canonical_annotation = annotation,
      canonical_abbreviation = abbreviation,
      canonical_class = class,
      canonical_level = level,
      canonical_n_rows = n_rows,
      canonical_score = score,
      canonical_total_score = total_score
    )

  ambiguous <- candidates %>%
    group_by(.abbreviation_clean) %>%
    mutate(
      n_candidate_metadata_rows = n(),
      best_total_score = max(total_score, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    filter(n_candidate_metadata_rows > 1) %>%
    arrange(.abbreviation_clean, desc(total_score), desc(n_rows))

  list(canonical = canonical, ambiguous = ambiguous)
}

normalize_region_metadata <- function(raw_df, canonical_ref) {
  raw_df %>%
    left_join(canonical_ref, by = ".abbreviation_clean") %>%
    mutate(
      .class_before = class,
      .annotation_before = annotation,
      .abbreviation_before = abbreviation,
      .level_before = level,

      # Replace Annotation only if missing or if it is only the abbreviation while
      # a fuller canonical annotation exists.
      annotation = case_when(
        is.na(annotation) ~ canonical_annotation,
        !is.na(canonical_annotation) & stringr::str_to_lower(annotation) == stringr::str_to_lower(abbreviation) ~ canonical_annotation,
        TRUE ~ annotation
      ),

      abbreviation = coalesce(abbreviation, canonical_abbreviation),

      # Replace only clearly bad/truncated classes, or fill missing classes.
      class = case_when(
        is_bad_class(class) & !is.na(canonical_class) ~ canonical_class,
        TRUE ~ class
      ),

      # Fill missing Level from the canonical reference. Do not overwrite a
      # non-missing Level unless you decide to do so manually after QC.
      level = coalesce(level, canonical_level),

      .metadata_score_after = metadata_quality_score(annotation, class, level),
      .metadata_changed =
        !identical(.annotation_before, annotation) |
        !identical(.abbreviation_before, abbreviation) |
        !identical(.class_before, class) |
        !identical(.level_before, level),

      RegionLabel = paste0(annotation, " | ", abbreviation),
      RegionKey = paste(class, RegionLabel, sep = " :: ")
    )
}

write_normalized_outputs <- function(cell_norm, intensity_norm, canonical_ref, ambiguous_ref) {
  correction_log <- bind_rows(cell_norm, intensity_norm) %>%
    filter(.metadata_changed | .metadata_score_after != .metadata_score_before) %>%
    transmute(
      source_file = .source_file,
      file_role = .file_role,
      abbreviation_clean = .abbreviation_clean,
      annotation_before = .annotation_before,
      annotation_after = annotation,
      abbreviation_before = .abbreviation_before,
      abbreviation_after = abbreviation,
      class_before = .class_before,
      class_after = class,
      level_before = .level_before,
      level_after = level,
      metadata_score_before = .metadata_score_before,
      metadata_score_after = .metadata_score_after
    ) %>%
    distinct()

  readr::write_csv(correction_log, file.path(out_dir, "QC_region_metadata_corrections.csv"))
  readr::write_csv(canonical_ref, file.path(out_dir, "QC_abbreviation_canonical_reference.csv"))
  readr::write_csv(ambiguous_ref, file.path(out_dir, "QC_abbreviation_ambiguous_candidates.csv"))

  drop_internal <- function(x) {
    x %>% select(-starts_with("."), -starts_with("canonical_"))
  }

  cell_out <- drop_internal(cell_norm)
  intensity_out <- drop_internal(intensity_norm)

  readr::write_csv(cell_out, file.path(out_dir, "CellCountNoOutliersFin_normalized.csv"))
  readr::write_csv(intensity_out, file.path(out_dir, "IntensityNoOutliersFin_normalized.csv"))
  openxlsx::write.xlsx(intensity_out, file.path(out_dir, "IntensityNoOutliersFin_normalized.xlsx"), overwrite = TRUE)

  if (overwrite_originals) {
    backup_dir <- file.path(input_dir, "raw_input_original_backup")
    dir.create(backup_dir, recursive = TRUE, showWarnings = FALSE)

    file.copy(cell_count_raw_file, file.path(backup_dir, fs::path_file(cell_count_raw_file)), overwrite = TRUE)
    file.copy(intensity_raw_file, file.path(backup_dir, fs::path_file(intensity_raw_file)), overwrite = TRUE)

    readr::write_csv(cell_out, cell_count_raw_file)
    openxlsx::write.xlsx(intensity_out, intensity_raw_file, overwrite = TRUE)
  }

  message("Done. Normalized files written to: ", out_dir)
  message("Corrections log: ", file.path(out_dir, "QC_region_metadata_corrections.csv"))
}

# -----------------------------
# 3. Run normalization
# -----------------------------
cell_raw <- read_tracing_raw_file(cell_count_raw_file, file_role = "cell_count")
intensity_raw <- read_tracing_raw_file(intensity_raw_file, file_role = "intensity")

combined_raw <- bind_rows(cell_raw, intensity_raw)
refs <- make_canonical_reference(combined_raw)

cell_norm <- normalize_region_metadata(cell_raw, refs$canonical)
intensity_norm <- normalize_region_metadata(intensity_raw, refs$canonical)

write_normalized_outputs(cell_norm, intensity_norm, refs$canonical, refs$ambiguous)

# -----------------------------
# 4. Quick console summary
# -----------------------------
summary_table <- bind_rows(cell_norm, intensity_norm) %>%
  group_by(.file_role) %>%
  summarise(
    n_rows = n(),
    n_regions = n_distinct(.abbreviation_clean),
    n_missing_class_after = sum(is.na(class) | class == ""),
    n_missing_level_after = sum(is.na(level)),
    n_changed_rows = sum(.metadata_changed | .metadata_score_after != .metadata_score_before),
    .groups = "drop"
  )

print(summary_table)
