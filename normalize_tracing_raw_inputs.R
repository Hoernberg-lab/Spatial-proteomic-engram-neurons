# ================================================================
# Normalize raw tracing/cFos input files before downstream analysis
# ================================================================
# Purpose:
#   Fix inconsistent upstream region metadata in the raw input files used by
#   tracing_cfos_correlation.r, especially duplicated/truncated Class/Level
#   metadata for the same Allen-style region abbreviation.
#
# Important logic:
#   1) Read CellCountNoOutliersFin.csv and IntensityNoOutliersFin.xlsx.
#   2) Optionally read FINAL.csv as a metadata-only reference if present.
#      This is useful because CellCountNoOutliersFin.csv can contain regions
#      that are not present in IntensityNoOutliersFin.xlsx, so those regions
#      cannot be repaired from intensity alone.
#   3) Build one canonical metadata row per abbreviation.
#   4) Fill missing/truncated Annotation, Class, and Level in the raw files.
#
# Output:
#   raw_input_normalized/CellCountNoOutliersFin_normalized.csv
#   raw_input_normalized/IntensityNoOutliersFin_normalized.xlsx
#   raw_input_normalized/IntensityNoOutliersFin_normalized.csv
#   raw_input_normalized/QC_region_metadata_corrections.csv
#   raw_input_normalized/QC_abbreviation_canonical_reference.csv
#   raw_input_normalized/QC_abbreviation_ambiguous_candidates.csv
#   raw_input_normalized/QC_missing_metadata_after_normalization.csv
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
final_reference_file <- file.path(input_dir, "FINAL.csv")

out_dir <- file.path(input_dir, "raw_input_normalized")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Safety default: never overwrite originals unless explicitly changed.
overwrite_originals <- FALSE

# If TRUE, missing Class values that cannot be repaired from any reference are
# set to "Unknown". Level is still left as NA because fabricating anatomical
# depth would be worse than explicit missingness.
fill_unresolved_class_with_unknown <- TRUE

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
    stringr::str_detect(x_low, "^(cortica|hypothala|striat|pallid|thalam|midbrain|hindbrain|cerebell)") |
    stringr::str_length(x_low) < 5
}

metadata_quality_score <- function(annotation, class, level) {
  annotation <- clean_chr(annotation)
  class <- clean_chr(class)

  100 * !is.na(level) +
    50 * !is_bad_class(class) +
    0.10 * stringr::str_length(dplyr::coalesce(class, "")) +
    0.01 * stringr::str_length(dplyr::coalesce(annotation, ""))
}

changed_chr <- function(before, after) {
  dplyr::coalesce(as.character(before), "<NA>") != dplyr::coalesce(as.character(after), "<NA>")
}

changed_int <- function(before, after) {
  dplyr::coalesce(as.integer(before), -999999L) != dplyr::coalesce(as.integer(after), -999999L)
}

read_table_flexible <- function(path) {
  ext <- tolower(fs::path_ext(path))
  switch(
    ext,
    csv  = readr::read_csv(path, show_col_types = FALSE, guess_max = 100000),
    txt  = readr::read_csv(path, show_col_types = FALSE, guess_max = 100000),
    xlsx = openxlsx::read.xlsx(path),
    stop("Unsupported file extension: ", ext, ". Expected csv, txt, or xlsx.")
  )
}

standardize_region_columns <- function(x, source_file, file_role) {
  x <- x %>% janitor::clean_names()

  # Common variants.
  names(x) <- stringr::str_replace_all(names(x), "^cellcount$", "cell_count")
  names(x) <- stringr::str_replace_all(names(x), "^abbrev$", "abbreviation")
  names(x) <- stringr::str_replace_all(names(x), "^abbr$", "abbreviation")

  if (!"annotation" %in% names(x)) {
    stop("Missing required column 'annotation' after clean_names() in ", source_file)
  }
  if (!"abbreviation" %in% names(x)) x <- x %>% mutate(abbreviation = annotation)
  if (!"class" %in% names(x)) x <- x %>% mutate(class = NA_character_)
  if (!"level" %in% names(x)) x <- x %>% mutate(level = NA_integer_)

  x %>%
    mutate(
      .source_file = source_file,
      .file_role = file_role,
      annotation = clean_chr(annotation),
      abbreviation = clean_chr(abbreviation),
      class = clean_chr(class),
      level = suppressWarnings(as.integer(level)),
      .annotation_code_clean = clean_abbrev(annotation),
      .abbreviation_clean = clean_abbrev(abbreviation),
      .annotation_clean = stringr::str_to_lower(clean_chr(annotation)),
      .metadata_score_before = metadata_quality_score(annotation, class, level)
    )
}

read_tracing_raw_file <- function(path, file_role = c("cell_count", "intensity")) {
  file_role <- match.arg(file_role)
  if (!fs::file_exists(path)) stop("Input file not found: ", path)

  message("Reading ", file_role, " file: ", path)
  read_table_flexible(path) %>%
    standardize_region_columns(source_file = fs::path_file(path), file_role = file_role)
}

read_metadata_reference_file <- function(path) {
  if (!fs::file_exists(path)) {
    message("No FINAL.csv metadata reference found at: ", path)
    return(tibble())
  }

  message("Reading metadata reference file: ", path)
  read_table_flexible(path) %>%
    standardize_region_columns(source_file = fs::path_file(path), file_role = "metadata_reference") %>%
    select(
      .source_file, .file_role,
      annotation, abbreviation, class, level,
      .annotation_code_clean, .abbreviation_clean, .annotation_clean, .metadata_score_before
    ) %>%
    distinct()
}

make_canonical_reference <- function(combined_metadata) {
  candidates <- combined_metadata %>%
    filter(!is.na(.annotation_code_clean), .annotation_code_clean != "") %>%
    group_by(.annotation_code_clean, annotation, abbreviation, class, level) %>%
    summarise(
      n_rows = n(),
      n_files = n_distinct(.source_file),
      source_files = paste(sort(unique(.source_file)), collapse = "; "),
      source_roles = paste(sort(unique(.file_role)), collapse = "; "),
      score = first(metadata_quality_score(annotation, class, level)),
      .groups = "drop"
    ) %>%
    mutate(
      # Strongly prefer complete ontology metadata, then repeated evidence.
      total_score = score + log1p(n_rows) + 0.5 * n_files
    )

  canonical <- candidates %>%
    arrange(.annotation_code_clean, desc(total_score), desc(n_rows), desc(score)) %>%
    group_by(.annotation_code_clean) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(
      .annotation_code_clean,
      canonical_annotation = annotation,
      canonical_abbreviation = abbreviation,
      canonical_class = class,
      canonical_level = level,
      canonical_source_files = source_files,
      canonical_source_roles = source_roles,
      canonical_n_rows = n_rows,
      canonical_score = score,
      canonical_total_score = total_score
    )

  ambiguous <- candidates %>%
    group_by(.annotation_code_clean) %>%
    mutate(
      n_candidate_metadata_rows = n(),
      best_total_score = max(total_score, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    filter(n_candidate_metadata_rows > 1) %>%
    arrange(.annotation_code_clean, desc(total_score), desc(n_rows))

  list(canonical = canonical, ambiguous = ambiguous)
}

normalize_region_metadata <- function(raw_df, canonical_ref) {
  out <- raw_df %>%
    left_join(canonical_ref, by = ".annotation_code_clean") %>%
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
      # non-missing Level unless manually reviewed.
      level = coalesce(level, canonical_level)
    )

  if (fill_unresolved_class_with_unknown) {
    out <- out %>% mutate(class = if_else(is.na(class) | class == "", "Unknown", class))
  }

  out %>%
    mutate(
      .metadata_score_after = metadata_quality_score(annotation, class, level),
      .metadata_changed =
        changed_chr(.annotation_before, annotation) |
        changed_chr(.abbreviation_before, abbreviation) |
        changed_chr(.class_before, class) |
        changed_int(.level_before, level),
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
      canonical_source_files,
      canonical_source_roles,
      metadata_score_before = .metadata_score_before,
      metadata_score_after = .metadata_score_after
    ) %>%
    distinct()

  missing_after <- bind_rows(cell_norm, intensity_norm) %>%
    filter(is.na(class) | class == "" | is.na(level)) %>%
    transmute(
      source_file = .source_file,
      file_role = .file_role,
      annotation_code_clean = .annotation_code_clean,
      abbreviation_clean = .abbreviation_clean,
      annotation,
      abbreviation,
      class,
      level,
      canonical_source_files,
      canonical_source_roles
    ) %>%
    distinct() %>%
    arrange(file_role, abbreviation_clean)

  readr::write_csv(correction_log, file.path(out_dir, "QC_region_metadata_corrections.csv"))
  readr::write_csv(canonical_ref, file.path(out_dir, "QC_abbreviation_canonical_reference.csv"))
  readr::write_csv(ambiguous_ref, file.path(out_dir, "QC_abbreviation_ambiguous_candidates.csv"))
  readr::write_csv(missing_after, file.path(out_dir, "QC_missing_metadata_after_normalization.csv"))
  readr::write_tsv(missing_after, file.path(out_dir, "QC_missing_metadata_after_normalization.tsv"))
  openxlsx::write.xlsx(missing_after, file.path(out_dir, "QC_missing_metadata_after_normalization.xlsx"), overwrite = TRUE)

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
  message("Remaining missing metadata QC: ", file.path(out_dir, "QC_missing_metadata_after_normalization.csv"))
}

# -----------------------------
# 3. Run normalization
# -----------------------------
cell_raw <- read_tracing_raw_file(cell_count_raw_file, file_role = "cell_count")
intensity_raw <- read_tracing_raw_file(intensity_raw_file, file_role = "intensity")
final_ref <- read_metadata_reference_file(final_reference_file)

combined_metadata <- bind_rows(cell_raw, intensity_raw, final_ref)
refs <- make_canonical_reference(combined_metadata)

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

if (any(summary_table$n_missing_level_after > 0)) {
  message(
    "Some Level values are still missing. Inspect QC_missing_metadata_after_normalization.csv. ",
    "These are likely abbreviations that are absent from both IntensityNoOutliersFin.xlsx and FINAL.csv metadata."
  )
}
