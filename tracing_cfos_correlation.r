# ================================================================
# Final tracing region analysis - individual-level FINAL.csv
# Paper framing:
#   VEH_paired   = associative learning condition
#   VEH_unpaired = stress / non-associative aversive exposure
#   CNO_paired   = CeM manipulation during learning
#   CNO_unpaired = CeM manipulation during stress / unpaired exposure
# ================================================================

# -----------------------------
# 0. Setup
# -----------------------------
packages <- c(
  "tidyverse", "readr", "janitor", "stringr", "fs", "openxlsx",
  "pheatmap", "ComplexHeatmap", "circlize", "RColorBrewer",
  "igraph", "ggraph", "tidygraph", "limma", "scales", "ggrepel",
  "uwot", "patchwork", "svglite"
)

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
input_file <- file.path(input_dir, "FINAL.csv")
intensity_raw_file <- file.path(input_dir, "raw_input_normalized", "IntensityNoOutliersFin_normalized.xlsx")
cell_count_raw_file <- file.path(input_dir, "raw_input_normalized", "CellCountNoOutliersFin_normalized.csv")

out_dir <- file.path(dirname(input_file), "region_analysis_outputs_individual_level")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

results_root <- file.path(out_dir, "results")
publication_fig4_fig_dir <- file.path(results_root, "01_publication", "figures", "fig4")
publication_fig4_panel_fig_dir <- file.path(publication_fig4_fig_dir, "individual_panels")
publication_fig4_composite_fig_dir <- file.path(publication_fig4_fig_dir, "composite")
publication_fig4_key_fig_dir <- file.path(publication_fig4_fig_dir, "keys")
publication_supp_fig_dir <- file.path(results_root, "01_publication", "figures", "supplementary")
publication_dashboard_fig_dir <- file.path(results_root, "01_publication", "figures", "dashboards")
publication_fig4_tab_dir <- file.path(results_root, "01_publication", "tables", "fig4")
publication_supp_tab_dir <- file.path(results_root, "01_publication", "tables", "supplementary")
publication_source_tab_dir <- file.path(results_root, "01_publication", "tables", "source_data")
publication_manifest_dir <- file.path(results_root, "01_publication", "manifests")
exploratory_effect_fig_dir <- file.path(results_root, "02_exploratory", "effect_size_maps")
exploratory_profile_fig_dir <- file.path(results_root, "02_exploratory", "regional_profiles")
exploratory_covariance_dir <- file.path(results_root, "02_exploratory", "covariance")
exploratory_network_dir <- file.path(results_root, "02_exploratory", "network_analysis")
exploratory_dimred_dir <- file.path(results_root, "02_exploratory", "dimensionality_reduction")
exploratory_sensitivity_dir <- file.path(results_root, "02_exploratory", "sensitivity_analysis")
qc_missingness_dir <- file.path(results_root, "03_qc", "missingness")
qc_normalization_dir <- file.path(results_root, "03_qc", "normalization")
qc_region_selection_dir <- file.path(results_root, "03_qc", "region_selection")
qc_covariance_dir <- file.path(results_root, "03_qc", "covariance_qc")
legacy_dir <- file.path(results_root, "04_legacy")
legacy_fig_dir <- file.path(legacy_dir, "figures")
legacy_tab_dir <- file.path(legacy_dir, "tables")
session_info_dir <- file.path(results_root, "session_info")
purrr::walk(
  c(
    publication_fig4_fig_dir, publication_fig4_panel_fig_dir, publication_fig4_composite_fig_dir,
    publication_fig4_key_fig_dir, publication_supp_fig_dir, publication_dashboard_fig_dir,
    publication_fig4_tab_dir, publication_supp_tab_dir, publication_source_tab_dir,
    publication_manifest_dir, exploratory_effect_fig_dir, exploratory_profile_fig_dir,
    exploratory_covariance_dir, exploratory_network_dir,
    exploratory_dimred_dir, exploratory_sensitivity_dir, qc_missingness_dir,
    qc_normalization_dir, qc_region_selection_dir, qc_covariance_dir,
    legacy_dir, legacy_fig_dir, legacy_tab_dir, session_info_dir
  ),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
)

metrics_to_analyse <- c("Cell_Count", "Intensity")

condition_lookup <- tibble::tribble(
  ~Group, ~Condition,
  1, "CNO_paired",
  2, "VEH_paired",
  3, "CNO_unpaired",
  4, "VEH_unpaired"
)

# Biologically prioritized contrasts for the paper.
contrast_definitions <- c(
  Learning_effect =
    "VEH_paired - VEH_unpaired",

  CeM_manipulation_during_learning =
    "CNO_paired - VEH_paired",

  CeM_manipulation_during_stress =
    "CNO_unpaired - VEH_unpaired",

  Learning_x_CeM_interaction =
    "(CNO_paired - VEH_paired) - (CNO_unpaired - VEH_unpaired)",

  Paired_vs_unpaired =
    "((VEH_paired + CNO_paired)/2) - ((VEH_unpaired + CNO_unpaired)/2)",

  CNO_vs_VEH =
    "((CNO_paired + CNO_unpaired)/2) - ((VEH_paired + VEH_unpaired)/2)"
)

central_contrasts <- c(
  "Learning_effect",
  "Learning_x_CeM_interaction",
  "CeM_manipulation_during_learning",
  "CeM_manipulation_during_stress"
)

min_pairwise_n <- 3
network_abs_r_cutoff <- 0.70
network_fdr_cutoff <- 0.10
covariance_mode <- "pooled_raw"
supported_covariance_modes <- c("pooled_raw", "residualized_by_condition", "per_condition")
if (!covariance_mode %in% supported_covariance_modes) {
  stop("Unsupported covariance_mode: ", covariance_mode)
}
# TODO: Add covariance_mode == "residualized_by_condition" by correlating
# per-region residuals from log1p(metric) ~ Condition.
# TODO: Add covariance_mode == "per_condition" by running the covariance/network
# workflow separately within each condition when pairwise n is sufficient.

# -----------------------------
# 2. Readers for merged FINAL.csv or raw no-outlier Excel files
# -----------------------------
read_final_tracing_csv <- function(path) {
  message("Reading merged individual-level file: ", path)

  x <- readr::read_csv(path, show_col_types = FALSE, guess_max = 100000) %>%
    janitor::clean_names()

  required <- c(
    "annotation", "cell_count", "intensity", "group",
    "abbreviation", "class", "level", "sample_id"
  )

  missing_cols <- setdiff(required, names(x))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in FINAL.csv: ", paste(missing_cols, collapse = ", "))
  }

  x %>%
    mutate(
      SourceFile = fs::path_file(path),
      FileID = "FINAL",
      Annotation = as.character(annotation),
      Cell_Count = suppressWarnings(as.numeric(cell_count)),
      Intensity = suppressWarnings(as.numeric(intensity)),
      Group = suppressWarnings(as.integer(group)),
      Abbreviation = stringr::str_squish(as.character(abbreviation)),
      Class = stringr::str_squish(stringr::str_to_title(as.character(class))),
      Level = suppressWarnings(as.integer(level)),
      Sample = as.character(sample_id),
      Animal = as.character(sample_id),
      SampleID = paste(Animal, paste0("G", Group), sep = "__")
    ) %>%
    select(
      SourceFile, FileID, SampleID, Sample, Animal, Group,
      Annotation, Abbreviation, Class, Level, Cell_Count, Intensity
    )
}

read_raw_metric_file <- function(path, metric_col) {
  message("Reading raw metric file: ", path)

  ext <- tolower(fs::path_ext(path))
  x <- switch(
    ext,
    xlsx = openxlsx::read.xlsx(path),
    csv = readr::read_csv(path, show_col_types = FALSE, guess_max = 100000),
    stop("Unsupported raw metric file type for ", path, ". Expected .xlsx or .csv.")
  ) %>%
    janitor::clean_names()

  names(x) <- stringr::str_replace_all(names(x), "^abbrev$", "abbreviation")
  names(x) <- stringr::str_replace_all(names(x), "^abbr$", "abbreviation")
  names(x) <- stringr::str_replace_all(names(x), "^cellcount$", "cell_count")

  required <- c("annotation", metric_col, "group", "abbreviation", "class", "sample_id")
  missing_cols <- setdiff(required, names(x))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in ",
      fs::path_file(path),
      ": ",
      paste(missing_cols, collapse = ", ")
    )
  }

  first_present <- function(x) {
    x <- x[!is.na(x) & x != ""]
    if (length(x) == 0) NA else x[[1]]
  }

  x %>%
    mutate(
      source_file = fs::path_file(path),
      annotation = as.character(annotation),
      group = suppressWarnings(as.integer(group)),
      sample_id = as.character(sample_id),
      condition_raw = if ("condition" %in% names(.)) as.character(condition) else NA_character_,
      treatment_raw = if ("treatment" %in% names(.)) as.character(treatment) else NA_character_,
      abbreviation = stringr::str_squish(as.character(abbreviation)),
      class = stringr::str_squish(stringr::str_to_title(as.character(class))),
      level = if ("level" %in% names(.)) suppressWarnings(as.integer(level)) else NA_integer_,
      value = suppressWarnings(as.numeric(.data[[metric_col]]))
    ) %>%
    select(
      source_file, annotation, group, sample_id, condition_raw, treatment_raw,
      abbreviation, class, level, value
    ) %>%
    group_by(annotation, group, sample_id, condition_raw, treatment_raw) %>%
    summarise(
      source_file = paste(sort(unique(source_file)), collapse = "; "),
      abbreviation = first_present(abbreviation),
      class = first_present(class),
      level = {
        level_values <- level[!is.na(level)]
        if (length(level_values) == 0) NA_integer_ else level_values[[1]]
      },
      value = mean(value, na.rm = TRUE),
      n_raw_rows_collapsed = n(),
      .groups = "drop"
    ) %>%
    mutate(
      value = if_else(is.nan(value), NA_real_, value)
    )
}

read_raw_tracing_metrics <- function(intensity_path, cell_count_path) {
  intensity <- read_raw_metric_file(intensity_path, "intensity") %>%
    rename(
      intensity_source_file = source_file,
      intensity_abbreviation = abbreviation,
      intensity_class = class,
      intensity_level = level,
      Intensity = value
    )

  cell_count <- read_raw_metric_file(cell_count_path, "cell_count") %>%
    rename(
      cell_count_source_file = source_file,
      cell_count_abbreviation = abbreviation,
      cell_count_class = class,
      cell_count_level = level,
      Cell_Count = value
    )

  join_cols <- c("annotation", "group", "sample_id", "condition_raw", "treatment_raw")

  full_join(intensity, cell_count, by = join_cols) %>%
    mutate(
      SourceFile = stringr::str_c(
        coalesce(intensity_source_file, ""),
        coalesce(cell_count_source_file, ""),
        sep = "; "
      ) %>%
        stringr::str_replace_all("(^;\\s*|;\\s*$)", "") %>%
        stringr::str_replace_all(";\\s*;", ";"),
      FileID = "RawNoOutliers",
      Annotation = as.character(annotation),
      Group = suppressWarnings(as.integer(group)),
      Abbreviation = coalesce(intensity_abbreviation, cell_count_abbreviation),
      Class = coalesce(intensity_class, cell_count_class),
      Level = coalesce(intensity_level, cell_count_level),
      Sample = as.character(sample_id),
      Animal = as.character(sample_id),
      SampleID = paste(Animal, paste0("G", Group), sep = "__")
    ) %>%
    select(
      SourceFile, FileID, SampleID, Sample, Animal, Group,
      Annotation, Abbreviation, Class, Level, Cell_Count, Intensity
    )
}

read_tracing_input <- function(final_path, intensity_path, cell_count_path) {
  if (fs::file_exists(intensity_path) && fs::file_exists(cell_count_path)) {
    message("Using raw no-outlier metric files as analysis input.")
    return(read_raw_tracing_metrics(intensity_path, cell_count_path))
  }

  warning(
    "Raw no-outlier metric files were not both found. Falling back to FINAL.csv: ",
    final_path
  )
  read_final_tracing_csv(final_path)
}

# -----------------------------
# 3. Read tracing data and prepare long table
# -----------------------------
raw_long <- read_tracing_input(input_file, intensity_raw_file, cell_count_raw_file) %>%
  left_join(condition_lookup, by = "Group") %>%
  mutate(
    Condition = factor(
      Condition,
      levels = c("CNO_paired", "VEH_paired", "CNO_unpaired", "VEH_unpaired")
    ),
    RegionLabel = paste0(Annotation, " | ", Abbreviation),
    Class = if_else(is.na(Class) | Class == "", "Unknown", Class),
    RegionKey = paste(Class, RegionLabel, sep = " :: ")
  )

unknown_group <- raw_long %>% filter(is.na(Condition))
if (nrow(unknown_group) > 0) {
  warning("Some rows have Group values outside 1:4. See QC_unknown_group_rows.xlsx")
  openxlsx::write.xlsx(
    unknown_group,
    file.path(out_dir, "tables", "QC_unknown_group_rows.xlsx"),
    overwrite = TRUE
  )
}

raw_long <- raw_long %>% filter(!is.na(Condition))

long <- raw_long %>%
  group_by(
    SampleID, Animal, Sample, Group, Condition,
    Annotation, Abbreviation, Class, Level, RegionLabel, RegionKey
  ) %>%
  summarise(
    SourceFile = paste(sort(unique(SourceFile)), collapse = "; "),
    Cell_Count = mean(Cell_Count, na.rm = TRUE),
    Intensity = mean(Intensity, na.rm = TRUE),
    n_rows_collapsed = n(),
    .groups = "drop"
  ) %>%
  mutate(
    Cell_Count = if_else(is.nan(Cell_Count), NA_real_, Cell_Count),
    Intensity = if_else(is.nan(Intensity), NA_real_, Intensity)
  ) %>%
  group_by(SampleID) %>%
  mutate(
    total_cell_count_per_animal = sum(Cell_Count, na.rm = TRUE),
    total_intensity_per_animal = sum(Intensity, na.rm = TRUE),
    Cell_Count_norm = if_else(total_cell_count_per_animal > 0, Cell_Count / total_cell_count_per_animal, NA_real_),
    Intensity_norm = if_else(total_intensity_per_animal > 0, Intensity / total_intensity_per_animal, NA_real_)
  ) %>%
  ungroup() %>%
  mutate(
    Cell_Count_norm = if_else(is.nan(Cell_Count_norm), NA_real_, Cell_Count_norm),
    Intensity_norm = if_else(is.nan(Intensity_norm), NA_real_, Intensity_norm)
  )

openxlsx::write.xlsx(long, file.path(out_dir, "tables", "merged_long_region_data.xlsx"), overwrite = TRUE)
readr::write_csv(long, file.path(out_dir, "tables", "merged_long_region_data.csv"))
openxlsx::write.xlsx(long, file.path(legacy_tab_dir, "merged_long_region_data.xlsx"), overwrite = TRUE)
readr::write_csv(long, file.path(legacy_tab_dir, "merged_long_region_data.csv"))

# Animal counts sanity check
animal_counts <- long %>%
  distinct(SampleID, Group, Condition) %>%
  count(Group, Condition, name = "n_animals")

openxlsx::write.xlsx(animal_counts, file.path(out_dir, "tables", "QC_animal_counts_by_condition.xlsx"), overwrite = TRUE)
openxlsx::write.xlsx(animal_counts, file.path(legacy_tab_dir, "QC_animal_counts_by_condition.xlsx"), overwrite = TRUE)

animal_level_qc <- long %>%
  group_by(SampleID, Animal, Sample, Group, Condition) %>%
  summarise(
    total_cell_count_per_animal = first(total_cell_count_per_animal),
    total_intensity_per_animal = first(total_intensity_per_animal),
    n_regions_detected_cell_count = sum(!is.na(Cell_Count)),
    n_regions_detected_intensity = sum(!is.na(Intensity)),
    n_regions_detected_any = n_distinct(Annotation[!is.na(Cell_Count) | !is.na(Intensity)]),
    missing_cell_count = sum(is.na(Cell_Count)),
    missing_intensity = sum(is.na(Intensity)),
    .groups = "drop"
  )

QC_animal_level_coverage_signal <- long %>%
  group_by(SampleID, Animal, Condition) %>%
  summarise(
    n_regions_detected = n_distinct(RegionKey[!is.na(Cell_Count) | !is.na(Intensity)]),
    n_regions_total = n_distinct(RegionKey),
    n_regions_missing = n_regions_total - n_regions_detected,
    percent_regions_missing = 100 * n_regions_missing / n_regions_total,
    total_Cell_Count = sum(Cell_Count, na.rm = TRUE),
    total_Intensity = sum(Intensity, na.rm = TRUE),
    median_Cell_Count = median(Cell_Count, na.rm = TRUE),
    median_Intensity = median(Intensity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    median_Cell_Count = if_else(is.finite(median_Cell_Count), median_Cell_Count, NA_real_),
    median_Intensity = if_else(is.finite(median_Intensity), median_Intensity, NA_real_)
  )

region_condition_missingness_qc <- long %>%
  group_by(Class, Annotation, Abbreviation, RegionLabel, RegionKey, Condition) %>%
  summarise(
    Level = {
      level_values <- Level[!is.na(Level)]
      if (length(level_values) == 0) NA_integer_ else level_values[[1]]
    },
    n_samples_total = n_distinct(SampleID),
    n_cell_count = n_distinct(SampleID[!is.na(Cell_Count)]),
    n_intensity = n_distinct(SampleID[!is.na(Intensity)]),
    missing_cell_count = n_samples_total - n_cell_count,
    missing_intensity = n_samples_total - n_intensity,
    pct_missing_cell_count = 100 * missing_cell_count / n_samples_total,
    pct_missing_intensity = 100 * missing_intensity / n_samples_total,
    .groups = "drop"
  )

QC_region_level_availability <- region_condition_missingness_qc %>%
  select(RegionKey, Annotation, Abbreviation, Class, Level, Condition, n_samples_total, n_cell_count, n_intensity) %>%
  mutate(n_present_any = pmax(n_cell_count, n_intensity)) %>%
  select(RegionKey, Annotation, Abbreviation, Class, Level, Condition, n_samples_total, n_present_any) %>%
  pivot_wider(
    names_from = Condition,
    values_from = c(n_samples_total, n_present_any),
    values_fill = 0
  ) %>%
  mutate(
    n_VEH_paired = n_present_any_VEH_paired,
    n_VEH_unpaired = n_present_any_VEH_unpaired,
    n_CNO_paired = n_present_any_CNO_paired,
    n_CNO_unpaired = n_present_any_CNO_unpaired,
    missing_VEH_paired = n_samples_total_VEH_paired - n_VEH_paired,
    missing_VEH_unpaired = n_samples_total_VEH_unpaired - n_VEH_unpaired,
    missing_CNO_paired = n_samples_total_CNO_paired - n_CNO_paired,
    missing_CNO_unpaired = n_samples_total_CNO_unpaired - n_CNO_unpaired,
    total_n_present = n_VEH_paired + n_VEH_unpaired + n_CNO_paired + n_CNO_unpaired,
    min_n_per_group = pmin(n_VEH_paired, n_VEH_unpaired, n_CNO_paired, n_CNO_unpaired),
    n_groups_present = (n_VEH_paired > 0) + (n_VEH_unpaired > 0) + (n_CNO_paired > 0) + (n_CNO_unpaired > 0),
    passes_learning_effect_filter = n_VEH_paired >= 2 & n_VEH_unpaired >= 2,
    passes_CNO_learning_filter = n_CNO_paired >= 2 & n_VEH_paired >= 2,
    passes_CNO_stress_filter = n_CNO_unpaired >= 2 & n_VEH_unpaired >= 2,
    passes_interaction_filter = n_VEH_paired >= 2 & n_VEH_unpaired >= 2 & n_CNO_paired >= 2 & n_CNO_unpaired >= 2,
    passes_main_display_filter = passes_learning_effect_filter | passes_CNO_learning_filter | passes_CNO_stress_filter | passes_interaction_filter,
    availability_note =
      "This table is any-metric region-level contrast/display availability using pmax(n_cell_count, n_intensity) and n >= 2. It does not determine correlation heatmap or network inclusion."
  ) %>%
  select(
    RegionKey, Annotation, Abbreviation, Class, Level,
    n_VEH_paired, n_VEH_unpaired, n_CNO_paired, n_CNO_unpaired,
    missing_VEH_paired, missing_VEH_unpaired, missing_CNO_paired, missing_CNO_unpaired,
    total_n_present, min_n_per_group, n_groups_present,
    passes_learning_effect_filter, passes_CNO_learning_filter, passes_CNO_stress_filter,
    passes_interaction_filter, passes_main_display_filter, availability_note
  )

QC_region_level_availability_by_metric <- region_condition_missingness_qc %>%
  select(
    RegionKey, Annotation, Abbreviation, Class, Level, Condition,
    n_samples_total, n_cell_count, n_intensity,
    missing_cell_count, missing_intensity
  ) %>%
  pivot_wider(
    names_from = Condition,
    values_from = c(n_samples_total, n_cell_count, n_intensity, missing_cell_count, missing_intensity),
    values_fill = 0
  ) %>%
  mutate(
    total_n_cell_count =
      n_cell_count_VEH_paired + n_cell_count_VEH_unpaired +
      n_cell_count_CNO_paired + n_cell_count_CNO_unpaired,
    total_n_intensity =
      n_intensity_VEH_paired + n_intensity_VEH_unpaired +
      n_intensity_CNO_paired + n_intensity_CNO_unpaired,
    min_n_cell_count_per_group = pmin(
      n_cell_count_VEH_paired, n_cell_count_VEH_unpaired,
      n_cell_count_CNO_paired, n_cell_count_CNO_unpaired
    ),
    min_n_intensity_per_group = pmin(
      n_intensity_VEH_paired, n_intensity_VEH_unpaired,
      n_intensity_CNO_paired, n_intensity_CNO_unpaired
    ),
    availability_note =
      "Metric-specific companion to QC_region_level_availability; the main availability table uses pmax(n_cell_count, n_intensity) for n_present_any."
  )

safe_write_csv <- function(x, path) {
  tryCatch(
    readr::write_csv(x, path),
    error = function(e) warning("Could not write CSV: ", path, " | ", conditionMessage(e), call. = FALSE)
  )
}

safe_write_tsv <- function(x, path) {
  tryCatch(
    readr::write_tsv(x, path),
    error = function(e) warning("Could not write TSV: ", path, " | ", conditionMessage(e), call. = FALSE)
  )
}

safe_write_xlsx <- function(x, path) {
  tryCatch(
    openxlsx::write.xlsx(x, path, overwrite = TRUE),
    error = function(e) warning("Could not write XLSX: ", path, " | ", conditionMessage(e), call. = FALSE)
  )
}

minimum_n_filter_qc <- region_condition_missingness_qc %>%
  transmute(
    Class, Annotation, Abbreviation, Level, RegionLabel, RegionKey, Condition,
    n_cell_count, n_intensity,
    Cell_Count_pass_min_n = n_cell_count >= 2,
    Intensity_pass_min_n = n_intensity >= 2
  )

openxlsx::write.xlsx(
  list(
    animal_counts = animal_counts,
    animal_level_qc = animal_level_qc,
    animal_coverage_signal = QC_animal_level_coverage_signal,
    region_condition_missingness = region_condition_missingness_qc,
    region_level_availability = QC_region_level_availability,
    region_avail_by_metric = QC_region_level_availability_by_metric,
    minimum_n_filter = minimum_n_filter_qc
  ),
  file.path(out_dir, "tables", "animal_region_qc_tables.xlsx"),
  overwrite = TRUE
)
openxlsx::write.xlsx(
  list(
    animal_counts = animal_counts,
    animal_level_qc = animal_level_qc,
    animal_coverage_signal = QC_animal_level_coverage_signal,
    region_condition_missingness = region_condition_missingness_qc,
    region_level_availability = QC_region_level_availability,
    region_avail_by_metric = QC_region_level_availability_by_metric,
    minimum_n_filter = minimum_n_filter_qc
  ),
  file.path(legacy_tab_dir, "animal_region_qc_tables.xlsx"),
  overwrite = TRUE
)
readr::write_csv(QC_animal_level_coverage_signal, file.path(out_dir, "tables", "QC_animal_level_coverage_signal.csv"))
safe_write_csv(QC_region_level_availability, file.path(out_dir, "tables", "QC_region_level_availability.csv"))
safe_write_tsv(QC_region_level_availability, file.path(out_dir, "tables", "QC_region_level_availability.tsv"))
safe_write_xlsx(QC_region_level_availability, file.path(out_dir, "tables", "QC_region_level_availability.xlsx"))
safe_write_csv(QC_region_level_availability_by_metric, file.path(out_dir, "tables", "QC_region_level_availability_by_metric.csv"))
safe_write_tsv(QC_region_level_availability_by_metric, file.path(out_dir, "tables", "QC_region_level_availability_by_metric.tsv"))
safe_write_xlsx(QC_region_level_availability_by_metric, file.path(out_dir, "tables", "QC_region_level_availability_by_metric.xlsx"))
readr::write_csv(QC_animal_level_coverage_signal, file.path(legacy_tab_dir, "QC_animal_level_coverage_signal.csv"))
safe_write_csv(QC_region_level_availability, file.path(legacy_tab_dir, "QC_region_level_availability.csv"))
safe_write_tsv(QC_region_level_availability, file.path(legacy_tab_dir, "QC_region_level_availability.tsv"))
safe_write_xlsx(QC_region_level_availability, file.path(legacy_tab_dir, "QC_region_level_availability.xlsx"))
safe_write_csv(QC_region_level_availability_by_metric, file.path(legacy_tab_dir, "QC_region_level_availability_by_metric.csv"))
safe_write_tsv(QC_region_level_availability_by_metric, file.path(legacy_tab_dir, "QC_region_level_availability_by_metric.tsv"))
safe_write_xlsx(QC_region_level_availability_by_metric, file.path(legacy_tab_dir, "QC_region_level_availability_by_metric.xlsx"))
readr::write_csv(QC_animal_level_coverage_signal, file.path(qc_missingness_dir, "QC_animal_level_coverage_signal.csv"))
safe_write_csv(QC_region_level_availability, file.path(qc_missingness_dir, "QC_region_level_availability.csv"))
safe_write_tsv(QC_region_level_availability, file.path(qc_missingness_dir, "QC_region_level_availability.tsv"))
safe_write_xlsx(QC_region_level_availability, file.path(qc_missingness_dir, "QC_region_level_availability.xlsx"))
safe_write_csv(QC_region_level_availability_by_metric, file.path(qc_missingness_dir, "QC_region_level_availability_by_metric.csv"))
safe_write_tsv(QC_region_level_availability_by_metric, file.path(qc_missingness_dir, "QC_region_level_availability_by_metric.tsv"))
safe_write_xlsx(QC_region_level_availability_by_metric, file.path(qc_missingness_dir, "QC_region_level_availability_by_metric.xlsx"))
readr::write_csv(region_condition_missingness_qc, file.path(qc_missingness_dir, "QC_region_condition_missingness_by_metric.csv"))

qc_condition_plot_labels <- c(
  VEH_paired = "VEH-L",
  VEH_unpaired = "VEH-S",
  CNO_paired = "CNO-L",
  CNO_unpaired = "CNO-S"
)

qc_condition_colors <- c(
  VEH_paired = "#0072B2",
  VEH_unpaired = "#6B6B6B",
  CNO_paired = "#D55E00",
  CNO_unpaired = "#CC79A7"
)

save_qc_figure <- function(plot, filename_base, width, height, dpi = 600, subdir = "qc") {
  target_dir <- file.path(out_dir, "figures", subdir)
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(legacy_fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(legacy_fig_dir, subdir), recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(target_dir, paste0(filename_base, ".pdf")), plot, width = width, height = height)
  ggsave(file.path(target_dir, paste0(filename_base, ".svg")), plot, width = width, height = height, device = svglite::svglite)
  ggsave(file.path(target_dir, paste0(filename_base, ".png")), plot, width = width, height = height, dpi = dpi, bg = "white")
  ggsave(file.path(legacy_fig_dir, paste0(filename_base, ".pdf")), plot, width = width, height = height)
  ggsave(file.path(legacy_fig_dir, paste0(filename_base, ".svg")), plot, width = width, height = height, device = svglite::svglite)
  ggsave(file.path(legacy_fig_dir, paste0(filename_base, ".png")), plot, width = width, height = height, dpi = dpi, bg = "white")
  ggsave(file.path(legacy_fig_dir, subdir, paste0(filename_base, ".pdf")), plot, width = width, height = height)
  ggsave(file.path(legacy_fig_dir, subdir, paste0(filename_base, ".svg")), plot, width = width, height = height, device = svglite::svglite)
  ggsave(file.path(legacy_fig_dir, subdir, paste0(filename_base, ".png")), plot, width = width, height = height, dpi = dpi, bg = "white")
}

qc_plot_animal_metric <- function(data, y_col, title, y_label, filename_base) {
  y_sym <- sym(y_col)
  p <- data %>%
    mutate(Condition = factor(as.character(Condition), levels = names(qc_condition_plot_labels))) %>%
    ggplot(aes(x = Condition, y = !!y_sym, fill = Condition)) +
    geom_point(shape = 21, size = 2.2, colour = "grey20", stroke = 0.2,
               position = position_jitter(width = 0.08, height = 0, seed = 1), alpha = 0.85) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 2.5, fill = "white", colour = "grey20") +
    scale_x_discrete(labels = qc_condition_plot_labels, drop = FALSE) +
    scale_fill_manual(values = qc_condition_colors, guide = "none", drop = FALSE) +
    theme_classic(base_size = 8) +
    labs(title = title, x = NULL, y = y_label)

  save_qc_figure(p, filename_base, width = 4.5, height = 3.5)
}

qc_plot_animal_metric(QC_animal_level_coverage_signal, "n_regions_detected", "Detected regions per animal", "Regions detected", "QC_detected_regions_per_animal")
qc_plot_animal_metric(QC_animal_level_coverage_signal, "total_Cell_Count", "Total cFos+ cell count per animal", "Total Cell_Count", "QC_total_cell_count_per_animal")
qc_plot_animal_metric(QC_animal_level_coverage_signal, "total_Intensity", "Total projection intensity per animal", "Total Intensity", "QC_total_intensity_per_animal")

animal_region_missingness <- long %>%
  mutate(
    missing_any = is.na(Cell_Count) & is.na(Intensity),
    SampleLabel = paste(SampleID, as.character(Condition), sep = " | ")
  ) %>%
  select(SampleID, SampleLabel, Condition, RegionKey, missing_any) %>%
  distinct() %>%
  left_join(long %>% distinct(RegionKey, Annotation), by = "RegionKey") %>%
  mutate(
    SampleLabel = factor(SampleLabel, levels = unique(SampleLabel[order(Condition, SampleID)])),
    RegionShort = str_trunc(coalesce(Annotation, RegionKey), 18)
  )

qc_missing_plot <- ggplot(animal_region_missingness, aes(x = RegionShort, y = SampleLabel, fill = missing_any)) +
  geom_tile(colour = NA) +
  scale_fill_manual(values = c(`FALSE` = "grey15", `TRUE` = "grey88"), name = "Missing") +
  theme_classic(base_size = 6) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 4.8),
    legend.position = "bottom"
  ) +
  labs(title = "Animal x region missingness", x = "Regions", y = NULL)
save_qc_figure(qc_missing_plot, "QC_animal_region_missingness_heatmap", width = 6.8, height = 5)

# -----------------------------
# 4. QC tables
# -----------------------------
qc_sample_summary <- long %>%
  group_by(SampleID, Animal, Sample, SourceFile, Group, Condition) %>%
  summarise(
    n_regions = n_distinct(Annotation),
    n_classes = n_distinct(Class),
    total_cell_count = sum(Cell_Count, na.rm = TRUE),
    total_intensity = sum(Intensity, na.rm = TRUE),
    missing_cell_count = sum(is.na(Cell_Count)),
    missing_intensity = sum(is.na(Intensity)),
    .groups = "drop"
  )

qc_region_availability <- region_condition_missingness_qc %>%
  transmute(
    Class, Annotation, Abbreviation, RegionLabel, RegionKey, Condition,
    n_samples_present = pmax(n_cell_count, n_intensity),
    n_cell_count_present = n_cell_count,
    n_intensity_present = n_intensity,
    availability_note =
      "Availability requires explicit non-missing Cell_Count or Intensity values; row presence alone is not counted."
  ) %>%
  pivot_wider(
    names_from = Condition,
    values_from = c(n_samples_present, n_cell_count_present, n_intensity_present),
    values_fill = 0
  ) %>%
  rename_with(~ str_remove(.x, "^n_samples_present_"), starts_with("n_samples_present_")) %>%
  mutate(
    n_conditions_present = rowSums(across(all_of(levels(long$Condition)), ~ .x > 0)),
    complete_all_conditions = n_conditions_present == length(levels(long$Condition))
  )

qc_class_availability <- long %>%
  group_by(Class, Condition) %>%
  summarise(
    n_regions = n_distinct(Annotation),
    n_samples = n_distinct(SampleID),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Condition,
    values_from = c(n_regions, n_samples),
    values_fill = 0
  )

openxlsx::write.xlsx(
  list(
    animal_counts = animal_counts,
    sample_summary = qc_sample_summary,
    region_availability = qc_region_availability,
    class_availability = qc_class_availability
  ),
  file.path(out_dir, "tables", "QC_region_tracing_tables.xlsx"),
  overwrite = TRUE
)

safe_pheatmap <- function(mat, ..., cluster_rows = TRUE, cluster_cols = TRUE) {
  mat <- as.matrix(mat)

  if (nrow(mat) < 2) cluster_rows <- FALSE
  if (ncol(mat) < 2) cluster_cols <- FALSE

  tryCatch(
    pheatmap::pheatmap(
      mat,
      ...,
      cluster_rows = cluster_rows,
      cluster_cols = cluster_cols
    ),
    error = function(e) {
      msg <- conditionMessage(e)
      if (grepl("non-conformable arrays|must have n >= 2 objects", msg, ignore.case = TRUE)) {
        warning(
          "pheatmap clustering failed due to matrix shape (",
          msg,
          "). Retrying with clustering disabled.",
          call. = FALSE
        )
        return(tryCatch(
          pheatmap::pheatmap(
            mat,
            ...,
            cluster_rows = FALSE,
            cluster_cols = FALSE
          ),
          error = function(e2) {
            warning(
              "pheatmap failed even after disabling clustering (",
              conditionMessage(e2),
              "). Skipping this heatmap.",
              call. = FALSE
            )
            NULL
          }
        ))
      }
      warning("pheatmap failed (", msg, "). Skipping this heatmap.", call. = FALSE)
      NULL
    }
  )
}

avail_mat <- qc_region_availability %>%
  select(RegionKey, all_of(levels(long$Condition))) %>%
  column_to_rownames("RegionKey") %>%
  as.matrix()

pdf(file.path(out_dir, "figures", "QC_region_availability_by_condition.pdf"), width = 8, height = 12)
safe_pheatmap(
  avail_mat,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  main = "Region availability by condition\nvalue = number of animals present",
  fontsize_row = 5
)
dev.off()

# -----------------------------
# 5. Helper functions
# -----------------------------
zscore_rows <- function(mat) {
  out <- t(scale(t(mat)))
  out[!is.finite(out)] <- 0
  out
}

zscore_cols <- function(mat) {
  out <- scale(mat)
  out[!is.finite(out)] <- 0
  out
}

safe_log1p <- function(x) log1p(pmax(x, 0))

clean_heatmap_matrix <- function(mat, fill = 0) {
  mat <- as.matrix(mat)
  mat[!is.finite(mat)] <- fill
  mat
}

condition_colors <- c(
  VEH_paired = "#0072B2",
  VEH_unpaired = "#6B6B6B",
  CNO_paired = "#D55E00",
  CNO_unpaired = "#CC79A7"
)

condition_display_labels <- c(
  VEH_paired = "VEH paired learning",
  VEH_unpaired = "VEH unpaired stress",
  CNO_paired = "CNO during learning",
  CNO_unpaired = "CNO during stress"
)

condition_short_labels <- c(
  VEH_paired = "VEH-L",
  VEH_unpaired = "VEH-S",
  CNO_paired = "CNO-L",
  CNO_unpaired = "CNO-S"
)

metric_display_labels <- c(
  Cell_Count = "cFos+ cell count",
  Intensity = "Projection intensity"
)

contrast_display_labels <- c(
  Learning_effect = "Paired learning effect",
  CeM_manipulation_during_learning = "CNO effect during learning",
  CeM_manipulation_during_stress = "CNO effect during stress",
  Learning_x_CeM_interaction = "Learning x CNO interaction",
  Paired_vs_unpaired = "Paired vs unpaired average",
  CNO_vs_VEH = "CNO vs VEH average"
)

contrast_short_labels <- c(
  Learning_effect = "Learning",
  CeM_manipulation_during_learning = "CNO-learning",
  CeM_manipulation_during_stress = "CNO-stress",
  Learning_x_CeM_interaction = "Interaction",
  Paired_vs_unpaired = "Paired avg.",
  CNO_vs_VEH = "CNO avg."
)

contrast_plain_language <- c(
  Learning_effect = "Positive values mean VEH paired learning is higher than VEH unpaired stress.",
  CeM_manipulation_during_learning = "Positive values mean CNO during learning is higher than VEH during learning.",
  CeM_manipulation_during_stress = "Positive values mean CNO during stress is higher than VEH during stress.",
  Learning_x_CeM_interaction = "Interaction asks whether the CNO-VEH difference changes between learning and stress: positive means the CNO effect is stronger, or more positive, during paired learning than during unpaired stress.",
  Paired_vs_unpaired = "Positive values mean paired groups are higher than unpaired groups on average.",
  CNO_vs_VEH = "Positive values mean CNO groups are higher than VEH groups on average."
)

display_metric <- function(metric) {
  dplyr::coalesce(unname(metric_display_labels[metric]), metric)
}

display_contrast <- function(contrast, short = FALSE) {
  labels <- if (isTRUE(short)) contrast_short_labels else contrast_display_labels
  dplyr::coalesce(unname(labels[contrast]), contrast)
}

effect_colors <- c(
  negative = "#2166AC",
  neutral = "#F7F7F7",
  positive = "#B2182B"
)

edge_class_colors <- c(
  retained = "#8C8C8C",
  gained = "#009E73",
  lost = "#D55E00",
  `sign switch` = "#0072B2"
)

fig4_system_group <- function(class, annotation, abbreviation, region_label) {
  text <- str_to_lower(str_c(class, annotation, abbreviation, region_label, sep = " "))
  case_when(
    str_detect(text, "amygdal|\\bcea\\b|\\bcea[clm]?\\b|\\bbla|\\bbma|\\bla\\b|\\bia\\b|\\bpa\\b|\\bpaa\\b|\\bcoa\\b|\\baaa\\b") ~ "Amygdala",
    str_detect(text, "hypothalam|preoptic|supraoptic|subfornical|subthalamic|retrochiasmatic|median eminence|\\bpv|\\bavp\\b|\\bdmh\\b|\\bme(p|po)?\\b|\\bstn\\b|\\bso\\b|\\bsfo\\b|\\brch\\b|\\bvmh|\\bvmpo\\b") ~ "Hypothalamus",
    str_detect(text, "thalam|habenula|geniculate|epithalam") ~ "Thalamus",
    str_detect(text, "cortical plate|cortex|cingulate|insular|retrosplenial|somatosensory|gustatory|ectorhinal|piriform|hippocamp|retrohippocamp|visual|auditory|orbital|prelimbic|infralimbic") ~ "Cortex / hippocampus",
    str_detect(text, "pallid|globus|bed nuclei|stria terminalis|triangular nucleus") ~ "Pallidum / BST",
    str_detect(text, "striat|accumbens|caudoputamen|septal|septofimbrial") ~ "Striatum / septum",
    str_detect(text, "midbrain|medulla|pons|perifornical") ~ "Brainstem",
    TRUE ~ coalesce(class, "Other")
  )
}

fig4_system_levels <- c(
  "Amygdala",
  "Hypothalamus",
  "Thalamus",
  "Cortex / hippocampus",
  "Striatum / septum",
  "Pallidum / BST",
  "Brainstem",
  "Other"
)

theme_nature <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(colour = "black"),
      plot.title = element_text(face = "bold", size = base_size + 0.8, hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.2, colour = "grey25", hjust = 0),
      axis.title = element_text(colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.line = element_line(linewidth = 0.25, colour = "black"),
      axis.ticks = element_line(linewidth = 0.25, colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size - 1),
      legend.title = element_text(size = base_size - 1),
      legend.text = element_text(size = base_size - 1),
      legend.key.size = grid::unit(3.4, "mm"),
      plot.margin = margin(4, 5, 4, 5)
    )
}

legacy_figure_subdir <- function(filename_base) {
  case_when(
    str_detect(filename_base, "^volcano_") ~ "volcano_fdr",
    str_detect(filename_base, "^top_region_profiles_") ~ "regional_profiles",
    str_detect(filename_base, "^class_level_profiles_") ~ "class_level_profiles",
    str_detect(filename_base, "^PCA_|^UMAP_") ~ "dimensionality_reduction",
    str_detect(filename_base, "^CeM_centered") ~ "covariance",
    str_detect(filename_base, "heatmap|effect_size") ~ "effect_size_maps",
    TRUE ~ "misc"
  )
}

figure_target_dir <- function(subdir) {
  case_when(
    str_starts(subdir, "volcano_fdr") ~ file.path(exploratory_profile_fig_dir, subdir),
    str_starts(subdir, "regional_profiles") ~ file.path(exploratory_profile_fig_dir, str_remove(subdir, "^regional_profiles/?")),
    str_starts(subdir, "class_level_profiles") ~ file.path(exploratory_profile_fig_dir, subdir),
    str_starts(subdir, "effect_size_maps") ~ file.path(exploratory_effect_fig_dir, str_remove(subdir, "^effect_size_maps/?")),
    str_starts(subdir, "covariance") ~ file.path(exploratory_covariance_dir, "figures", str_remove(subdir, "^covariance/?")),
    str_starts(subdir, "dimensionality_reduction") ~ file.path(exploratory_dimred_dir, "figures", str_remove(subdir, "^dimensionality_reduction/?")),
    str_starts(subdir, "qc") ~ file.path(results_root, "03_qc", "figures", str_remove(subdir, "^qc/?")),
    TRUE ~ file.path(out_dir, "figures", subdir)
  ) %>%
    str_replace("/$", "")
}

save_figure <- function(plot, filename_base, width, height, dpi = 600, subdir = NULL) {
  if (is.null(plot)) return(invisible(NULL))
  subdir <- subdir %||% legacy_figure_subdir(filename_base)
  target_dir <- figure_target_dir(subdir)
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "figures", subdir), recursive = TRUE, showWarnings = FALSE)
  dir.create(legacy_fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(legacy_fig_dir, subdir), recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(target_dir, paste0(filename_base, ".pdf")), plot, width = width, height = height)
  ggsave(file.path(target_dir, paste0(filename_base, ".svg")), plot, width = width, height = height, device = svglite::svglite)
  ggsave(file.path(target_dir, paste0(filename_base, ".png")), plot, width = width, height = height, dpi = dpi, bg = "white")
  ggsave(file.path(out_dir, "figures", subdir, paste0(filename_base, ".pdf")), plot, width = width, height = height)
  ggsave(file.path(out_dir, "figures", subdir, paste0(filename_base, ".svg")), plot, width = width, height = height, device = svglite::svglite)
  ggsave(file.path(out_dir, "figures", subdir, paste0(filename_base, ".png")), plot, width = width, height = height, dpi = dpi, bg = "white")
  ggsave(file.path(legacy_fig_dir, paste0(filename_base, ".pdf")), plot, width = width, height = height)
  ggsave(file.path(legacy_fig_dir, paste0(filename_base, ".svg")), plot, width = width, height = height, device = svglite::svglite)
  ggsave(file.path(legacy_fig_dir, paste0(filename_base, ".png")), plot, width = width, height = height, dpi = dpi, bg = "white")
  ggsave(file.path(legacy_fig_dir, subdir, paste0(filename_base, ".pdf")), plot, width = width, height = height)
  ggsave(file.path(legacy_fig_dir, subdir, paste0(filename_base, ".svg")), plot, width = width, height = height, device = svglite::svglite)
  ggsave(file.path(legacy_fig_dir, subdir, paste0(filename_base, ".png")), plot, width = width, height = height, dpi = dpi, bg = "white")
  invisible(plot)
}

make_group_summary_matrix <- function(data, metric) {
  metric_sym <- sym(metric)

  data %>%
    group_by(Class, Annotation, Abbreviation, RegionLabel, RegionKey, Condition) %>%
    summarise(
      mean_value = mean(!!metric_sym, na.rm = TRUE),
      median_value = median(!!metric_sym, na.rm = TRUE),
      sd_value = sd(!!metric_sym, na.rm = TRUE),
      n = sum(!is.na(!!metric_sym)),
      .groups = "drop"
    ) %>%
    mutate(mean_log1p = safe_log1p(mean_value)) %>%
    select(Class, RegionKey, RegionLabel, Condition, mean_log1p) %>%
    pivot_wider(names_from = Condition, values_from = mean_log1p) %>%
    arrange(Class, RegionLabel)
}

plot_group_heatmap <- function(summary_df, metric, class_filter = NULL) {
  df <- summary_df
  title_suffix <- "all_classes"

  if (!is.null(class_filter)) {
    df <- df %>% filter(Class == class_filter)
    title_suffix <- str_replace_all(class_filter, "[^A-Za-z0-9]+", "_")
  }

  if (nrow(df) < 2) return(NULL)

  mat <- df %>%
    select(RegionKey, all_of(levels(long$Condition))) %>%
    column_to_rownames("RegionKey") %>%
    as.matrix()

  mat_z <- zscore_rows(mat)

  finite_values <- mat_z[is.finite(mat_z)]
  if (length(unique(finite_values)) < 2) return(NULL)

  pdf(file.path(out_dir, "figures", paste0("heatmap_group_", metric, "_", title_suffix, ".pdf")),
      width = 7, height = max(5, min(18, nrow(mat_z) * 0.12)))
  safe_pheatmap(
    mat_z,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
    main = paste0(metric, " group mean log1p, row z-score: ", title_suffix),
    fontsize_row = ifelse(nrow(mat_z) > 80, 4, 6),
    border_color = NA
  )
  dev.off()
}

make_sample_matrix <- function(data, metric, class_filter = NULL, complete_region_only = FALSE) {
  metric_sym <- sym(metric)

  df <- data
  if (!is.null(class_filter)) df <- df %>% filter(Class == class_filter)

  if (complete_region_only) {
    keep_regions <- df %>%
      group_by(RegionKey) %>%
      summarise(n_conditions = n_distinct(Condition), .groups = "drop") %>%
      filter(n_conditions == length(levels(long$Condition))) %>%
      pull(RegionKey)
    df <- df %>% filter(RegionKey %in% keep_regions)
  }

  mat <- df %>%
    select(SampleID, Condition, RegionKey, value = !!metric_sym) %>%
    group_by(SampleID, Condition, RegionKey) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
    mutate(value = if_else(is.nan(value), NA_real_, value)) %>%
    pivot_wider(names_from = RegionKey, values_from = value) %>%
    arrange(Condition, SampleID)

  sample_anno <- mat %>% select(SampleID, Condition)
  mat_numeric <- mat %>% select(-SampleID, -Condition) %>% as.data.frame()
  rownames(mat_numeric) <- mat$SampleID

  list(mat = as.matrix(mat_numeric), annotation = sample_anno)
}

# Helper for supplementary condition-specific covariance displays below. This is
# not used by the active pooled heatmap/network workflow.
prune_regions_by_pairwise_overlap <- function(data, metric, region_keys,
                                              min_pairwise_n = 3,
                                              conditions = nature_condition_levels) {
  metric_sym <- rlang::sym(metric)

  current_regions <- region_keys
  dropped <- tibble()

  repeat {
    if (length(current_regions) < 2) break

    pair_qc <- purrr::map_dfr(conditions, function(cond) {
      wide <- data %>%
        filter(Condition == cond, RegionKey %in% current_regions) %>%
        select(SampleID, RegionKey, RawValue = !!metric_sym) %>%
        mutate(Value = safe_log1p(RawValue)) %>%
        group_by(SampleID, RegionKey) %>%
        summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
        mutate(Value = if_else(is.nan(Value), NA_real_, Value)) %>%
        pivot_wider(names_from = RegionKey, values_from = Value)

      mat <- wide %>%
        select(any_of(current_regions)) %>%
        as.matrix()

      missing_cols <- setdiff(current_regions, colnames(mat))
      if (length(missing_cols) > 0) {
        mat <- cbind(
          mat,
          matrix(
            NA_real_,
            nrow = nrow(mat),
            ncol = length(missing_cols),
            dimnames = list(NULL, missing_cols)
          )
        )
      }

      mat <- mat[, current_regions, drop = FALSE]

      expand.grid(
        Region1 = current_regions,
        Region2 = current_regions,
        stringsAsFactors = FALSE
      ) %>%
        as_tibble() %>%
        filter(Region1 < Region2) %>%
        rowwise() %>%
        mutate(
          Condition = cond,
          n_pair = sum(is.finite(mat[, Region1]) & is.finite(mat[, Region2])),
          var1_ok = sd(mat[, Region1], na.rm = TRUE) > 0,
          var2_ok = sd(mat[, Region2], na.rm = TRUE) > 0,
          pair_ok = n_pair >= min_pairwise_n & var1_ok & var2_ok
        ) %>%
        ungroup()
    })

    failing_pairs <- pair_qc %>% filter(!pair_ok)

    if (nrow(failing_pairs) == 0) {
      return(list(
        kept_regions = current_regions,
        dropped_regions = dropped,
        pair_qc = pair_qc
      ))
    }

    region_fail_counts <- failing_pairs %>%
      select(Condition, Region1, Region2, n_pair, var1_ok, var2_ok) %>%
      pivot_longer(
        cols = c(Region1, Region2),
        names_to = "side",
        values_to = "RegionKey"
      ) %>%
      group_by(RegionKey) %>%
      summarise(
        n_failing_pairs = n(),
        min_n_pair = min(n_pair, na.rm = TRUE),
        n_zero_variance_flags = sum(!var1_ok | !var2_ok, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(n_failing_pairs), min_n_pair)

    drop_region <- region_fail_counts$RegionKey[[1]]

    dropped <- bind_rows(
      dropped,
      region_fail_counts %>%
        filter(RegionKey == drop_region) %>%
        mutate(drop_step = nrow(dropped) + 1)
    )

    current_regions <- setdiff(current_regions, drop_region)
  }

  list(
    kept_regions = current_regions,
    dropped_regions = dropped,
    pair_qc = tibble()
  )
}

# -----------------------------
# 6. Group heatmaps, sample heatmaps, correlations, networks
# -----------------------------
cor_with_p <- function(mat, min_n = 3) {
  regions <- colnames(mat)

  expand.grid(region1 = regions, region2 = regions, stringsAsFactors = FALSE) %>%
    as_tibble() %>%
    rowwise() %>%
    mutate(
      n_pair = sum(complete.cases(mat[, region1], mat[, region2])),
      r = ifelse(
        n_pair >= min_n,
        suppressWarnings(tryCatch(
          cor(mat[, region1], mat[, region2], use = "pairwise.complete.obs", method = "spearman"),
          error = function(e) NA_real_
        )),
        NA_real_
      ),
      p = ifelse(
        n_pair >= min_n && !is.na(r) && region1 != region2,
        suppressWarnings(tryCatch(
          cor.test(mat[, region1], mat[, region2], method = "spearman", exact = FALSE)$p.value,
          error = function(e) NA_real_
        )),
        NA_real_
      )
    ) %>%
    ungroup() %>%
    mutate(fdr = p.adjust(p, method = "BH"))
}

pairwise_complete_n_matrix <- function(mat) {
  observed <- is.finite(mat)
  storage.mode(observed) <- "integer"
  crossprod(observed)
}

spearman_cor_matrix <- function(mat) {
  mat <- as.matrix(mat)
  if (ncol(mat) == 0) {
    return(matrix(numeric(), nrow = 0, ncol = 0))
  }
  if (ncol(mat) == 1) {
    return(matrix(1, nrow = 1, ncol = 1, dimnames = list(colnames(mat), colnames(mat))))
  }
  suppressWarnings(cor(mat, use = "pairwise.complete.obs", method = "spearman"))
}

get_region_metadata <- function(region_keys) {
  long %>%
    filter(RegionKey %in% region_keys) %>%
    distinct(RegionKey, Annotation, Abbreviation, Class, Level) %>%
    group_by(RegionKey) %>%
    summarise(
      Annotation = dplyr::first(Annotation),
      Abbreviation = dplyr::first(Abbreviation),
      Class = dplyr::first(Class),
      Level = {
        level_values <- Level[!is.na(Level)]
        if (length(level_values) == 0) NA_integer_ else level_values[[1]]
      },
      .groups = "drop"
    )
}

class_filter_label <- function(class_filter) {
  if (is.null(class_filter)) "all_classes" else class_filter
}

make_covariance_pairwise_qc <- function(mat_filtered, metric, class_filter, title_suffix,
                                        min_pairwise_n = 3,
                                        abs_r_cutoff = network_abs_r_cutoff,
                                        fdr_cutoff = network_fdr_cutoff) {
  empty <- tibble(
    Metric = character(),
    covariance_mode = character(),
    class_filter = character(),
    title_suffix = character(),
    region1 = character(),
    region2 = character(),
    n_pair = integer(),
    r = double(),
    p = double(),
    fdr = double(),
    passes_pairwise_n = logical(),
    passes_abs_r = logical(),
    passes_fdr = logical(),
    retained_edge = logical()
  )

  if (ncol(mat_filtered) < 2) return(empty)

  cor_with_p(mat_filtered, min_n = min_pairwise_n) %>%
    filter(region1 < region2) %>%
    mutate(
      Metric = metric,
      covariance_mode = covariance_mode,
      class_filter = class_filter_label(class_filter),
      title_suffix = title_suffix,
      passes_pairwise_n = n_pair >= min_pairwise_n,
      passes_abs_r = !is.na(r) & abs(r) >= abs_r_cutoff,
      passes_fdr = !is.na(fdr) & fdr <= fdr_cutoff,
      retained_edge = passes_pairwise_n & passes_abs_r & passes_fdr
    ) %>%
    select(
      Metric, covariance_mode, class_filter, title_suffix,
      region1, region2, n_pair, r, p, fdr,
      passes_pairwise_n, passes_abs_r, passes_fdr, retained_edge
    )
}

write_covariance_network_summary <- function(region_qc, pairwise_qc, metric, class_filter, title_suffix) {
  summary_df <- tibble(
    Metric = metric,
    covariance_mode = covariance_mode,
    class_filter = class_filter_label(class_filter),
    title_suffix = title_suffix,
    min_pairwise_n = min_pairwise_n,
    network_abs_r_cutoff = network_abs_r_cutoff,
    network_fdr_cutoff = network_fdr_cutoff,
    n_regions_raw = nrow(region_qc),
    n_regions_after_nonmissing_variance_filter = sum(region_qc$passes_covariance_region_filter, na.rm = TRUE),
    n_possible_pairs = if (nrow(pairwise_qc) == 0) 0L else nrow(pairwise_qc),
    n_pairs_with_pairwise_n_ge_min = sum(pairwise_qc$passes_pairwise_n, na.rm = TRUE),
    n_pairs_with_abs_r_ge_cutoff = sum(pairwise_qc$passes_abs_r, na.rm = TRUE),
    n_pairs_with_fdr_le_cutoff = sum(pairwise_qc$passes_fdr, na.rm = TRUE),
    n_edges_retained = sum(pairwise_qc$retained_edge, na.rm = TRUE)
  )

  readr::write_csv(
    summary_df,
    file.path(qc_covariance_dir, paste0("Covariance_network_summary_", metric, "_", title_suffix, ".csv"))
  )
  readr::write_csv(
    summary_df,
    file.path(exploratory_network_dir, paste0("Network_QC_summary_", metric, "_", title_suffix, ".csv"))
  )

  invisible(summary_df)
}

write_correlation_heatmap_qc <- function(mat_before, mat_after, cor_mat, n_pair_mat,
                                         metric, class_filter, title_suffix,
                                         min_pairwise_n = 3) {
  region_qc <- tibble(
    RegionKey = colnames(mat_before),
    n_nonmissing_overall = colSums(is.finite(mat_before)),
    sd_log1p = apply(mat_before, 2, function(x) stats::sd(x, na.rm = TRUE))
  ) %>%
    mutate(
      Metric = metric,
      covariance_mode = covariance_mode,
      class_filter = class_filter_label(class_filter),
      title_suffix = title_suffix,
      sd_log1p = if_else(is.finite(sd_log1p), sd_log1p, NA_real_),
      passes_n = n_nonmissing_overall >= min_pairwise_n,
      passes_variance = !is.na(sd_log1p) & sd_log1p > 0,
      passes_covariance_region_filter = passes_n & passes_variance,
      passes_heatmap_filter = passes_covariance_region_filter,
      n_regions_before_filter = ncol(mat_before),
      n_regions_after_filter = ncol(mat_after)
    )

  if (ncol(mat_after) > 0) {
    filtered_regions <- colnames(mat_after)
    n_valid_possible <- pmax(length(filtered_regions) - 1L, 0L)

    pairwise_summary <- tibble(
      RegionKey = filtered_regions,
      n_pairwise_regions_with_n_ge_min = rowSums(n_pair_mat >= min_pairwise_n) - 1L,
      prop_pairwise_regions_with_n_ge_min = if (n_valid_possible > 0) {
        (rowSums(n_pair_mat >= min_pairwise_n) - 1L) / n_valid_possible
      } else {
        rep(NA_real_, length(filtered_regions))
      },
      n_estimable_correlations = rowSums(!is.na(cor_mat)) - 1L,
      prop_estimable_correlations = if (n_valid_possible > 0) {
        (rowSums(!is.na(cor_mat)) - 1L) / n_valid_possible
      } else {
        rep(NA_real_, length(filtered_regions))
      },
      min_pairwise_n_observed = purrr::map_int(filtered_regions, function(region) {
        off_diag <- n_pair_mat[region, setdiff(filtered_regions, region)]
        if (length(off_diag) == 0) NA_integer_ else min(off_diag, na.rm = TRUE)
      })
    )

    n_pair_out <- as_tibble(n_pair_mat, rownames = "RegionKey")
    cor_out <- as_tibble(cor_mat, rownames = "RegionKey")
  } else {
    pairwise_summary <- tibble(
      RegionKey = character(),
      n_pairwise_regions_with_n_ge_min = integer(),
      prop_pairwise_regions_with_n_ge_min = double(),
      n_estimable_correlations = integer(),
      prop_estimable_correlations = double(),
      min_pairwise_n_observed = integer()
    )
    n_pair_out <- tibble()
    cor_out <- tibble()
  }

  region_qc <- region_qc %>%
    left_join(get_region_metadata(colnames(mat_before)), by = "RegionKey") %>%
    left_join(pairwise_summary, by = "RegionKey") %>%
    select(
      Metric, covariance_mode, class_filter, title_suffix,
      RegionKey, Annotation, Abbreviation, Class, Level,
      n_nonmissing_overall, sd_log1p, passes_covariance_region_filter,
      n_pairwise_regions_with_n_ge_min, prop_pairwise_regions_with_n_ge_min,
      min_pairwise_n_observed, n_estimable_correlations,
      prop_estimable_correlations, passes_n, passes_variance,
      n_regions_before_filter, n_regions_after_filter, passes_heatmap_filter
    ) %>%
    arrange(desc(passes_covariance_region_filter), RegionKey)

  readr::write_csv(
    region_qc,
    file.path(qc_covariance_dir, paste0("Correlation_heatmap_region_qc_", metric, "_", title_suffix, ".csv"))
  )
  readr::write_csv(
    region_qc,
    file.path(qc_covariance_dir, paste0("Covariance_region_qc_", metric, "_", title_suffix, ".csv"))
  )
  readr::write_csv(
    n_pair_out,
    file.path(qc_covariance_dir, paste0("Correlation_heatmap_pairwise_n_matrix_", metric, "_", title_suffix, ".csv"))
  )
  readr::write_csv(
    cor_out,
    file.path(qc_covariance_dir, paste0("Correlation_heatmap_masked_spearman_matrix_", metric, "_", title_suffix, ".csv"))
  )

  invisible(region_qc)
}

plot_correlation_heatmap <- function(mat, metric, class_filter = NULL) {
  title_suffix <- ifelse(is.null(class_filter), "all_classes", str_replace_all(class_filter, "[^A-Za-z0-9]+", "_"))

  mat_before <- as.matrix(mat)

  keep <- colSums(is.finite(mat_before)) >= min_pairwise_n &
    apply(mat_before, 2, function(x) {
      x_sd <- stats::sd(x, na.rm = TRUE)
      is.finite(x_sd) && x_sd > 0
    })

  mat <- mat_before[, keep, drop = FALSE]

  if (ncol(mat) > 0) {
    n_pair_mat <- pairwise_complete_n_matrix(mat)
    cor_mat <- spearman_cor_matrix(mat)
    cor_mat[n_pair_mat < min_pairwise_n] <- NA_real_
    cor_mat[!is.finite(cor_mat)] <- NA_real_
    diag(cor_mat) <- 1
  } else {
    n_pair_mat <- matrix(integer(), nrow = 0, ncol = 0)
    cor_mat <- matrix(numeric(), nrow = 0, ncol = 0)
  }

  write_correlation_heatmap_qc(
    mat_before = mat_before,
    mat_after = mat,
    cor_mat = cor_mat,
    n_pair_mat = n_pair_mat,
    metric = metric,
    class_filter = class_filter,
    title_suffix = title_suffix,
    min_pairwise_n = min_pairwise_n
  )

  if (ncol(mat) < 3) return(NULL)

  # Keep correlation heatmaps robust: avoid clustering path that can fail
  # on sparse/degenerate covariance structures.
  heatmap_cluster <- FALSE

  pdf(file.path(out_dir, "figures", paste0("correlation_heatmap_", metric, "_", title_suffix, ".pdf")), width = 6.8, height = 6.8)
  safe_pheatmap(
    cor_mat,
    cluster_rows = heatmap_cluster,
    cluster_cols = heatmap_cluster,
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
    breaks = seq(-1, 1, length.out = 101),
    na_col = "grey88",
    main = paste0("Spearman region-region correlation: ", metric, " | ", title_suffix, " | ", covariance_mode, " pooled across all conditions (grey = not estimable)"),
    fontsize_row = ifelse(ncol(cor_mat) > 80, 3, 6),
    fontsize_col = ifelse(ncol(cor_mat) > 80, 3, 6),
    border_color = NA
  )
  dev.off()

  pdf(file.path(qc_covariance_dir, paste0("Correlation_heatmap_pairwise_n_", metric, "_", title_suffix, ".pdf")), width = 6.8, height = 6.8)
  safe_pheatmap(
    n_pair_mat,
    color = colorRampPalette(c("grey95", "#56B4E9", "#0072B2"))(100),
    main = paste0("Pairwise complete n: ", metric, " | ", title_suffix, " | ", covariance_mode, " pooled across all conditions"),
    fontsize_row = ifelse(ncol(n_pair_mat) > 80, 3, 6),
    fontsize_col = ifelse(ncol(n_pair_mat) > 80, 3, 6),
    border_color = NA
  )
  dev.off()
}

plot_network <- function(mat, metric, class_filter = NULL) {
  title_suffix <- ifelse(is.null(class_filter), "all_classes", str_replace_all(class_filter, "[^A-Za-z0-9]+", "_"))

  mat_before <- as.matrix(mat)

  keep <- colSums(is.finite(mat_before)) >= min_pairwise_n &
    apply(mat_before, 2, function(x) {
      x_sd <- stats::sd(x, na.rm = TRUE)
      is.finite(x_sd) && x_sd > 0
    })

  mat <- mat_before[, keep, drop = FALSE]

  if (ncol(mat) > 0) {
    n_pair_mat <- pairwise_complete_n_matrix(mat)
    cor_mat <- spearman_cor_matrix(mat)
    cor_mat[n_pair_mat < min_pairwise_n] <- NA_real_
    cor_mat[!is.finite(cor_mat)] <- NA_real_
    diag(cor_mat) <- 1
  } else {
    n_pair_mat <- matrix(integer(), nrow = 0, ncol = 0)
    cor_mat <- matrix(numeric(), nrow = 0, ncol = 0)
  }

  region_qc <- write_correlation_heatmap_qc(
    mat_before = mat_before,
    mat_after = mat,
    cor_mat = cor_mat,
    n_pair_mat = n_pair_mat,
    metric = metric,
    class_filter = class_filter,
    title_suffix = title_suffix,
    min_pairwise_n = min_pairwise_n
  )

  pairwise_qc <- make_covariance_pairwise_qc(
    mat_filtered = mat,
    metric = metric,
    class_filter = class_filter,
    title_suffix = title_suffix,
    min_pairwise_n = min_pairwise_n,
    abs_r_cutoff = network_abs_r_cutoff,
    fdr_cutoff = network_fdr_cutoff
  )

  readr::write_csv(
    pairwise_qc,
    file.path(qc_covariance_dir, paste0("Network_pairwise_edge_qc_", metric, "_", title_suffix, ".csv"))
  )
  readr::write_csv(
    pairwise_qc,
    file.path(exploratory_network_dir, paste0("Network_pairwise_edge_qc_", metric, "_", title_suffix, ".csv"))
  )

  write_covariance_network_summary(
    region_qc = region_qc,
    pairwise_qc = pairwise_qc,
    metric = metric,
    class_filter = class_filter,
    title_suffix = title_suffix
  )

  cors <- pairwise_qc %>%
    filter(retained_edge) %>%
    mutate(fdr_scope = "BH across all pairwise region tests in this metric/class-filtered exploratory network") %>%
    select(region1, region2, n_pair, r, p, fdr, fdr_scope)

  readr::write_csv(cors, file.path(out_dir, "tables", paste0("network_edges_", metric, "_", title_suffix, ".csv")))
  readr::write_csv(cors, file.path(legacy_tab_dir, paste0("network_edges_", metric, "_", title_suffix, ".csv")))
  readr::write_csv(cors, file.path(exploratory_network_dir, paste0("Network_edges_", metric, "_", title_suffix, ".csv")))

  if (ncol(mat) < 4) return(NULL)
  if (nrow(cors) < 2) return(NULL)

  graph <- igraph::graph_from_data_frame(
    d = cors %>%
      transmute(
        from = region1,
        to = region2,
        r,
        fdr,
        n_pair,
        sign = if_else(r > 0, "positive", "negative")
      ),
    directed = FALSE
  )

  set.seed(1)
  p <- ggraph(graph, layout = "fr") +
    geom_edge_link(aes(width = abs(r), linetype = sign), alpha = 0.45) +
    geom_node_point(size = 3) +
    geom_node_text(
      aes(label = str_remove(name, "^.* :: ") %>% str_extract("^[^|]+")),
      repel = TRUE,
      size = 2.5
    ) +
    scale_edge_width(range = c(0.2, 2.2)) +
    theme_void(base_size = 9) +
    labs(
      title = paste0("Region correlation network: ", metric, " | ", title_suffix, " | ", covariance_mode),
      subtitle = paste0("Pooled across all conditions; Spearman |r| >= ", network_abs_r_cutoff, ", FDR <= ", network_fdr_cutoff)
    )

  ggsave(file.path(out_dir, "figures", paste0("network_", metric, "_", title_suffix, ".pdf")),
         p, width = 6.8, height = 5.8)
  ggsave(file.path(legacy_fig_dir, paste0("network_", metric, "_", title_suffix, ".pdf")),
         p, width = 6.8, height = 5.8)
  ggsave(file.path(exploratory_network_dir, paste0("Network_", metric, "_", title_suffix, ".pdf")),
         p, width = 6.8, height = 5.8)
  ggsave(file.path(exploratory_network_dir, paste0("Network_", metric, "_", title_suffix, ".png")),
         p, width = 6.8, height = 5.8, dpi = 600, bg = "white")
}

run_step_safely <- function(step_label, fn) {
  tryCatch(
    fn(),
    error = function(e) {
      warning(step_label, " failed: ", conditionMessage(e), call. = FALSE)
      NULL
    }
  )
}

all_classes <- sort(unique(long$Class))

for (metric in metrics_to_analyse) {
  message("Analysing metric: ", metric)

  summary_df <- make_group_summary_matrix(long, metric)
  readr::write_csv(summary_df, file.path(out_dir, "tables", paste0("group_summary_", metric, ".csv")))
  readr::write_csv(summary_df, file.path(legacy_tab_dir, paste0("group_summary_", metric, ".csv")))

  run_step_safely(
    paste0("plot_group_heatmap[", metric, ",all_classes]"),
    function() plot_group_heatmap(summary_df, metric, class_filter = NULL)
  )
  purrr::walk(
    all_classes,
    ~ run_step_safely(
      paste0("plot_group_heatmap[", metric, ",", .x, "]"),
      function() plot_group_heatmap(summary_df, metric, class_filter = .x)
    )
  )

  sm <- make_sample_matrix(long, metric, class_filter = NULL)
  mat_log <- safe_log1p(sm$mat)
  keep <- colSums(!is.na(mat_log)) >= min_pairwise_n
  mat_log <- mat_log[, keep, drop = FALSE]

  if (nrow(mat_log) >= 2 && ncol(mat_log) >= 2) {
    row_anno <- sm$annotation %>% as.data.frame()
    rownames(row_anno) <- sm$annotation$SampleID
    row_anno <- row_anno[rownames(mat_log), "Condition", drop = FALSE]
    mat_log_z <- zscore_rows(mat_log)

    finite_values <- mat_log_z[is.finite(mat_log_z)]
    if (length(unique(finite_values)) >= 2) {
      pdf(file.path(out_dir, "figures", paste0("sample_region_heatmap_", metric, "_all_classes.pdf")),
          width = 6.8, height = max(5, nrow(mat_log) * 0.20))
      safe_pheatmap(
        mat_log_z,
        annotation_row = row_anno,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        show_colnames = FALSE,
        color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
        main = paste0("Sample x region heatmap: log1p ", metric, ", row z-score"),
        border_color = NA
      )
      dev.off()
    }
  }

  run_step_safely(
    paste0("plot_correlation_heatmap[", metric, ",all_classes]"),
    function() plot_correlation_heatmap(safe_log1p(sm$mat), metric, class_filter = NULL)
  )
  run_step_safely(
    paste0("plot_network[", metric, ",all_classes]"),
    function() plot_network(safe_log1p(sm$mat), metric, class_filter = NULL)
  )

  for (cl in all_classes) {
    sm_class <- make_sample_matrix(long, metric, class_filter = cl)
    run_step_safely(
      paste0("plot_correlation_heatmap[", metric, ",", cl, "]"),
      function() plot_correlation_heatmap(safe_log1p(sm_class$mat), metric, class_filter = cl)
    )
    run_step_safely(
      paste0("plot_network[", metric, ",", cl, "]"),
      function() plot_network(safe_log1p(sm_class$mat), metric, class_filter = cl)
    )
  }
}

# -----------------------------
# 7. Region-level limma contrasts
# -----------------------------
run_limma_region_contrasts <- function(data, metric) {
  metric_sym <- sym(metric)

  df <- data %>%
    select(SampleID, Condition, Class, Annotation, Abbreviation, Level, RegionLabel, RegionKey, value = !!metric_sym) %>%
    mutate(value = safe_log1p(value)) %>%
    filter(!is.na(Condition))

  sample_meta <- df %>%
    distinct(SampleID, Condition) %>%
    arrange(factor(as.character(Condition), levels = levels(long$Condition)), SampleID) %>%
    mutate(Condition = factor(as.character(Condition), levels = levels(long$Condition)))

  first_present <- function(x) {
    x <- x[!is.na(x) & x != ""]
    if (length(x) == 0) NA else x[[1]]
  }

  region_meta <- df %>%
    group_by(RegionKey) %>%
    summarise(
      Class = first_present(Class),
      Annotation = first_present(Annotation),
      Abbreviation = first_present(Abbreviation),
      Level = {
        level_values <- Level[!is.na(Level)]
        if (length(level_values) == 0) NA_integer_ else level_values[[1]]
      },
      RegionLabel = first_present(RegionLabel),
      .groups = "drop"
    ) %>%
    arrange(RegionKey)

  wide <- df %>%
    group_by(RegionKey, SampleID) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
    mutate(value = if_else(is.nan(value), NA_real_, value)) %>%
    pivot_wider(names_from = SampleID, values_from = value) %>%
    right_join(region_meta %>% distinct(RegionKey), by = "RegionKey") %>%
    arrange(RegionKey)

  mat <- wide %>%
    select(all_of(sample_meta$SampleID)) %>%
    as.matrix()
  rownames(mat) <- wide$RegionKey

  n_by_condition <- df %>%
    filter(!is.na(value)) %>%
    distinct(RegionKey, SampleID, Condition) %>%
    count(RegionKey, Condition, name = "n") %>%
    complete(RegionKey = region_meta$RegionKey, Condition = levels(long$Condition), fill = list(n = 0)) %>%
    pivot_wider(names_from = Condition, values_from = n, names_prefix = "n_")

  n_samples <- rowSums(!is.na(mat))
  variable_region <- apply(mat, 1, function(x) sum(!is.na(x)) >= 4 && sd(x, na.rm = TRUE) > 0)

  result_template <- tidyr::expand_grid(
    RegionKey = region_meta$RegionKey,
    contrast = names(contrast_definitions)
  ) %>%
    left_join(region_meta, by = "RegionKey") %>%
    left_join(n_by_condition, by = "RegionKey") %>%
    mutate(
      logFC = NA_real_,
      AveExpr = NA_real_,
      t = NA_real_,
      P.Value = NA_real_,
      adj.P.Val = NA_real_,
      B = NA_real_,
      n_samples = n_samples[RegionKey],
      note = if_else(variable_region[RegionKey], NA_character_, "insufficient non-missing or non-variable data")
    )

  if (sum(variable_region) < 2 || n_distinct(sample_meta$Condition) < 2) {
    return(
      result_template %>%
        group_by(contrast) %>%
        mutate(adj.P.Val.global = p.adjust(P.Value, method = "BH")) %>%
        ungroup() %>%
        arrange(contrast, adj.P.Val.global, P.Value)
    )
  }

  mat_fit <- mat[variable_region, , drop = FALSE]
  for (j in seq_len(ncol(mat_fit))) {
    missing_rows <- is.na(mat_fit[, j])
    if (any(missing_rows)) {
      row_medians <- apply(mat_fit[missing_rows, , drop = FALSE], 1, median, na.rm = TRUE)
      mat_fit[missing_rows, j] <- row_medians
    }
  }

  design <- model.matrix(~ 0 + Condition, data = sample_meta)
  colnames(design) <- str_remove(colnames(design), "^Condition")

  possible_contrasts <- contrast_definitions[
    vapply(contrast_definitions, function(expr) {
      all(str_extract_all(expr, "[A-Za-z]+_[A-Za-z]+")[[1]] %in% colnames(design))
    }, logical(1))
  ]

  if (length(possible_contrasts) == 0) {
    return(
      result_template %>%
        mutate(note = coalesce(note, "no possible contrast")) %>%
        group_by(contrast) %>%
        mutate(adj.P.Val.global = p.adjust(P.Value, method = "BH")) %>%
        ungroup() %>%
        arrange(contrast, adj.P.Val.global, P.Value)
    )
  }

  fit <- limma::lmFit(mat_fit, design)
  cm <- limma::makeContrasts(contrasts = possible_contrasts, levels = design)
  fit2 <- limma::contrasts.fit(fit, cm)
  fit2 <- limma::eBayes(fit2)

  fitted_results <- purrr::map_dfr(colnames(cm), function(co) {
    limma::topTable(fit2, coef = co, number = Inf, sort.by = "none") %>%
      rownames_to_column("RegionKey") %>%
      as_tibble() %>%
      mutate(contrast = co, note = NA_character_) %>%
      select(RegionKey, contrast, logFC, AveExpr, t, P.Value, adj.P.Val, B, note)
  })

  result_template %>%
    select(-logFC, -AveExpr, -t, -P.Value, -adj.P.Val, -B, -note) %>%
    left_join(fitted_results, by = c("RegionKey", "contrast")) %>%
    mutate(
      n_samples = n_samples[RegionKey],
      note = case_when(
        !variable_region[RegionKey] ~ "insufficient non-missing or non-variable data",
        contrast %in% names(possible_contrasts) ~ note,
        TRUE ~ "contrast not estimated"
      )
    ) %>%
    group_by(contrast) %>%
    mutate(adj.P.Val.global = p.adjust(P.Value, method = "BH")) %>%
    ungroup() %>%
    arrange(contrast, adj.P.Val.global, P.Value)
}

contrast_tables <- list()
for (metric in metrics_to_analyse) {
  contrast_tables[[metric]] <- run_limma_region_contrasts(long, metric)
}

normalized_metrics_to_analyse <- c("Cell_Count_norm", "Intensity_norm")
normalized_metric_lookup <- c(Cell_Count_norm = "Cell_Count", Intensity_norm = "Intensity")

contrast_tables_normalized <- list()
for (metric in normalized_metrics_to_analyse) {
  contrast_tables_normalized[[metric]] <- run_limma_region_contrasts(long, metric)
}

normalization_comparison <- purrr::imap_dfr(normalized_metric_lookup, function(raw_metric, norm_metric) {
  raw_res <- contrast_tables[[raw_metric]] %>%
    filter(contrast %in% central_contrasts, !is.na(logFC)) %>%
    group_by(contrast) %>%
    arrange(desc(abs(logFC)), P.Value, .by_group = TRUE) %>%
    mutate(raw_abs_effect_rank = row_number()) %>%
    ungroup() %>%
    select(RegionKey, contrast, raw_logFC = logFC, raw_P.Value = P.Value, raw_adj.P.Val.global = adj.P.Val.global, raw_abs_effect_rank)

  norm_res <- contrast_tables_normalized[[norm_metric]] %>%
    filter(contrast %in% central_contrasts, !is.na(logFC)) %>%
    group_by(contrast) %>%
    arrange(desc(abs(logFC)), P.Value, .by_group = TRUE) %>%
    mutate(norm_abs_effect_rank = row_number()) %>%
    ungroup() %>%
    select(RegionKey, contrast, norm_logFC = logFC, norm_P.Value = P.Value, norm_adj.P.Val.global = adj.P.Val.global, norm_abs_effect_rank)

  full_join(raw_res, norm_res, by = c("RegionKey", "contrast")) %>%
    mutate(
      Metric = raw_metric,
      NormalizedMetric = norm_metric,
      rank_shift = norm_abs_effect_rank - raw_abs_effect_rank,
      same_direction = sign(raw_logFC) == sign(norm_logFC)
    ) %>%
    left_join(long %>% distinct(RegionKey, Class, Annotation, Abbreviation, RegionLabel), by = "RegionKey") %>%
    select(Metric, NormalizedMetric, contrast, Class, Annotation, Abbreviation, RegionLabel, RegionKey,
           raw_logFC, norm_logFC, same_direction, raw_abs_effect_rank, norm_abs_effect_rank, rank_shift,
           raw_P.Value, norm_P.Value, raw_adj.P.Val.global, norm_adj.P.Val.global)
})

readr::write_csv(normalization_comparison, file.path(out_dir, "tables", "normalization_comparison_raw_vs_normalized_log1p.csv"))
readr::write_csv(normalization_comparison, file.path(legacy_tab_dir, "normalization_comparison_raw_vs_normalized_log1p.csv"))
readr::write_csv(normalization_comparison, file.path(qc_normalization_dir, "normalization_comparison_raw_vs_normalized_log1p.csv"))
openxlsx::write.xlsx(
  c(contrast_tables_normalized, list(normalization_comparison = normalization_comparison)),
  file.path(out_dir, "tables", "region_group_contrasts_limma_normalized_log1p.xlsx"),
  overwrite = TRUE
)
openxlsx::write.xlsx(
  c(contrast_tables_normalized, list(normalization_comparison = normalization_comparison)),
  file.path(legacy_tab_dir, "region_group_contrasts_limma_normalized_log1p.xlsx"),
  overwrite = TRUE
)
openxlsx::write.xlsx(
  c(contrast_tables_normalized, list(normalization_comparison = normalization_comparison)),
  file.path(qc_normalization_dir, "region_group_contrasts_limma_normalized_log1p.xlsx"),
  overwrite = TRUE
)

interpret_region_effects <- function(res, metric, top_n_per_contrast = 25) {
  metric_label <- case_when(
    metric == "Cell_Count" ~ "labelled-cell burden",
    metric == "Intensity" ~ "signal intensity",
    TRUE ~ metric
  )

  res %>%
    filter(!is.na(P.Value), !is.na(logFC)) %>%
    mutate(
      metric = metric,
      direction = case_when(
        logFC > 0 ~ "higher in first condition",
        logFC < 0 ~ "higher in second condition",
        TRUE ~ "no directional change"
      ),
      evidence = case_when(
        adj.P.Val.global < 0.05 ~ "FDR q < 0.05",
        adj.P.Val.global < 0.10 ~ "FDR q < 0.10",
        TRUE ~ "ranked exploratory"
      ),
      biological_read = paste0(
        Annotation, " shows ", direction, " for ", metric_label,
        " in ", contrast, " (logFC = ", scales::number(logFC, accuracy = 0.01),
        ", FDR q = ", scales::pvalue(adj.P.Val.global), ")."
      )
    ) %>%
    group_by(contrast) %>%
    arrange(adj.P.Val.global, P.Value, .by_group = TRUE) %>%
    slice_head(n = top_n_per_contrast) %>%
    ungroup() %>%
    select(
      metric, contrast, evidence, direction,
      Class, Annotation, Abbreviation, RegionLabel, RegionKey,
      logFC, AveExpr, t, P.Value, adj.P.Val.global, n_samples, biological_read
    )
}

candidate_tables <- purrr::imap(
  contrast_tables,
  ~ interpret_region_effects(.x, metric = .y, top_n_per_contrast = 25)
)

openxlsx::write.xlsx(
  c(
    contrast_tables,
    set_names(candidate_tables, paste0(names(candidate_tables), "_ranked_candidates"))
  ),
  file.path(out_dir, "tables", "region_group_contrasts_limma_log1p.xlsx"),
  overwrite = TRUE
)
openxlsx::write.xlsx(
  c(
    contrast_tables,
    set_names(candidate_tables, paste0(names(candidate_tables), "_ranked_candidates"))
  ),
  file.path(legacy_tab_dir, "region_group_contrasts_limma_log1p.xlsx"),
  overwrite = TRUE
)

for (metric in names(candidate_tables)) {
  readr::write_csv(
    candidate_tables[[metric]],
    file.path(out_dir, "tables", paste0("ranked_candidate_regions_", metric, ".csv"))
  )
  readr::write_csv(
    candidate_tables[[metric]],
    file.path(legacy_tab_dir, paste0("ranked_candidate_regions_", metric, ".csv"))
  )
  readr::write_csv(
    candidate_tables[[metric]],
    file.path(exploratory_sensitivity_dir, paste0("Ranked_candidate_regions_", metric, ".csv"))
  )
}

# -----------------------------
# 8. Effect-size heatmaps
# -----------------------------
plot_effect_size_heatmap <- function(res, metric) {
  d <- res %>%
    filter(!is.na(logFC), contrast %in% names(contrast_definitions)) %>%
    select(RegionKey, contrast, logFC) %>%
    pivot_wider(names_from = contrast, values_from = logFC)

  if (nrow(d) < 2) return(NULL)

  mat <- d %>%
    column_to_rownames("RegionKey") %>%
    as.matrix()

  mat <- clean_heatmap_matrix(mat, fill = 0)

  pdf(file.path(exploratory_effect_fig_dir, paste0("effect_size_heatmap_logFC_", metric, ".pdf")),
      width = 5.8, height = max(4.5, min(12, nrow(mat) * 0.10)))
  safe_pheatmap(
    mat,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
    main = paste0("Regional effect-size map: ", display_metric(metric), " logFC"),
    fontsize_row = ifelse(nrow(mat) > 80, 4, 6),
    border_color = NA
  )
  dev.off()
}

purrr::walk2(contrast_tables, names(contrast_tables), ~ plot_effect_size_heatmap(.x, .y))

# -----------------------------
# 9. Volcano plots for biologically central contrasts
# -----------------------------
plot_contrast_volcano <- function(res, metric, contrast_name) {
  d <- res %>%
    filter(contrast == contrast_name, !is.na(adj.P.Val.global), !is.na(logFC)) %>%
    mutate(
      neglog10q = -log10(pmax(adj.P.Val.global, .Machine$double.xmin)),
      fdr_pass = adj.P.Val.global < 0.10,
      label = if_else(fdr_pass, Annotation, NA_character_)
    )

  if (nrow(d) < 3) return(NULL)

  p <- ggplot(d, aes(x = logFC, y = neglog10q)) +
    geom_hline(yintercept = -log10(0.10), linewidth = 0.25, linetype = "dashed", colour = "grey45") +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey45") +
    geom_point(aes(fill = Class), shape = 21, alpha = 0.85, size = 2.0, colour = "grey20", stroke = 0.2) +
    geom_point(data = d %>% filter(fdr_pass), shape = 21, size = 2.8, colour = "black", fill = NA, stroke = 0.35) +
    ggrepel::geom_text_repel(aes(label = label), size = 2.2, max.overlaps = 30) +
    theme_nature(base_size = 8) +
    labs(
      title = paste0(display_metric(metric), ": ", display_contrast(contrast_name)),
      subtitle = paste0("Circle and labels mark FDR q < 0.10. ", contrast_plain_language[contrast_name]),
      x = "Effect size, logFC on log1p scale",
      y = "-log10 FDR q"
    )

  save_figure(p, paste0("volcano_", metric, "_", contrast_name), width = 4.8, height = 4.2, subdir = "volcano_fdr")
}

for (metric in names(contrast_tables)) {
  for (co in central_contrasts) {
    plot_contrast_volcano(contrast_tables[[metric]], metric, co)
  }
}

# -----------------------------
# 10. Top region profile plots
# -----------------------------
plot_top_region_profiles <- function(data, res, metric, contrast_name, top_n = 20) {
  metric_sym <- sym(metric)

  top_regions <- res %>%
    filter(contrast == contrast_name, !is.na(P.Value)) %>%
    arrange(P.Value) %>%
    slice_head(n = top_n) %>%
    pull(RegionKey)

  if (length(top_regions) == 0) return(NULL)

  d <- data %>%
    filter(RegionKey %in% top_regions) %>%
    mutate(value = safe_log1p(!!metric_sym))

  p <- ggplot(d, aes(x = Condition, y = value)) +
    geom_boxplot(aes(colour = Condition), outlier.shape = NA, width = 0.55, linewidth = 0.25) +
    geom_jitter(aes(fill = Condition), width = 0.12, size = 1.2, alpha = 0.7, shape = 21, colour = "white", stroke = 0.15) +
    facet_wrap(~ RegionLabel, scales = "free_y", ncol = 3) +
    scale_x_discrete(labels = condition_short_labels, drop = FALSE) +
    scale_colour_manual(values = condition_colors, guide = "none", drop = FALSE) +
    scale_fill_manual(values = condition_colors, labels = condition_display_labels, name = "Condition", drop = FALSE) +
    theme_nature(base_size = 8) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      strip.text = element_text(size = 6),
      legend.position = "bottom"
    ) +
    labs(
      title = paste0("Regional profiles: ", display_metric(metric)),
      subtitle = display_contrast(contrast_name),
      x = NULL,
      y = paste0("log1p ", display_metric(metric))
    )

  save_figure(p, paste0("top_region_profiles_", metric, "_", contrast_name), width = 6.8, height = 7.2, subdir = "regional_profiles")
}

for (metric in names(contrast_tables)) {
  for (co in central_contrasts) {
    plot_top_region_profiles(long, contrast_tables[[metric]], metric, co, top_n = 20)
  }
}

# -----------------------------
# 11. Class-level aggregation and contrasts
# -----------------------------
class_level_long <- long %>%
  pivot_longer(
    cols = all_of(metrics_to_analyse),
    names_to = "Metric",
    values_to = "RawValue"
  ) %>%
  mutate(Value = safe_log1p(RawValue)) %>%
  group_by(SampleID, Animal, Group, Condition, Class, Metric) %>%
  summarise(
    ClassMean = mean(Value, na.rm = TRUE),
    ClassMedian = median(Value, na.rm = TRUE),
    n_regions = n_distinct(Annotation),
    .groups = "drop"
  )

openxlsx::write.xlsx(class_level_long, file.path(out_dir, "tables", "class_level_long.xlsx"), overwrite = TRUE)
openxlsx::write.xlsx(class_level_long, file.path(legacy_tab_dir, "class_level_long.xlsx"), overwrite = TRUE)

plot_class_profiles <- function(metric) {
  d <- class_level_long %>% filter(Metric == metric)

  p <- ggplot(d, aes(x = Condition, y = ClassMean)) +
    geom_boxplot(aes(colour = Condition), outlier.shape = NA, width = 0.55, linewidth = 0.25) +
    geom_jitter(aes(fill = Condition), width = 0.12, size = 1.1, alpha = 0.75, shape = 21, colour = "white", stroke = 0.15) +
    facet_wrap(~ Class, scales = "free_y", ncol = 4) +
    scale_x_discrete(labels = condition_short_labels, drop = FALSE) +
    scale_colour_manual(values = condition_colors, guide = "none", drop = FALSE) +
    scale_fill_manual(values = condition_colors, labels = condition_display_labels, name = "Condition", drop = FALSE) +
    theme_nature(base_size = 8) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom") +
    labs(
      title = paste0("Class-level profiles: ", display_metric(metric)),
      x = NULL,
      y = paste0("Mean log1p ", display_metric(metric), " across regions")
    )

  save_figure(p, paste0("class_level_profiles_", metric), width = 6.8, height = 6.2, subdir = "regional_profiles/class_level")
}

purrr::walk(metrics_to_analyse, plot_class_profiles)

# -----------------------------
# 12. PCA and UMAP: animal-level systems structure
# -----------------------------
plot_dimensionality <- function(data, metric) {
  sm <- make_sample_matrix(data, metric)
  mat <- safe_log1p(sm$mat)

  keep <- colSums(!is.na(mat)) >= min_pairwise_n &
    apply(mat, 2, function(x) sd(x, na.rm = TRUE) > 0)

  mat <- mat[, keep, drop = FALSE]
  if (nrow(mat) < 4 || ncol(mat) < 3) return(NULL)

  # Impute missing values by region median.
  for (j in seq_len(ncol(mat))) {
    mat[is.na(mat[, j]), j] <- median(mat[, j], na.rm = TRUE)
  }

  mat_z <- zscore_cols(mat)

  pca <- prcomp(mat_z, center = FALSE, scale. = FALSE)
  pca_df <- as_tibble(pca$x[, 1:2], rownames = "SampleID") %>%
    left_join(sm$annotation, by = "SampleID")

  var_exp <- round(100 * summary(pca)$importance[2, 1:2], 1)

  p1 <- ggplot(pca_df, aes(x = PC1, y = PC2, fill = Condition)) +
    geom_point(shape = 21, size = 3.2, colour = "grey20", stroke = 0.25, alpha = 0.9) +
    scale_fill_manual(values = condition_colors, drop = FALSE) +
    theme_nature(base_size = 8) +
    labs(
      title = paste0("PCA of regional tracing profiles: ", metric),
      x = paste0("PC1 (", var_exp[1], "%)"),
      y = paste0("PC2 (", var_exp[2], "%)")
    )

  save_figure(p1, paste0("PCA_systems_structure_", metric), width = 5.5, height = 4.5)
  ggsave(file.path(exploratory_dimred_dir, paste0("PCA_descriptive_animal_level_", metric, ".pdf")), p1, width = 5.5, height = 4.5)
  ggsave(file.path(exploratory_dimred_dir, paste0("PCA_descriptive_animal_level_", metric, ".svg")), p1, width = 5.5, height = 4.5, device = svglite::svglite)
  ggsave(file.path(exploratory_dimred_dir, paste0("PCA_descriptive_animal_level_", metric, ".png")), p1, width = 5.5, height = 4.5, dpi = 600, bg = "white")

  set.seed(1)
  umap_mat <- uwot::umap(mat_z, n_neighbors = min(6, nrow(mat_z) - 1), min_dist = 0.2, metric = "euclidean")

  umap_df <- tibble(
    SampleID = rownames(mat_z),
    UMAP1 = umap_mat[, 1],
    UMAP2 = umap_mat[, 2]
  ) %>%
    left_join(sm$annotation, by = "SampleID")

  p2 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, fill = Condition)) +
    geom_point(shape = 21, size = 3.2, colour = "grey20", stroke = 0.25, alpha = 0.9) +
    scale_fill_manual(values = condition_colors, drop = FALSE) +
    theme_nature(base_size = 8) +
    labs(
      title = paste0("UMAP of regional tracing profiles: ", metric),
      x = "UMAP1",
      y = "UMAP2"
    )

  ggsave(file.path(exploratory_dimred_dir, paste0("UMAP_exploratory_animal_level_", metric, ".pdf")), p2, width = 5.5, height = 4.5)
  ggsave(file.path(exploratory_dimred_dir, paste0("UMAP_exploratory_animal_level_", metric, ".svg")), p2, width = 5.5, height = 4.5, device = svglite::svglite)
  ggsave(file.path(exploratory_dimred_dir, paste0("UMAP_exploratory_animal_level_", metric, ".png")), p2, width = 5.5, height = 4.5, dpi = 600, bg = "white")

  openxlsx::write.xlsx(
    list(PCA = pca_df, UMAP = umap_df),
    file.path(out_dir, "tables", paste0("dimensionality_outputs_", metric, ".xlsx")),
    overwrite = TRUE
  )
  openxlsx::write.xlsx(
    list(PCA = pca_df, UMAP = umap_df),
    file.path(legacy_tab_dir, paste0("dimensionality_outputs_", metric, ".xlsx")),
    overwrite = TRUE
  )
  openxlsx::write.xlsx(
    list(PCA_descriptive = pca_df, UMAP_exploratory = umap_df),
    file.path(exploratory_dimred_dir, paste0("Dimensionality_reduction_descriptive_exploratory_", metric, ".xlsx")),
    overwrite = TRUE
  )
}

purrr::walk(metrics_to_analyse, ~ plot_dimensionality(long, .x))

# -----------------------------
# 13. CeM-centered connectivity
# -----------------------------
find_cem_region <- function(data) {
  candidates <- data %>%
    distinct(RegionKey, Annotation, Abbreviation, Class) %>%
    filter(
      str_detect(str_to_lower(Annotation), "central amygdalar nucleus, medial") |
        str_to_lower(Annotation) %in% c("cem", "ceam") |
        str_to_lower(Abbreviation) %in% c("cem", "ceam")
    )

  if (nrow(candidates) == 0) {
    warning("No CeM/CeAm region detected. Check Annotation/Abbreviation labels.")
    return(NA_character_)
  }

  candidates$RegionKey[1]
}

run_cem_connectivity <- function(data, metric) {
  cem_region <- find_cem_region(data)
  if (is.na(cem_region)) return(NULL)

  sm <- make_sample_matrix(data, metric)
  mat <- safe_log1p(sm$mat)

  if (!cem_region %in% colnames(mat)) {
    warning("Detected CeM region is not present in matrix columns: ", cem_region)
    return(NULL)
  }

  keep <- colSums(!is.na(mat)) >= min_pairwise_n &
    apply(mat, 2, function(x) sd(x, na.rm = TRUE) > 0)

  mat <- mat[, keep, drop = FALSE]
  if (!cem_region %in% colnames(mat)) return(NULL)

  anno <- sm$annotation
  conditions <- levels(long$Condition)

  cem_all <- purrr::map_dfr(conditions, function(cond) {
    ids <- anno %>%
      filter(Condition == cond) %>%
      pull(SampleID)

    submat <- mat[rownames(mat) %in% ids, , drop = FALSE]
    if (nrow(submat) < min_pairwise_n) return(NULL)

    tibble(region = colnames(submat)) %>%
      rowwise() %>%
      mutate(
        Condition = cond,
        n_pair = sum(complete.cases(submat[, cem_region], submat[, region])),
        r = ifelse(
          region != cem_region && n_pair >= min_pairwise_n,
          suppressWarnings(cor(submat[, cem_region], submat[, region], method = "spearman", use = "pairwise.complete.obs")),
          NA_real_
        ),
        p = ifelse(
          region != cem_region && n_pair >= min_pairwise_n && !is.na(r),
          suppressWarnings(cor.test(submat[, cem_region], submat[, region], method = "spearman", exact = FALSE)$p.value),
          NA_real_
        )
      ) %>%
      ungroup()
  }) %>%
    mutate(fdr = p.adjust(p, method = "BH"))

  readr::write_csv(cem_all, file.path(out_dir, "tables", paste0("CeM_centered_connectivity_", metric, ".csv")))

  top_cem <- cem_all %>%
    filter(!is.na(r), region != cem_region) %>%
    group_by(Condition) %>%
    arrange(p, .by_group = TRUE) %>%
    slice_head(n = 20) %>%
    ungroup() %>%
    mutate(
      RegionShort = str_remove(region, "^.* :: "),
      RegionShort = str_trunc(RegionShort, 45)
    )

  if (nrow(top_cem) < 3) return(NULL)

  p <- ggplot(top_cem, aes(x = r, y = reorder(RegionShort, r))) +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey55") +
    geom_point(aes(size = -log10(p), fill = Condition), shape = 21, colour = "grey20", stroke = 0.2, alpha = 0.9) +
    facet_wrap(~ Condition, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = condition_colors, drop = FALSE) +
    scale_size_continuous(range = c(1.2, 4.0), name = "-log10 P") +
    theme_nature(base_size = 8) +
    theme(axis.text.y = element_text(size = 6)) +
    labs(
      title = paste0("CeM-centered regional connectivity: ", metric),
      subtitle = paste0("Seed region: ", cem_region),
      x = "Spearman correlation with CeM",
      y = NULL
    )

  save_figure(p, paste0("CeM_centered_connectivity_", metric), width = 6.4, height = 6.2, subdir = "covariance")
}

purrr::walk(metrics_to_analyse, ~ run_cem_connectivity(long, .x))

# -----------------------------
# 14. Legacy exploratory main figure
# -----------------------------
main_fig_dir <- file.path(legacy_dir, "figures", "main_figure")
main_tab_dir <- file.path(legacy_dir, "tables", "main_figure")
dir.create(main_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(main_tab_dir, recursive = TRUE, showWarnings = FALSE)

nature_condition_levels <- c("VEH_paired", "VEH_unpaired", "CNO_paired", "CNO_unpaired")
nature_condition_colors <- condition_colors

nature_condition_labels <- tibble::tribble(
  ~Condition, ~Biological_interpretation,
  "VEH_paired", "Associative learning",
  "VEH_unpaired", "Stress / non-associative aversive exposure",
  "CNO_paired", "CeM manipulation during learning",
  "CNO_unpaired", "CeM manipulation during stress"
) %>%
  mutate(Condition = factor(Condition, levels = nature_condition_levels))

readr::write_csv(nature_condition_labels, file.path(main_tab_dir, "main_figure_condition_key.csv"))

theme_nature_main <- function(base_size = 7.5) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(colour = "black"),
      plot.title = element_text(face = "bold", size = base_size + 0.5),
      plot.subtitle = element_text(size = base_size - 0.2, colour = "grey25"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.5, colour = "black"),
      axis.line = element_line(linewidth = 0.25, colour = "black"),
      axis.ticks = element_line(linewidth = 0.25, colour = "black"),
      legend.title = element_text(size = base_size - 0.5),
      legend.text = element_text(size = base_size - 0.5),
      legend.key.size = grid::unit(3.2, "mm"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size - 0.5),
      plot.margin = margin(5, 6, 5, 6)
    )
}

write_main_warning <- function(filename, text) {
  writeLines(text, con = file.path(main_tab_dir, filename))
  warning(text)
}

safe_main_ggsave <- function(filename, plot, ...) {
  tryCatch(
    {
      ggsave(filename, plot, ...)
      invisible(filename)
    },
    error = function(e) {
      fallback <- file.path(
        dirname(filename),
        paste0(
          tools::file_path_sans_ext(basename(filename)),
          "_", format(Sys.time(), "%Y%m%d_%H%M%S"),
          ".", tools::file_ext(filename)
        )
      )
      ggsave(fallback, plot, ...)
      write_main_warning(
        "main_figure_export_warning.txt",
        paste0("Could not overwrite ", filename, ". Wrote timestamped fallback instead: ", fallback,
               ". Original error: ", conditionMessage(e))
      )
      invisible(fallback)
    }
  )
}

find_cem_region_main <- function(data) {
  candidates <- data %>%
    distinct(RegionKey, Annotation, Abbreviation, Class) %>%
    mutate(
      ann_lower = str_to_lower(coalesce(Annotation, "")),
      abbr_lower = str_to_lower(coalesce(Abbreviation, "")),
      score = case_when(
        abbr_lower %in% c("cem", "ceam", "cea-m", "cea m") ~ 1L,
        ann_lower == "cem" | ann_lower == "ceam" ~ 2L,
        str_detect(ann_lower, "central amygdalar nucleus, medial") ~ 3L,
        str_detect(ann_lower, "central amygdalar nucleus medial") ~ 4L,
        str_detect(abbr_lower, "cem|ceam|cea") ~ 5L,
        str_detect(ann_lower, "cem|ceam|central amygdalar") ~ 6L,
        TRUE ~ NA_integer_
      )
    ) %>%
    filter(!is.na(score)) %>%
    arrange(score, RegionKey)

  if (nrow(candidates) == 0) return(NULL)
  candidates %>% slice(1)
}

make_panel_a <- function() {
  d <- nature_condition_labels %>%
    mutate(
      y = rev(seq_len(n())),
      Biological_interpretation = str_wrap(Biological_interpretation, width = 46)
    )

  ggplot(d, aes(y = y)) +
    geom_tile(aes(x = 0.04, fill = Condition), width = 0.025, height = 0.55, colour = NA) +
    geom_text(aes(x = 0.08, label = Condition), hjust = 0, size = 2.25, fontface = "bold") +
    geom_text(aes(x = 0.34, label = Biological_interpretation), hjust = 0, size = 2.15, lineheight = 0.9) +
    scale_fill_manual(values = nature_condition_colors, guide = "none", drop = FALSE) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0.45, 4.55), clip = "off") +
    theme_void(base_size = 7.2) +
    theme(
      plot.title = element_text(face = "bold", size = 8),
      plot.margin = margin(3, 6, 3, 6)
    ) +
    labs(title = "Condition framework")
}

get_cem_values <- function(data, cem_region, metric) {
  metric_sym <- sym(metric)

  data %>%
    filter(RegionKey == cem_region) %>%
    transmute(
      SampleID, Animal, Condition,
      Class, Annotation, Abbreviation, RegionLabel, RegionKey,
      Metric = metric,
      RawValue = !!metric_sym,
      Value = safe_log1p(!!metric_sym)
    ) %>%
    mutate(Condition = factor(as.character(Condition), levels = nature_condition_levels))
}

get_cem_stats <- function(cem_region) {
  target_contrasts <- c(
    "Learning_effect",
    "CeM_manipulation_during_learning",
    "CeM_manipulation_during_stress"
  )

  purrr::imap_dfr(contrast_tables, function(res, metric) {
    res %>%
      filter(RegionKey == cem_region, contrast %in% target_contrasts) %>%
      mutate(
        Metric = metric,
        DisplayContrast = recode(
          contrast,
          Learning_effect = "VEH_paired - VEH_unpaired",
          CeM_manipulation_during_learning = "CNO_paired - VEH_paired",
          CeM_manipulation_during_stress = "CNO_unpaired - VEH_unpaired"
        )
      ) %>%
      select(
        Metric, contrast, DisplayContrast, Class, Annotation, Abbreviation,
        RegionLabel, RegionKey, logFC, P.Value, adj.P.Val.global, n_samples, note
      )
  })
}

make_cem_metric_plot <- function(cem_values, cem_stats, metric, y_label) {
  d <- cem_values %>% filter(Metric == metric)
  stats_text <- cem_stats %>%
    filter(Metric == metric) %>%
    mutate(label = paste0(
      DisplayContrast,
      ": logFC ", scales::number(logFC, accuracy = 0.01),
      ", P ", scales::pvalue(P.Value, accuracy = 0.001)
    )) %>%
    pull(label) %>%
    paste(collapse = "\n")

  summary_d <- d %>%
    group_by(Condition) %>%
    summarise(
      mean_value = mean(Value, na.rm = TRUE),
      n = sum(!is.na(Value)),
      se = sd(Value, na.rm = TRUE) / sqrt(n),
      ci = qt(0.975, pmax(n - 1, 1)) * se,
      .groups = "drop"
    ) %>%
    mutate(
      ci = if_else(is.finite(ci), ci, 0),
      ymin = mean_value - ci,
      ymax = mean_value + ci
    )

  y_top <- max(c(d$Value, summary_d$ymax), na.rm = TRUE)
  y_bottom <- min(c(d$Value, summary_d$ymin), na.rm = TRUE)
  y_range <- max(y_top - y_bottom, 1)

  ggplot(d, aes(x = Condition, y = Value)) +
    geom_jitter(aes(fill = Condition), width = 0.10, height = 0, shape = 21, size = 1.8,
                colour = "grey15", stroke = 0.2, alpha = 0.86) +
    geom_errorbar(data = summary_d, aes(x = Condition, y = mean_value, ymin = ymin, ymax = ymax, colour = Condition),
                  width = 0.12, linewidth = 0.35, inherit.aes = FALSE) +
    geom_point(data = summary_d, aes(y = mean_value, fill = Condition),
               shape = 23, size = 2.3, colour = "white", stroke = 0.25) +
    annotate("text", x = 0.65, y = y_top + 0.13 * y_range, label = stats_text,
             hjust = 0, vjust = 1, size = 1.85, lineheight = 0.93) +
    scale_fill_manual(values = nature_condition_colors, drop = FALSE) +
    scale_colour_manual(values = nature_condition_colors, guide = "none", drop = FALSE) +
    coord_cartesian(ylim = c(y_bottom, y_top + 0.18 * y_range), clip = "off") +
    theme_nature_main(base_size = 7.5) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      legend.position = "none"
    ) +
    labs(x = NULL, y = y_label, title = paste0("CeM ", if_else(metric == "Cell_Count", "cFos+ cell count", "projection signal")))
}

make_cfos_projection_scatter <- function() {
  d <- full_join(
    contrast_tables$Cell_Count %>%
      filter(contrast == "Learning_effect") %>%
      select(Class, Annotation, Abbreviation, RegionLabel, RegionKey, Cell_Count_Learning_effect = logFC,
             Cell_Count_P.Value = P.Value),
    contrast_tables$Intensity %>%
      filter(contrast == "Learning_effect") %>%
      select(RegionKey, Intensity_Learning_effect = logFC, Intensity_P.Value = P.Value),
    by = "RegionKey"
  ) %>%
    filter(!is.na(Cell_Count_Learning_effect), !is.na(Intensity_Learning_effect)) %>%
    mutate(
      combined_abs_effect = sqrt(Cell_Count_Learning_effect^2 + Intensity_Learning_effect^2),
      label = if_else(rank(-combined_abs_effect, ties.method = "first") <= 8, Annotation, NA_character_)
    ) %>%
    arrange(desc(combined_abs_effect))

  readr::write_csv(d, file.path(main_tab_dir, "main_figure_cfos_cell_count_projection_scatter.csv"))

  if (nrow(d) < 3) {
    write_main_warning("main_figure_scatter_warning.txt", "cFos cell-count/projection scatter skipped: fewer than 3 regions with paired learning effects.")
    return(NULL)
  }

  ggplot(d, aes(x = Cell_Count_Learning_effect, y = Intensity_Learning_effect)) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey55") +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey55") +
    geom_point(aes(size = combined_abs_effect), shape = 21, fill = "grey70",
               colour = "grey18", stroke = 0.2, alpha = 0.86) +
    ggrepel::geom_text_repel(
      data = d %>% filter(!is.na(label)),
      aes(label = label),
      size = 2.0, min.segment.length = 0,
      segment.size = 0.18, max.overlaps = Inf, seed = 1
    ) +
    scale_size_continuous(range = c(1.2, 3.6), guide = "none") +
    theme_nature_main(base_size = 7.5) +
    labs(
      title = "Activity-projection dissociation",
      x = "Learning effect, cFos+ cell count logFC",
      y = "Learning effect, projection intensity logFC"
    )
}

make_learning_stress_heatmap <- function(top_n = 12) {
  central_heatmap_contrasts <- c(
    "Learning_effect",
    "CeM_manipulation_during_learning",
    "CeM_manipulation_during_stress"
  )

  d <- purrr::imap_dfr(contrast_tables, function(res, metric) {
    res %>%
      filter(contrast %in% central_heatmap_contrasts) %>%
      mutate(
        MetricContrast = paste(metric, contrast, sep = "__"),
        ColumnLabel = recode(
          MetricContrast,
          "Cell_Count__Learning_effect" = "Cell_Count Learning_effect",
          "Intensity__Learning_effect" = "Intensity Learning_effect",
          "Cell_Count__CeM_manipulation_during_learning" = "Cell_Count CeM_manipulation_during_learning",
          "Intensity__CeM_manipulation_during_learning" = "Intensity CeM_manipulation_during_learning",
          "Cell_Count__CeM_manipulation_during_stress" = "Cell_Count CeM_manipulation_during_stress",
          "Intensity__CeM_manipulation_during_stress" = "Intensity CeM_manipulation_during_stress"
        )
      ) %>%
      select(Class, Annotation, Abbreviation, RegionLabel, RegionKey, ColumnLabel, logFC)
  })

  if (nrow(d) == 0) {
    write_main_warning("main_figure_heatmap_warning.txt", "Learning vs stress heatmap skipped: no central contrast logFC values available.")
    return(NULL)
  }

  heatmap_columns <- c(
    "Cell_Count Learning_effect",
    "Intensity Learning_effect",
    "Cell_Count CeM_manipulation_during_learning",
    "Intensity CeM_manipulation_during_learning",
    "Cell_Count CeM_manipulation_during_stress",
    "Intensity CeM_manipulation_during_stress"
  )

  wide_all <- d %>%
    mutate(ColumnLabel = factor(ColumnLabel, levels = heatmap_columns)) %>%
    filter(!is.na(ColumnLabel)) %>%
    pivot_wider(names_from = ColumnLabel, values_from = logFC)

  missing_heatmap_columns <- setdiff(heatmap_columns, names(wide_all))
  if (length(missing_heatmap_columns) > 0) {
    wide_all[missing_heatmap_columns] <- NA_real_
  }

  wide_ranked <- wide_all %>%
    mutate(
      SystemGroup = factor(
        fig4_system_group(Class, Annotation, Abbreviation, RegionLabel),
        levels = fig4_system_levels
      )
    ) %>%
    rowwise() %>%
    mutate(
      max_abs_effect = {
        row_values <- c_across(all_of(heatmap_columns))
        if (all(is.na(row_values))) NA_real_ else max(abs(row_values), na.rm = TRUE)
      }
    ) %>%
    ungroup() %>%
    mutate(max_abs_effect = if_else(is.finite(max_abs_effect), max_abs_effect, NA_real_)) %>%
    filter(!is.na(max_abs_effect), max_abs_effect >= 0.25) %>%
    arrange(SystemGroup, desc(abs(max_abs_effect)), Annotation)

  wide <- wide_ranked %>%
    slice_head(n = min(top_n, nrow(wide_ranked))) %>%
    mutate(RowLabel = str_trunc(paste0(Annotation, " (", Abbreviation, ")"), 45))

  if (nrow(wide) == 0) {
    write_main_warning("main_figure_heatmap_warning.txt", "Learning vs stress heatmap skipped: no regions available after filtering central contrasts.")
    return(NULL)
  }

  readr::write_csv(wide, file.path(main_tab_dir, "main_figure_heatmap_matrix.csv"))

  plot_d <- wide %>%
    select(SystemGroup, RowLabel, all_of(heatmap_columns)) %>%
    pivot_longer(cols = all_of(heatmap_columns), names_to = "ContrastMetric", values_to = "logFC") %>%
    mutate(
      RowLabel = factor(RowLabel, levels = rev(wide$RowLabel)),
      ContrastMetric = factor(ContrastMetric, levels = heatmap_columns),
      Metric = if_else(str_starts(as.character(ContrastMetric), "Cell_Count"), "Cell count", "Intensity"),
      ContrastShort = case_when(
        str_detect(as.character(ContrastMetric), "Learning_effect") ~ "Learning",
        str_detect(as.character(ContrastMetric), "during_learning") ~ "CeM x learning",
        str_detect(as.character(ContrastMetric), "during_stress") ~ "CeM x stress",
        TRUE ~ as.character(ContrastMetric)
      ),
      ColumnShort = paste(ContrastShort, Metric, sep = "\n")
    )

  max_abs <- max(abs(plot_d$logFC), na.rm = TRUE)
  max_abs <- ifelse(is.finite(max_abs) && max_abs > 0, max_abs, 1)
  column_levels <- plot_d %>%
    distinct(ContrastMetric, ColumnShort) %>%
    arrange(ContrastMetric) %>%
    pull(ColumnShort)

  ggplot(plot_d, aes(x = factor(ColumnShort, levels = column_levels), y = RowLabel, fill = logFC)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    facet_grid(SystemGroup ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient2(low = effect_colors[["negative"]], mid = effect_colors[["neutral"]], high = effect_colors[["positive"]],
                         midpoint = 0, limits = c(-max_abs, max_abs), oob = scales::squish,
                         name = "Effect size\n(logFC)") +
    theme_nature_main(base_size = 7.2) +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "right"
    ) +
    labs(title = "Learning vs stress regional effect map", x = NULL, y = NULL)
}

make_main_pca <- function(data) {
  pca_warning <- file.path(main_tab_dir, "main_figure_pca_warning.txt")
  if (file.exists(pca_warning)) file.remove(pca_warning)

  feature_long <- data %>%
    select(SampleID, Animal, Condition, RegionKey, Annotation, Abbreviation, all_of(metrics_to_analyse)) %>%
    pivot_longer(cols = all_of(metrics_to_analyse), names_to = "Metric", values_to = "RawValue") %>%
    mutate(
      Value = safe_log1p(RawValue),
      Feature = paste(Metric, RegionKey, sep = "__")
    ) %>%
    group_by(SampleID, Animal, Condition, Feature) %>%
    summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
    mutate(Value = if_else(is.nan(Value), NA_real_, Value))

  wide <- feature_long %>%
    pivot_wider(names_from = Feature, values_from = Value) %>%
    arrange(factor(as.character(Condition), levels = nature_condition_levels), SampleID)

  if (nrow(wide) < 4) {
    write_main_warning("main_figure_pca_warning.txt", "PCA skipped: fewer than 4 animal-level samples available.")
    return(NULL)
  }

  meta <- wide %>% select(SampleID, Animal, Condition)
  mat <- wide %>% select(-SampleID, -Animal, -Condition) %>% as.data.frame()
  rownames(mat) <- meta$SampleID
  mat <- as.matrix(mat)

  keep <- colSums(!is.na(mat)) >= min_pairwise_n &
    apply(mat, 2, function(x) sd(x, na.rm = TRUE) > 0)
  mat <- mat[, keep, drop = FALSE]

  if (nrow(mat) < 4 || ncol(mat) < 3) {
    write_main_warning(
      "main_figure_pca_warning.txt",
      paste0("PCA skipped: need at least 4 samples and 3 variable features after filtering; found ",
             nrow(mat), " samples and ", ncol(mat), " features.")
    )
    return(NULL)
  }

  for (j in seq_len(ncol(mat))) {
    mat[is.na(mat[, j]), j] <- median(mat[, j], na.rm = TRUE)
  }

  mat_z <- zscore_cols(mat)
  pca <- prcomp(mat_z, center = FALSE, scale. = FALSE)
  var_exp <- 100 * summary(pca)$importance[2, 1:2]

  pca_df <- as_tibble(pca$x[, 1:2], rownames = "SampleID") %>%
    left_join(meta, by = "SampleID") %>%
    mutate(
      Condition = factor(as.character(Condition), levels = nature_condition_levels),
      PC1_variance_percent = var_exp[1],
      PC2_variance_percent = var_exp[2],
      n_features = ncol(mat)
    )

  readr::write_csv(pca_df, file.path(main_tab_dir, "main_figure_pca_coordinates.csv"))
  readr::write_csv(
    tibble(PC = c("PC1", "PC2"), variance_explained_percent = var_exp),
    file.path(main_tab_dir, "main_figure_pca_variance_explained.csv")
  )

  ggplot(pca_df, aes(x = PC1, y = PC2, fill = Condition)) +
    geom_point(shape = 21, size = 2.8, colour = "grey15", stroke = 0.25, alpha = 0.9) +
    scale_fill_manual(values = nature_condition_colors, drop = FALSE) +
    theme_nature_main(base_size = 7.5) +
    theme(legend.position = "right") +
    labs(
      title = "Animal-level systems separation",
      x = paste0("PC1 (", scales::number(var_exp[1], accuracy = 0.1), "%)"),
      y = paste0("PC2 (", scales::number(var_exp[2], accuracy = 0.1), "%)")
    )
}

make_network_region_set <- function(top_n = 14) {
  central_network_contrasts <- c(
    "Learning_effect",
    "CeM_manipulation_during_learning",
    "CeM_manipulation_during_stress"
  )

  network_columns <- c(
    "Cell_Count Learning_effect",
    "Intensity Learning_effect",
    "Cell_Count CeM_manipulation_during_learning",
    "Intensity CeM_manipulation_during_learning",
    "Cell_Count CeM_manipulation_during_stress",
    "Intensity CeM_manipulation_during_stress"
  )

  d <- purrr::imap_dfr(contrast_tables, function(res, metric) {
    res %>%
      filter(contrast %in% central_network_contrasts) %>%
      mutate(ColumnLabel = paste(metric, contrast)) %>%
      select(Class, Annotation, Abbreviation, RegionLabel, RegionKey, ColumnLabel, logFC, P.Value)
  }) %>%
    mutate(
      interpretable = !str_detect(
        str_to_lower(paste(Annotation, Abbreviation, Class)),
        "tract|fiber|commissure|ventricle|root|peduncle|bundle|nerve|arbor vitae|corpus callosum"
      )
    ) %>%
    filter(interpretable)

  if (nrow(d) == 0) {
    write_main_warning("main_figure_network_warning.txt", "Network analysis skipped: no interpretable regions available from central contrasts.")
    return(character())
  }

  wide <- d %>%
    select(Class, Annotation, Abbreviation, RegionLabel, RegionKey, ColumnLabel, logFC) %>%
    pivot_wider(names_from = ColumnLabel, values_from = logFC)

  missing_network_columns <- setdiff(network_columns, names(wide))
  if (length(missing_network_columns) > 0) {
    wide[missing_network_columns] <- NA_real_
  }

  ranked_all <- wide %>%
    mutate(
      SystemGroup = factor(
        fig4_system_group(Class, Annotation, Abbreviation, RegionLabel),
        levels = fig4_system_levels
      )
    ) %>%
    rowwise() %>%
    mutate(
      max_abs_effect = {
        row_values <- c_across(all_of(network_columns))
        if (all(is.na(row_values))) NA_real_ else max(abs(row_values), na.rm = TRUE)
      }
    ) %>%
    ungroup() %>%
    mutate(
      max_abs_effect = if_else(is.finite(max_abs_effect), max_abs_effect, NA_real_),
      RegionShort = str_trunc(paste0(Annotation, " (", Abbreviation, ")"), 42)
    ) %>%
    filter(!is.na(max_abs_effect), max_abs_effect >= 0.25) %>%
    arrange(SystemGroup, desc(max_abs_effect), RegionShort)

  ranked <- ranked_all %>%
    slice_head(n = min(top_n, nrow(ranked_all)))

  readr::write_csv(ranked, file.path(main_tab_dir, "main_figure_network_region_set.csv"))
  ranked$RegionKey
}

make_network_feature_matrix <- function(data, region_keys, metric) {
  metric_sym <- sym(metric)

  data %>%
    filter(RegionKey %in% region_keys) %>%
    mutate(
      NetworkCondition = factor(as.character(Condition), levels = nature_condition_levels)
    ) %>%
    filter(!is.na(NetworkCondition)) %>%
    select(SampleID, Animal, NetworkCondition, RegionKey, RawValue = !!metric_sym) %>%
    mutate(
      Value = safe_log1p(RawValue),
      Metric = metric,
      Feature = paste(metric, RegionKey, sep = "__")
    ) %>%
    group_by(SampleID, Animal, NetworkCondition, Feature) %>%
    summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
    mutate(Value = if_else(is.nan(Value), NA_real_, Value)) %>%
    pivot_wider(names_from = Feature, values_from = Value) %>%
    arrange(NetworkCondition, SampleID)
}

build_condition_network <- function(mat,
                                    abs_r_cutoff = 0.60,
                                    min_complete_n = 4,
                                    p_cutoff = 0.05,
                                    fdr_cutoff = 0.10,
                                    edge_rule = c("fdr", "effect_size_only", "nominal_p")) {
  edge_rule <- match.arg(edge_rule)
  features <- colnames(mat)
  empty_edges <- tibble(
    feature1 = character(),
    feature2 = character(),
    n_pair = integer(),
    edge_eligible = logical(),
    r = numeric(),
    p_value = numeric(),
    abs_r = numeric(),
    sign = character(),
    fdr = numeric(),
    edge_effect_size = logical(),
    edge_nominal = logical(),
    edge_fdr = logical(),
    edge_present = logical(),
    edge_rule_used = character()
  )

  if (nrow(mat) < min_complete_n || length(features) < 3) {
    return(list(
      edges = empty_edges,
      metrics = tibble(),
      hubs = tibble(),
      graph = igraph::make_empty_graph()
    ))
  }

  edges <- combn(features, 2, simplify = FALSE) %>%
    purrr::map_dfr(function(pair) {
      x <- mat[, pair[1]]
      y <- mat[, pair[2]]

      keep <- is.finite(x) & is.finite(y)
      n_pair <- sum(keep)

      x_ok <- n_pair >= min_complete_n && stats::sd(x[keep], na.rm = TRUE) > 0
      y_ok <- n_pair >= min_complete_n && stats::sd(y[keep], na.rm = TRUE) > 0
      edge_eligible <- x_ok && y_ok

      if (edge_eligible) {
        test <- suppressWarnings(
          tryCatch(
            stats::cor.test(
              x[keep],
              y[keep],
              method = "spearman",
              exact = FALSE
            ),
            error = function(e) NULL
          )
        )

        r <- if (is.null(test)) NA_real_ else unname(test$estimate)
        p_value <- if (is.null(test)) NA_real_ else test$p.value
      } else {
        r <- NA_real_
        p_value <- NA_real_
      }

      tibble(
        feature1 = pair[1],
        feature2 = pair[2],
        n_pair = n_pair,
        edge_eligible = edge_eligible,
        r = r,
        p_value = p_value,
        abs_r = abs(r),
        sign = case_when(
          r > 0 ~ "positive",
          r < 0 ~ "negative",
          TRUE ~ NA_character_
        )
      )
    }) %>%
    mutate(
      fdr = p.adjust(p_value, method = "BH"),
      edge_effect_size = edge_eligible & !is.na(r) & abs_r >= abs_r_cutoff,
      edge_nominal = edge_effect_size & !is.na(p_value) & p_value <= p_cutoff,
      edge_fdr = edge_effect_size & !is.na(fdr) & fdr <= fdr_cutoff,
      edge_present = if (edge_rule == "effect_size_only") {
        edge_effect_size
      } else if (edge_rule == "nominal_p") {
        edge_nominal
      } else if (edge_rule == "fdr") {
        edge_fdr
      } else {
        FALSE
      },
      edge_rule_used = case_when(
        edge_rule == "effect_size_only" ~ paste0("|rho| >= ", abs_r_cutoff, ", n_pair >= ", min_complete_n),
        edge_rule == "nominal_p" ~ paste0("|rho| >= ", abs_r_cutoff, ", P <= ", p_cutoff, ", n_pair >= ", min_complete_n),
        edge_rule == "fdr" ~ paste0("|rho| >= ", abs_r_cutoff, ", FDR <= ", fdr_cutoff, ", n_pair >= ", min_complete_n)
      )
    )

  present_edges <- edges %>%
    filter(edge_present) %>%
    transmute(
      from = feature1,
      to = feature2,
      r,
      p_value,
      fdr,
      abs_r,
      sign,
      n_pair,
      edge_eligible,
      edge_effect_size,
      edge_nominal,
      edge_fdr,
      weight = abs_r
    )

  graph <- igraph::graph_from_data_frame(
    present_edges,
    directed = FALSE,
    vertices = tibble(name = features)
  )

  graph_density <- igraph::edge_density(graph, loops = FALSE)

  mean_abs_r <- if (nrow(present_edges) > 0) {
    mean(present_edges$abs_r, na.rm = TRUE)
  } else {
    NA_real_
  }

  modularity_value <- NA_real_
  if (igraph::ecount(graph) > 0 && igraph::vcount(graph) > 2) {
    community <- suppressWarnings(
      igraph::cluster_louvain(graph, weights = igraph::E(graph)$weight)
    )
    modularity_value <- igraph::modularity(community)
  }

  metrics <- tibble(
    n_samples = nrow(mat),
    n_features = length(features),
    n_edges = igraph::ecount(graph),
    possible_edges = length(features) * (length(features) - 1) / 2,
    n_eligible_edges = sum(edges$edge_eligible, na.rm = TRUE),
    n_effect_size_edges = sum(edges$edge_effect_size, na.rm = TRUE),
    n_nominal_edges = sum(edges$edge_nominal, na.rm = TRUE),
    n_fdr_edges = sum(edges$edge_fdr, na.rm = TRUE),
    density = graph_density,
    mean_abs_r = mean_abs_r,
    modularity = modularity_value,
    edge_rule = unique(edges$edge_rule_used)[1]
  )

  hubs <- tibble(
    Feature = features,
    degree = igraph::degree(graph, v = features),
    strength = igraph::strength(graph, v = features, weights = igraph::E(graph)$weight),
    betweenness = igraph::betweenness(graph, v = features, directed = FALSE, weights = NA),
    eigen_centrality = igraph::eigen_centrality(
      graph,
      directed = FALSE,
      weights = igraph::E(graph)$weight
    )$vector[features]
  ) %>%
    mutate(
      Metric = str_extract(Feature, "^[^_]+(?:_[^_]+)?"),
      RegionKey = str_remove(Feature, "^[^_]+(?:_[^_]+)?__")
    )

  list(edges = edges, metrics = metrics, hubs = hubs, graph = graph)
}

compare_network_edges <- function(edge_tables, reference_condition = "VEH_unpaired") {
  d <- edge_tables %>%
    select(
      NetworkCondition,
      feature1,
      feature2,
      r,
      p_value,
      fdr,
      abs_r,
      sign,
      n_pair,
      edge_eligible,
      edge_effect_size,
      edge_nominal,
      edge_fdr,
      edge_present
    ) %>%
    mutate(
      edge_id = paste(pmin(feature1, feature2), pmax(feature1, feature2), sep = " -- ")
    ) %>%
    select(
      NetworkCondition,
      edge_id,
      r,
      p_value,
      fdr,
      abs_r,
      sign,
      n_pair,
      edge_eligible,
      edge_effect_size,
      edge_nominal,
      edge_fdr,
      edge_present
    ) %>%
    pivot_wider(
      names_from = NetworkCondition,
      values_from = c(
        r,
        p_value,
        fdr,
        abs_r,
        sign,
        n_pair,
        edge_eligible,
        edge_effect_size,
        edge_nominal,
        edge_fdr,
        edge_present
      ),
      values_fill = list(
        edge_eligible = FALSE,
        edge_effect_size = FALSE,
        edge_nominal = FALSE,
        edge_fdr = FALSE,
        edge_present = FALSE
      )
    )

  compare_conditions <- setdiff(nature_condition_levels, reference_condition)

  for (cond in compare_conditions) {
    ref_present <- paste0("edge_present_", reference_condition)
    cond_present <- paste0("edge_present_", cond)

    ref_eligible <- paste0("edge_eligible_", reference_condition)
    cond_eligible <- paste0("edge_eligible_", cond)

    ref_sign <- paste0("sign_", reference_condition)
    cond_sign <- paste0("sign_", cond)

    ref_abs <- paste0("abs_r_", reference_condition)
    cond_abs <- paste0("abs_r_", cond)

    ref_n <- paste0("n_pair_", reference_condition)
    cond_n <- paste0("n_pair_", cond)

    rewiring_col <- paste0("rewiring_", cond, "_vs_", reference_condition)
    comparable_col <- paste0("comparable_", cond, "_vs_", reference_condition)
    delta_col <- paste0("delta_abs_r_", cond, "_vs_", reference_condition)

    if (!ref_present %in% names(d)) d[[ref_present]] <- FALSE
    if (!cond_present %in% names(d)) d[[cond_present]] <- FALSE

    if (!ref_eligible %in% names(d)) d[[ref_eligible]] <- FALSE
    if (!cond_eligible %in% names(d)) d[[cond_eligible]] <- FALSE

    if (!ref_sign %in% names(d)) d[[ref_sign]] <- NA_character_
    if (!cond_sign %in% names(d)) d[[cond_sign]] <- NA_character_

    if (!ref_abs %in% names(d)) d[[ref_abs]] <- NA_real_
    if (!cond_abs %in% names(d)) d[[cond_abs]] <- NA_real_

    if (!ref_n %in% names(d)) d[[ref_n]] <- NA_integer_
    if (!cond_n %in% names(d)) d[[cond_n]] <- NA_integer_

    d[[comparable_col]] <- d[[ref_eligible]] & d[[cond_eligible]]

    d[[rewiring_col]] <- case_when(
      !d[[comparable_col]] ~ "not comparable",
      d[[ref_present]] & !d[[cond_present]] ~ "lost",
      !d[[ref_present]] & d[[cond_present]] ~ "gained",
      d[[ref_present]] & d[[cond_present]] & d[[ref_sign]] != d[[cond_sign]] ~ "sign switch",
      d[[ref_present]] & d[[cond_present]] ~ "retained",
      TRUE ~ "absent"
    )

    d[[delta_col]] <- if_else(
      d[[comparable_col]],
      d[[cond_abs]] - d[[ref_abs]],
      NA_real_
    )
  }

  d
}

format_network_feature_label <- function(feature) {
  metric <- case_when(
    str_starts(feature, "Cell_Count__") ~ "Cell",
    str_starts(feature, "Intensity__") ~ "Int",
    TRUE ~ "Signal"
  )

  region_key <- str_remove(feature, "^(Cell_Count|Intensity)__")
  region_label <- str_remove(region_key, "^.* :: ")
  annotation <- str_squish(str_extract(region_label, "^[^|]+"))
  abbreviation <- str_squish(str_replace(region_label, "^.*\\|", ""))

  region_short <- case_when(
    !is.na(annotation) & annotation != "" & nchar(annotation) <= 14 ~ annotation,
    !is.na(abbreviation) & abbreviation != "" & abbreviation != region_label & nchar(abbreviation) <= 14 ~ abbreviation,
    !is.na(annotation) & annotation != "" ~ annotation,
    TRUE ~ abbreviation
  )
  paste(metric, str_trunc(region_short, 18))
}

format_network_edge_label <- function(edge_id) {
  pieces <- str_split(edge_id, " -- ", simplify = TRUE)
  if (ncol(pieces) < 2) return(str_trunc(edge_id, 70))

  paste(
    format_network_feature_label(pieces[, 1]),
    format_network_feature_label(pieces[, 2]),
    sep = " - "
  )
}

make_network_rewiring_figure <- function(data,
                                         region_top_n = 14,
                                         abs_r_cutoff = 0.60,
                                         edge_rule = c("fdr", "effect_size_only", "nominal_p"),
                                         reference_condition = "VEH_unpaired") {
  edge_rule <- match.arg(edge_rule)
  edge_rule_suffix <- paste0("_", edge_rule)
  region_keys <- make_network_region_set(top_n = region_top_n)
  if (length(region_keys) < 4) {
    write_main_warning("main_figure_network_warning.txt", "Network analysis skipped: fewer than 4 interpretable regions selected.")
    return(NULL)
  }

  network_levels <- nature_condition_levels
  network_metric_levels <- c("Cell_Count", "Intensity")
  network_outputs <- purrr::map(network_metric_levels, function(network_metric) {
    wide <- make_network_feature_matrix(data, region_keys, metric = network_metric)
    if (nrow(wide) < 6) {
      write_main_warning(
        paste0("main_figure_network_warning_", network_metric, ".txt"),
        paste0("Network analysis skipped for ", network_metric, ": fewer than 6 samples across conditions.")
      )
      return(NULL)
    }

    meta <- wide %>% select(SampleID, Animal, NetworkCondition)
    mat_all <- wide %>% select(-SampleID, -Animal, -NetworkCondition) %>% as.data.frame()
    rownames(mat_all) <- meta$SampleID

    keep <- colSums(!is.na(mat_all)) >= min_pairwise_n &
      apply(mat_all, 2, function(x) sd(x, na.rm = TRUE) > 0)
    mat_all <- as.matrix(mat_all[, keep, drop = FALSE])

    if (ncol(mat_all) < 4) {
      write_main_warning(
        paste0("main_figure_network_warning_", network_metric, ".txt"),
        paste0("Network analysis skipped for ", network_metric, ": fewer than 4 variable region features after filtering.")
      )
      return(NULL)
    }

    network_results <- purrr::map(network_levels, function(cond) {
      ids <- meta %>% filter(NetworkCondition == cond) %>% pull(SampleID)
      mat <- mat_all[rownames(mat_all) %in% ids, , drop = FALSE]
      build_condition_network(
        mat,
        abs_r_cutoff = abs_r_cutoff,
        min_complete_n = min_pairwise_n,
        p_cutoff = 0.05,
        fdr_cutoff = 0.10,
        edge_rule = edge_rule
      )
    })
    names(network_results) <- network_levels

    metrics <- purrr::imap_dfr(network_results, ~ .x$metrics %>% mutate(NetworkCondition = .y, .before = 1)) %>%
      mutate(NetworkMetric = network_metric, .before = 1)
    edges <- purrr::imap_dfr(network_results, ~ .x$edges %>% mutate(NetworkCondition = .y, .before = 1)) %>%
      mutate(NetworkMetric = network_metric, .before = 1)
    hubs <- purrr::imap_dfr(network_results, ~ .x$hubs %>% mutate(NetworkCondition = .y, .before = 1)) %>%
      mutate(NetworkMetric = network_metric, .before = 1)
    edge_rewiring <- compare_network_edges(edges, reference_condition = reference_condition) %>%
      mutate(NetworkMetric = network_metric, .before = 1)

    list(metrics = metrics, edges = edges, hubs = hubs, edge_rewiring = edge_rewiring)
  }) %>%
    purrr::compact()

  if (length(network_outputs) == 0) {
    write_main_warning("main_figure_network_warning.txt", "Network analysis skipped: no metric-specific networks could be constructed.")
    return(NULL)
  }

  metrics <- purrr::map_dfr(network_outputs, "metrics")
  edges <- purrr::map_dfr(network_outputs, "edges")
  hubs <- purrr::map_dfr(network_outputs, "hubs")
  edge_rewiring <- purrr::map_dfr(network_outputs, "edge_rewiring")

  region_lookup <- data %>%
    distinct(RegionKey, Annotation, Abbreviation, RegionLabel, Class) %>%
    mutate(
      SystemGroup = factor(fig4_system_group(Class, Annotation, Abbreviation, RegionLabel), levels = fig4_system_levels),
      RegionShort = str_trunc(paste0(Annotation, " (", Abbreviation, ")"), 38)
    )

  hubs <- hubs %>%
    left_join(region_lookup, by = "RegionKey") %>%
    group_by(NetworkMetric, NetworkCondition) %>%
    mutate(hub_rank = dense_rank(desc(strength))) %>%
    ungroup()

  readr::write_csv(
  metrics,
  file.path(main_tab_dir, paste0("main_figure_network_metrics_by_metric", edge_rule_suffix, ".csv"))
  )
  readr::write_csv(hubs, file.path(main_tab_dir, paste0("main_figure_network_hub_centrality_by_metric", edge_rule_suffix, ".csv")))
  readr::write_csv(edges, file.path(main_tab_dir, paste0("main_figure_network_edges_by_metric", edge_rule_suffix, ".csv")))
  readr::write_csv(edge_rewiring, file.path(main_tab_dir, paste0("main_figure_network_edge_rewiring_by_metric", edge_rule_suffix, ".csv")))
  readr::write_csv(metrics %>% mutate(fdr_scope = "within condition x metric network edge tests"),
                   file.path(exploratory_network_dir, paste0("Network_metrics", edge_rule_suffix, ".csv")))
  readr::write_csv(hubs %>% mutate(fdr_scope = "within condition x metric network edge tests"),
                   file.path(exploratory_network_dir, paste0("Network_hub_centrality", edge_rule_suffix, ".csv")))
  readr::write_csv(edges %>% mutate(fdr_scope = "within condition x metric network edge tests"),
                   file.path(exploratory_network_dir, paste0("Network_edges", edge_rule_suffix, ".csv")))
  readr::write_csv(edge_rewiring %>% mutate(fdr_scope = "within condition x metric network edge tests"),
                   file.path(exploratory_network_dir, paste0("Network_rewiring", edge_rule_suffix, ".csv")))

  reference_suffix <- paste0("_vs_", reference_condition)
  rewiring_columns <- grep(paste0("^rewiring_.*", reference_suffix, "$"), names(edge_rewiring), value = TRUE)
  edge_rewiring_changed <- edge_rewiring %>%
    filter(if_any(all_of(rewiring_columns), ~ .x %in% c("gained", "lost", "sign switch", "retained"))) %>%
    arrange(NetworkMetric, edge_id)
  readr::write_csv(
    edge_rewiring_changed,
    file.path(main_tab_dir, paste0("main_figure_network_edge_rewiring_changed_or_retained", edge_rule_suffix, ".csv"))
  )

  # Backward-compatible copies now contain the separated, metric-labelled analyses.
  readr::write_csv(metrics, file.path(main_tab_dir, paste0("main_figure_network_metrics", edge_rule_suffix, ".csv")))
  readr::write_csv(hubs, file.path(main_tab_dir, paste0("main_figure_network_hub_centrality", edge_rule_suffix, ".csv")))
  readr::write_csv(edges, file.path(main_tab_dir, paste0("main_figure_network_edges", edge_rule_suffix, ".csv")))
  readr::write_csv(edge_rewiring, file.path(main_tab_dir, paste0("main_figure_network_edge_rewiring", edge_rule_suffix, ".csv")))

  edge_overlap <- edges %>%
    filter(edge_present) %>%
    mutate(
      edge_id = paste(pmin(feature1, feature2), pmax(feature1, feature2), sep = " -- "),
      EdgeShort = format_network_edge_label(edge_id)
    ) %>%
    select(
      NetworkCondition, NetworkMetric, edge_id, EdgeShort,
      r, p_value, fdr, abs_r, sign, n_pair
    ) %>%
    pivot_wider(
      names_from = NetworkMetric,
      values_from = c(r, p_value, fdr, abs_r, sign, n_pair),
      names_sep = "__"
    ) %>%
    {
      d <- .
      numeric_cols <- c(
        "r__Cell_Count", "r__Intensity",
        "p_value__Cell_Count", "p_value__Intensity",
        "fdr__Cell_Count", "fdr__Intensity",
        "abs_r__Cell_Count", "abs_r__Intensity",
        "n_pair__Cell_Count", "n_pair__Intensity"
      )
      character_cols <- c("sign__Cell_Count", "sign__Intensity")
      for (col in setdiff(numeric_cols, names(d))) d[[col]] <- NA_real_
      for (col in setdiff(character_cols, names(d))) d[[col]] <- NA_character_
      d
    } %>%
    mutate(
      present_Cell_Count = !is.na(r__Cell_Count),
      present_Intensity = !is.na(r__Intensity),
      overlap_class = case_when(
        present_Cell_Count & present_Intensity ~ "overlap",
        present_Cell_Count ~ "cell-count only",
        present_Intensity ~ "intensity only",
        TRUE ~ "not present"
      ),
      same_direction = present_Cell_Count & present_Intensity & sign__Cell_Count == sign__Intensity,
      rho_difference_intensity_minus_cell = r__Intensity - r__Cell_Count,
      mean_abs_rho = rowMeans(cbind(abs_r__Cell_Count, abs_r__Intensity), na.rm = TRUE),
      edge_rule = edge_rule,
      reference_condition = reference_condition,
      interpretation_note = "Edges are thresholded within each readout and condition. Overlap means the same region pair passes the selected edge rule in both cFos+ cell-count and projection-intensity networks."
    ) %>%
    arrange(
      factor(NetworkCondition, levels = network_levels),
      desc(overlap_class == "overlap"),
      desc(mean_abs_rho),
      EdgeShort
    )

  edge_overlap_summary <- edge_overlap %>%
    count(NetworkCondition, overlap_class, name = "n_edges") %>%
    mutate(
      NetworkCondition = factor(NetworkCondition, levels = network_levels),
      overlap_class = factor(overlap_class, levels = c("overlap", "cell-count only", "intensity only"))
    ) %>%
    arrange(NetworkCondition, overlap_class)

  readr::write_csv(edge_overlap, file.path(main_tab_dir, paste0("main_figure_network_cell_intensity_edge_overlap", edge_rule_suffix, ".csv")))
  readr::write_csv(edge_overlap_summary, file.path(main_tab_dir, paste0("main_figure_network_cell_intensity_edge_overlap_summary", edge_rule_suffix, ".csv")))
  readr::write_csv(edge_overlap %>% mutate(fdr_scope = "within condition x metric network edge tests"),
                   file.path(exploratory_network_dir, paste0("Network_overlap_cell_count_intensity", edge_rule_suffix, ".csv")))
  readr::write_csv(edge_overlap_summary %>% mutate(fdr_scope = "within condition x metric network edge tests"),
                   file.path(exploratory_network_dir, paste0("Network_overlap_cell_count_intensity_summary", edge_rule_suffix, ".csv")))

  metric_plot <- metrics %>%
    select(NetworkMetric, NetworkCondition, density, modularity, mean_abs_r) %>%
    pivot_longer(cols = c(density, modularity, mean_abs_r), names_to = "Measure", values_to = "Value") %>%
    mutate(
      NetworkMetric = factor(NetworkMetric, levels = network_metric_levels),
      NetworkCondition = factor(NetworkCondition, levels = network_levels),
      Measure = recode(Measure, density = "Density", modularity = "Modularity", mean_abs_r = "Mean |r|")
    ) %>%
    ggplot(aes(x = NetworkCondition, y = Value, fill = NetworkCondition)) +
    geom_col(width = 0.62, colour = "grey20", linewidth = 0.18) +
    facet_grid(NetworkMetric ~ Measure, scales = "free_y") +
    scale_fill_manual(values = nature_condition_colors, guide = "none", drop = FALSE) +
    theme_nature_main(base_size = 7.2) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 5.8)) +
    labs(title = "Network architecture by signal type", x = NULL, y = NULL)

  region_hubs <- hubs %>%
    group_by(NetworkMetric, NetworkCondition, RegionKey, Annotation, Abbreviation, RegionShort) %>%
    summarise(
      hub_strength = max(strength, na.rm = TRUE),
      hub_degree = max(degree, na.rm = TRUE),
      hub_betweenness = max(betweenness, na.rm = TRUE),
      hub_eigen_centrality = max(eigen_centrality, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(across(starts_with("hub_"), ~ if_else(is.finite(.x), .x, NA_real_)))

  readr::write_csv(region_hubs, file.path(main_tab_dir, "main_figure_network_region_hub_centrality_by_metric.csv"))
  readr::write_csv(region_hubs, file.path(main_tab_dir, "main_figure_network_region_hub_centrality.csv"))

  hub_ranked_regions <- region_hubs %>%
    group_by(NetworkMetric, RegionKey, Annotation, Abbreviation, RegionShort) %>%
    summarise(max_strength = max(hub_strength, na.rm = TRUE), .groups = "drop") %>%
    mutate(max_strength = if_else(is.finite(max_strength), max_strength, NA_real_)) %>%
    group_by(NetworkMetric) %>%
    arrange(desc(max_strength), .by_group = TRUE) %>%
    slice_head(n = 8) %>%
    ungroup()

  hub_plot_data <- hub_ranked_regions %>%
    select(NetworkMetric, RegionKey, RegionShort) %>%
    left_join(region_hubs, by = c("NetworkMetric", "RegionKey")) %>%
    mutate(
      RegionShort = coalesce(RegionShort.x, RegionShort.y),
      RegionShort = factor(RegionShort, levels = rev(unique(RegionShort))),
      NetworkMetric = factor(NetworkMetric, levels = network_metric_levels),
      NetworkCondition = factor(NetworkCondition, levels = network_levels)
    )

  hub_plot <- ggplot(hub_plot_data, aes(x = NetworkCondition, y = RegionShort, fill = hub_strength)) +
    geom_tile(colour = "white", linewidth = 0.22) +
    facet_wrap(~ NetworkMetric, scales = "free_y", ncol = 1) +
    scale_fill_gradient(low = effect_colors[["neutral"]], high = "#2B2B2B", name = "Hub\nstrength") +
    theme_nature_main(base_size = 7.0) +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1, size = 5.8),
      axis.text.y = element_text(size = 5.7)
    ) +
    labs(title = "Top hubs, separated by signal type", x = NULL, y = NULL)

  rewiring_summary <- edge_rewiring %>%
    pivot_longer(cols = all_of(rewiring_columns), names_to = "Comparison", values_to = "Class") %>%
    group_by(NetworkMetric, Comparison, Class) %>%
    summarise(n = n(), .groups = "drop") %>%
    filter(!Class %in% c("absent", "not comparable")) %>%
    mutate(
      Comparison = str_remove(Comparison, "^rewiring_"),
      Comparison = str_replace(Comparison, paste0(reference_suffix, "$"), ""),
      Comparison = factor(Comparison, levels = setdiff(network_levels, reference_condition)),
      Class = factor(Class, levels = c("retained", "gained", "lost", "sign switch")),
      NetworkMetric = factor(NetworkMetric, levels = network_metric_levels)
    )

  rewiring_plot <- ggplot(rewiring_summary, aes(x = Class, y = n, fill = Class)) +
    geom_col(width = 0.65, colour = "grey20", linewidth = 0.18) +
    facet_grid(NetworkMetric ~ Comparison, scales = "free_y") +
    scale_fill_manual(values = edge_class_colors, guide = "none") +
    theme_nature_main(base_size = 7.0) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, size = 5.7),
      strip.text = element_text(size = 6.1, face = "bold")
    ) +
    labs(title = paste0("Comparable covariance pairs vs ", reference_condition), x = NULL, y = "Region pairs")

  delta_columns <- grep(paste0("^delta_abs_r_.*", reference_suffix, "$"), names(edge_rewiring), value = TRUE)
  edge_ranked <- edge_rewiring %>%
    filter(if_any(all_of(rewiring_columns), ~ .x %in% c("gained", "lost", "sign switch"))) %>%
    rowwise() %>%
    mutate(
      max_delta = {
        delta_values <- c_across(all_of(delta_columns))
        if (all(is.na(delta_values))) NA_real_ else max(abs(delta_values), na.rm = TRUE)
      }
    ) %>%
    ungroup() %>%
    mutate(max_delta = if_else(is.finite(max_delta), max_delta, NA_real_)) %>%
    group_by(NetworkMetric) %>%
    arrange(desc(max_delta), .by_group = TRUE) %>%
    slice_head(n = 8) %>%
    mutate(edge_rank = row_number()) %>%
    ungroup()

  edge_plot_data <- edge_ranked %>%
    mutate(
      EdgeShort = format_network_edge_label(edge_id),
      EdgeShort = paste0(edge_rank, ". ", str_trunc(EdgeShort, 64)),
      EdgeShort = factor(EdgeShort, levels = rev(unique(EdgeShort))),
      NetworkMetric = factor(NetworkMetric, levels = network_metric_levels)
    ) %>%
    select(NetworkMetric, EdgeShort, all_of(paste0("r_", network_levels))) %>%
    pivot_longer(cols = starts_with("r_"), names_to = "NetworkCondition", values_to = "r") %>%
    mutate(
      NetworkCondition = str_remove(NetworkCondition, "^r_"),
      NetworkCondition = factor(NetworkCondition, levels = network_levels)
    )

  if (nrow(edge_plot_data) > 0) {
    edge_plot <- ggplot(edge_plot_data, aes(x = NetworkCondition, y = EdgeShort, fill = r)) +
      geom_tile(colour = "white", linewidth = 0.22) +
      facet_wrap(~ NetworkMetric, scales = "free_y", ncol = 1) +
      scale_fill_gradient2(low = effect_colors[["negative"]], mid = effect_colors[["neutral"]], high = effect_colors[["positive"]], midpoint = 0,
                           limits = c(-1, 1), oob = scales::squish, name = "Spearman r") +
      theme_nature_main(base_size = 6.8) +
      theme(
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 35, hjust = 1, size = 5.6),
        axis.text.y = element_text(size = 5.4)
      ) +
      labs(title = "Top changing covariance pairs", x = NULL, y = NULL)
  } else {
    edge_plot <- ggplot() +
      annotate("text", x = 0, y = 0, label = "No thresholded edge rewiring detected", size = 2.3) +
      theme_void(base_size = 7.5) +
      labs(title = "Top changing covariance pairs")
  }

  overlap_plot_data <- edge_overlap_summary %>%
    filter(overlap_class %in% c("overlap", "cell-count only", "intensity only"))

  overlap_bar_plot <- ggplot(overlap_plot_data, aes(x = NetworkCondition, y = n_edges, fill = overlap_class)) +
    geom_col(width = 0.68, colour = "grey20", linewidth = 0.15) +
    scale_fill_manual(
      values = c(
        "overlap" = "#2B2B2B",
        "cell-count only" = condition_colors[["VEH_paired"]],
        "intensity only" = condition_colors[["CNO_paired"]]
      ),
      name = NULL,
      drop = FALSE
    ) +
    theme_nature_main(base_size = 7.0) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 5.7)) +
    labs(title = "Thresholded edge overlap by readout", x = NULL, y = "Edges")

  overlap_top_edges <- edge_overlap %>%
    filter(overlap_class == "overlap") %>%
    group_by(NetworkCondition) %>%
    slice_max(order_by = mean_abs_rho, n = 5, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      EdgeShort = str_trunc(EdgeShort, 58),
      EdgeShort = factor(EdgeShort, levels = rev(unique(EdgeShort))),
      NetworkCondition = factor(NetworkCondition, levels = network_levels)
    ) %>%
    select(NetworkCondition, EdgeShort, r__Cell_Count, r__Intensity) %>%
    pivot_longer(cols = starts_with("r__"), names_to = "Readout", values_to = "rho") %>%
    mutate(Readout = recode(Readout, r__Cell_Count = "Cell count", r__Intensity = "Intensity"))

  if (nrow(overlap_top_edges) > 0) {
    overlap_edge_plot <- ggplot(overlap_top_edges, aes(x = Readout, y = EdgeShort, fill = rho)) +
      geom_tile(colour = "white", linewidth = 0.22) +
      facet_wrap(~ NetworkCondition, scales = "free_y", ncol = 2) +
      scale_fill_gradient2(low = effect_colors[["negative"]], mid = effect_colors[["neutral"]], high = effect_colors[["positive"]], midpoint = 0,
                           limits = c(-1, 1), oob = scales::squish, name = "Spearman\nrho") +
      fig4_theme(base_size = 6.4) +
      theme(
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 25, hjust = 1),
        axis.text.y = element_text(size = 5.0)
      ) +
      labs(title = "Shared cell-count/intensity covariance edges", x = NULL, y = NULL)
  } else {
    overlap_edge_plot <- ggplot() +
      annotate("text", x = 0, y = 0, label = "No cell-count/intensity edge overlap", size = 2.3) +
      theme_void(base_size = 7.5) +
      labs(title = "Shared cell-count/intensity covariance edges")
  }

  overlap_plot <- (overlap_bar_plot / overlap_edge_plot) +
    patchwork::plot_layout(heights = c(0.9, 1.8), guides = "collect") &
    theme(legend.position = "bottom")

  safe_main_ggsave(
    file.path(main_fig_dir, paste0("main_figure_network_cell_intensity_edge_overlap", edge_rule_suffix, ".pdf")),
    overlap_plot,
    width = 180 / 25.4,
    height = 150 / 25.4,
    units = "in"
  )
  safe_main_ggsave(
    file.path(main_fig_dir, paste0("main_figure_network_cell_intensity_edge_overlap", edge_rule_suffix, ".svg")),
    overlap_plot,
    width = 180 / 25.4,
    height = 150 / 25.4,
    units = "in",
    device = svglite::svglite
  )
  safe_main_ggsave(
    file.path(main_fig_dir, paste0("main_figure_network_cell_intensity_edge_overlap", edge_rule_suffix, ".png")),
    overlap_plot,
    width = 180 / 25.4,
    height = 150 / 25.4,
    units = "in",
    dpi = 600,
    bg = "white"
  )

  network_plot <- (metric_plot / hub_plot / rewiring_plot / edge_plot) +
    patchwork::plot_layout(heights = c(1.05, 1.65, 1.25, 1.55), guides = "collect") +
    patchwork::plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 9), legend.position = "bottom")

  safe_main_ggsave(
    file.path(main_fig_dir, paste0("main_figure_network_rewiring_separate_cell_count_intensity", edge_rule_suffix, ".pdf")),
    network_plot,
    width = 180 / 25.4,
    height = 245 / 25.4,
    units = "in"
  )
  safe_main_ggsave(
    file.path(main_fig_dir, paste0("main_figure_network_rewiring_separate_cell_count_intensity", edge_rule_suffix, ".svg")),
    network_plot,
    width = 180 / 25.4,
    height = 245 / 25.4,
    units = "in",
    device = svglite::svglite
  )
  safe_main_ggsave(
    file.path(main_fig_dir, paste0("main_figure_network_rewiring_separate_cell_count_intensity", edge_rule_suffix, ".png")),
    network_plot,
    width = 180 / 25.4,
    height = 245 / 25.4,
    units = "in"
  )

  # Deprecated compatibility filenames are written only for the nominal-P
  # display rule to avoid silently overwriting them with stricter/looser rules.
  if (edge_rule == "nominal_p") {
    safe_main_ggsave(
      file.path(main_fig_dir, "main_figure_network_rewiring_four_conditions.pdf"),
      network_plot,
      width = 180 / 25.4,
      height = 245 / 25.4,
      units = "in"
    )
    safe_main_ggsave(
      file.path(main_fig_dir, "main_figure_network_rewiring_four_conditions.svg"),
      network_plot,
      width = 180 / 25.4,
      height = 245 / 25.4,
      units = "in",
      device = svglite::svglite
    )
    safe_main_ggsave(
      file.path(main_fig_dir, "main_figure_network_rewiring_four_conditions.png"),
      network_plot,
      width = 180 / 25.4,
      height = 245 / 25.4,
      units = "in"
    )
  }

  network_plot
}

panel_a <- make_panel_a()

cem_candidate <- find_cem_region_main(long)
panel_b <- NULL
panel_c <- NULL

if (is.null(cem_candidate)) {
  write_main_warning(
    "main_figure_CeM_warning.txt",
    "CeM panels skipped: no Annotation/Abbreviation match for CeM, CeA, CeAm, or central amygdalar nucleus medial."
  )
} else {
  cem_region <- cem_candidate$RegionKey[[1]]
  cem_values <- bind_rows(
    get_cem_values(long, cem_region, "Cell_Count"),
    get_cem_values(long, cem_region, "Intensity")
  )
  cem_stats <- get_cem_stats(cem_region)

  readr::write_csv(cem_values, file.path(main_tab_dir, "main_figure_CeM_values.csv"))
  readr::write_csv(cem_stats, file.path(main_tab_dir, "main_figure_CeM_statistics.csv"))
  readr::write_csv(cem_candidate, file.path(main_tab_dir, "main_figure_CeM_detected_region.csv"))

  panel_b <- make_cem_metric_plot(cem_values, cem_stats, "Cell_Count", "cFos+ cell count, log1p")
  panel_c <- make_cem_metric_plot(cem_values, cem_stats, "Intensity", "Projection intensity, log1p")
}

panel_d <- make_cfos_projection_scatter()
panel_e <- make_learning_stress_heatmap(top_n = 12)
panel_f <- make_main_pca(long)
panel_network <- make_network_rewiring_figure(
  long,
  region_top_n = 14,
  abs_r_cutoff = 0.60,
  edge_rule = "fdr"
)

panel_network_effect_size_only <- make_network_rewiring_figure(
  long,
  region_top_n = 14,
  abs_r_cutoff = 0.60,
  edge_rule = "effect_size_only"
)

panel_network_fdr <- make_network_rewiring_figure(
  long,
  region_top_n = 14,
  abs_r_cutoff = 0.60,
  edge_rule = "fdr"
)

main_panels <- list(panel_a, panel_b, panel_c, panel_d, panel_e, panel_f)
main_panels <- main_panels[!vapply(main_panels, is.null, logical(1))]

if (length(main_panels) >= 1) {
  if (!is.null(panel_f)) {
    panel_f <- panel_f + theme(legend.position = "none")
  }

  main_plot <- (
    panel_a /
      (panel_b | panel_c) /
      (panel_d | panel_f) /
      panel_e
  ) +
    patchwork::plot_layout(heights = c(0.45, 1.0, 1.0, 1.18), guides = "keep") +
    patchwork::plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 9))

  safe_main_ggsave(
    file.path(main_fig_dir, "main_figure_learning_vs_stress_cfos_cell_count_projection.pdf"),
    main_plot,
    width = 180 / 25.4,
    height = 255 / 25.4,
    units = "in"
  )
  safe_main_ggsave(
    file.path(main_fig_dir, "main_figure_learning_vs_stress_cfos_cell_count_projection.svg"),
    main_plot,
    width = 180 / 25.4,
    height = 255 / 25.4,
    units = "in",
    device = svglite::svglite
  )
  safe_main_ggsave(
    file.path(main_fig_dir, "main_figure_learning_vs_stress_cfos_cell_count_projection.png"),
    main_plot,
    width = 180 / 25.4,
    height = 255 / 25.4,
    units = "in",
    dpi = 600
  )

} else {
  write_main_warning("main_figure_warning.txt", "Legacy exploratory main figure was not generated because every panel was skipped.")
}

# -----------------------------
# 15. Manuscript-ready Fig. 4 and supplementary exploratory outputs
# -----------------------------
fig4_fig_dir <- publication_fig4_fig_dir
fig4_tab_dir <- publication_fig4_tab_dir
fig4_supp_fig_dir <- publication_supp_fig_dir
fig4_supp_tab_dir <- publication_supp_tab_dir
dir.create(fig4_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig4_tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig4_supp_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig4_supp_tab_dir, recursive = TRUE, showWarnings = FALSE)

fig4_deprecated_main_bases <- c(
  "fig4_main_manuscript_matched_combined",
  "fig4_main_learning_stress_cfos_cell_count_projection",
  "fig4F_seed_based_cfos_covariation",
  "fig4G_network_rewiring_summary"
)
fig4_deprecated_main_files <- as.vector(outer(fig4_deprecated_main_bases, c("svg", "pdf", "png"), paste, sep = "."))
invisible(file.remove(file.path(fig4_fig_dir, fig4_deprecated_main_files[file.exists(file.path(fig4_fig_dir, fig4_deprecated_main_files))])))

fig4_deprecated_table_bases <- c(
  "fig4BC_effect_map_source",
  "fig4BC_region_condition_availability",
  "fig4D_key_region_condition_profiles_source",
  "fig4D_key_region_condition_profiles_summary",
  "fig4D_key_region_condition_profiles_availability",
  "fig4E_cfos_cell_count_projection_dissociation_source",
  "fig4F_seed_based_cfos_covariation_source",
  "fig4F_seed_based_cfos_covariation_warnings",
  "fig4F_seed_region_matching",
  "fig4G_network_rewiring_edges_source",
  "fig4G_network_rewiring_summary",
  "fig4_region_prior_matching",
  "fig4_region_prior_match_display_decisions"
)
fig4_deprecated_table_files <- as.vector(outer(fig4_deprecated_table_bases, c("csv", "xlsx"), paste, sep = "."))
invisible(file.remove(file.path(fig4_tab_dir, fig4_deprecated_table_files[file.exists(file.path(fig4_tab_dir, fig4_deprecated_table_files))])))

fig4_main_contrasts <- c(
  "Learning_effect",
  "CeM_manipulation_during_learning",
  "CeM_manipulation_during_stress",
  "Learning_x_CeM_interaction"
)
fig4_contrast_labels <- c(
  Learning_effect = "Learning",
  CeM_manipulation_during_learning = "CNO during\nlearning",
  CeM_manipulation_during_stress = "CNO during\nstress",
  Learning_x_CeM_interaction = "CNO effect\nchange"
)
fig4_condition_labels_short <- c(
  VEH_paired = "VEH-L",
  VEH_unpaired = "VEH-S",
  CNO_paired = "CNO-L",
  CNO_unpaired = "CNO-S"
)
fig4_condition_labels <- c(
  VEH_paired = "VEH paired",
  VEH_unpaired = "VEH unpaired",
  CNO_paired = "CNO paired",
  CNO_unpaired = "CNO unpaired"
)
fig4_max_heatmap_regions <- 10
fig4_min_main_non_na_display <- 2
fig4_profile_region_n <- 3
fig4_min_profile_n_per_condition <- 2
fig4_cov_abs_r_cutoff <- 0.70
fig4_cov_fdr_cutoff <- 0.10
fig4_cov_min_n_pair <- 3
fig4_skipped_panels <- character()

fig4_theme <- function(base_size = 7.0) {
  theme_nature_main(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = base_size + 0.7),
      plot.subtitle = element_text(size = base_size - 0.2, colour = "grey35"),
      strip.text = element_text(face = "bold", size = base_size - 0.5)
    )
}

fig4_note_skip <- function(panel, reason) {
  fig4_skipped_panels <<- c(fig4_skipped_panels, paste0(panel, ": ", reason))
  invisible(NULL)
}

save_fig4_plot <- function(plot, filename_base, width, height, dir = publication_fig4_panel_fig_dir, dpi = 600, mirror_root = TRUE) {
  if (is.null(plot)) return(invisible(NULL))
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  for (ext in c("svg", "pdf", "png")) {
    filename <- file.path(dir, paste0(filename_base, ".", ext))
    if (ext == "svg") {
      ggsave(filename, plot, width = width, height = height, units = "in", device = svglite::svglite)
    } else if (ext == "png") {
      ggsave(filename, plot, width = width, height = height, units = "in", dpi = dpi, bg = "white")
    } else {
      ggsave(filename, plot, width = width, height = height, units = "in")
    }
  }
  if (isTRUE(mirror_root) && normalizePath(dir, winslash = "/", mustWork = FALSE) != normalizePath(fig4_fig_dir, winslash = "/", mustWork = FALSE)) {
    dir.create(fig4_fig_dir, recursive = TRUE, showWarnings = FALSE)
    for (ext in c("svg", "pdf", "png")) {
      filename <- file.path(fig4_fig_dir, paste0(filename_base, ".", ext))
      if (ext == "svg") {
        ggsave(filename, plot, width = width, height = height, units = "in", device = svglite::svglite)
      } else if (ext == "png") {
        ggsave(filename, plot, width = width, height = height, units = "in", dpi = dpi, bg = "white")
      } else {
        ggsave(filename, plot, width = width, height = height, units = "in")
      }
    }
  }
  invisible(plot)
}

write_fig4_table <- function(x, filename_base, dir = fig4_tab_dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(x, file.path(dir, paste0(filename_base, ".csv")))
  openxlsx::write.xlsx(x, file.path(dir, paste0(filename_base, ".xlsx")), overwrite = TRUE)
  invisible(x)
}

fig4_short_region <- function(annotation, abbreviation) {
  label <- coalesce(annotation, abbreviation, "Unknown")
  label <- if_else(!is.na(annotation) & annotation != "", annotation, label)
  str_trunc(label, 16)
}

fig4_parent_abbreviation_candidate <- function(abbreviation) {
  abbreviation <- coalesce(abbreviation, "")
  candidate <- str_replace(abbreviation, "[a-z]+$", "")
  if_else(candidate != "" & candidate != abbreviation, candidate, abbreviation)
}

fig4_first_present <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) NA else x[[1]]
}

fig4_region_catalog <- long %>%
  group_by(RegionKey) %>%
  summarise(
    Annotation = fig4_first_present(Annotation),
    Abbreviation = fig4_first_present(Abbreviation),
    Class = fig4_first_present(Class),
    Level = {
      level_values <- Level[!is.na(Level)]
      if (length(level_values) == 0) NA_integer_ else level_values[[1]]
    },
    RegionLabel = fig4_first_present(RegionLabel),
    .groups = "drop"
  ) %>%
  mutate(
    Annotation_lower = str_to_lower(coalesce(Annotation, "")),
    Abbreviation_lower = str_to_lower(coalesce(Abbreviation, "")),
    RegionLabel_lower = str_to_lower(coalesce(RegionLabel, "")),
    RegionShort = fig4_short_region(Annotation, Abbreviation),
    ParentAbbreviationCandidate = fig4_parent_abbreviation_candidate(Abbreviation)
  )

fig4_parent_lookup <- fig4_region_catalog %>%
  distinct(ParentAbbreviationCandidate = Abbreviation, MainRegionAnnotation = Annotation) %>%
  filter(
    !is.na(ParentAbbreviationCandidate),
    ParentAbbreviationCandidate != "",
    ParentAbbreviationCandidate %in% fig4_region_catalog$Abbreviation
  ) %>%
  group_by(ParentAbbreviationCandidate) %>%
  slice(1) %>%
  ungroup()

fig4_region_catalog <- fig4_region_catalog %>%
  left_join(fig4_parent_lookup, by = "ParentAbbreviationCandidate") %>%
  mutate(
    parent_collapsing_note = if_else(
      !is.na(MainRegionAnnotation) & ParentAbbreviationCandidate != Abbreviation,
      "parent label inferred from abbreviation stem present in data",
      "original region label retained; parent inference unavailable or unnecessary"
    ),
    MainRegionAbbreviation = if_else(!is.na(MainRegionAnnotation), ParentAbbreviationCandidate, Abbreviation),
    MainRegion = coalesce(MainRegionAnnotation, Annotation, Abbreviation, "Unknown"),
    MainRegionShort = str_trunc(coalesce(MainRegionAbbreviation, MainRegion), 16),
    SystemGroup = fig4_system_group(Class, Annotation, Abbreviation, RegionLabel),
    SystemGroup = factor(SystemGroup, levels = fig4_system_levels),
    MainRegionGroup = paste(coalesce(Class, "Unknown"), MainRegionShort, sep = " / "),
    RegionDisplay = if_else(
      !is.na(MainRegionAbbreviation) & !is.na(Annotation) & MainRegionAbbreviation != Annotation,
      paste0(MainRegionShort, ": ", RegionShort),
      RegionShort
    )
  )

fig4_region_condition_n <- QC_region_level_availability %>%
  select(
    RegionKey,
    n_VEH_paired, n_VEH_unpaired, n_CNO_paired, n_CNO_unpaired,
    missing_VEH_paired, missing_VEH_unpaired, missing_CNO_paired, missing_CNO_unpaired
  )

fig4_metric_region_condition_n <- region_condition_missingness_qc %>%
  select(
    RegionKey, Condition, n_cell_count, n_intensity,
    missing_cell_count, missing_intensity
  ) %>%
  pivot_longer(
    cols = c(n_cell_count, n_intensity, missing_cell_count, missing_intensity),
    names_to = c(".value", "metric_stub"),
    names_pattern = "^(n|missing)_(cell_count|intensity)$"
  ) %>%
  mutate(
    Metric = recode(metric_stub, cell_count = "Cell_Count", intensity = "Intensity")
  ) %>%
  select(RegionKey, Metric, Condition, n, missing) %>%
  pivot_wider(
    names_from = Condition,
    values_from = c(n, missing),
    values_fill = 0
  )

fig4_add_contrast_qc <- function(data) {
  data %>%
    left_join(fig4_metric_region_condition_n, by = c("RegionKey", "Metric")) %>%
    rowwise() %>%
    mutate(
      min_n_per_relevant_group = case_when(
        contrast == "Learning_effect" ~ min(c(n_VEH_paired, n_VEH_unpaired), na.rm = TRUE),
        contrast == "CeM_manipulation_during_learning" ~ min(c(n_CNO_paired, n_VEH_paired), na.rm = TRUE),
        contrast == "CeM_manipulation_during_stress" ~ min(c(n_CNO_unpaired, n_VEH_unpaired), na.rm = TRUE),
        contrast == "Learning_x_CeM_interaction" ~ min(c(n_VEH_paired, n_VEH_unpaired, n_CNO_paired, n_CNO_unpaired), na.rm = TRUE),
        TRUE ~ NA_real_
      ),
      qc_flag = case_when(
        is.na(min_n_per_relevant_group) | min_n_per_relevant_group < 2 ~ "excluded_insufficient_n",
        min_n_per_relevant_group == 2 ~ "low_n_exploratory",
        min_n_per_relevant_group >= 3 ~ "OK",
        TRUE ~ "excluded_insufficient_n"
      ),
      logFC_main_display = if_else(qc_flag == "excluded_insufficient_n", NA_real_, logFC),
      missingness_relevant = case_when(
        contrast == "Learning_effect" ~ paste0("VEH_paired=", missing_VEH_paired, "; VEH_unpaired=", missing_VEH_unpaired),
        contrast == "CeM_manipulation_during_learning" ~ paste0("CNO_paired=", missing_CNO_paired, "; VEH_paired=", missing_VEH_paired),
        contrast == "CeM_manipulation_during_stress" ~ paste0("CNO_unpaired=", missing_CNO_unpaired, "; VEH_unpaired=", missing_VEH_unpaired),
        contrast == "Learning_x_CeM_interaction" ~ paste0("VEH_paired=", missing_VEH_paired, "; VEH_unpaired=", missing_VEH_unpaired, "; CNO_paired=", missing_CNO_paired, "; CNO_unpaired=", missing_CNO_unpaired),
        TRUE ~ NA_character_
      )
    ) %>%
    ungroup()
}

fig4_sign_key <- tibble::tribble(
  ~contrast, ~short_label, ~contrast_formula, ~positive_logFC_means, ~plain_language,
  "Learning_effect", "Learning", "VEH_paired - VEH_unpaired", "VEH_paired > VEH_unpaired", "How paired learning differs from unpaired stress in VEH animals.",
  "CeM_manipulation_during_learning", "CNO during learning", "CNO_paired - VEH_paired", "CNO_paired > VEH_paired", "How CNO changes the paired-learning condition.",
  "CeM_manipulation_during_stress", "CNO during stress", "CNO_unpaired - VEH_unpaired", "CNO_unpaired > VEH_unpaired", "How CNO changes the unpaired-stress condition.",
  "Learning_x_CeM_interaction", "CNO effect change", "(CNO_paired - VEH_paired) - (CNO_unpaired - VEH_unpaired)", "CNO effect is stronger / more positive during learning than during stress", "Interaction means the CNO-VEH difference is compared across contexts; it asks whether the CNO effect depends on learning versus stress."
)
write_fig4_table(fig4_sign_key, "fig4_contrast_sign_key")

fig4_all_effects <- purrr::imap_dfr(contrast_tables, function(res, metric) {
  res %>%
    filter(contrast %in% fig4_main_contrasts) %>%
    mutate(Metric = metric) %>%
    select(Metric, contrast, Class, Annotation, Abbreviation, Level, RegionLabel, RegionKey,
           logFC, AveExpr, t, P.Value, adj.P.Val.global, n_samples, note)
}) %>%
  fig4_add_contrast_qc()
write_fig4_table(fig4_all_effects, "Fig4_all_region_central_contrast_effects", dir = publication_source_tab_dir)
write_fig4_table(fig4_all_effects, "Fig4_all_region_central_contrast_effects", dir = fig4_supp_tab_dir)

run_fig4_complete_case_contrasts <- function(data, metric) {
  metric_sym <- sym(metric)
  region_meta <- data %>%
    distinct(RegionKey, Class, Annotation, Abbreviation, Level, RegionLabel)

  purrr::map_dfr(fig4_main_contrasts, function(co) {
    needed_conditions <- switch(
      co,
      Learning_effect = c("VEH_paired", "VEH_unpaired"),
      CeM_manipulation_during_learning = c("CNO_paired", "VEH_paired"),
      CeM_manipulation_during_stress = c("CNO_unpaired", "VEH_unpaired"),
      Learning_x_CeM_interaction = c("VEH_paired", "VEH_unpaired", "CNO_paired", "CNO_unpaired")
    )

    data %>%
      filter(Condition %in% needed_conditions) %>%
      select(SampleID, Condition, RegionKey, RawValue = !!metric_sym) %>%
      mutate(Value = safe_log1p(RawValue)) %>%
      group_by(RegionKey) %>%
      group_modify(function(d, key) {
        d <- d %>% filter(!is.na(Value), !is.na(Condition))
        n_by_condition <- d %>%
          count(Condition, name = "n") %>%
          complete(Condition = needed_conditions, fill = list(n = 0))
        min_n <- min(n_by_condition$n)
        region_sd <- stats::sd(d$Value, na.rm = TRUE)
        if (min_n < 2 || !is.finite(region_sd) || region_sd == 0) {
          return(tibble(
            contrast = co,
            complete_case_logFC = NA_real_,
            complete_case_P.Value = NA_real_,
            complete_case_n = nrow(d),
            complete_case_min_n_per_group = min_n,
            complete_case_note = "not estimable: fewer than 2 complete observations per relevant group or zero variance"
          ))
        }

        d <- d %>%
          mutate(
            Condition = factor(as.character(Condition), levels = nature_condition_levels),
            Paired = if_else(Condition %in% c("VEH_paired", "CNO_paired"), "paired", "unpaired"),
            Drug = if_else(Condition %in% c("CNO_paired", "CNO_unpaired"), "CNO", "VEH")
          )

        if (co == "Learning_x_CeM_interaction") {
          fit <- tryCatch(stats::lm(Value ~ Paired * Drug, data = d), error = function(e) NULL)
          if (is.null(fit)) {
            estimate <- NA_real_
            p_value <- NA_real_
          } else {
            means <- d %>% group_by(Condition) %>% summarise(m = mean(Value), .groups = "drop") %>% deframe()
            estimate <- (means[["CNO_paired"]] - means[["VEH_paired"]]) -
              (means[["CNO_unpaired"]] - means[["VEH_unpaired"]])
            coef_tab <- summary(fit)$coefficients
            interaction_row <- grep("Paired.*Drug|Drug.*Paired", rownames(coef_tab), value = TRUE)
            p_value <- if (length(interaction_row) > 0) coef_tab[interaction_row[[1]], "Pr(>|t|)"] else NA_real_
          }
        } else {
          first_condition <- needed_conditions[[1]]
          second_condition <- needed_conditions[[2]]
          estimate <- mean(d$Value[d$Condition == first_condition], na.rm = TRUE) -
            mean(d$Value[d$Condition == second_condition], na.rm = TRUE)
          p_value <- suppressWarnings(
            tryCatch(stats::t.test(Value ~ Condition, data = d)$p.value, error = function(e) NA_real_)
          )
        }

        tibble(
          contrast = co,
          complete_case_logFC = unname(estimate),
          complete_case_P.Value = p_value,
          complete_case_n = nrow(d),
          complete_case_min_n_per_group = min_n,
          complete_case_note = "complete-case raw log1p contrast"
        )
      }) %>%
      ungroup()
  }) %>%
    mutate(Metric = metric, .before = 1) %>%
    left_join(region_meta, by = "RegionKey") %>%
    group_by(Metric, contrast) %>%
    mutate(complete_case_adj.P.Val.global = p.adjust(complete_case_P.Value, method = "BH")) %>%
    ungroup()
}

fig4_complete_case_effects <- purrr::map_dfr(metrics_to_analyse, ~ run_fig4_complete_case_contrasts(long, .x)) %>%
  mutate(
    logFC = complete_case_logFC,
    P.Value = complete_case_P.Value,
    adj.P.Val.global = complete_case_adj.P.Val.global
  ) %>%
  fig4_add_contrast_qc()
write_fig4_table(fig4_complete_case_effects, "Fig4_complete_case_central_contrast_effects", dir = exploratory_sensitivity_dir)

fig4_stability_inputs <- bind_rows(
  purrr::imap_dfr(contrast_tables, function(res, metric) {
    res %>%
      filter(contrast %in% fig4_main_contrasts) %>%
      transmute(
        analysis_variant = "raw_log1p_limma_median_imputed",
        Metric = metric,
        contrast,
        RegionKey,
        logFC,
        P.Value,
        adj.P.Val.global,
        normalization = "none",
        imputation = "region median before limma"
      )
  }),
  purrr::imap_dfr(normalized_metric_lookup, function(raw_metric, norm_metric) {
    contrast_tables_normalized[[norm_metric]] %>%
      filter(contrast %in% fig4_main_contrasts) %>%
      transmute(
        analysis_variant = "normalized_log1p_limma_median_imputed",
        Metric = raw_metric,
        contrast,
        RegionKey,
        logFC,
        P.Value,
        adj.P.Val.global,
        normalization = "animal-level normalized metric",
        imputation = "region median before limma"
      )
  }),
  fig4_complete_case_effects %>%
    transmute(
      analysis_variant = "raw_log1p_complete_case",
      Metric,
      contrast,
      RegionKey,
      logFC = complete_case_logFC,
      P.Value = complete_case_P.Value,
      adj.P.Val.global = complete_case_adj.P.Val.global,
      normalization = "none",
      imputation = "none; complete cases only"
    )
) %>%
  mutate(
    logFC = suppressWarnings(as.numeric(logFC)),
    P.Value = suppressWarnings(as.numeric(P.Value)),
    adj.P.Val.global = suppressWarnings(as.numeric(adj.P.Val.global))
  ) %>%
  filter(!is.na(logFC)) %>%
  group_by(analysis_variant, Metric, contrast) %>%
  arrange(desc(abs(logFC)), P.Value, .by_group = TRUE) %>%
  mutate(abs_effect_rank = row_number()) %>%
  ungroup() %>%
  group_by(analysis_variant, Metric, contrast, RegionKey) %>%
  summarise(
    logFC = first(logFC),
    P.Value = first(P.Value),
    adj.P.Val.global = first(adj.P.Val.global),
    abs_effect_rank = first(abs_effect_rank),
    normalization = first(normalization),
    imputation = first(imputation),
    .groups = "drop"
  )

fig4_stability_comparison <- fig4_stability_inputs %>%
  select(analysis_variant, Metric, contrast, RegionKey, logFC, P.Value, adj.P.Val.global, abs_effect_rank, normalization, imputation) %>%
  pivot_wider(
    names_from = analysis_variant,
    values_from = c(logFC, P.Value, adj.P.Val.global, abs_effect_rank, normalization, imputation),
    names_sep = "__"
  ) %>%
  {
    d <- .
    expected_numeric_cols <- c(
      "logFC__raw_log1p_limma_median_imputed",
      "logFC__normalized_log1p_limma_median_imputed",
      "logFC__raw_log1p_complete_case",
      "P.Value__raw_log1p_limma_median_imputed",
      "P.Value__normalized_log1p_limma_median_imputed",
      "P.Value__raw_log1p_complete_case",
      "adj.P.Val.global__raw_log1p_limma_median_imputed",
      "adj.P.Val.global__normalized_log1p_limma_median_imputed",
      "adj.P.Val.global__raw_log1p_complete_case",
      "abs_effect_rank__raw_log1p_limma_median_imputed",
      "abs_effect_rank__normalized_log1p_limma_median_imputed",
      "abs_effect_rank__raw_log1p_complete_case"
    )
    for (col in setdiff(expected_numeric_cols, names(d))) d[[col]] <- NA_real_
    d
  } %>%
  mutate(
    across(starts_with("logFC__"), ~ suppressWarnings(as.numeric(.x))),
    across(starts_with("P.Value__"), ~ suppressWarnings(as.numeric(.x))),
    across(starts_with("adj.P.Val.global__"), ~ suppressWarnings(as.numeric(.x))),
    across(starts_with("abs_effect_rank__"), ~ suppressWarnings(as.numeric(.x)))
  ) %>%
  mutate(
    direction_stable_raw_vs_normalized = sign(logFC__raw_log1p_limma_median_imputed) ==
      sign(logFC__normalized_log1p_limma_median_imputed),
    direction_stable_raw_vs_complete_case = sign(logFC__raw_log1p_limma_median_imputed) ==
      sign(logFC__raw_log1p_complete_case),
    rank_shift_raw_vs_normalized = abs_effect_rank__normalized_log1p_limma_median_imputed -
      abs_effect_rank__raw_log1p_limma_median_imputed,
    rank_shift_raw_vs_complete_case = abs_effect_rank__raw_log1p_complete_case -
      abs_effect_rank__raw_log1p_limma_median_imputed
  ) %>%
  left_join(fig4_region_catalog %>% select(RegionKey, Class, Annotation, Abbreviation, RegionLabel), by = "RegionKey") %>%
  arrange(Metric, contrast, abs_effect_rank__raw_log1p_limma_median_imputed)
write_fig4_table(fig4_stability_comparison, "Sensitivity_stability_raw_normalized_complete_case_imputed", dir = exploratory_sensitivity_dir)

fig4_effect_rank <- fig4_all_effects %>%
  group_by(Metric, RegionKey, Class, Annotation, Abbreviation, Level, RegionLabel) %>%
  summarise(
    n_estimable_contrasts_total = sum(qc_flag != "excluded_insufficient_n" & !is.na(logFC_main_display)),
    max_abs_effect = if (any(!is.na(logFC_main_display))) max(abs(logFC_main_display), na.rm = TRUE) else NA_real_,
    min_raw_p = if (any(!is.na(P.Value))) min(P.Value, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(max_abs_effect = coalesce(max_abs_effect, 0), min_raw_p = coalesce(min_raw_p, Inf))

fig4_effect_rank_region_summary <- fig4_effect_rank %>%
  group_by(RegionKey, Class, Annotation, Abbreviation, Level, RegionLabel) %>%
  summarise(
    n_estimable_contrasts_total = max(n_estimable_contrasts_total, na.rm = TRUE),
    max_abs_effect = max(max_abs_effect, na.rm = TRUE),
    min_raw_p = min(min_raw_p, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    n_estimable_contrasts_total = if_else(is.finite(n_estimable_contrasts_total), n_estimable_contrasts_total, 0),
    max_abs_effect = if_else(is.finite(max_abs_effect), max_abs_effect, 0),
    min_raw_p = if_else(is.finite(min_raw_p), min_raw_p, Inf)
  )

fig4_metric_region_selection <- purrr::map_dfr(metrics_to_analyse, function(metric_name) {
  metric_rank <- fig4_effect_rank %>%
    filter(Metric == metric_name) %>%
    filter(
      n_estimable_contrasts_total >= fig4_min_main_non_na_display
    ) %>%
    arrange(desc(max_abs_effect), min_raw_p, RegionKey) %>%
    slice_head(n = fig4_max_heatmap_regions)

  metric_rank %>%
    mutate(
      Metric = metric_name,
      selection_basis = "data-driven effect-size ranking after n/estimability filters",
      inclusion_reason = "kept: >=2 displayed central effects after n filter; ranked by max absolute displayed logFC"
    )
})

fig4_region_selection <- fig4_metric_region_selection %>%
  group_by(RegionKey) %>%
  summarise(
    metric_specific_inclusion = paste(sort(unique(Metric)), collapse = "; "),
    n_estimable_contrasts_total = max(n_estimable_contrasts_total, na.rm = TRUE),
    max_abs_effect = max(max_abs_effect, na.rm = TRUE),
    min_raw_p = min(min_raw_p, na.rm = TRUE),
    inclusion_reason = paste(sort(unique(inclusion_reason)), collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(
    n_estimable_contrasts_total = if_else(is.finite(n_estimable_contrasts_total), n_estimable_contrasts_total, 0),
    max_abs_effect = if_else(is.finite(max_abs_effect), max_abs_effect, 0),
    min_raw_p = if_else(is.finite(min_raw_p), min_raw_p, Inf)
) %>%
  distinct(RegionKey, .keep_all = TRUE) %>%
  left_join(
    fig4_region_catalog %>%
      select(
        RegionKey, Annotation, Abbreviation, Class, Level, RegionLabel,
        RegionShort, SystemGroup, MainRegion, MainRegionAbbreviation, MainRegionShort,
        MainRegionGroup, RegionDisplay, parent_collapsing_note
      ),
    by = "RegionKey"
  ) %>%
  filter(n_estimable_contrasts_total > 0) %>%
  mutate(
    main_display_filter = "kept: >=2 displayed central effects; top ranked by max absolute displayed logFC",
    selection_basis = "data-driven effect-size ranking after n/estimability filters",
    RegionShort = make.unique(RegionShort, sep = " "),
    RegionDisplay = make.unique(coalesce(RegionDisplay, RegionShort), sep = " "),
    Class = coalesce(Class, "Unknown")
  ) %>%
  arrange(SystemGroup, MainRegionShort, desc(max_abs_effect), RegionShort) %>%
  mutate(plot_order = row_number())

fig4_region_ontology_qc <- fig4_region_catalog %>%
  select(
    Annotation, Abbreviation, Class, SystemGroup, MainRegion, MainRegionAbbreviation,
    MainRegionShort, MainRegionGroup, parent_collapsing_note, Level, RegionKey, RegionLabel
  ) %>%
  left_join(fig4_region_selection %>% select(RegionKey, inclusion_reason), by = "RegionKey") %>%
  left_join(
    region_condition_missingness_qc %>%
      group_by(RegionKey) %>%
      summarise(
        n_CNO_paired_cell = n_cell_count[Condition == "CNO_paired"][1],
        n_VEH_paired_cell = n_cell_count[Condition == "VEH_paired"][1],
        n_CNO_unpaired_cell = n_cell_count[Condition == "CNO_unpaired"][1],
        n_VEH_unpaired_cell = n_cell_count[Condition == "VEH_unpaired"][1],
        n_CNO_paired_intensity = n_intensity[Condition == "CNO_paired"][1],
        n_VEH_paired_intensity = n_intensity[Condition == "VEH_paired"][1],
        n_CNO_unpaired_intensity = n_intensity[Condition == "CNO_unpaired"][1],
        n_VEH_unpaired_intensity = n_intensity[Condition == "VEH_unpaired"][1],
        mean_missing_cell_count = mean(pct_missing_cell_count, na.rm = TRUE),
        mean_missing_intensity = mean(pct_missing_intensity, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "RegionKey"
  ) %>%
  mutate(
    inclusion_reason = coalesce(inclusion_reason, "not included in main Fig. 4")
  )
write_fig4_table(fig4_region_ontology_qc, "fig4_region_ontology_qc")
write_fig4_table(fig4_region_ontology_qc, "Fig4_region_ontology_qc", dir = qc_region_selection_dir)
write_fig4_table(normalization_comparison, "fig4_normalization_comparison")
write_fig4_table(normalization_comparison, "Normalization_comparison_raw_vs_normalized_log1p", dir = qc_normalization_dir)

fig4_region_selection_audit <- fig4_effect_rank %>%
  left_join(
    fig4_metric_region_selection %>%
      select(Metric, RegionKey, selected_in_metric_heatmap = inclusion_reason),
    by = c("Metric", "RegionKey")
  ) %>%
  mutate(
    selected_in_metric_heatmap = !is.na(selected_in_metric_heatmap),
    inclusion_reason = case_when(
      selected_in_metric_heatmap ~ "included: >=2 displayed central effects; selected among top max absolute displayed logFC",
      n_estimable_contrasts_total < fig4_min_main_non_na_display ~ "excluded: fewer than 2 displayed central effects after n filter",
      TRUE ~ "excluded: outside top ranked display limit"
    ),
    selection_rule = paste0(
      "Main heatmap inclusion is data-driven: min relevant group n >= 2 for displayed contrasts, at least ",
      fig4_min_main_non_na_display,
      " displayed central effects per metric-region, then top ",
      fig4_max_heatmap_regions,
      " by max absolute displayed logFC. Raw P is retained for audit/exploratory ranking, not as an inclusion rule."
    )
  ) %>%
  left_join(
    fig4_region_catalog %>% select(RegionKey, SystemGroup, MainRegion, MainRegionAbbreviation, RegionDisplay, parent_collapsing_note),
    by = "RegionKey"
  ) %>%
  arrange(Metric, desc(selected_in_metric_heatmap), desc(max_abs_effect), min_raw_p)
write_fig4_table(fig4_region_selection_audit, "Fig4_region_selection_inclusion_exclusion", dir = qc_region_selection_dir)
write_fig4_table(fig4_region_selection_audit, "Fig4_region_selection_inclusion_exclusion", dir = publication_source_tab_dir)

fig4_region_label_map <- fig4_region_selection %>%
  select(
    RegionDisplay, RegionShort, Annotation, Abbreviation, Class,
    SystemGroup, MainRegion, MainRegionAbbreviation, MainRegionShort, MainRegionGroup,
    Level, RegionLabel, RegionKey, inclusion_reason, selection_basis,
    metric_specific_inclusion, main_display_filter, parent_collapsing_note,
    n_estimable_contrasts_total, max_abs_effect, min_raw_p
  )
write_fig4_table(fig4_region_label_map, "fig4_region_abbreviation_key")

fig4_heatmap_source <- fig4_all_effects %>%
  inner_join(
    fig4_metric_region_selection %>%
      select(
        Metric, RegionKey, selection_basis, inclusion_reason,
        metric_n_estimable_contrasts = n_estimable_contrasts_total,
        metric_max_abs_effect = max_abs_effect,
        metric_min_raw_p = min_raw_p
      ),
    by = c("Metric", "RegionKey")
  ) %>%
  left_join(
    fig4_region_selection %>%
      select(
        RegionKey, RegionShort, RegionDisplay, MainRegion, MainRegionAbbreviation,
        SystemGroup, MainRegionShort, MainRegionGroup, plot_order
      ),
    by = "RegionKey"
  ) %>%
  mutate(
    contrast = factor(contrast, levels = fig4_main_contrasts),
    contrast_estimable = !is.na(logFC),
    low_n_exploratory = qc_flag == "low_n_exploratory",
    not_estimable_reason = case_when(
      !contrast_estimable ~ "grey: contrast not estimable from available non-missing values",
      is.na(logFC_main_display) ~ "grey: insufficient n, missing, non-estimable, or excluded from main display",
      TRUE ~ NA_character_
    )
  ) %>%
  select(
    Metric, contrast, RegionKey, Annotation, Abbreviation, RegionDisplay, SystemGroup,
    logFC, logFC_main_display, P.Value, adj.P.Val.global, min_n_per_relevant_group,
    qc_flag, missingness_relevant, selection_basis, inclusion_reason,
    low_n_exploratory, not_estimable_reason, contrast_estimable,
    MainRegion, MainRegionAbbreviation, MainRegionShort, MainRegionGroup,
    metric_n_estimable_contrasts, metric_max_abs_effect, metric_min_raw_p,
    Class, Level, RegionLabel, RegionShort, plot_order
  )

write_fig4_table(fig4_heatmap_source %>% filter(Metric == "Intensity"), "Fig4B_projection_effect_sizes", dir = publication_source_tab_dir)
write_fig4_table(fig4_heatmap_source %>% filter(Metric == "Cell_Count"), "Fig4C_cfos_cell_count_effect_sizes", dir = publication_source_tab_dir)
write_fig4_table(fig4_heatmap_source %>% filter(Metric == "Intensity"), "Fig4B_projection_effect_sizes")
write_fig4_table(fig4_heatmap_source %>% filter(Metric == "Cell_Count"), "Fig4C_cfos_cell_count_effect_sizes")

fig4_heatmap_limit <- max(abs(fig4_heatmap_source$logFC_main_display), na.rm = TRUE)
fig4_heatmap_limit <- ifelse(is.finite(fig4_heatmap_limit) && fig4_heatmap_limit > 0, fig4_heatmap_limit, 1)

make_fig4_heatmap <- function(metric, title_label) {
  metric_region_order <- fig4_heatmap_source %>%
    filter(Metric == metric) %>%
    group_by(RegionKey, RegionDisplay, SystemGroup, MainRegionShort, metric_max_abs_effect, plot_order) %>%
    summarise(n_display_values = sum(!is.na(logFC_main_display)), .groups = "drop") %>%
    filter(
      n_display_values >= fig4_min_main_non_na_display
    ) %>%
    arrange(desc(metric_max_abs_effect), plot_order, RegionDisplay) %>%
    slice_head(n = fig4_max_heatmap_regions) %>%
    arrange(SystemGroup, MainRegionShort, desc(metric_max_abs_effect), RegionDisplay)

  d <- fig4_heatmap_source %>%
    filter(Metric == metric, RegionKey %in% metric_region_order$RegionKey) %>%
    mutate(
      x = as.integer(contrast),
      y = factor(RegionDisplay, levels = rev(metric_region_order$RegionDisplay))
    )
  if (nrow(d) == 0) return(NULL)

  ggplot(d, aes(x = x, y = y)) +
    geom_tile(aes(fill = logFC_main_display), colour = "white", linewidth = 0.18, width = 0.92, height = 0.92) +
    geom_point(
      data = d %>% filter(low_n_exploratory, !is.na(logFC_main_display)),
      aes(x = x, y = y),
      inherit.aes = FALSE,
      shape = 21, size = 0.8, stroke = 0.18, colour = "grey20", fill = NA
    ) +
    facet_grid(SystemGroup ~ ., scales = "free_y", space = "free_y") +
    scale_x_continuous(
      breaks = seq_along(fig4_main_contrasts),
      labels = fig4_contrast_labels,
      limits = c(0.1, length(fig4_main_contrasts) + 0.55),
      expand = expansion(mult = c(0, 0.02))
    ) +
    scale_fill_gradient2(
      low = effect_colors[["negative"]], mid = effect_colors[["neutral"]], high = effect_colors[["positive"]],
      midpoint = 0, limits = c(-fig4_heatmap_limit, fig4_heatmap_limit),
      oob = scales::squish, na.value = "grey88", name = "Effect size\n(logFC)"
    ) +
    fig4_theme(base_size = 6.4) +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1, size = 5.6),
      axis.text.y = element_text(size = 6.2),
      strip.text.y = element_text(angle = 0, size = 5.3),
      legend.position = "right"
    ) +
    labs(
      title = title_label,
      subtitle = "Colour shows logFC effect size. Open circle = low-n exploratory. Grey = missing or not estimable.",
      x = NULL,
      y = NULL
    )
}

fig4A <- {
  condition_d <- tibble(
    Condition = factor(names(fig4_condition_labels), levels = nature_condition_levels),
    Label = fig4_condition_labels[names(fig4_condition_labels)],
    Meaning = c(
      "Associative fear learning",
      "Stress / non-associative exposure",
      "CeM manipulation during learning",
      "CeM manipulation during stress"
    )
  ) %>%
    mutate(y = rev(row_number()))

  contrast_d <- fig4_sign_key %>%
    mutate(
      label = paste0(short_label, ": positive logFC = ", positive_logFC_means, ". ", plain_language),
      y = rev(row_number())
    )

  ggplot() +
    geom_point(data = condition_d, aes(x = 0.03, y = y, fill = Condition), shape = 21, size = 2.6, colour = "grey20", stroke = 0.2) +
    geom_text(data = condition_d, aes(x = 0.08, y = y, label = Label), hjust = 0, size = 2.2, fontface = "bold") +
    geom_text(data = condition_d, aes(x = 0.29, y = y, label = Meaning), hjust = 0, size = 2.05) +
    geom_text(data = contrast_d, aes(x = 0.03, y = y - 4.7, label = str_wrap(label, 78)), hjust = 0, size = 1.9, lineheight = 0.9) +
    scale_fill_manual(values = nature_condition_colors, guide = "none", drop = FALSE) +
    coord_cartesian(xlim = c(0, 1), ylim = c(-4.2, 4.45), clip = "off") +
    theme_void(base_size = 7) +
    labs(
      title = "Conditions and effect directions",
      subtitle = "Interaction means the CNO-VEH difference is tested for whether it changes between learning and stress."
    )
}

fig4B <- make_fig4_heatmap("Intensity", "Projection intensity")
fig4C <- make_fig4_heatmap("Cell_Count", "cFos+ cell count")

profile_candidates <- fig4_region_selection %>%
  arrange(desc(max_abs_effect), min_raw_p)

profile_qc <- long %>%
  filter(RegionKey %in% profile_candidates$RegionKey) %>%
  select(SampleID, Condition, RegionKey, Cell_Count, Intensity) %>%
  pivot_longer(cols = c(Cell_Count, Intensity), names_to = "Metric", values_to = "RawValue") %>%
  group_by(RegionKey, Metric, Condition) %>%
  summarise(n_non_missing = n_distinct(SampleID[!is.na(RawValue)]), .groups = "drop") %>%
  group_by(RegionKey) %>%
  summarise(
    min_cell_count_n = min(n_non_missing[Metric == "Cell_Count"], na.rm = TRUE),
    min_intensity_n = min(n_non_missing[Metric == "Intensity"], na.rm = TRUE),
    profile_pass_min_n = min(min_cell_count_n, min_intensity_n, na.rm = TRUE) >= fig4_min_profile_n_per_condition,
    .groups = "drop"
  ) %>%
  mutate(across(starts_with("min_"), ~ if_else(is.finite(.x), .x, NA_real_))) %>%
  left_join(profile_candidates, by = "RegionKey")

fig4_profile_regions <- profile_qc %>%
  filter(profile_pass_min_n) %>%
  arrange(desc(max_abs_effect), min_raw_p) %>%
  slice_head(n = fig4_profile_region_n)
write_fig4_table(profile_qc, "fig4_panelD_key_region_profiles_qc")
write_fig4_table(profile_qc, "Fig4D_key_region_profiles_qc", dir = publication_source_tab_dir)

fig4D_source <- long %>%
  filter(RegionKey %in% fig4_profile_regions$RegionKey) %>%
  select(SampleID, Animal, Group, Condition, Class, Annotation, Abbreviation, RegionLabel, RegionKey, Cell_Count, Intensity) %>%
  pivot_longer(cols = c(Cell_Count, Intensity), names_to = "Metric", values_to = "RawValue") %>%
  mutate(
    Value = safe_log1p(RawValue),
    MetricLabel = recode(Metric, Cell_Count = "cFos+ cell count", Intensity = "Projection intensity"),
    Condition = factor(as.character(Condition), levels = nature_condition_levels),
    ConditionLabel = factor(fig4_condition_labels_short[as.character(Condition)], levels = fig4_condition_labels_short),
    RegionShort = fig4_short_region(Annotation, Abbreviation)
  ) %>%
  left_join(
    fig4_region_selection %>%
      select(RegionKey, MainRegion, MainRegionAbbreviation, MainRegionShort, MainRegionGroup, RegionDisplay),
    by = "RegionKey"
  ) %>%
  left_join(
    profile_qc %>%
      select(RegionKey, min_cell_count_n, min_intensity_n, profile_pass_min_n) %>%
      mutate(profile_qc_flag = if_else(profile_pass_min_n, "OK", "excluded_insufficient_n")),
    by = "RegionKey"
  ) %>%
  left_join(
    fig4_region_condition_n,
    by = "RegionKey"
  ) %>%
  mutate(
    min_n_per_relevant_group = if_else(Metric == "Cell_Count", min_cell_count_n, min_intensity_n),
    qc_flag = profile_qc_flag
  )

fig4D_summary <- fig4D_source %>%
  group_by(Metric, MetricLabel, RegionKey, RegionDisplay, RegionShort, Condition, ConditionLabel) %>%
  summarise(
    mean_value = mean(Value, na.rm = TRUE),
    n = sum(!is.na(Value)),
    ci95 = qt(0.975, pmax(n - 1, 1)) * sd(Value, na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(ci95 = if_else(is.finite(ci95), ci95, 0), ymin = mean_value - ci95, ymax = mean_value + ci95)
write_fig4_table(fig4D_source, "fig4_panelD_key_region_profiles")
write_fig4_table(fig4D_summary, "fig4_panelD_key_region_profiles_summary")
write_fig4_table(fig4D_source, "Fig4D_key_region_profiles", dir = publication_source_tab_dir)
write_fig4_table(fig4D_summary, "Fig4D_key_region_profiles_summary", dir = publication_source_tab_dir)

fig4D <- if (nrow(fig4D_source) > 0) {
  profile_region_levels <- fig4_profile_regions %>%
    arrange(SystemGroup, MainRegionShort, desc(max_abs_effect), RegionShort) %>%
    pull(RegionDisplay)

  ggplot(fig4D_source %>% mutate(RegionDisplay = factor(RegionDisplay, levels = profile_region_levels)), aes(x = ConditionLabel, y = Value)) +
    geom_point(aes(fill = Condition), shape = 21, size = 1.15, colour = "grey20", stroke = 0.15,
               alpha = 0.8, position = position_jitter(width = 0.08, height = 0, seed = 1)) +
    geom_errorbar(data = fig4D_summary, aes(x = ConditionLabel, ymin = ymin, ymax = ymax, colour = Condition),
                  inherit.aes = FALSE, width = 0.08, linewidth = 0.25) +
    geom_point(data = fig4D_summary, aes(x = ConditionLabel, y = mean_value, fill = Condition),
               inherit.aes = FALSE, shape = 23, size = 1.6, colour = "white", stroke = 0.18) +
    facet_grid(MetricLabel ~ RegionDisplay, scales = "free_y") +
    scale_fill_manual(values = nature_condition_colors, drop = FALSE) +
    scale_colour_manual(values = nature_condition_colors, guide = "none", drop = FALSE) +
    fig4_theme(base_size = 6.0) +
    theme(axis.text.x = element_text(size = 5.2), legend.position = "none") +
    labs(title = "Selected regional profiles", subtitle = "L = paired learning; S = unpaired stress", x = NULL, y = "log1p raw signal")
} else {
  fig4_note_skip("Panel D", "No profile regions passed the minimum-n filter.")
  NULL
}

fig4E_source <- full_join(
  contrast_tables$Cell_Count %>%
    filter(contrast == "Learning_effect") %>%
    select(Class, Annotation, Abbreviation, RegionLabel, RegionKey,
           Cell_Count_logFC = logFC, Cell_Count_P.Value = P.Value, Cell_Count_adj.P.Val.global = adj.P.Val.global),
  contrast_tables$Intensity %>%
    filter(contrast == "Learning_effect") %>%
    select(RegionKey, Intensity_logFC = logFC, Intensity_P.Value = P.Value, Intensity_adj.P.Val.global = adj.P.Val.global),
  by = "RegionKey"
) %>%
  left_join(fig4_region_selection %>% select(RegionKey, selection_basis), by = "RegionKey") %>%
  left_join(
    fig4_metric_region_condition_n %>%
      filter(Metric == "Cell_Count") %>%
      transmute(
        RegionKey,
        Cell_Count_min_learning_n = pmin(n_VEH_paired, n_VEH_unpaired),
        Cell_Count_missingness_learning = paste0("VEH_paired=", missing_VEH_paired, "; VEH_unpaired=", missing_VEH_unpaired)
      ),
    by = "RegionKey"
  ) %>%
  left_join(
    fig4_metric_region_condition_n %>%
      filter(Metric == "Intensity") %>%
      transmute(
        RegionKey,
        Intensity_min_learning_n = pmin(n_VEH_paired, n_VEH_unpaired),
        Intensity_missingness_learning = paste0("VEH_paired=", missing_VEH_paired, "; VEH_unpaired=", missing_VEH_unpaired)
      ),
    by = "RegionKey"
  ) %>%
  mutate(
    selection_basis = coalesce(selection_basis, "not selected for Fig. 4 heatmap"),
    min_n_per_relevant_group = pmin(Cell_Count_min_learning_n, Intensity_min_learning_n, na.rm = TRUE),
    min_n_per_relevant_group = if_else(is.finite(min_n_per_relevant_group), min_n_per_relevant_group, NA_real_),
    qc_flag = case_when(
      is.na(min_n_per_relevant_group) | min_n_per_relevant_group < 2 ~ "excluded_insufficient_n",
      min_n_per_relevant_group == 2 ~ "low_n_exploratory",
      min_n_per_relevant_group >= 3 ~ "OK",
      TRUE ~ "excluded_insufficient_n"
    ),
    plotted = !is.na(Cell_Count_logFC) & !is.na(Intensity_logFC),
    not_plotted_reason = case_when(
      plotted ~ NA_character_,
      is.na(Cell_Count_logFC) & is.na(Intensity_logFC) ~ "missing both metrics",
      is.na(Cell_Count_logFC) ~ "missing regional cFos recruitment effect",
      is.na(Intensity_logFC) ~ "missing projection intensity effect",
      TRUE ~ NA_character_
    ),
    combined_abs_effect = sqrt(Cell_Count_logFC^2 + Intensity_logFC^2),
    RegionShort = fig4_short_region(Annotation, Abbreviation),
    label = if_else(plotted & qc_flag != "excluded_insufficient_n" & rank(-combined_abs_effect, ties.method = "first", na.last = "keep") <= 8, RegionShort, NA_character_),
    quadrant = case_when(
      !plotted | abs(Cell_Count_logFC) < 0.05 | abs(Intensity_logFC) < 0.05 ~ "mixed or near-zero",
      Cell_Count_logFC > 0 & Intensity_logFC > 0 ~ "coupled increase",
      Cell_Count_logFC > 0 & Intensity_logFC < 0 ~ "cFos-only increase",
      Cell_Count_logFC < 0 & Intensity_logFC > 0 ~ "projection-only increase",
      Cell_Count_logFC < 0 & Intensity_logFC < 0 ~ "coupled decrease / reduced signal",
      TRUE ~ NA_character_
    )
  )
write_fig4_table(fig4E_source, "fig4_panelE_cfos_cell_count_projection_scatter")
write_fig4_table(fig4E_source, "Fig4E_cfos_cell_count_projection_scatter", dir = publication_source_tab_dir)

fig4E_qc <- fig4E_source %>%
  summarise(
    n_regions_total = n(),
    n_regions_with_both_metrics = sum(plotted, na.rm = TRUE),
    n_missing_cell_count_effect = sum(is.na(Cell_Count_logFC)),
    n_missing_intensity_effect = sum(is.na(Intensity_logFC))
  )
write_fig4_table(fig4E_qc, "fig4_panelE_cfos_cell_count_projection_matching_qc")
write_fig4_table(fig4E_qc, "Fig4E_cfos_cell_count_projection_matching_qc", dir = publication_source_tab_dir)

fig4E <- if (sum(fig4E_source$plotted, na.rm = TRUE) >= 3) {
  ggplot(fig4E_source %>% filter(plotted, qc_flag != "excluded_insufficient_n"), aes(x = Cell_Count_logFC, y = Intensity_logFC)) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey55") +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey55") +
    geom_point(aes(size = combined_abs_effect), shape = 21, fill = "grey72", colour = "grey15", stroke = 0.18, alpha = 0.85) +
    ggrepel::geom_text_repel(data = fig4E_source %>% filter(!is.na(label)), aes(label = label),
                             size = 1.9, segment.size = 0.16, min.segment.length = 0, max.overlaps = 18, force = 1.8, seed = 4) +
    scale_size_continuous(range = c(0.7, 3.2), guide = "none") +
    fig4_theme(base_size = 6.8) +
    theme(legend.position = "none") +
    labs(title = "cFos cell-count / projection dissociation", x = "Learning logFC, cFos+ cell count", y = "Learning logFC, projection intensity")
} else {
  fig4_note_skip("Panel E", "Fewer than 3 regions had matched Cell_Count and Intensity effects.")
  NULL
}

make_fig4_pca <- function(data) {
  feature_long <- data %>%
    filter(RegionKey %in% fig4_region_selection$RegionKey) %>%
    select(SampleID, Animal, Condition, RegionKey, Cell_Count, Intensity) %>%
    pivot_longer(cols = c(Cell_Count, Intensity), names_to = "Metric", values_to = "RawValue") %>%
    mutate(Value = safe_log1p(RawValue), Feature = paste(Metric, RegionKey, sep = "__")) %>%
    group_by(SampleID, Animal, Condition, Feature) %>%
    summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
    mutate(Value = if_else(is.nan(Value), NA_real_, Value))
  wide <- feature_long %>%
    pivot_wider(names_from = Feature, values_from = Value) %>%
    arrange(factor(as.character(Condition), levels = nature_condition_levels), SampleID)
  if (nrow(wide) < 4) return(NULL)
  meta <- wide %>% select(SampleID, Animal, Condition)
  mat <- wide %>% select(-SampleID, -Animal, -Condition) %>% as.data.frame()
  rownames(mat) <- meta$SampleID
  mat <- as.matrix(mat)
  keep <- colSums(!is.na(mat)) >= 4 & apply(mat, 2, function(x) sd(x, na.rm = TRUE) > 0)
  mat <- mat[, keep, drop = FALSE]
  if (ncol(mat) < 3) return(NULL)
  for (j in seq_len(ncol(mat))) mat[is.na(mat[, j]), j] <- median(mat[, j], na.rm = TRUE)
  mat_z <- zscore_cols(mat)
  pca <- prcomp(mat_z, center = FALSE, scale. = FALSE)
  var_exp <- 100 * summary(pca)$importance[2, 1:2]
  pca_df <- as_tibble(pca$x[, 1:2], rownames = "SampleID") %>%
    left_join(meta, by = "SampleID") %>%
    mutate(Condition = factor(as.character(Condition), levels = nature_condition_levels), n_features = ncol(mat)) %>%
    group_by(Condition) %>%
    mutate(n_animals_condition = n_distinct(SampleID)) %>%
    ungroup()
  write_fig4_table(pca_df, "fig4_panelF_pca_coordinates")
  write_fig4_table(pca_df, "Fig4F_pca_coordinates", dir = publication_source_tab_dir)
  ggplot(pca_df, aes(x = PC1, y = PC2, fill = Condition)) +
    geom_point(shape = 21, size = 2.4, colour = "grey15", stroke = 0.22, alpha = 0.9) +
    scale_fill_manual(values = nature_condition_colors, labels = fig4_condition_labels, drop = FALSE) +
    fig4_theme(base_size = 6.8) +
    theme(legend.position = "bottom") +
    labs(title = "Descriptive animal-level profiles", subtitle = paste0(ncol(mat), " features"), x = paste0("PC1 (", scales::number(var_exp[1], accuracy = 0.1), "%)"), y = paste0("PC2 (", scales::number(var_exp[2], accuracy = 0.1), "%)"))
}
fig4F <- make_fig4_pca(long)
if (is.null(fig4F)) fig4_note_skip("Panel F", "PCA skipped because too few samples or variable features were available.")

save_fig4_plot(fig4A, "fig4A_condition_and_contrast_key", width = 6.8, height = 2.6, dir = publication_fig4_key_fig_dir)
save_fig4_plot(fig4B, "fig4B_projection_intensity_effect_size_map", width = 4.3, height = 5.8)
save_fig4_plot(fig4C, "fig4C_cfos_cell_count_effect_size_map", width = 4.3, height = 5.8)
save_fig4_plot(fig4D, "fig4D_key_region_condition_profiles", width = 6.8, height = 3.6)
save_fig4_plot(fig4E, "fig4E_cfos_cell_count_projection_dissociation", width = 3.8, height = 3.6)
save_fig4_plot(fig4F, "fig4F_systems_level_pca", width = 3.8, height = 3.5)

fig4_main <- (
  fig4A /
    (fig4B | fig4C) /
    fig4D /
    (fig4E | fig4F)
) +
  patchwork::plot_layout(heights = c(0.55, 1.35, 1.1, 1.05), guides = "collect") +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 9), legend.position = "bottom")

save_fig4_plot(fig4_main, "Fig4_learning_stress_cfos_cell_count_projection", width = 7.2, height = 10.4, dir = publication_fig4_composite_fig_dir)
save_fig4_plot(fig4_main, "Fig4_dashboard_learning_stress_cfos_cell_count_projection", width = 7.2, height = 10.4, dir = publication_dashboard_fig_dir, mirror_root = FALSE)

make_full_effect_heatmap <- function(metric) {
  d <- fig4_all_effects %>%
    filter(Metric == metric, !is.na(logFC)) %>%
    left_join(
      fig4_region_catalog %>%
        select(RegionKey, RegionDisplay, SystemGroup, MainRegionShort, MainRegionGroup),
      by = "RegionKey"
    ) %>%
    mutate(
      RegionShort = fig4_short_region(Annotation, Abbreviation),
      RegionDisplay = coalesce(RegionDisplay, RegionShort),
      contrast = factor(contrast, levels = fig4_main_contrasts)
    )
  if (nrow(d) == 0) return(NULL)
  top_regions <- d %>%
    group_by(RegionKey, RegionDisplay, SystemGroup, MainRegionShort) %>%
    summarise(max_abs = max(abs(logFC), na.rm = TRUE), .groups = "drop") %>%
    arrange(SystemGroup, MainRegionShort, desc(max_abs), RegionDisplay) %>%
    slice_head(n = 80)
  d <- d %>%
    filter(RegionKey %in% top_regions$RegionKey) %>%
    mutate(RegionDisplay = factor(RegionDisplay, levels = rev(top_regions$RegionDisplay)))
  lim <- max(abs(d$logFC), na.rm = TRUE)
  ggplot(d, aes(x = contrast, y = RegionDisplay, fill = logFC)) +
    geom_tile(colour = "white", linewidth = 0.12) +
    facet_grid(SystemGroup ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient2(
      low = effect_colors[["negative"]], mid = effect_colors[["neutral"]], high = effect_colors[["positive"]],
      midpoint = 0, limits = c(-lim, lim), oob = scales::squish
    ) +
    scale_x_discrete(labels = fig4_contrast_labels, drop = FALSE) +
    fig4_theme(base_size = 5.6) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 4.5),
      strip.text.y = element_text(angle = 0, size = 4.9)
    ) +
    labs(
      title = paste0("Supplementary all-region central effects: ", metric),
      subtitle = "Regions are grouped by anatomical class and ordered by parent-region stem",
      x = NULL,
      y = NULL
    )
}
save_fig4_plot(make_full_effect_heatmap("Intensity"), "supp_exploratory_projection_all_region_effect_heatmap", width = 7.0, height = 10.0, dir = fig4_supp_fig_dir)
save_fig4_plot(make_full_effect_heatmap("Cell_Count"), "supp_exploratory_cfos_all_region_effect_heatmap", width = 7.0, height = 10.0, dir = fig4_supp_fig_dir)
save_fig4_plot(make_full_effect_heatmap("Intensity"), "SuppFig4_projection_all_region_effect_heatmap", width = 7.0, height = 10.0, dir = exploratory_sensitivity_dir)
save_fig4_plot(make_full_effect_heatmap("Cell_Count"), "SuppFig4_cfos_cell_count_all_region_effect_heatmap", width = 7.0, height = 10.0, dir = exploratory_sensitivity_dir)

fig4_covariance_note <- "Exploratory Spearman covariance across animals; not interpreted as anatomical connectivity, functional connectivity, causality, or rewiring."

make_condition_correlation_outputs <- function(metric, value_label, use_normalized = FALSE) {
  metric_sym <- sym(metric)
  suffix <- if (use_normalized) paste0(metric, "_normalized") else metric

  min_region_n_all_conditions <- 3
  min_pairwise_n_covariance <- 3
  max_covariance_regions <- 20

  candidate_regions <- fig4_region_selection$RegionKey

  label_lookup <- fig4_region_selection %>%
    distinct(RegionKey, RegionDisplay) %>%
    deframe()

  covariance_caption_note <- paste(
    "Only regions with data from at least",
    min_region_n_all_conditions,
    "animals in every condition are shown.",
    "Remaining grey cells indicate region pairs with fewer than",
    min_pairwise_n_covariance,
    "animals with paired measurements in that condition."
  )

  covariance_region_meta <- long %>%
    filter(RegionKey %in% candidate_regions) %>%
    group_by(RegionKey) %>%
    summarise(
      Annotation = first(Annotation[!is.na(Annotation)]),
      Abbreviation = first(Abbreviation[!is.na(Abbreviation)]),
      Class = first(Class[!is.na(Class)]),
      Level = {
        level_values <- Level[!is.na(Level)]
        if (length(level_values) == 0) NA_integer_ else level_values[[1]]
      },
      RegionShort = label_lookup[RegionKey[[1]]],
      .groups = "drop"
    )

  covariance_region_coverage <- long %>%
    filter(RegionKey %in% candidate_regions) %>%
    select(SampleID, RegionKey, Condition, RawValue = !!metric_sym) %>%
    group_by(RegionKey, Condition) %>%
    summarise(
      n_non_missing = n_distinct(SampleID[!is.na(RawValue)]),
      .groups = "drop"
    ) %>%
    complete(
      RegionKey = candidate_regions,
      Condition = factor(nature_condition_levels, levels = nature_condition_levels),
      fill = list(n_non_missing = 0)
    ) %>%
    pivot_wider(
      names_from = Condition,
      values_from = n_non_missing,
      values_fill = 0,
      names_prefix = "n_"
    ) %>%
    left_join(covariance_region_meta, by = "RegionKey") %>%
    mutate(
      min_n_all_conditions = pmin(
        n_VEH_paired,
        n_VEH_unpaired,
        n_CNO_paired,
        n_CNO_unpaired
      ),
      passes_complete_covariance_filter =
        min_n_all_conditions >= min_region_n_all_conditions,
      covariance_filter_note = if_else(
        passes_complete_covariance_filter,
        paste0("included: n >= ", min_region_n_all_conditions, " in all four conditions"),
        paste0("excluded: fewer than ", min_region_n_all_conditions, " animals in at least one condition")
      ),
      covariance_note = fig4_covariance_note
    ) %>%
    arrange(
      desc(passes_complete_covariance_filter),
      match(RegionKey, candidate_regions)
    )

  write_fig4_table(
    covariance_region_coverage,
    paste0("Condition_covariance_region_coverage_", suffix),
    dir = qc_covariance_dir
  )

  selected_regions <- covariance_region_coverage %>%
    filter(passes_complete_covariance_filter) %>%
    arrange(match(RegionKey, candidate_regions)) %>%
    pull(RegionKey)

  pairwise_pruned <- prune_regions_by_pairwise_overlap(
    data = long,
    metric = metric,
    region_keys = selected_regions,
    min_pairwise_n = min_pairwise_n_covariance,
    conditions = nature_condition_levels
  )

  selected_regions <- pairwise_pruned$kept_regions

  write_fig4_table(
    pairwise_pruned$dropped_regions,
    paste0("Condition_covariance_pairwise_pruned_regions_", suffix),
    dir = qc_covariance_dir
  )

  write_fig4_table(
    pairwise_pruned$pair_qc,
    paste0("Condition_covariance_pairwise_qc_", suffix),
    dir = qc_covariance_dir
  )

  if (length(selected_regions) > max_covariance_regions) {
    selected_regions <- selected_regions[seq_len(max_covariance_regions)]
  }

  if (length(selected_regions) < 3) {
    fig4_note_skip(
      paste0("supplementary covariance heatmap ", suffix),
      paste0(
        "fewer than 3 regions passed n >= ",
        min_region_n_all_conditions,
        " in all four conditions"
      )
    )
    return(NULL)
  }

  selected_labels <- label_lookup[selected_regions]

  covariance_region_selection <- covariance_region_coverage %>%
    mutate(
      selected_for_complete_region_covariance = RegionKey %in% selected_regions
    ) %>%
    arrange(desc(selected_for_complete_region_covariance), match(RegionKey, candidate_regions))

  write_fig4_table(
    covariance_region_selection,
    paste0("Condition_covariance_region_selection_", suffix),
    dir = qc_covariance_dir
  )

  wide <- long %>%
    filter(RegionKey %in% selected_regions) %>%
    select(SampleID, Condition, RegionKey, RawValue = !!metric_sym) %>%
    mutate(Value = safe_log1p(RawValue)) %>%
    group_by(SampleID, Condition, RegionKey) %>%
    summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
    mutate(Value = if_else(is.nan(Value), NA_real_, Value)) %>%
    pivot_wider(names_from = RegionKey, values_from = Value)

  plot_rows <- purrr::map_dfr(nature_condition_levels, function(cond) {
    d <- wide %>% filter(Condition == cond)

    mat <- d %>%
      select(any_of(selected_regions)) %>%
      as.matrix()

    rownames(mat) <- d$SampleID

    missing_cols <- setdiff(selected_regions, colnames(mat))
    if (length(missing_cols) > 0) {
      mat <- cbind(
        mat,
        matrix(
          NA_real_,
          nrow = nrow(mat),
          ncol = length(missing_cols),
          dimnames = list(rownames(mat), missing_cols)
        )
      )
    }

    mat <- mat[, selected_regions, drop = FALSE]
    regs <- selected_regions

    cor_mat <- matrix(
      NA_real_,
      nrow = length(regs),
      ncol = length(regs),
      dimnames = list(regs, regs)
    )

    p_mat <- matrix(
      NA_real_,
      nrow = length(regs),
      ncol = length(regs),
      dimnames = list(regs, regs)
    )

    npair_mat <- matrix(
      0L,
      nrow = length(regs),
      ncol = length(regs),
      dimnames = list(regs, regs)
    )

    for (i in seq_along(regs)) {
      for (j in seq_along(regs)) {
        x <- mat[, i]
        y <- mat[, j]
        keep_pair <- is.finite(x) & is.finite(y)

        n_pair <- sum(keep_pair)
        npair_mat[i, j] <- n_pair

        if (
          n_pair >= min_pairwise_n_covariance &&
          stats::sd(x[keep_pair], na.rm = TRUE) > 0 &&
          stats::sd(y[keep_pair], na.rm = TRUE) > 0
        ) {
          test <- suppressWarnings(
            tryCatch(
              stats::cor.test(
                x[keep_pair],
                y[keep_pair],
                method = "spearman",
                exact = FALSE
              ),
              error = function(e) NULL
            )
          )

          if (!is.null(test)) {
            cor_mat[i, j] <- unname(test$estimate)
            p_mat[i, j] <- test$p.value
          }
        }
      }
    }

    diag(cor_mat) <- 1
    diag(p_mat) <- NA_real_

    cor_out <- as.data.frame(cor_mat) %>%
      rownames_to_column("RegionKey")

    p_out <- as.data.frame(p_mat) %>%
      rownames_to_column("RegionKey")

    npair_out <- as.data.frame(npair_mat) %>%
      rownames_to_column("RegionKey")

    readr::write_csv(
      cor_out,
      file.path(
        fig4_supp_tab_dir,
        paste0("supp_fig4_condition_correlation_matrix_", suffix, "_", cond, ".csv")
      )
    )
    readr::write_csv(
      cor_out,
      file.path(
        exploratory_covariance_dir,
        paste0("Condition_covariance_matrix_", suffix, "_", cond, ".csv")
      )
    )

    readr::write_csv(
      p_out,
      file.path(
        fig4_supp_tab_dir,
        paste0("supp_fig4_condition_pvalue_matrix_", suffix, "_", cond, ".csv")
      )
    )
    readr::write_csv(
      p_out,
      file.path(
        exploratory_covariance_dir,
        paste0("Condition_covariance_pvalue_matrix_", suffix, "_", cond, ".csv")
      )
    )

    readr::write_csv(
      npair_out,
      file.path(
        fig4_supp_tab_dir,
        paste0("supp_fig4_condition_npair_matrix_", suffix, "_", cond, ".csv")
      )
    )
    readr::write_csv(
      npair_out,
      file.path(
        exploratory_covariance_dir,
        paste0("Condition_covariance_npair_matrix_", suffix, "_", cond, ".csv")
      )
    )

    cor_long <- as_tibble(cor_mat, rownames = "Region1") %>%
      pivot_longer(cols = -Region1, names_to = "Region2", values_to = "rho") %>%
      left_join(
        as_tibble(npair_mat, rownames = "Region1") %>%
          pivot_longer(cols = -Region1, names_to = "Region2", values_to = "n_pair"),
        by = c("Region1", "Region2")
      ) %>%
      left_join(
        as_tibble(p_mat, rownames = "Region1") %>%
          pivot_longer(cols = -Region1, names_to = "Region2", values_to = "p_value"),
        by = c("Region1", "Region2")
      ) %>%
      left_join(
        tibble(
          Region1 = selected_regions,
          Region1Short = selected_labels[selected_regions]
        ),
        by = "Region1"
      ) %>%
      left_join(
        tibble(
          Region2 = selected_regions,
          Region2Short = selected_labels[selected_regions]
        ),
        by = "Region2"
      ) %>%
      mutate(
        Condition = cond,
        fdr = p.adjust(
          if_else(Region1 == Region2, NA_real_, p_value),
          method = "BH"
        ),
        fdr_sig = !is.na(fdr) & fdr < 0.10 & Region1 != Region2,
        significance_marker = if_else(fdr_sig, "FDR q < 0.10", NA_character_),
        Region1Short = factor(
          Region1Short,
        levels = rev(selected_labels[selected_regions])
        ),
        Region2Short = factor(
          Region2Short,
        levels = selected_labels[selected_regions]
        ),
        covariance_qc_flag = case_when(
          !is.na(rho) ~ "OK",
          n_pair < min_pairwise_n_covariance ~ "excluded_insufficient_pair_n",
          TRUE ~ "excluded_nonvariable_or_missing"
        ),
        insufficient_pair_note = if_else(
          n_pair < min_pairwise_n_covariance,
          paste0(
            "grey: fewer than ",
            min_pairwise_n_covariance,
            " animals with paired measurements"
          ),
          NA_character_
        ),
        complete_region_filter_note = covariance_caption_note,
        covariance_note = fig4_covariance_note,
        fdr_scope = "BH within each condition-specific covariance matrix; exploratory"
      )

    cor_long
  })

  write_fig4_table(
    plot_rows,
    paste0("SuppFig4_condition_covariance_", suffix),
    dir = exploratory_covariance_dir
  )

  system_group_colors <- c(
    "Amygdala" = "#D55E00",
    "Hypothalamus" = "#0072B2",
    "Thalamus" = "#CC79A7",
    "Cortex / hippocampus" = "#009E73",
    "Striatum / septum" = "#E69F00",
    "Pallidum / BST" = "#56B4E9",
    "Brainstem" = "#F0E442",
    "Other" = "#7F7F7F"
  )

  selected_region_meta <- covariance_region_selection %>%
    filter(selected_for_complete_region_covariance) %>%
    transmute(
      RegionKey,
      RegionShort = selected_labels[RegionKey],
      SystemGroup = factor(as.character(fig4_system_group(Class, Annotation, Abbreviation, RegionShort)),
                           levels = fig4_system_levels)
    ) %>%
    distinct(RegionKey, RegionShort, SystemGroup)

  x_levels <- selected_labels[selected_regions]
  y_levels <- rev(x_levels)

  plot_rows <- plot_rows %>%
    mutate(
      x_idx = as.integer(factor(as.character(Region2Short), levels = x_levels)),
      y_idx = as.integer(factor(as.character(Region1Short), levels = y_levels))
    )

  x_indicator <- tibble(
    RegionShort = x_levels,
    x_idx = seq_along(x_levels),
    y_idx = 0.35
  ) %>%
    left_join(selected_region_meta %>% select(RegionShort, SystemGroup), by = "RegionShort")

  y_indicator <- tibble(
    RegionShort = y_levels,
    y_idx = seq_along(y_levels),
    x_idx = 0.35
  ) %>%
    left_join(selected_region_meta %>% select(RegionShort, SystemGroup), by = "RegionShort")

  x_indicator_blocks <- x_indicator %>%
    mutate(
      SystemGroup = as.character(SystemGroup),
      block_id = cumsum(dplyr::coalesce(SystemGroup != dplyr::lag(SystemGroup), FALSE))
    ) %>%
    group_by(SystemGroup, block_id) %>%
    summarise(
      x_start = min(x_idx),
      x_end = max(x_idx),
      .groups = "drop"
    )

  y_indicator_blocks <- y_indicator %>%
    mutate(
      SystemGroup = as.character(SystemGroup),
      block_id = cumsum(dplyr::coalesce(SystemGroup != dplyr::lag(SystemGroup), FALSE))
    ) %>%
    group_by(SystemGroup, block_id) %>%
    summarise(
      y_start = min(y_idx),
      y_end = max(y_idx),
      .groups = "drop"
    )

  p <- ggplot(plot_rows, aes(x = x_idx, y = y_idx, fill = rho)) +
    geom_tile(colour = "white", linewidth = 0.08) +
    geom_point(
      data = plot_rows %>%
        filter(fdr_sig, Region1 != Region2),
      aes(x = x_idx, y = y_idx),
      inherit.aes = FALSE,
      shape = 21,
      size = 0.85,
      stroke = 0.25,
      colour = "black",
      fill = NA
    ) +
    geom_segment(
      data = x_indicator_blocks,
      aes(
        x = x_start - 0.5,
        xend = x_end + 0.5,
        y = 0.25,
        yend = 0.25,
        colour = SystemGroup
      ),
      inherit.aes = FALSE,
      linewidth = 2.6,
      lineend = "butt"
    ) +
    geom_segment(
      data = y_indicator_blocks,
      aes(
        x = 0.25,
        xend = 0.25,
        y = y_start - 0.5,
        yend = y_end + 0.5,
        colour = SystemGroup
      ),
      inherit.aes = FALSE,
      linewidth = 2.6,
      lineend = "butt"
    ) +
    facet_wrap(~ Condition, ncol = 2) +
    scale_x_continuous(
      breaks = seq_along(x_levels),
      labels = x_levels,
      limits = c(0, length(x_levels) + 0.5),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      breaks = seq_along(y_levels),
      labels = y_levels,
      limits = c(0, length(y_levels) + 0.5),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_fill_gradient2(
      low = effect_colors[["negative"]],
      mid = effect_colors[["neutral"]],
      high = effect_colors[["positive"]],
      midpoint = 0,
      limits = c(-1, 1),
      oob = scales::squish,
      na.value = "grey88",
      name = "Spearman\nrho"
    ) +
    scale_colour_manual(
      values = system_group_colors,
      breaks = fig4_system_levels,
      drop = FALSE,
      na.value = "grey70",
      name = "Region system"
    ) +
    fig4_theme(base_size = 5.8) +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_text(angle = 60, hjust = 1, size = 4.7),
      axis.text.y = element_text(size = 4.7),
      legend.position = "right"
    ) +
    labs(
      title = paste0("Exploratory inter-regional covariance: ", value_label),
      subtitle = paste0(
        "Fill = Spearman rho; open circle = FDR q < 0.10; ",
        "colored axis-side strips = contiguous region system groups. ",
        "Grey = insufficient paired n or zero variance."
      ),
      x = NULL,
      y = NULL
    )

  save_fig4_plot(
    p,
    paste0("SuppFig4_condition_covariance_", suffix),
    width = 6.8,
    height = 7.0,
    dir = exploratory_covariance_dir
  )

  invisible(plot_rows)
}

make_condition_correlation_outputs("Cell_Count", "cFos+ cell count")
make_condition_correlation_outputs("Intensity", "projection intensity")

make_cem_seed_covariance <- function(metric, value_label, max_targets_to_plot = 24, min_seed_pair_n = 3) {
  metric_sym <- sym(metric)
  cem_key <- find_cem_region(long)
  if (is.na(cem_key) || length(cem_key) == 0) return(NULL)

  selected_regions <- setdiff(fig4_region_catalog$RegionKey, cem_key)
  wide <- long %>%
    filter(RegionKey %in% c(cem_key, selected_regions)) %>%
    select(SampleID, Condition, RegionKey, RawValue = !!metric_sym) %>%
    mutate(Value = safe_log1p(RawValue)) %>%
    group_by(SampleID, Condition, RegionKey) %>%
    summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
    mutate(Value = if_else(is.nan(Value), NA_real_, Value)) %>%
    pivot_wider(names_from = RegionKey, values_from = Value)

  out <- purrr::map_dfr(nature_condition_levels, function(cond) {
    d <- wide %>% filter(Condition == cond)
    if (!cem_key %in% names(d)) return(tibble())
    purrr::map_dfr(selected_regions[selected_regions %in% names(d)], function(target) {
      n_pair <- sum(is.finite(d[[cem_key]]) & is.finite(d[[target]]))
      seed_sd <- if (n_pair >= min_seed_pair_n) stats::sd(d[[cem_key]][is.finite(d[[cem_key]]) & is.finite(d[[target]])], na.rm = TRUE) else NA_real_
      target_sd <- if (n_pair >= min_seed_pair_n) stats::sd(d[[target]][is.finite(d[[cem_key]]) & is.finite(d[[target]])], na.rm = TRUE) else NA_real_
      estimable <- n_pair >= min_seed_pair_n && is.finite(seed_sd) && is.finite(target_sd) && seed_sd > 0 && target_sd > 0
      rho <- if (estimable) suppressWarnings(cor(d[[cem_key]], d[[target]], method = "spearman", use = "pairwise.complete.obs")) else NA_real_
      p_value <- if (estimable && !is.na(rho)) suppressWarnings(tryCatch(cor.test(d[[cem_key]], d[[target]], method = "spearman", exact = FALSE)$p.value, error = function(e) NA_real_)) else NA_real_
      tibble(
        Condition = cond,
        SeedRegionKey = cem_key,
        TargetRegionKey = target,
        n_pair = n_pair,
        seed_sd = seed_sd,
        target_sd = target_sd,
        rho = rho,
        p_value = p_value,
        covariance_qc_flag = case_when(
          estimable & n_pair >= 4 ~ "OK",
          estimable & n_pair == 3 ~ "low_n_exploratory",
          n_pair >= min_seed_pair_n ~ "excluded_nonvariable",
          TRUE ~ "excluded_insufficient_n"
        )
      )
    }) %>%
      mutate(fdr = p.adjust(p_value, method = "BH"))
  }) %>%
    left_join(
      fig4_region_catalog %>%
        select(
          TargetRegionKey = RegionKey,
          TargetAnnotation = Annotation,
          TargetAbbreviation = Abbreviation,
          TargetClass = Class,
          TargetSystemGroup = SystemGroup,
          TargetMainRegion = MainRegion,
          TargetMainRegionAbbreviation = MainRegionAbbreviation,
          TargetMainRegionShort = MainRegionShort,
          TargetMainRegionGroup = MainRegionGroup,
          TargetShort = RegionShort,
          TargetDisplay = RegionDisplay
        ),
      by = "TargetRegionKey"
    ) %>%
    mutate(
      Metric = metric,
      fdr_sig = !is.na(fdr) & fdr < 0.10,
      covariance_note = paste0(
        "CeM-seed covariance tests whether inter-animal variation in CeM signal is associated with variation in other regions. ",
        "All detected regions are tested; the plot shows the strongest estimable targets for readability. ",
        "It is not interpreted as causal or anatomical connectivity."
      ),
      fdr_scope = "BH within each condition for CeM seed versus all estimable target regions; exploratory"
    )

  cem_target_ranking <- out %>%
    filter(!is.na(rho)) %>%
    group_by(
      TargetRegionKey, TargetAnnotation, TargetAbbreviation, TargetClass, TargetSystemGroup,
      TargetMainRegion, TargetMainRegionAbbreviation, TargetMainRegionShort,
      TargetMainRegionGroup, TargetShort, TargetDisplay
    ) %>%
    summarise(
      max_abs_rho = max(abs(rho), na.rm = TRUE),
      min_p_value = min(p_value, na.rm = TRUE),
      min_fdr = min(fdr, na.rm = TRUE),
      n_estimable_conditions = n_distinct(Condition[!is.na(rho)]),
      any_fdr_sig = any(fdr_sig, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      min_p_value = if_else(is.finite(min_p_value), min_p_value, NA_real_),
      min_fdr = if_else(is.finite(min_fdr), min_fdr, NA_real_)
    ) %>%
    arrange(desc(any_fdr_sig), desc(max_abs_rho), min_p_value)

  out <- out %>%
    left_join(cem_target_ranking %>% select(TargetRegionKey, cem_target_rank = max_abs_rho), by = "TargetRegionKey") %>%
    arrange(desc(!is.na(rho)), TargetSystemGroup, TargetMainRegionShort, TargetAnnotation, factor(Condition, levels = nature_condition_levels))

  readr::write_csv(out, file.path(exploratory_covariance_dir, paste0("CeM_seed_covariance_", metric, ".csv")))
  openxlsx::write.xlsx(out, file.path(exploratory_covariance_dir, paste0("CeM_seed_covariance_", metric, ".xlsx")), overwrite = TRUE)
  write_fig4_table(cem_target_ranking %>% mutate(fdr_scope = "BH within each condition for CeM seed versus all estimable target regions; exploratory"),
                   paste0("CeM_seed_covariance_target_ranking_", metric), dir = exploratory_covariance_dir)

  plot_targets <- cem_target_ranking %>%
    slice_head(n = max_targets_to_plot) %>%
    pull(TargetRegionKey)

  plot_target_levels <- cem_target_ranking %>%
    filter(TargetRegionKey %in% plot_targets) %>%
    arrange(TargetSystemGroup, TargetMainRegionShort, desc(any_fdr_sig), desc(max_abs_rho), TargetDisplay) %>%
    pull(TargetDisplay) %>%
    unique()

  p <- out %>%
    filter(TargetRegionKey %in% plot_targets) %>%
    mutate(
      TargetDisplay = factor(
        TargetDisplay,
        levels = rev(plot_target_levels)
      ),
      Condition = factor(Condition, levels = nature_condition_levels),
      fdr_marker = if_else(fdr_sig, "FDR q < 0.10", NA_character_)
    ) %>%
    ggplot(aes(x = Condition, y = TargetDisplay, fill = rho)) +
    geom_tile(colour = "white", linewidth = 0.12) +
    geom_point(aes(size = n_pair), shape = 21, colour = "grey35", stroke = 0.2, fill = NA) +
    geom_point(data = ~ dplyr::filter(.x, fdr_sig), shape = 21, size = 1.9, colour = "black", stroke = 0.35, fill = NA) +
    facet_grid(TargetSystemGroup ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient2(low = effect_colors[["negative"]], mid = effect_colors[["neutral"]], high = effect_colors[["positive"]], midpoint = 0,
                         limits = c(-1, 1), oob = scales::squish, na.value = "grey88", name = "Spearman\nrho") +
    scale_size_continuous(range = c(0.4, 1.5), breaks = sort(unique(out$n_pair)), name = "n") +
    fig4_theme(base_size = 6.0) +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1),
      axis.text.y = element_text(size = 5.0),
      strip.text.y = element_text(angle = 0, size = 4.9),
      legend.position = "right"
    ) +
    labs(
      title = paste0("Exploratory CeM-seed covariance: ", value_label),
      subtitle = paste0("Fill = Spearman rho; open black circle = FDR q < 0.10. Top ", max_targets_to_plot, " estimable targets shown."),
      x = NULL,
      y = NULL
    )

  save_fig4_plot(p, paste0("CeM_seed_covariance_", metric), width = 5.6, height = 6.5, dir = exploratory_covariance_dir)
  invisible(out)
}

make_cem_seed_covariance("Cell_Count", "cFos+ cell count")
make_cem_seed_covariance("Intensity", "projection intensity")

fig4_seed_terms <- c("CEAl", "CeM", "LA", "PVH")
fig4_seed_candidates <- purrr::map_dfr(fig4_seed_terms, function(seed_term) {
  term_lower <- str_to_lower(seed_term)
  candidates <- fig4_region_catalog %>%
    filter(
      Annotation_lower == term_lower |
        Abbreviation_lower == term_lower |
        str_detect(Annotation_lower, fixed(term_lower)) |
        str_detect(Abbreviation_lower, fixed(term_lower))
    ) %>%
    left_join(fig4_effect_rank_region_summary %>% select(RegionKey, n_estimable_contrasts_total, max_abs_effect, min_raw_p), by = "RegionKey") %>%
    mutate(
      n_estimable_contrasts_total = coalesce(n_estimable_contrasts_total, 0),
      max_abs_effect = coalesce(max_abs_effect, 0),
      min_raw_p = coalesce(min_raw_p, Inf),
      SeedLabel = seed_term
    ) %>%
    arrange(desc(n_estimable_contrasts_total), desc(max_abs_effect), min_raw_p, RegionKey)
  if (nrow(candidates) == 0) {
    tibble(SeedLabel = seed_term, SeedRegionKey = NA_character_, seed_selection_note = "not found")
  } else {
    candidates %>%
      slice(1) %>%
      transmute(SeedLabel, SeedRegionKey = RegionKey, seed_selection_note = "first data-matched seed candidate")
  }
}) %>%
  filter(!is.na(SeedRegionKey))
write_fig4_table(fig4_seed_candidates, "supp_exploratory_covariance_seed_matching", dir = fig4_supp_tab_dir)
write_fig4_table(fig4_seed_candidates, "Covariance_seed_matching", dir = exploratory_covariance_dir)

fig4_covariation_results <- tibble()
if (nrow(fig4_seed_candidates) > 0) {
  sm_cell_fig4 <- make_sample_matrix(long, "Cell_Count")
  mat_cell_fig4 <- safe_log1p(sm_cell_fig4$mat)
  meta_cell_fig4 <- sm_cell_fig4$annotation %>% mutate(Condition = factor(as.character(Condition), levels = nature_condition_levels))
  fig4_covariation_results <- purrr::pmap_dfr(
    fig4_seed_candidates %>% select(SeedLabel, SeedRegionKey),
    function(SeedLabel, SeedRegionKey) {
    purrr::map_dfr(nature_condition_levels, function(cond) {
      ids <- meta_cell_fig4 %>% filter(Condition == cond) %>% pull(SampleID)
      submat <- mat_cell_fig4[rownames(mat_cell_fig4) %in% ids, , drop = FALSE]
      if (!SeedRegionKey %in% colnames(submat) || nrow(submat) < min_pairwise_n) return(tibble())
      tibble(TargetRegionKey = setdiff(colnames(submat), SeedRegionKey)) %>%
        rowwise() %>%
        mutate(
          SeedLabel = SeedLabel,
          SeedRegionKey = SeedRegionKey,
          Condition = cond,
          n_pair = sum(is.finite(submat[, SeedRegionKey]) & is.finite(submat[, TargetRegionKey])),
          rho = if (n_pair >= min_pairwise_n) suppressWarnings(cor(submat[, SeedRegionKey], submat[, TargetRegionKey], method = "spearman", use = "pairwise.complete.obs")) else NA_real_,
          p_value = if (n_pair >= min_pairwise_n && !is.na(rho)) suppressWarnings(tryCatch(cor.test(submat[, SeedRegionKey], submat[, TargetRegionKey], method = "spearman", exact = FALSE)$p.value, error = function(e) NA_real_)) else NA_real_
        ) %>%
        ungroup() %>%
        group_by(SeedLabel, SeedRegionKey, Condition) %>%
        mutate(fdr = p.adjust(p_value, method = "BH")) %>%
        ungroup()
    })
  }) %>%
    left_join(fig4_seed_candidates, by = c("SeedLabel", "SeedRegionKey")) %>%
    left_join(fig4_region_catalog %>% select(TargetRegionKey = RegionKey, TargetAnnotation = Annotation, TargetAbbreviation = Abbreviation, TargetClass = Class), by = "TargetRegionKey") %>%
    mutate(
      covariance_display = n_pair >= fig4_cov_min_n_pair & !is.na(rho),
      thresholded_covariance = covariance_display & abs(rho) >= fig4_cov_abs_r_cutoff & fdr <= fig4_cov_fdr_cutoff,
      terminology_note = fig4_covariance_note,
      fdr_scope = "BH within each seed and condition across tested target regions; exploratory"
    )
}
write_fig4_table(fig4_covariation_results, "supp_exploratory_seed_based_cfos_covariance", dir = fig4_supp_tab_dir)
write_fig4_table(fig4_covariation_results, "Seed_based_cfos_cell_count_covariance", dir = exploratory_covariance_dir)

fig4_cov_summary <- fig4_covariation_results %>%
  filter(!is.na(rho)) %>%
  group_by(SeedLabel, Condition) %>%
  summarise(
    n_estimated_pairs = n(),
    n_display_pairs = sum(covariance_display, na.rm = TRUE),
    n_thresholded_positive = sum(thresholded_covariance & rho > 0, na.rm = TRUE),
    n_thresholded_negative = sum(thresholded_covariance & rho < 0, na.rm = TRUE),
    mean_abs_rho = mean(abs(rho[covariance_display]), na.rm = TRUE),
    terminology_note = fig4_covariance_note,
    .groups = "drop"
  )
write_fig4_table(fig4_cov_summary, "supp_exploratory_cfos_covariance_summary", dir = fig4_supp_tab_dir)
write_fig4_table(fig4_cov_summary %>% mutate(fdr_scope = "BH within each seed and condition across tested target regions; exploratory"),
                 "Seed_based_cfos_cell_count_covariance_summary", dir = exploratory_covariance_dir)

publication_output_index <- tibble(
  output_group = c(
    "Main manuscript figure",
    "Main manuscript panel figures",
    "Main manuscript keys",
    "Main manuscript tables",
    "Supplementary exploratory figures",
    "Supplementary exploratory tables",
    "Exploratory regional profile figures",
    "Exploratory effect-size figures",
    "Legacy exploratory figures",
    "Legacy exploratory tables"
  ),
  directory = c(
    publication_fig4_composite_fig_dir,
    publication_fig4_panel_fig_dir,
    publication_fig4_key_fig_dir,
    fig4_tab_dir,
    fig4_supp_fig_dir,
    fig4_supp_tab_dir,
    exploratory_profile_fig_dir,
    exploratory_effect_fig_dir,
    legacy_fig_dir,
    legacy_tab_dir
  ),
  recommended_use = c(
    "Use first for the assembled Fig. 4 composite; SVG/PDF are vector outputs and PNG is 600 dpi.",
    "Use for editable standalone Fig. 4 panels.",
    "Use for condition and contrast explanation panels.",
    "Use for source data, region selection documentation, contrast sign key, and QC tables tied to Fig. 4.",
    "Use as supplementary/exploratory displays for full effects, covariance heatmaps, CeM-seed covariance, and overlap summaries.",
    "Use as supplementary source tables; includes full edge lists, cell-count/intensity overlap, CeM rankings, and pairwise-n QC.",
    "Use for exploratory volcano, regional profile, and class-level profile figures.",
    "Use for exploratory all-region effect-size maps.",
    "Earlier exploratory figures retained for audit trail; not the primary publication set.",
    "Earlier exploratory tables retained for audit trail; not the primary publication set."
  )
)
write_fig4_table(publication_output_index, "publication_output_index", dir = fig4_tab_dir)

fig4_output_manifest <- tibble::tribble(
  ~output_file, ~figure_panel, ~metric, ~analysis_type, ~normalization, ~imputation, ~contrast_definition, ~source_function, ~source_table, ~intended_use, ~exploratory_flag, ~fdr_scope,
  file.path(fig4_fig_dir, "fig4B_projection_intensity_effect_size_map.svg"), "Fig4B", "Intensity", "regional effect-size heatmap", "none", "none for displayed effects; limma model used median imputation internally", paste(fig4_main_contrasts, collapse = "; "), "make_fig4_heatmap", file.path(publication_source_tab_dir, "Fig4B_projection_effect_sizes.csv"), "main publication effect-size map", FALSE, "FDR values are global within contrast/metric tables; heatmap inclusion is not significance based",
  file.path(fig4_fig_dir, "fig4C_cfos_cell_count_effect_size_map.svg"), "Fig4C", "Cell_Count", "regional effect-size heatmap", "none", "none for displayed effects; limma model used median imputation internally", paste(fig4_main_contrasts, collapse = "; "), "make_fig4_heatmap", file.path(publication_source_tab_dir, "Fig4C_cfos_cell_count_effect_sizes.csv"), "main publication effect-size map", FALSE, "FDR values are global within contrast/metric tables; heatmap inclusion is not significance based",
  file.path(fig4_fig_dir, "fig4D_key_region_condition_profiles.svg"), "Fig4D", "Cell_Count; Intensity", "raw condition profiles", "none", "none", "condition-level descriptive profiles", "fig4D_source block", file.path(fig4_tab_dir, "fig4_panelD_key_region_profiles.csv"), "main publication descriptive profiles", FALSE, "not applicable",
  file.path(fig4_fig_dir, "fig4E_cfos_cell_count_projection_dissociation.svg"), "Fig4E", "Cell_Count; Intensity", "paired effect-size scatter", "none", "none for displayed effects; limma model used median imputation internally", "Learning_effect", "fig4E_source block", file.path(fig4_tab_dir, "fig4_panelE_cfos_cell_count_projection_scatter.csv"), "main publication descriptive effect comparison", FALSE, "FDR values are global within contrast/metric tables; panel is descriptive",
  file.path(fig4_fig_dir, "fig4F_systems_level_pca.svg"), "Fig4F", "Cell_Count; Intensity", "descriptive PCA", "none", "region-median feature imputation for PCA only", "animal-level structure, no contrast test", "make_fig4_pca", file.path(fig4_tab_dir, "fig4_panelF_pca_coordinates.csv"), "descriptive animal-level structure only", FALSE, "not applicable",
  file.path(exploratory_covariance_dir, "SuppFig4_condition_covariance_Cell_Count.csv"), "Supplementary", "Cell_Count", "condition covariance", "none", "none", "within-condition Spearman covariance", "make_condition_correlation_outputs", file.path(exploratory_covariance_dir, "SuppFig4_condition_covariance_Cell_Count.csv"), "supplementary exploratory covariance", TRUE, "BH within each condition-specific covariance matrix; exploratory",
  file.path(exploratory_covariance_dir, "CeM_seed_covariance_Cell_Count.csv"), "Supplementary", "Cell_Count", "CeM seed covariance", "none", "none", "within-condition Spearman covariance with CeM seed", "make_cem_seed_covariance", file.path(exploratory_covariance_dir, "CeM_seed_covariance_Cell_Count.csv"), "supplementary exploratory covariance", TRUE, "BH within each condition for CeM seed versus all estimable target regions; exploratory",
  file.path(exploratory_sensitivity_dir, "Sensitivity_stability_raw_normalized_complete_case_imputed.csv"), "QC/Sensitivity", "Cell_Count; Intensity", "stability comparison", "raw and normalized log1p", "complete-case and median-imputed limma variants", paste(fig4_main_contrasts, collapse = "; "), "run_fig4_complete_case_contrasts", file.path(exploratory_sensitivity_dir, "Sensitivity_stability_raw_normalized_complete_case_imputed.csv"), "sensitivity and robustness audit", TRUE, "FDR values inherited from each variant's contrast/metric scope"
)
readr::write_csv(fig4_output_manifest, file.path(publication_manifest_dir, "output_manifest.csv"))

fig4_readme <- c(
  "Manuscript-ready Fig. 4 generation",
  "",
  paste0("Input data used: ", paste(sort(unique(long$SourceFile)), collapse = "; ")),
  "Main Fig. 4 output: results/01_publication/figures/fig4/composite/Fig4_learning_stress_cfos_cell_count_projection.svg/pdf/png.",
  "Publication outputs are organized under results/01_publication; exploratory covariance, network, dimensionality-reduction, and sensitivity outputs are under results/02_exploratory; QC outputs are under results/03_qc.",
  "Main panels: A condition/contrast key, B projection intensity effects, C cFos+ cell-count effects, D raw key-region profiles, E cFos cell-count/projection dissociation, F descriptive animal-level PCA.",
  "Central contrasts visualized:",
  paste0("- ", fig4_sign_key$short_label, ": ", fig4_sign_key$contrast_formula, "; positive = ", fig4_sign_key$positive_logFC_means, ". ", fig4_sign_key$plain_language),
  "",
  "Transformation: log1p. Main effect maps use raw, non-normalized values by default. Normalized log1p analyses are exported for comparison.",
  "Imputation: main heatmap source tables retain NA displayed effects; grey heatmap tiles mean missing, non-estimable, insufficient n, or excluded from main display. Median-imputed limma results are retained as one model variant and are not the only sensitivity result.",
  paste0("Region inclusion: no manuscript-prioritized anatomy is used for main Fig. 4 selection. Regions require min relevant group n >= 2 for displayed contrasts and >= ", fig4_min_main_non_na_display, " displayed central effects within a metric; each heatmap then shows the top ", fig4_max_heatmap_regions, " regions ranked by maximum absolute displayed logFC. Raw P is not used as a main-display inclusion rule."),
  "Main Fig. 4 heatmaps are regional effect-size maps, not whole-brain discovery significance maps.",
  "Visible significance markings use one statistic: FDR q < 0.10 within the stated correction scope. Low-n effects are marked separately as exploratory.",
  "Projection intensity and regional cFos recruitment are treated as related but distinct readouts; no causal relation is implied.",
  "Network edge comparisons use VEH_unpaired as the explicit baseline condition and are exploratory.",
  "Supplementary cell-count/intensity edge-overlap tables and plots identify thresholded region pairs shared by both readouts.",
  "CeM-seed covariance tables test CeM against all detected target regions; plots show the strongest estimable targets for readability.",
  "Supplementary condition-specific inter-regional covariance heatmaps, CeM-seed covariance, and network analyses are exploratory Spearman covariance analyses. They are not anatomical connectivity, functional connectivity, causality, or literal rewiring.",
  paste0("FDR scope: see explicit fdr_scope columns and manifest at ", file.path(publication_manifest_dir, "output_manifest.csv"), "."),
  paste0("Session info: ", file.path(session_info_dir, "session_info.txt"), "."),
  "",
  "Skipped panels or warnings:",
  if (length(fig4_skipped_panels) == 0) "- None" else paste0("- ", fig4_skipped_panels)
)
writeLines(fig4_readme, con = file.path(fig4_tab_dir, "README_Fig4_publication.txt"))
writeLines(fig4_readme, con = file.path(fig4_fig_dir, "README_Fig4_publication.txt"))

sync_legacy_outputs <- function() {
  copy_nonrecursive <- function(from_dir, to_dir) {
    if (!dir.exists(from_dir)) return(invisible(NULL))
    dir.create(to_dir, recursive = TRUE, showWarnings = FALSE)
    files <- list.files(from_dir, full.names = TRUE, recursive = FALSE)
    files <- files[file.exists(files) & !dir.exists(files)]
    if (length(files) > 0) {
      file.copy(files, file.path(to_dir, basename(files)), overwrite = TRUE)
    }
    invisible(NULL)
  }
  copy_nonrecursive(file.path(out_dir, "figures"), legacy_fig_dir)
  copy_nonrecursive(file.path(out_dir, "tables"), legacy_tab_dir)
}
sync_legacy_outputs()

output_directory_audit <- tibble(
  output_directory = c(
    publication_fig4_fig_dir,
    publication_fig4_panel_fig_dir,
    publication_fig4_composite_fig_dir,
    publication_fig4_key_fig_dir,
    publication_supp_fig_dir,
    publication_dashboard_fig_dir,
    publication_fig4_tab_dir,
    publication_supp_tab_dir,
    publication_source_tab_dir,
    publication_manifest_dir,
    exploratory_effect_fig_dir,
    exploratory_profile_fig_dir,
    exploratory_covariance_dir,
    exploratory_network_dir,
    exploratory_dimred_dir,
    exploratory_sensitivity_dir,
    qc_missingness_dir,
    qc_normalization_dir,
    qc_region_selection_dir,
    qc_covariance_dir,
    legacy_fig_dir,
    legacy_tab_dir,
    session_info_dir
  )
) %>%
  mutate(
    n_files = purrr::map_int(output_directory, ~ if (dir.exists(.x)) length(list.files(.x, recursive = TRUE, all.files = FALSE)) else 0L),
    note = case_when(
      n_files > 0 ~ "populated",
      output_directory == publication_dashboard_fig_dir ~ "empty unless dashboard exports are enabled/generated",
      TRUE ~ "empty; upstream analysis may have skipped this output family"
    )
  )
readr::write_csv(output_directory_audit, file.path(publication_manifest_dir, "output_directory_audit.csv"))

# -----------------------------
# 16. Save session info
# -----------------------------
sink(file.path(out_dir, "session_info.txt"))
print(sessionInfo())
sink()
sink(file.path(session_info_dir, "session_info.txt"))
print(sessionInfo())
sink()

message("Done. Outputs written to: ", out_dir)
