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
intensity_raw_file <- file.path(input_dir, "IntensityNoOutliersFin.xlsx")
cell_count_raw_file <- file.path(input_dir, "CellCountNoOutliersFin.csv")

out_dir <- file.path(dirname(input_file), "region_analysis_outputs_individual_level")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

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

min_pairwise_n <- 4
network_abs_r_cutoff <- 0.70
network_fdr_cutoff <- 0.10

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
  )

openxlsx::write.xlsx(long, file.path(out_dir, "tables", "merged_long_region_data.xlsx"), overwrite = TRUE)
readr::write_csv(long, file.path(out_dir, "tables", "merged_long_region_data.csv"))

# Animal counts sanity check
animal_counts <- long %>%
  distinct(SampleID, Group, Condition) %>%
  count(Group, Condition, name = "n_animals")

openxlsx::write.xlsx(animal_counts, file.path(out_dir, "tables", "QC_animal_counts_by_condition.xlsx"), overwrite = TRUE)

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

qc_region_availability <- long %>%
  group_by(Class, Annotation, Abbreviation, RegionLabel, RegionKey, Condition) %>%
  summarise(n_samples_present = n_distinct(SampleID), .groups = "drop") %>%
  pivot_wider(names_from = Condition, values_from = n_samples_present, values_fill = 0) %>%
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

avail_mat <- qc_region_availability %>%
  select(RegionKey, all_of(levels(long$Condition))) %>%
  column_to_rownames("RegionKey") %>%
  as.matrix()

pdf(file.path(out_dir, "figures", "QC_region_availability_by_condition.pdf"), width = 8, height = 12)
pheatmap::pheatmap(
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
  CNO_paired = "#3B6EA8",
  VEH_paired = "#111111",
  CNO_unpaired = "#C44E52",
  VEH_unpaired = "#4C956C"
)

theme_nature <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1),
      plot.subtitle = element_text(size = base_size, colour = "grey30"),
      axis.line = element_line(linewidth = 0.25),
      axis.ticks = element_line(linewidth = 0.25),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size - 1),
      legend.title = element_text(size = base_size - 1),
      legend.text = element_text(size = base_size - 1)
    )
}

save_figure <- function(plot, filename_base, width, height, dpi = 600) {
  ggsave(file.path(out_dir, "figures", paste0(filename_base, ".pdf")), plot, width = width, height = height)
  ggsave(file.path(out_dir, "figures", paste0(filename_base, ".png")), plot, width = width, height = height, dpi = dpi, bg = "white")
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
  pheatmap::pheatmap(
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

# -----------------------------
# 6. Group heatmaps, sample heatmaps, correlations, networks
# -----------------------------
cor_with_p <- function(mat, min_n = 4) {
  regions <- colnames(mat)

  expand.grid(region1 = regions, region2 = regions, stringsAsFactors = FALSE) %>%
    as_tibble() %>%
    rowwise() %>%
    mutate(
      n_pair = sum(complete.cases(mat[, region1], mat[, region2])),
      r = ifelse(
        n_pair >= min_n,
        suppressWarnings(cor(mat[, region1], mat[, region2], use = "pairwise.complete.obs", method = "spearman")),
        NA_real_
      ),
      p = ifelse(
        n_pair >= min_n && !is.na(r) && region1 != region2,
        suppressWarnings(cor.test(mat[, region1], mat[, region2], method = "spearman", exact = FALSE)$p.value),
        NA_real_
      )
    ) %>%
    ungroup() %>%
    mutate(fdr = p.adjust(p, method = "BH"))
}

plot_correlation_heatmap <- function(mat, metric, class_filter = NULL) {
  title_suffix <- ifelse(is.null(class_filter), "all_classes", str_replace_all(class_filter, "[^A-Za-z0-9]+", "_"))

  keep <- colSums(!is.na(mat)) >= min_pairwise_n &
    apply(mat, 2, function(x) sd(x, na.rm = TRUE) > 0)

  mat <- mat[, keep, drop = FALSE]
  if (ncol(mat) < 3) return(NULL)

  cor_mat <- suppressWarnings(cor(mat, use = "pairwise.complete.obs", method = "spearman"))
  cor_mat <- clean_heatmap_matrix(cor_mat, fill = 0)
  diag(cor_mat) <- 1

  pdf(file.path(out_dir, "figures", paste0("correlation_heatmap_", metric, "_", title_suffix, ".pdf")), width = 9, height = 9)
  pheatmap::pheatmap(
    cor_mat,
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
    breaks = seq(-1, 1, length.out = 101),
    main = paste0("Spearman region-region correlation: ", metric, " | ", title_suffix),
    fontsize_row = ifelse(ncol(cor_mat) > 80, 3, 6),
    fontsize_col = ifelse(ncol(cor_mat) > 80, 3, 6),
    border_color = NA
  )
  dev.off()
}

plot_network <- function(mat, metric, class_filter = NULL) {
  title_suffix <- ifelse(is.null(class_filter), "all_classes", str_replace_all(class_filter, "[^A-Za-z0-9]+", "_"))

  keep <- colSums(!is.na(mat)) >= min_pairwise_n &
    apply(mat, 2, function(x) sd(x, na.rm = TRUE) > 0)

  mat <- mat[, keep, drop = FALSE]
  if (ncol(mat) < 4) return(NULL)

  cors <- cor_with_p(mat, min_n = min_pairwise_n) %>%
    filter(region1 < region2, !is.na(r), abs(r) >= network_abs_r_cutoff, fdr <= network_fdr_cutoff)

  readr::write_csv(cors, file.path(out_dir, "tables", paste0("network_edges_", metric, "_", title_suffix, ".csv")))

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
      title = paste0("Region correlation network: ", metric, " | ", title_suffix),
      subtitle = paste0("Spearman |r| >= ", network_abs_r_cutoff, ", FDR <= ", network_fdr_cutoff)
    )

  ggsave(file.path(out_dir, "figures", paste0("network_", metric, "_", title_suffix, ".pdf")),
         p, width = 9, height = 7)
}

all_classes <- sort(unique(long$Class))

for (metric in metrics_to_analyse) {
  message("Analysing metric: ", metric)

  summary_df <- make_group_summary_matrix(long, metric)
  readr::write_csv(summary_df, file.path(out_dir, "tables", paste0("group_summary_", metric, ".csv")))

  plot_group_heatmap(summary_df, metric, class_filter = NULL)
  purrr::walk(all_classes, ~ plot_group_heatmap(summary_df, metric, class_filter = .x))

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
          width = 10, height = max(5, nrow(mat_log) * 0.20))
      pheatmap::pheatmap(
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

  plot_correlation_heatmap(safe_log1p(sm$mat), metric, class_filter = NULL)
  plot_network(safe_log1p(sm$mat), metric, class_filter = NULL)

  for (cl in all_classes) {
    sm_class <- make_sample_matrix(long, metric, class_filter = cl)
    plot_correlation_heatmap(safe_log1p(sm_class$mat), metric, class_filter = cl)
    plot_network(safe_log1p(sm_class$mat), metric, class_filter = cl)
  }
}

# -----------------------------
# 7. Region-level limma contrasts
# -----------------------------
run_limma_region_contrasts <- function(data, metric) {
  metric_sym <- sym(metric)

  df <- data %>%
    select(SampleID, Condition, Class, Annotation, Abbreviation, RegionLabel, RegionKey, value = !!metric_sym) %>%
    mutate(value = safe_log1p(value)) %>%
    filter(!is.na(value), !is.na(Condition))

  region_results <- df %>%
    group_by(Class, Annotation, Abbreviation, RegionLabel, RegionKey) %>%
    group_modify(~ {
      d <- .x

      cond_counts <- table(d$Condition)
      valid_conds <- names(cond_counts)[cond_counts >= 2]
      d <- d %>% filter(as.character(Condition) %in% valid_conds)

      if (n_distinct(d$Condition) < 2 || nrow(d) < 4) {
        return(tibble(
          contrast = names(contrast_definitions),
          logFC = NA_real_,
          AveExpr = NA_real_,
          t = NA_real_,
          P.Value = NA_real_,
          adj.P.Val = NA_real_,
          B = NA_real_,
          n_samples = n_distinct(.x$SampleID),
          note = "insufficient data"
        ))
      }

      d <- d %>% mutate(Condition = droplevels(factor(Condition)))
      design <- model.matrix(~ 0 + Condition, data = d)
      colnames(design) <- levels(d$Condition)

      y <- matrix(d$value, nrow = 1)
      colnames(y) <- d$SampleID
      rownames(y) <- unique(.y$RegionKey)

      fit <- limma::lmFit(y, design)

      possible_contrasts <- contrast_definitions[
        vapply(contrast_definitions, function(expr) {
          all(str_extract_all(expr, "[A-Za-z]+_[A-Za-z]+")[[1]] %in% colnames(design))
        }, logical(1))
      ]

      if (length(possible_contrasts) == 0) {
        return(tibble(
          contrast = names(contrast_definitions),
          logFC = NA_real_,
          AveExpr = NA_real_,
          t = NA_real_,
          P.Value = NA_real_,
          adj.P.Val = NA_real_,
          B = NA_real_,
          n_samples = n_distinct(d$SampleID),
          note = "no possible contrast"
        ))
      }

      cm <- limma::makeContrasts(contrasts = possible_contrasts, levels = design)
      fit2 <- limma::contrasts.fit(fit, cm)
      fit2 <- limma::eBayes(fit2)

      purrr::map_dfr(colnames(cm), function(co) {
        limma::topTable(fit2, coef = co, number = Inf, sort.by = "none") %>%
          rownames_to_column("region_tmp") %>%
          as_tibble() %>%
          mutate(
            contrast = co,
            n_samples = n_distinct(d$SampleID),
            note = NA_character_
          ) %>%
          select(contrast, logFC, AveExpr, t, P.Value, adj.P.Val, B, n_samples, note)
      })
    }) %>%
    ungroup() %>%
    group_by(contrast) %>%
    mutate(adj.P.Val.global = p.adjust(P.Value, method = "BH")) %>%
    ungroup() %>%
    arrange(contrast, adj.P.Val.global, P.Value)

  region_results
}

contrast_tables <- list()
for (metric in metrics_to_analyse) {
  contrast_tables[[metric]] <- run_limma_region_contrasts(long, metric)
}

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
        adj.P.Val.global < 0.05 ~ "FDR < 0.05",
        adj.P.Val.global < 0.10 ~ "FDR < 0.10",
        P.Value < 0.05 ~ "nominal P < 0.05",
        TRUE ~ "ranked exploratory"
      ),
      biological_read = paste0(
        Annotation, " shows ", direction, " for ", metric_label,
        " in ", contrast, " (logFC = ", scales::number(logFC, accuracy = 0.01),
        ", raw P = ", scales::pvalue(P.Value), ")."
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

for (metric in names(candidate_tables)) {
  readr::write_csv(
    candidate_tables[[metric]],
    file.path(out_dir, "tables", paste0("ranked_candidate_regions_", metric, ".csv"))
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

  pdf(file.path(out_dir, "figures", paste0("effect_size_heatmap_logFC_", metric, ".pdf")),
      width = 8, height = max(5, min(18, nrow(mat) * 0.12)))
  pheatmap::pheatmap(
    mat,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
    main = paste0("Regional effect-size map: ", metric, " logFC"),
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
    filter(contrast == contrast_name, !is.na(P.Value), !is.na(logFC)) %>%
    mutate(
      neglog10p = -log10(P.Value),
      significant = adj.P.Val.global < 0.10,
      label = if_else(significant | rank(P.Value) <= 10, Annotation, NA_character_)
    )

  if (nrow(d) < 3) return(NULL)

  p <- ggplot(d, aes(x = logFC, y = neglog10p)) +
    geom_hline(yintercept = -log10(0.05), linewidth = 0.25, linetype = "dashed") +
    geom_vline(xintercept = 0, linewidth = 0.25) +
    geom_point(aes(fill = Class, shape = significant), alpha = 0.85, size = 2.2, colour = "grey20", stroke = 0.2) +
    ggrepel::geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 30) +
    scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 24)) +
    theme_nature(base_size = 8) +
    labs(
      title = paste0(metric, ": ", contrast_name),
      subtitle = "Labelled regions are FDR < 0.10 or among the 10 lowest raw P values",
      x = "logFC on log1p scale",
      y = "-log10 raw P"
    )

  save_figure(p, paste0("volcano_", metric, "_", contrast_name), width = 6, height = 5)
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
    facet_wrap(~ RegionLabel, scales = "free_y", ncol = 4) +
    scale_colour_manual(values = condition_colors, drop = FALSE) +
    scale_fill_manual(values = condition_colors, drop = FALSE) +
    theme_nature(base_size = 8) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 6)
    ) +
    labs(
      title = paste0("Top region profiles: ", metric, " | ", contrast_name),
      x = NULL,
      y = paste0("log1p ", metric)
    )

  save_figure(p, paste0("top_region_profiles_", metric, "_", contrast_name), width = 10, height = 8)
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

plot_class_profiles <- function(metric) {
  d <- class_level_long %>% filter(Metric == metric)

  p <- ggplot(d, aes(x = Condition, y = ClassMean)) +
    geom_boxplot(aes(colour = Condition), outlier.shape = NA, width = 0.55, linewidth = 0.25) +
    geom_jitter(aes(fill = Condition), width = 0.12, size = 1.1, alpha = 0.75, shape = 21, colour = "white", stroke = 0.15) +
    facet_wrap(~ Class, scales = "free_y", ncol = 4) +
    scale_colour_manual(values = condition_colors, drop = FALSE) +
    scale_fill_manual(values = condition_colors, drop = FALSE) +
    theme_nature(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = paste0("Class-level tracing profiles: ", metric),
      x = NULL,
      y = paste0("Mean log1p ", metric, " across regions")
    )

  save_figure(p, paste0("class_level_profiles_", metric), width = 10, height = 7)
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

  save_figure(p2, paste0("UMAP_systems_structure_", metric), width = 5.5, height = 4.5)

  openxlsx::write.xlsx(
    list(PCA = pca_df, UMAP = umap_df),
    file.path(out_dir, "tables", paste0("dimensionality_outputs_", metric, ".xlsx")),
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

  save_figure(p, paste0("CeM_centered_connectivity_", metric), width = 9, height = 7)
}

purrr::walk(metrics_to_analyse, ~ run_cem_connectivity(long, .x))

# -----------------------------
# 14. Nature-style main figure
# -----------------------------
main_fig_dir <- file.path(out_dir, "figures", "main_figure_nature")
main_tab_dir <- file.path(out_dir, "tables", "main_figure_nature")
dir.create(main_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(main_tab_dir, recursive = TRUE, showWarnings = FALSE)

nature_condition_levels <- c("VEH_paired", "VEH_unpaired", "CNO_paired", "CNO_unpaired")
nature_condition_colors <- c(
  VEH_paired = "#1F1F1F",
  VEH_unpaired = "#4C956C",
  CNO_paired = "#3B6EA8",
  CNO_unpaired = "#C44E52"
)

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
      plot.title = element_text(face = "bold", size = base_size + 0.5),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.5, colour = "grey15"),
      axis.line = element_line(linewidth = 0.25),
      axis.ticks = element_line(linewidth = 0.25),
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
    labs(x = NULL, y = y_label, title = paste0("CeM ", if_else(metric == "Cell_Count", "activity", "projection signal")))
}

make_activity_projection_scatter <- function() {
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

  readr::write_csv(d, file.path(main_tab_dir, "main_figure_activity_projection_scatter.csv"))

  if (nrow(d) < 3) {
    write_main_warning("main_figure_scatter_warning.txt", "Activity-projection scatter skipped: fewer than 3 regions with paired learning effects.")
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
    rowwise() %>%
    mutate(
      max_abs_effect = {
        row_values <- c_across(all_of(heatmap_columns))
        if (all(is.na(row_values))) NA_real_ else max(abs(row_values), na.rm = TRUE)
      }
    ) %>%
    ungroup() %>%
    mutate(max_abs_effect = if_else(is.finite(max_abs_effect), max_abs_effect, NA_real_)) %>%
    arrange(desc(abs(max_abs_effect)))

  wide <- wide_ranked %>%
    slice_head(n = min(top_n, nrow(wide_ranked))) %>%
    mutate(RowLabel = str_trunc(paste0(Annotation, " (", Abbreviation, ")"), 45))

  if (nrow(wide) == 0) {
    write_main_warning("main_figure_heatmap_warning.txt", "Learning vs stress heatmap skipped: no regions available after filtering central contrasts.")
    return(NULL)
  }

  readr::write_csv(wide, file.path(main_tab_dir, "main_figure_heatmap_matrix.csv"))

  plot_d <- wide %>%
    select(RowLabel, all_of(heatmap_columns)) %>%
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
    scale_fill_gradient2(low = "#3B6EA8", mid = "white", high = "#C44E52",
                         midpoint = 0, limits = c(-max_abs, max_abs), oob = scales::squish,
                         name = "logFC") +
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
    arrange(desc(max_abs_effect))

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

build_condition_network <- function(mat, abs_r_cutoff = 0.60, min_complete_n = 4) {
  features <- colnames(mat)
  if (nrow(mat) < min_complete_n || length(features) < 3) {
    return(list(edges = tibble(), metrics = tibble(), hubs = tibble(), graph = igraph::make_empty_graph()))
  }

  cor_mat <- suppressWarnings(cor(mat, method = "spearman", use = "pairwise.complete.obs"))
  cor_mat[!is.finite(cor_mat)] <- NA_real_
  diag(cor_mat) <- NA_real_

  edges <- combn(features, 2, simplify = FALSE) %>%
    purrr::map_dfr(function(pair) {
      x <- mat[, pair[1]]
      y <- mat[, pair[2]]
      n_pair <- sum(complete.cases(x, y))
      r <- cor_mat[pair[1], pair[2]]
      tibble(
        feature1 = pair[1],
        feature2 = pair[2],
        n_pair = n_pair,
        r = r,
        abs_r = abs(r),
        sign = case_when(r > 0 ~ "positive", r < 0 ~ "negative", TRUE ~ NA_character_),
        edge_present = !is.na(r) & n_pair >= min_complete_n & abs(r) >= abs_r_cutoff
      )
    })

  present_edges <- edges %>%
    filter(edge_present) %>%
    transmute(from = feature1, to = feature2, r, abs_r, sign, n_pair, weight = abs_r)

  graph <- igraph::graph_from_data_frame(present_edges, directed = FALSE, vertices = tibble(name = features))
  graph_density <- igraph::edge_density(graph, loops = FALSE)
  mean_abs_r <- if (nrow(present_edges) > 0) mean(present_edges$abs_r, na.rm = TRUE) else NA_real_

  modularity_value <- NA_real_
  if (igraph::ecount(graph) > 0 && igraph::vcount(graph) > 2) {
    community <- suppressWarnings(igraph::cluster_louvain(graph, weights = igraph::E(graph)$weight))
    modularity_value <- igraph::modularity(community)
  }

  metrics <- tibble(
    n_samples = nrow(mat),
    n_features = length(features),
    n_edges = igraph::ecount(graph),
    possible_edges = length(features) * (length(features) - 1) / 2,
    density = graph_density,
    mean_abs_r = mean_abs_r,
    modularity = modularity_value
  )

  hubs <- tibble(
    Feature = features,
    degree = igraph::degree(graph, v = features),
    strength = igraph::strength(graph, v = features, weights = igraph::E(graph)$weight),
    betweenness = igraph::betweenness(graph, v = features, directed = FALSE, weights = NA),
    eigen_centrality = igraph::eigen_centrality(graph, directed = FALSE, weights = igraph::E(graph)$weight)$vector[features]
  ) %>%
    mutate(
      Metric = str_extract(Feature, "^[^_]+(?:_[^_]+)?"),
      RegionKey = str_remove(Feature, "^[^_]+(?:_[^_]+)?__")
    )

  list(edges = edges, metrics = metrics, hubs = hubs, graph = graph)
}

compare_network_edges <- function(edge_tables, reference_condition = "VEH_paired") {
  d <- edge_tables %>%
    select(NetworkCondition, feature1, feature2, r, abs_r, sign, edge_present) %>%
    mutate(edge_id = paste(pmin(feature1, feature2), pmax(feature1, feature2), sep = " -- ")) %>%
    select(NetworkCondition, edge_id, r, abs_r, sign, edge_present) %>%
    pivot_wider(
      names_from = NetworkCondition,
      values_from = c(r, abs_r, sign, edge_present),
      values_fill = list(edge_present = FALSE)
    )

  compare_conditions <- setdiff(nature_condition_levels, reference_condition)

  for (cond in compare_conditions) {
    ref_present <- paste0("edge_present_", reference_condition)
    cond_present <- paste0("edge_present_", cond)
    ref_sign <- paste0("sign_", reference_condition)
    cond_sign <- paste0("sign_", cond)
    ref_abs <- paste0("abs_r_", reference_condition)
    cond_abs <- paste0("abs_r_", cond)
    rewiring_col <- paste0("rewiring_", cond, "_vs_", reference_condition)
    delta_col <- paste0("delta_abs_r_", cond, "_vs_", reference_condition)

    if (!cond_present %in% names(d)) d[[cond_present]] <- FALSE
    if (!cond_sign %in% names(d)) d[[cond_sign]] <- NA_character_
    if (!cond_abs %in% names(d)) d[[cond_abs]] <- NA_real_

    d[[rewiring_col]] <- case_when(
      d[[ref_present]] & !d[[cond_present]] ~ "lost",
      !d[[ref_present]] & d[[cond_present]] ~ "gained",
      d[[ref_present]] & d[[cond_present]] & d[[ref_sign]] != d[[cond_sign]] ~ "sign switch",
      d[[ref_present]] & d[[cond_present]] ~ "retained",
      TRUE ~ "absent"
    )

    d[[delta_col]] <- d[[cond_abs]] - d[[ref_abs]]
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

make_network_rewiring_figure <- function(data, region_top_n = 14, abs_r_cutoff = 0.60) {
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
      build_condition_network(mat, abs_r_cutoff = abs_r_cutoff, min_complete_n = min_pairwise_n)
    })
    names(network_results) <- network_levels

    metrics <- purrr::imap_dfr(network_results, ~ .x$metrics %>% mutate(NetworkCondition = .y, .before = 1)) %>%
      mutate(NetworkMetric = network_metric, .before = 1)
    edges <- purrr::imap_dfr(network_results, ~ .x$edges %>% mutate(NetworkCondition = .y, .before = 1)) %>%
      mutate(NetworkMetric = network_metric, .before = 1)
    hubs <- purrr::imap_dfr(network_results, ~ .x$hubs %>% mutate(NetworkCondition = .y, .before = 1)) %>%
      mutate(NetworkMetric = network_metric, .before = 1)
    edge_rewiring <- compare_network_edges(edges, reference_condition = "VEH_paired") %>%
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
    mutate(RegionShort = str_trunc(paste0(Annotation, " (", Abbreviation, ")"), 38))

  hubs <- hubs %>%
    left_join(region_lookup, by = "RegionKey") %>%
    group_by(NetworkMetric, NetworkCondition) %>%
    mutate(hub_rank = dense_rank(desc(strength))) %>%
    ungroup()

  readr::write_csv(metrics, file.path(main_tab_dir, "main_figure_network_metrics_by_metric.csv"))
  readr::write_csv(hubs, file.path(main_tab_dir, "main_figure_network_hub_centrality_by_metric.csv"))
  readr::write_csv(edges, file.path(main_tab_dir, "main_figure_network_edges_by_metric.csv"))
  readr::write_csv(edge_rewiring, file.path(main_tab_dir, "main_figure_network_edge_rewiring_by_metric.csv"))

  # Backward-compatible copies now contain the separated, metric-labelled analyses.
  readr::write_csv(metrics, file.path(main_tab_dir, "main_figure_network_metrics.csv"))
  readr::write_csv(hubs, file.path(main_tab_dir, "main_figure_network_hub_centrality.csv"))
  readr::write_csv(edges, file.path(main_tab_dir, "main_figure_network_edges.csv"))
  readr::write_csv(edge_rewiring, file.path(main_tab_dir, "main_figure_network_edge_rewiring.csv"))

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
    scale_fill_gradient(low = "white", high = "#1F1F1F", name = "Hub\nstrength") +
    theme_nature_main(base_size = 7.0) +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1, size = 5.8),
      axis.text.y = element_text(size = 5.7)
    ) +
    labs(title = "Top hubs, separated by signal type", x = NULL, y = NULL)

  rewiring_columns <- grep("^rewiring_.*_vs_VEH_paired$", names(edge_rewiring), value = TRUE)
  rewiring_summary <- edge_rewiring %>%
    pivot_longer(cols = all_of(rewiring_columns), names_to = "Comparison", values_to = "Class") %>%
    group_by(NetworkMetric, Comparison, Class) %>%
    summarise(n = n(), .groups = "drop") %>%
    filter(Class != "absent") %>%
    mutate(
      Comparison = str_remove(Comparison, "^rewiring_"),
      Comparison = str_replace(Comparison, "_vs_VEH_paired$", ""),
      Comparison = factor(Comparison, levels = setdiff(network_levels, "VEH_paired")),
      Class = factor(Class, levels = c("retained", "gained", "lost", "sign switch")),
      NetworkMetric = factor(NetworkMetric, levels = network_metric_levels)
    )

  rewiring_plot <- ggplot(rewiring_summary, aes(x = Class, y = n, fill = Class)) +
    geom_col(width = 0.65, colour = "grey20", linewidth = 0.18) +
    facet_grid(NetworkMetric ~ Comparison, scales = "free_y") +
    scale_fill_manual(values = c(retained = "grey70", gained = "#4C956C", lost = "#C44E52", `sign switch` = "#3B6EA8"), guide = "none") +
    theme_nature_main(base_size = 7.0) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, size = 5.7),
      strip.text = element_text(size = 6.1, face = "bold")
    ) +
    labs(title = "Edge rewiring vs VEH_paired", x = NULL, y = "Edges")

  delta_columns <- grep("^delta_abs_r_.*_vs_VEH_paired$", names(edge_rewiring), value = TRUE)
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
      scale_fill_gradient2(low = "#3B6EA8", mid = "white", high = "#C44E52", midpoint = 0,
                           limits = c(-1, 1), oob = scales::squish, name = "Spearman r") +
      theme_nature_main(base_size = 6.8) +
      theme(
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 35, hjust = 1, size = 5.6),
        axis.text.y = element_text(size = 5.4)
      ) +
      labs(title = "Top rewired edges, separated by signal type", x = NULL, y = NULL)
  } else {
    edge_plot <- ggplot() +
      annotate("text", x = 0, y = 0, label = "No thresholded edge rewiring detected", size = 2.3) +
      theme_void(base_size = 7.5) +
      labs(title = "Top rewired edges")
  }

  network_plot <- (metric_plot / hub_plot / rewiring_plot / edge_plot) +
    patchwork::plot_layout(heights = c(1.05, 1.65, 1.25, 1.55), guides = "collect") +
    patchwork::plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 9), legend.position = "bottom")

  safe_main_ggsave(
    file.path(main_fig_dir, "main_figure_network_rewiring_separate_cell_count_intensity.pdf"),
    network_plot,
    width = 180 / 25.4,
    height = 285 / 25.4,
    units = "in"
  )
  safe_main_ggsave(
    file.path(main_fig_dir, "main_figure_network_rewiring_separate_cell_count_intensity.svg"),
    network_plot,
    width = 180 / 25.4,
    height = 285 / 25.4,
    units = "in",
    device = svglite::svglite
  )
  safe_main_ggsave(
    file.path(main_fig_dir, "main_figure_network_rewiring_separate_cell_count_intensity.png"),
    network_plot,
    width = 180 / 25.4,
    height = 285 / 25.4,
    units = "in"
  )

  # Deprecated compatibility filenames, now with separated cell-count/intensity networks.
  safe_main_ggsave(
    file.path(main_fig_dir, "main_figure_network_rewiring_four_conditions.pdf"),
    network_plot,
    width = 180 / 25.4,
    height = 285 / 25.4,
    units = "in"
  )
  safe_main_ggsave(
    file.path(main_fig_dir, "main_figure_network_rewiring_four_conditions.svg"),
    network_plot,
    width = 180 / 25.4,
    height = 285 / 25.4,
    units = "in",
    device = svglite::svglite
  )
  safe_main_ggsave(
    file.path(main_fig_dir, "main_figure_network_rewiring_four_conditions.png"),
    network_plot,
    width = 180 / 25.4,
    height = 285 / 25.4,
    units = "in"
  )

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

panel_d <- make_activity_projection_scatter()
panel_e <- make_learning_stress_heatmap(top_n = 12)
panel_f <- make_main_pca(long)
panel_network <- make_network_rewiring_figure(long, region_top_n = 14, abs_r_cutoff = 0.60)

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
    file.path(main_fig_dir, "main_figure_learning_vs_stress_activity_projection.pdf"),
    main_plot,
    width = 180 / 25.4,
    height = 255 / 25.4,
    units = "in"
  )
  safe_main_ggsave(
    file.path(main_fig_dir, "main_figure_learning_vs_stress_activity_projection.svg"),
    main_plot,
    width = 180 / 25.4,
    height = 255 / 25.4,
    units = "in",
    device = svglite::svglite
  )
  safe_main_ggsave(
    file.path(main_fig_dir, "main_figure_learning_vs_stress_activity_projection.png"),
    main_plot,
    width = 180 / 25.4,
    height = 255 / 25.4,
    units = "in",
    dpi = 600
  )

} else {
  write_main_warning("main_figure_warning.txt", "Nature-style main figure was not generated because every panel was skipped.")
}

# -----------------------------
# 15. Manuscript-matched Fig. 4 panels
# -----------------------------
fig4_fig_dir <- file.path(out_dir, "figures", "fig4_manuscript_matched")
fig4_tab_dir <- file.path(out_dir, "tables", "fig4_manuscript_matched")
dir.create(fig4_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig4_tab_dir, recursive = TRUE, showWarnings = FALSE)

fig4_main_contrasts <- c(
  "Learning_effect",
  "CeM_manipulation_during_learning",
  "CeM_manipulation_during_stress",
  "Learning_x_CeM_interaction"
)

fig4_cov_abs_r_cutoff <- 0.70
fig4_cov_fdr_cutoff <- 0.10
fig4_cov_min_n_pair <- 5
fig4_max_heatmap_regions <- 28
fig4_min_estimable_contrasts_total <- 2
fig4_profile_region_n <- 8
fig4_condition_labels <- c(
  VEH_paired = "VEH_paired",
  VEH_unpaired = "VEH_unpaired",
  CNO_paired = "CNO_paired",
  CNO_unpaired = "CNO_unpaired"
)
fig4_contrast_labels <- c(
  Learning_effect = "VEH_paired -\nVEH_unpaired",
  CeM_manipulation_during_learning = "CNO_paired -\nVEH_paired",
  CeM_manipulation_during_stress = "CNO_unpaired -\nVEH_unpaired",
  Learning_x_CeM_interaction = "(CNO_paired - VEH_paired) -\n(CNO_unpaired - VEH_unpaired)"
)

fig4_condition_colors <- nature_condition_colors
fig4_skipped_panels <- character()

fig4_theme <- function(base_size = 7.2) {
  theme_nature_main(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = base_size + 0.7),
      plot.subtitle = element_text(size = base_size - 0.3, colour = "grey35"),
      strip.text = element_text(face = "bold", size = base_size - 0.6),
      legend.key.height = grid::unit(3.0, "mm"),
      legend.key.width = grid::unit(3.6, "mm")
    )
}

save_fig4_plot <- function(plot, filename_base, width, height, dpi = 600) {
  if (is.null(plot)) return(invisible(NULL))
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
  invisible(plot)
}

write_fig4_table <- function(x, filename_base) {
  readr::write_csv(x, file.path(fig4_tab_dir, paste0(filename_base, ".csv")))
  openxlsx::write.xlsx(x, file.path(fig4_tab_dir, paste0(filename_base, ".xlsx")), overwrite = TRUE)
  invisible(x)
}

fig4_note_skip <- function(panel, reason) {
  fig4_skipped_panels <<- c(fig4_skipped_panels, paste0(panel, ": ", reason))
  invisible(NULL)
}

fig4_clean_region_label <- function(annotation, abbreviation, max_width = 34) {
  label <- if_else(
    !is.na(annotation) & annotation != "" & !is.na(abbreviation) & abbreviation != "",
    paste0(annotation, " (", abbreviation, ")"),
    coalesce(annotation, abbreviation, "Unknown")
  )
  str_trunc(label, max_width)
}

fig4_region_catalog <- long %>%
  distinct(RegionKey, Class, Annotation, Abbreviation, RegionLabel) %>%
  mutate(
    Annotation_lower = str_to_lower(coalesce(Annotation, "")),
    Abbreviation_lower = str_to_lower(coalesce(Abbreviation, "")),
    RegionLabel_lower = str_to_lower(coalesce(RegionLabel, "")),
    RegionShort = fig4_clean_region_label(Annotation, Abbreviation)
  )

fig4_prior_terms <- c(
  "ECT", "LA", "PVi", "PVa", "PVH", "AV", "SS", "RSP", "AVP", "GU",
  "EPv", "GPe", "GPi", "CEAl", "CeM", "PIR", "SF", "MEPO", "MEPA",
  "PA", "VMPO", "DMH", "RT", "SFO", "STN", "PALd"
)

fig4_match_prior_term <- function(term, catalog = fig4_region_catalog) {
  term_lower <- str_to_lower(term)
  if (term_lower == "ss") {
    matched <- catalog %>% filter(str_starts(Annotation_lower, "ss") | str_detect(Abbreviation_lower, "somatosensory"))
  } else if (term_lower == "cem") {
    matched <- catalog %>% filter(
      Annotation_lower %in% c("cem", "ceam") |
        Abbreviation_lower %in% c("cem", "ceam", "cea-m", "cea m") |
        str_detect(RegionLabel_lower, "central amygdalar nucleus, medial|central amygdalar nucleus medial")
    )
  } else {
    matched <- catalog %>%
      filter(
        Annotation_lower == term_lower |
          str_starts(Annotation_lower, paste0(term_lower, "-")) |
          str_starts(Annotation_lower, paste0(term_lower, "_")) |
          Abbreviation_lower == term_lower |
          str_detect(Abbreviation_lower, fixed(term_lower))
      )
  }

  if (nrow(matched) == 0) {
    tibble(
      prior_term = term,
      prior_found = FALSE,
      RegionKey = NA_character_,
      Class = NA_character_,
      Annotation = NA_character_,
      Abbreviation = NA_character_,
      RegionLabel = NA_character_
    )
  } else {
    matched %>%
      arrange(Annotation, RegionKey) %>%
      mutate(prior_term = term, prior_found = TRUE, .before = 1) %>%
      select(prior_term, prior_found, RegionKey, Class, Annotation, Abbreviation, RegionLabel)
  }
}

fig4_prior_matches <- purrr::map_dfr(fig4_prior_terms, fig4_match_prior_term)
write_fig4_table(fig4_prior_matches, "fig4_region_prior_matching")

fig4_effect_rank <- purrr::imap_dfr(contrast_tables, function(res, metric) {
  res %>%
    filter(contrast %in% fig4_main_contrasts) %>%
    mutate(Metric = metric) %>%
    select(Metric, contrast, Class, Annotation, Abbreviation, RegionLabel, RegionKey, logFC, P.Value, adj.P.Val.global)
}) %>%
  group_by(RegionKey, Class, Annotation, Abbreviation, RegionLabel) %>%
  summarise(
    n_estimable_contrasts_total = sum(!is.na(logFC)),
    n_estimable_cell_count = sum(Metric == "Cell_Count" & !is.na(logFC)),
    n_estimable_intensity = sum(Metric == "Intensity" & !is.na(logFC)),
    max_abs_effect = {
      if (any(!is.na(logFC))) max(abs(logFC), na.rm = TRUE) else NA_real_
    },
    min_raw_p = {
      if (any(!is.na(P.Value))) suppressWarnings(min(P.Value, na.rm = TRUE)) else NA_real_
    },
    .groups = "drop"
  ) %>%
  mutate(
    max_abs_effect = if_else(is.finite(max_abs_effect), max_abs_effect, NA_real_),
    min_raw_p = if_else(is.finite(min_raw_p), min_raw_p, NA_real_)
  ) %>%
  arrange(desc(n_estimable_contrasts_total), desc(max_abs_effect), min_raw_p)

fig4_prior_representatives <- fig4_prior_matches %>%
  filter(prior_found, !is.na(RegionKey)) %>%
  left_join(
    fig4_effect_rank %>%
      select(RegionKey, n_estimable_contrasts_total, n_estimable_cell_count, n_estimable_intensity, max_abs_effect, min_raw_p),
    by = "RegionKey"
  ) %>%
  mutate(
    n_estimable_contrasts_total = coalesce(n_estimable_contrasts_total, 0L),
    max_abs_effect = coalesce(max_abs_effect, 0)
  ) %>%
  group_by(prior_term) %>%
  arrange(desc(n_estimable_contrasts_total), desc(max_abs_effect), RegionKey, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

fig4_prior_region_keys <- fig4_prior_representatives %>%
  pull(RegionKey) %>%
  unique()

fig4_prior_match_display_decisions <- fig4_prior_matches %>%
  left_join(
    fig4_prior_representatives %>%
      transmute(prior_term, RegionKey, selected_as_representative = TRUE),
    by = c("prior_term", "RegionKey")
  ) %>%
  mutate(
    selected_as_representative = coalesce(selected_as_representative, FALSE),
    display_decision = case_when(
      !prior_found ~ "not found in data",
      selected_as_representative ~ "plotted representative",
      TRUE ~ "not plotted; alternative match for same manuscript term"
    )
  )
write_fig4_table(fig4_prior_match_display_decisions, "fig4_region_prior_match_display_decisions")

fig4_ranked_region_keys <- fig4_effect_rank %>%
  filter(
    !RegionKey %in% fig4_prior_region_keys,
    n_estimable_contrasts_total >= fig4_min_estimable_contrasts_total
  ) %>%
  slice_head(n = max(fig4_max_heatmap_regions - length(fig4_prior_region_keys), 0)) %>%
  pull(RegionKey)

fig4_region_selection <- bind_rows(
  tibble(RegionKey = fig4_prior_region_keys, manuscript_prioritized = TRUE),
  tibble(RegionKey = fig4_ranked_region_keys, manuscript_prioritized = FALSE)
) %>%
  distinct(RegionKey, .keep_all = TRUE) %>%
  left_join(fig4_region_catalog, by = "RegionKey") %>%
  left_join(
    fig4_effect_rank %>%
      select(RegionKey, n_estimable_contrasts_total, n_estimable_cell_count, n_estimable_intensity, max_abs_effect, min_raw_p),
    by = "RegionKey"
  ) %>%
  mutate(
    effect_ranked = RegionKey %in% fig4_ranked_region_keys,
    inclusion_reason = case_when(
      manuscript_prioritized & effect_ranked ~ "manuscript-prioritized; effect-ranked",
      manuscript_prioritized ~ "manuscript-prioritized",
      effect_ranked ~ "effect-ranked",
      TRUE ~ "selected"
    ),
    Class = if_else(is.na(Class) | Class == "", "Unknown", Class),
    RegionShort = fig4_clean_region_label(Annotation, Abbreviation),
    plot_order = row_number()
  ) %>%
  arrange(Class, desc(manuscript_prioritized), desc(max_abs_effect), RegionShort) %>%
  mutate(
    RegionShort = make.unique(RegionShort, sep = " "),
    plot_order = row_number()
  )

write_fig4_table(fig4_region_selection, "fig4_region_selection_documentation")

fig4_heatmap_source <- purrr::imap_dfr(contrast_tables, function(res, metric) {
  res %>%
    filter(contrast %in% fig4_main_contrasts, RegionKey %in% fig4_region_selection$RegionKey) %>%
    mutate(Metric = metric) %>%
    select(Metric, contrast, Class, Annotation, Abbreviation, RegionLabel, RegionKey, logFC, P.Value, adj.P.Val.global)
}) %>%
  left_join(fig4_region_selection %>% select(RegionKey, RegionShort, inclusion_reason, plot_order), by = "RegionKey") %>%
  mutate(
    contrast_estimable = !is.na(logFC),
    not_estimable_reason = if_else(contrast_estimable, NA_character_, "Contrast not estimable from available non-missing values")
  )

fig4_metric_availability <- long %>%
  filter(RegionKey %in% fig4_region_selection$RegionKey) %>%
  select(SampleID, Condition, RegionKey, all_of(metrics_to_analyse)) %>%
  pivot_longer(cols = all_of(metrics_to_analyse), names_to = "Metric", values_to = "RawValue") %>%
  group_by(Metric, RegionKey, Condition) %>%
  summarise(
    n_samples_total = n_distinct(SampleID),
    n_non_missing = n_distinct(SampleID[!is.na(RawValue)]),
    n_missing = n_samples_total - n_non_missing,
    pct_missing = 100 * n_missing / n_samples_total,
    .groups = "drop"
  ) %>%
  left_join(fig4_region_selection %>% select(RegionKey, RegionShort, inclusion_reason), by = "RegionKey") %>%
  mutate(Condition = factor(as.character(Condition), levels = nature_condition_levels)) %>%
  arrange(Metric, RegionShort, Condition)

write_fig4_table(fig4_heatmap_source, "fig4BC_effect_map_source")
write_fig4_table(fig4_metric_availability, "fig4BC_region_condition_availability")

make_fig4_effect_heatmap <- function(metric, title_label) {
  d <- fig4_heatmap_source %>%
    filter(Metric == metric) %>%
    mutate(
      contrast = factor(contrast, levels = fig4_main_contrasts),
      RowLabel = factor(RegionShort, levels = rev(fig4_region_selection$RegionShort))
    )

  if (nrow(d) == 0) {
    fig4_note_skip(paste0("Panel ", if_else(metric == "Intensity", "B", "C")), paste0("No ", metric, " contrast rows available."))
    return(NULL)
  }

  max_abs <- max(abs(d$logFC), na.rm = TRUE)
  max_abs <- ifelse(is.finite(max_abs) && max_abs > 0, max_abs, 1)

  ggplot(d, aes(x = contrast, y = RowLabel, fill = logFC)) +
    geom_tile(colour = "white", linewidth = 0.18) +
    facet_grid(Class ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient2(
      low = "#3B6EA8", mid = "white", high = "#C44E52",
      midpoint = 0, limits = c(-max_abs, max_abs), oob = scales::squish,
      name = "logFC",
      na.value = "grey88"
    ) +
    scale_x_discrete(
      labels = fig4_contrast_labels,
      drop = FALSE
    ) +
    fig4_theme(base_size = 6.8) +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 5.6),
      strip.text.y = element_text(angle = 0, size = 5.7),
      legend.position = "right"
    ) +
    labs(title = title_label, x = NULL, y = NULL)
}

fig4B <- make_fig4_effect_heatmap("Intensity", "Projection intensity")
fig4C <- make_fig4_effect_heatmap("Cell_Count", "cFos+ cell count")
save_fig4_plot(fig4B, "fig4B_projection_intensity_effect_map", width = 4.8, height = 7.2)
save_fig4_plot(fig4C, "fig4C_cfos_activity_effect_map", width = 4.8, height = 7.2)

fig4_profile_regions <- fig4_region_selection %>%
  filter(str_to_lower(Annotation) %in% str_to_lower(c("ECT", "LA", "PVi", "PVa", "PVH", "AV", "SS", "RSP", "GU", "EPv", "GPe", "GPi")) |
           str_starts(str_to_lower(Annotation), "ss")) %>%
  arrange(desc(max_abs_effect)) %>%
  slice_head(n = fig4_profile_region_n)

if (nrow(fig4_profile_regions) < fig4_profile_region_n) {
  fig4_profile_regions <- bind_rows(
    fig4_profile_regions,
    fig4_region_selection %>%
      filter(!RegionKey %in% fig4_profile_regions$RegionKey) %>%
      arrange(desc(max_abs_effect)) %>%
      slice_head(n = fig4_profile_region_n - nrow(fig4_profile_regions))
  ) %>%
    distinct(RegionKey, .keep_all = TRUE)
}

fig4D_source <- long %>%
  filter(RegionKey %in% fig4_profile_regions$RegionKey) %>%
  select(SampleID, Animal, Group, Condition, Class, Annotation, Abbreviation, RegionLabel, RegionKey, all_of(metrics_to_analyse)) %>%
  pivot_longer(cols = all_of(metrics_to_analyse), names_to = "Metric", values_to = "RawValue") %>%
  mutate(
    Value = safe_log1p(RawValue),
    MetricLabel = recode(Metric, Cell_Count = "cFos+ cells", Intensity = "Projection"),
    Condition = factor(as.character(Condition), levels = nature_condition_levels),
    ConditionLabel = factor(fig4_condition_labels[as.character(Condition)], levels = fig4_condition_labels),
    RegionShort = fig4_clean_region_label(Annotation, Abbreviation, max_width = 24)
  )

fig4D_availability <- fig4D_source %>%
  group_by(Metric, MetricLabel, RegionKey, RegionShort, Condition, ConditionLabel) %>%
  summarise(
    n_samples_total = n_distinct(SampleID),
    n_non_missing = n_distinct(SampleID[!is.na(RawValue)]),
    n_missing = n_samples_total - n_non_missing,
    pct_missing = 100 * n_missing / n_samples_total,
    .groups = "drop"
  )

fig4D_summary <- fig4D_source %>%
  group_by(Metric, MetricLabel, RegionKey, RegionShort, Condition, ConditionLabel) %>%
  summarise(
    mean_value = mean(Value, na.rm = TRUE),
    n = sum(!is.na(Value)),
    sem = sd(Value, na.rm = TRUE) / sqrt(n),
    ci95 = qt(0.975, pmax(n - 1, 1)) * sem,
    .groups = "drop"
  ) %>%
  mutate(
    sem = if_else(is.finite(sem), sem, 0),
    ci95 = if_else(is.finite(ci95), ci95, 0),
    ymin = mean_value - ci95,
    ymax = mean_value + ci95
  )

write_fig4_table(fig4D_source, "fig4D_key_region_condition_profiles_source")
write_fig4_table(fig4D_summary, "fig4D_key_region_condition_profiles_summary")
write_fig4_table(fig4D_availability, "fig4D_key_region_condition_profiles_availability")

fig4D <- if (nrow(fig4D_source) > 0) {
  ggplot(fig4D_source, aes(x = ConditionLabel, y = Value)) +
    geom_point(aes(fill = Condition), shape = 21, size = 1.25, colour = "grey20", stroke = 0.15,
               alpha = 0.78, position = position_jitter(width = 0.10, height = 0, seed = 1)) +
    geom_errorbar(
      data = fig4D_summary,
      aes(x = ConditionLabel, y = mean_value, ymin = ymin, ymax = ymax, colour = Condition),
      width = 0.10, linewidth = 0.25, inherit.aes = FALSE
    ) +
    geom_point(
      data = fig4D_summary,
      aes(y = mean_value, fill = Condition),
      shape = 23, size = 1.75, colour = "white", stroke = 0.18
    ) +
    facet_grid(MetricLabel ~ RegionShort, scales = "free_y") +
    scale_fill_manual(values = fig4_condition_colors, drop = FALSE) +
    scale_colour_manual(values = fig4_condition_colors, guide = "none", drop = FALSE) +
    fig4_theme(base_size = 6.6) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 5.3),
      strip.text.x = element_text(size = 5.6),
      legend.position = "bottom"
    ) +
    labs(title = "Key region condition profiles", x = NULL, y = "log1p signal")
} else {
  fig4_note_skip("Panel D", "No key profile regions were available.")
  NULL
}
save_fig4_plot(fig4D, "fig4D_key_region_condition_profiles", width = 9.0, height = 4.4)

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
  left_join(fig4_region_selection %>% select(RegionKey, manuscript_prioritized, inclusion_reason), by = "RegionKey") %>%
  mutate(
    manuscript_prioritized = coalesce(manuscript_prioritized, FALSE),
    combined_abs_effect = sqrt(Cell_Count_logFC^2 + Intensity_logFC^2),
    plotted = !is.na(Cell_Count_logFC) & !is.na(Intensity_logFC),
    not_plotted_reason = case_when(
      plotted ~ NA_character_,
      is.na(Cell_Count_logFC) & is.na(Intensity_logFC) ~ "Missing both Cell_Count and Intensity VEH_paired - VEH_unpaired logFC",
      is.na(Cell_Count_logFC) ~ "Missing Cell_Count VEH_paired - VEH_unpaired logFC",
      is.na(Intensity_logFC) ~ "Missing Intensity VEH_paired - VEH_unpaired logFC",
      TRUE ~ NA_character_
    ),
    RegionShort = fig4_clean_region_label(Annotation, Abbreviation, max_width = 24),
    label = if_else(
      plotted & (manuscript_prioritized | rank(-combined_abs_effect, ties.method = "first", na.last = "keep") <= 12),
      RegionShort,
      NA_character_
    )
  )

write_fig4_table(fig4E_source, "fig4E_activity_projection_dissociation_source")

fig4E <- if (sum(complete.cases(fig4E_source$Cell_Count_logFC, fig4E_source$Intensity_logFC)) >= 3) {
  ggplot(fig4E_source %>% filter(plotted), aes(x = Cell_Count_logFC, y = Intensity_logFC)) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey55") +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey55") +
    geom_point(aes(size = combined_abs_effect, fill = manuscript_prioritized),
               shape = 21, colour = "grey15", stroke = 0.18, alpha = 0.84) +
    ggrepel::geom_text_repel(
      data = fig4E_source %>% filter(!is.na(label)),
      aes(label = label),
      size = 1.9, min.segment.length = 0, segment.size = 0.16,
      max.overlaps = Inf, seed = 4
    ) +
    scale_fill_manual(values = c(`FALSE` = "grey72", `TRUE` = "#C44E52"), name = "Prior", labels = c("No", "Yes")) +
    scale_size_continuous(range = c(0.8, 3.4), guide = "none") +
    fig4_theme(base_size = 7.0) +
    theme(legend.position = "bottom") +
    labs(
      title = "Activity-projection dissociation",
      x = "VEH_paired - VEH_unpaired logFC, cFos+ cell count",
      y = "VEH_paired - VEH_unpaired logFC, projection intensity"
    )
} else {
  fig4_note_skip("Panel E", "Fewer than 3 regions had paired Cell_Count and Intensity VEH_paired - VEH_unpaired effects.")
  NULL
}
save_fig4_plot(fig4E, "fig4E_activity_projection_dissociation", width = 4.7, height = 4.2)

fig4_seed_terms <- c("CEAl", "CeM", "LA", "PVH")
fig4_seed_candidates <- purrr::map_dfr(fig4_seed_terms, fig4_match_prior_term) %>%
  filter(prior_found, !is.na(RegionKey)) %>%
  group_by(prior_term) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(SeedLabel = prior_term)
write_fig4_table(fig4_seed_candidates, "fig4F_seed_region_matching")

fig4_covariation_results <- tibble()
fig4_covariation_warnings <- tibble(seed = character(), condition = character(), reason = character())

if (nrow(fig4_seed_candidates) > 0) {
  sm_cell_fig4 <- make_sample_matrix(long, "Cell_Count")
  mat_cell_fig4 <- safe_log1p(sm_cell_fig4$mat)
  meta_cell_fig4 <- sm_cell_fig4$annotation %>%
    mutate(Condition = factor(as.character(Condition), levels = nature_condition_levels))

  candidate_targets <- unique(c(fig4_region_selection$RegionKey, colnames(mat_cell_fig4)))

  fig4_covariation_results <- purrr::pmap_dfr(
    fig4_seed_candidates %>% select(SeedLabel, SeedRegionKey = RegionKey),
    function(SeedLabel, SeedRegionKey) {
      purrr::map_dfr(nature_condition_levels, function(cond) {
        if (!SeedRegionKey %in% colnames(mat_cell_fig4)) {
          fig4_covariation_warnings <<- bind_rows(
            fig4_covariation_warnings,
            tibble(seed = SeedLabel, condition = cond, reason = "Seed region not present in Cell_Count matrix.")
          )
          return(tibble())
        }

        ids <- meta_cell_fig4 %>% filter(Condition == cond) %>% pull(SampleID)
        submat <- mat_cell_fig4[rownames(mat_cell_fig4) %in% ids, , drop = FALSE]
        if (nrow(submat) < min_pairwise_n) {
          fig4_covariation_warnings <<- bind_rows(
            fig4_covariation_warnings,
            tibble(seed = SeedLabel, condition = cond, reason = paste0("Fewer than ", min_pairwise_n, " animals."))
          )
          return(tibble())
        }

        tibble(TargetRegionKey = intersect(candidate_targets, colnames(submat))) %>%
          filter(TargetRegionKey != SeedRegionKey) %>%
          rowwise() %>%
          mutate(
            SeedLabel = SeedLabel,
            SeedRegionKey = SeedRegionKey,
            Condition = cond,
            n_pair = sum(is.finite(submat[, SeedRegionKey]) & is.finite(submat[, TargetRegionKey])),
            rho = {
              if (n_pair >= min_pairwise_n) {
                suppressWarnings(cor(submat[, SeedRegionKey], submat[, TargetRegionKey], method = "spearman", use = "pairwise.complete.obs"))
              } else {
                NA_real_
              }
            },
            p_value = {
              if (n_pair >= min_pairwise_n && !is.na(rho)) {
                suppressWarnings(
                  tryCatch(
                    cor.test(submat[, SeedRegionKey], submat[, TargetRegionKey], method = "spearman", exact = FALSE)$p.value,
                    error = function(e) NA_real_
                  )
                )
              } else {
                NA_real_
              }
            }
          ) %>%
          ungroup() %>%
          group_by(SeedLabel, SeedRegionKey, Condition) %>%
          mutate(fdr = p.adjust(p_value, method = "BH")) %>%
          ungroup()
      })
    }
  ) %>%
    left_join(
      fig4_region_catalog %>%
        select(TargetRegionKey = RegionKey, TargetClass = Class, TargetAnnotation = Annotation,
               TargetAbbreviation = Abbreviation, TargetRegionLabel = RegionLabel),
      by = "TargetRegionKey"
    ) %>%
    left_join(
      fig4_region_catalog %>%
        select(SeedRegionKey = RegionKey, SeedClass = Class, SeedAnnotation = Annotation,
               SeedAbbreviation = Abbreviation, SeedRegionLabel = RegionLabel),
      by = "SeedRegionKey"
    ) %>%
    mutate(
      TargetShort = fig4_clean_region_label(TargetAnnotation, TargetAbbreviation, max_width = 28),
      SeedCondition = paste(SeedLabel, fig4_condition_labels[Condition], sep = " | "),
      target_is_fig4_selected = TargetRegionKey %in% fig4_region_selection$RegionKey,
      covariation_display = n_pair >= fig4_cov_min_n_pair & !is.na(rho),
      covariation_display_reason = case_when(
        covariation_display ~ NA_character_,
        is.na(rho) & n_pair < fig4_cov_min_n_pair ~ paste0("n_pair below display threshold (", fig4_cov_min_n_pair, ")"),
        is.na(rho) ~ "rho not estimable",
        n_pair < fig4_cov_min_n_pair ~ paste0("n_pair below display threshold (", fig4_cov_min_n_pair, ")"),
        TRUE ~ NA_character_
      )
    )
}

if (!all(c("SeedLabel", "SeedRegionKey", "Condition", "TargetRegionKey", "rho", "p_value", "fdr", "n_pair") %in% names(fig4_covariation_results))) {
  fig4_note_skip("Panel F", "No CEAl/CeM, LA, or PVH seed regions were available after flexible matching.")
  fig4_covariation_results <- tibble(
    SeedLabel = character(),
    SeedRegionKey = character(),
    Condition = character(),
    TargetRegionKey = character(),
    n_pair = integer(),
    rho = numeric(),
    p_value = numeric(),
    fdr = numeric(),
    TargetClass = character(),
    TargetAnnotation = character(),
    TargetAbbreviation = character(),
    TargetRegionLabel = character(),
    SeedClass = character(),
    SeedAnnotation = character(),
    SeedAbbreviation = character(),
    SeedRegionLabel = character(),
    TargetShort = character(),
    SeedCondition = character(),
    target_is_fig4_selected = logical(),
    covariation_display = logical(),
    covariation_display_reason = character()
  )
}

write_fig4_table(fig4_covariation_results, "fig4F_seed_based_cfos_covariation_source")
write_fig4_table(fig4_covariation_warnings, "fig4F_seed_based_cfos_covariation_warnings")

fig4_cov_targets <- fig4_covariation_results %>%
  filter(covariation_display) %>%
  group_by(TargetRegionKey, TargetShort) %>%
  summarise(
    target_is_fig4_selected = any(target_is_fig4_selected),
    max_abs_rho = max(abs(rho), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(target_is_fig4_selected), desc(max_abs_rho)) %>%
  slice_head(n = 35) %>%
  mutate(TargetPlotShort = make.unique(TargetShort, sep = " "))

fig4F_source_plot <- fig4_covariation_results %>%
  filter(TargetRegionKey %in% fig4_cov_targets$TargetRegionKey) %>%
  left_join(fig4_cov_targets %>% select(TargetRegionKey, TargetPlotShort), by = "TargetRegionKey") %>%
  mutate(
    rho_display = if_else(covariation_display, rho, NA_real_),
    TargetPlotShort = factor(TargetPlotShort, levels = rev(fig4_cov_targets$TargetPlotShort)),
    SeedCondition = factor(SeedCondition, levels = unique(SeedCondition))
  )

fig4F <- if (nrow(fig4F_source_plot) > 0 && n_distinct(fig4F_source_plot$SeedCondition) > 0) {
  ggplot(fig4F_source_plot, aes(x = SeedCondition, y = TargetPlotShort, fill = rho_display)) +
    geom_tile(colour = "white", linewidth = 0.16) +
      scale_fill_gradient2(low = "#3B6EA8", mid = "white", high = "#C44E52",
                         midpoint = 0, limits = c(-1, 1), oob = scales::squish,
                         name = "rho", na.value = "grey88") +
    fig4_theme(base_size = 6.2) +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 4.9),
      axis.text.y = element_text(size = 5.0),
      legend.position = "right"
    ) +
    labs(title = "Seed-based cFos covariation", x = NULL, y = NULL)
} else {
  fig4_note_skip("Panel F", "No seed-condition covariation results passed the minimum-n requirement.")
  NULL
}
save_fig4_plot(fig4F, "fig4F_seed_based_cfos_covariation", width = 7.4, height = 6.6)

fig4G_source <- fig4_covariation_results %>%
  filter(!is.na(rho), n_pair >= min_pairwise_n) %>%
  mutate(
    edge_pass = n_pair >= fig4_cov_min_n_pair & abs(rho) >= fig4_cov_abs_r_cutoff & fdr <= fig4_cov_fdr_cutoff,
    edge_pass_reason = case_when(
      edge_pass ~ NA_character_,
      n_pair < fig4_cov_min_n_pair ~ paste0("n_pair below edge threshold (", fig4_cov_min_n_pair, ")"),
      abs(rho) < fig4_cov_abs_r_cutoff ~ paste0("abs(rho) below ", fig4_cov_abs_r_cutoff),
      is.na(fdr) | fdr > fig4_cov_fdr_cutoff ~ paste0("FDR above ", fig4_cov_fdr_cutoff),
      TRUE ~ NA_character_
    ),
    edge_sign = case_when(rho > 0 ~ "Positive", rho < 0 ~ "Negative", TRUE ~ "Zero")
  )

fig4G_summary <- fig4G_source %>%
  filter(edge_pass, edge_sign %in% c("Positive", "Negative")) %>%
  group_by(SeedLabel, Condition, edge_sign) %>%
  summarise(
    n_edges = n(),
    mean_rho = mean(rho, na.rm = TRUE),
    signed_covariation_strength = sum(rho, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  complete(
    SeedLabel,
    Condition = nature_condition_levels,
    edge_sign = c("Positive", "Negative"),
    fill = list(n_edges = 0, mean_rho = NA_real_, signed_covariation_strength = 0)
  ) %>%
  mutate(
    Condition = factor(Condition, levels = nature_condition_levels),
    ConditionLabel = factor(fig4_condition_labels[as.character(Condition)], levels = fig4_condition_labels)
  )

fig4G_strength <- fig4G_summary %>%
  group_by(SeedLabel, Condition, ConditionLabel) %>%
  summarise(signed_covariation_strength = sum(signed_covariation_strength, na.rm = TRUE), .groups = "drop")

write_fig4_table(fig4G_source, "fig4G_network_rewiring_edges_source")
write_fig4_table(fig4G_summary, "fig4G_network_rewiring_summary")

fig4G <- if (nrow(fig4G_summary) > 0 && any(fig4G_summary$n_edges > 0)) {
  edge_count_plot <- ggplot(fig4G_summary, aes(x = ConditionLabel, y = n_edges, fill = edge_sign)) +
    geom_col(width = 0.62, colour = "grey20", linewidth = 0.15, position = position_dodge(width = 0.68)) +
    facet_wrap(~ SeedLabel, nrow = 1) +
    scale_fill_manual(values = c(Positive = "#C44E52", Negative = "#3B6EA8"), name = NULL) +
    fig4_theme(base_size = 6.5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.2), legend.position = "bottom") +
    labs(title = "Condition-specific co-activation structure", x = NULL, y = "Thresholded edges")

  strength_plot <- ggplot(fig4G_strength, aes(x = ConditionLabel, y = signed_covariation_strength, colour = SeedLabel, group = SeedLabel)) +
    geom_hline(yintercept = 0, linewidth = 0.22, colour = "grey55") +
    geom_line(linewidth = 0.3) +
    geom_point(size = 1.7) +
    fig4_theme(base_size = 6.5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.2), legend.position = "bottom") +
    labs(title = "Signed covariation strength", x = NULL, y = "sum rho")

  edge_count_plot / strength_plot + patchwork::plot_layout(heights = c(1.15, 1))
} else {
  fig4_note_skip("Panel G", "No seed-based cFos covariation edges passed the rho/FDR threshold.")
  NULL
}
save_fig4_plot(fig4G, "fig4G_network_rewiring_summary", width = 7.4, height = 5.2)

fig4_panel_list <- list(fig4B, fig4C, fig4D, fig4E, fig4F, fig4G)
fig4_panel_list <- fig4_panel_list[!vapply(fig4_panel_list, is.null, logical(1))]

fig4_combined <- NULL
if (length(fig4_panel_list) > 0) {
  fig4_combined <- patchwork::wrap_plots(fig4_panel_list, ncol = 2, guides = "collect") +
    patchwork::plot_annotation(tag_levels = list(c("B", "C", "D", "E", "F", "G"))) &
    theme(plot.tag = element_text(face = "bold", size = 9), legend.position = "bottom")

  save_fig4_plot(fig4_combined, "fig4_main_manuscript_matched_combined", width = 13.5, height = 14.5)
} else {
  fig4_note_skip("Combined Fig. 4", "All manuscript-matched panels were skipped.")
}

fig4_readme <- c(
  "Manuscript-matched Fig. 4 generation",
  "",
  paste0("Input data used: ", paste(sort(unique(long$SourceFile)), collapse = "; ")),
  "",
  "Contrasts visualized:",
  paste0("- ", fig4_main_contrasts, ": ", fig4_contrast_labels[fig4_main_contrasts]),
  "",
  paste0("Region selection: manuscript-prioritized terms were matched flexibly against Annotation and Abbreviation; one best-covered representative per term was plotted, then additional top regions by maximum absolute logFC and at least ", fig4_min_estimable_contrasts_total, " estimable contrasts were added up to approximately ", fig4_max_heatmap_regions, " regions."),
  "Alternative matches for broad manuscript terms are documented but not plotted in the main figure to avoid turning Fig. 4 into a missingness map.",
  paste0("Manuscript-prioritized terms: ", paste(fig4_prior_terms, collapse = ", ")),
  "Missing data handling: missing/non-estimable effect-map cells are kept as NA and rendered grey; no missing region values are imputed for manuscript-facing panels.",
  "Availability tables report n_samples_total, n_non_missing, n_missing, and pct_missing by region, metric, and condition.",
  "",
  paste0("Covariation method: Spearman correlations across animals using Cell_Count only; computed minimum n_pair = ", min_pairwise_n, "."),
  paste0("Covariation heatmap display threshold: n_pair >= ", fig4_cov_min_n_pair, "; lower-n or non-estimable cells are rendered grey and documented in the source table."),
  paste0("Network/co-activation summary threshold: n_pair >= ", fig4_cov_min_n_pair, ", abs(rho) >= ", fig4_cov_abs_r_cutoff, " and FDR <= ", fig4_cov_fdr_cutoff, "."),
  "",
  "Terminology: panels describe regional cFos covariation and condition-specific co-activation structure, not causal connectivity.",
  "",
  "Skipped panels or warnings:",
  if (length(fig4_skipped_panels) == 0) "- None" else paste0("- ", fig4_skipped_panels)
)
writeLines(fig4_readme, con = file.path(fig4_tab_dir, "README_fig4_manuscript_matched.txt"))
writeLines(fig4_readme, con = file.path(fig4_fig_dir, "README_fig4_manuscript_matched.txt"))

# -----------------------------
# 16. Save session info
# -----------------------------
sink(file.path(out_dir, "session_info.txt"))
print(sessionInfo())
sink()

message("Done. Outputs written to: ", out_dir)
