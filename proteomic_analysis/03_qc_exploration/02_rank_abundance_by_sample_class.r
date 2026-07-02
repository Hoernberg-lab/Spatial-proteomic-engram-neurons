# ==============================================================================
# Rank-abundance and marker-validation QC by Neha sample class
# ==============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(tidyverse)
  library(scales)
  library(ggrepel)
  library(writexl)
  library(svglite)
})

input_dir <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Datasets/gct/data/imputed"
saving_dir <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Results/QC/rank_abundance"
dir.create(saving_dir, recursive = TRUE, showWarnings = FALSE)

sample_classes <- c("mcherry", "neuropil", "cfos", "neuron")

sample_class_files <- list(
  mcherry = file.path(input_dir, "pgmatrix_imputed_mcherry.xlsx"),
  neuropil = file.path(input_dir, "pgmatrix_imputed_neuropil.xlsx"),
  cfos = file.path(input_dir, "pgmatrix_imputed_cfos.xlsx"),
  neuron = file.path(input_dir, "pgmatrix_imputed_neuron.xlsx")
)

resolve_sample_class_file <- function(sample_class, configured_path) {
  if (file.exists(configured_path)) return(configured_path)

  candidates <- list.files(
    input_dir,
    pattern = paste0(sample_class, ".*\\.xlsx$"),
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (length(candidates) == 0) {
    stop("No input file found for sample_class '", sample_class, "'. Expected: ", configured_path)
  }

  candidates[order(file.info(candidates)$mtime, decreasing = TRUE)][1]
}

sample_class_files <- purrr::imap(sample_class_files, ~ resolve_sample_class_file(.y, .x))

sample_class_colors <- c(
  mcherry = "#4C78A8",
  neuropil = "#72B7B2",
  cfos = "#E45756",
  neuron = "#54A24B"
)

marker_sets <- list(
  neuropil = c("Vcan", "Kcna1", "Mog", "Cntnap1", "Cnp"),
  neuron = c("Rps9", "Rpl22", "Rps16", "Brd4", "Rpl35", "H1-1", "H1-3", "Rps18")
)

marker_set_labels <- c(
  mcherry = "mCherry marker-defined sample class",
  neuropil = "Neuropil marker-defined sample class",
  cfos = "cFos marker-defined sample class",
  neuron = "Neuron marker-defined sample class"
)

marker_colors <- c(sample_class_colors, unassigned = "grey40")

standardize_expression_table <- function(df) {
  gene_col <- intersect(c("Genes", "Gene", "gene_symbol", "Protein.Names", "T: Protein.Names", "id"), names(df))[1]
  if (is.na(gene_col)) stop("Could not identify a gene/protein identifier column.")

  df %>%
    dplyr::rename(Genes = dplyr::all_of(gene_col)) %>%
    dplyr::mutate(Genes = as.character(Genes))
}

read_sample_class_long <- function(file_path, sample_class) {
  readxl::read_excel(file_path) %>%
    standardize_expression_table() %>%
    tidyr::pivot_longer(
      cols = where(is.numeric),
      names_to = "Sample",
      values_to = "Log2Intensity"
    ) %>%
    dplyr::filter(!is.na(Log2Intensity), !is.na(Genes), Genes != "") %>%
    dplyr::mutate(sample_class = sample_class)
}

make_rank_data <- function(long_df) {
  long_df %>%
    dplyr::group_by(sample_class, Genes) %>%
    dplyr::summarise(
      MeanLog2 = mean(Log2Intensity, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(LinearValue = 2^MeanLog2) %>%
    dplyr::group_by(sample_class) %>%
    dplyr::arrange(dplyr::desc(LinearValue), .by_group = TRUE) %>%
    dplyr::mutate(Rank = dplyr::row_number()) %>%
    dplyr::ungroup()
}

annotate_markers <- function(rank_data) {
  marker_tbl <- purrr::imap_dfr(marker_sets, function(markers, marker_set) {
    tibble::tibble(Genes = markers, MarkerSet = marker_set)
  }) %>%
    dplyr::distinct()

  rank_data %>%
    dplyr::left_join(marker_tbl, by = "Genes", relationship = "many-to-many") %>%
    dplyr::mutate(
      MarkerSet = tidyr::replace_na(MarkerSet, "None"),
      MarkerColor = dplyr::if_else(MarkerSet %in% sample_classes, MarkerSet, "unassigned")
    )
}

save_rank_excel <- function(plot_data, file_name) {
  writexl::write_xlsx(
    list(
      ranks = plot_data %>%
        dplyr::select(sample_class, Rank, Genes, MeanLog2, LinearValue, MarkerSet) %>%
        dplyr::arrange(sample_class, Rank)
    ),
    path = file.path(saving_dir, paste0(file_name, "_processed_ranks.xlsx"))
  )
}

plot_rank_abundance <- function(plot_data, sample_class) {
  label_df <- plot_data %>% dplyr::filter(MarkerSet != "None")

  rank_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Rank, y = LinearValue)) +
    ggplot2::geom_point(alpha = 0.08, size = 0.15, color = "grey82") +
    ggplot2::geom_line(color = sample_class_colors[[sample_class]], alpha = 0.85, linewidth = 0.35) +
    ggplot2::geom_point(
      data = label_df,
      ggplot2::aes(color = MarkerColor),
      size = 1.45,
      alpha = 1
    ) +
    ggrepel::geom_label_repel(
      data = label_df,
      ggplot2::aes(label = Genes, fill = MarkerColor),
      color = "white",
      size = 2.35,
      fontface = "bold",
      family = "sans",
      label.padding = grid::unit(0.15, "lines"),
      label.r = grid::unit(0, "lines"),
      label.size = 0,
      segment.color = "grey20",
      segment.linewidth = 0.25,
      segment.alpha = 0.85,
      box.padding = 0.35,
      point.padding = 0.15,
      force = 10,
      max.overlaps = Inf
    ) +
    ggplot2::scale_y_log10(
      expand = ggplot2::expansion(mult = c(0.05, 0.22)),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.025, 0.02)),
      labels = scales::label_comma()
    ) +
    ggplot2::scale_color_manual(values = marker_colors) +
    ggplot2::scale_fill_manual(values = marker_colors) +
    ggplot2::labs(x = "Protein rank", y = "Intensity (log10)", title = sample_class) +
    ggplot2::theme_minimal(base_size = 8) +
    ggplot2::theme(
      aspect.ratio = 1,
      text = ggplot2::element_text(family = "sans", color = "black"),
      plot.title = ggplot2::element_text(face = "bold", size = 9),
      axis.title = ggplot2::element_text(size = 8, face = "bold"),
      axis.text = ggplot2::element_text(size = 6, color = "black"),
      axis.line = ggplot2::element_line(linewidth = 0.35, color = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.35, color = "black"),
      panel.grid = ggplot2::element_blank(),
      legend.position = "none",
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA)
    )

  file_name <- paste0("rank_abundance_", sample_class)
  svg_path <- file.path(saving_dir, paste0(file_name, ".svg"))

  ggplot2::ggsave(
    svg_path,
    rank_plot,
    width = 120,
    height = 120,
    units = "mm",
    device = svglite::svglite
  )

  save_rank_excel(plot_data, file_name)
  message("Saved plot: ", svg_path)

  rank_plot
}

long_by_class <- purrr::imap_dfr(sample_class_files, read_sample_class_long)
rank_data <- make_rank_data(long_by_class) %>% annotate_markers()

rank_plots <- purrr::imap(sample_classes, function(sample_class, .idx) {
  plot_rank_abundance(
    rank_data %>% dplyr::filter(.data$sample_class == sample_class),
    sample_class
  )
})

marker_summary <- rank_data %>%
  dplyr::filter(MarkerSet != "None") %>%
  dplyr::group_by(sample_class, Genes, MarkerSet) %>%
  dplyr::summarise(
    MeanLog2 = mean(MeanLog2, na.rm = TRUE),
    MedianRank = median(Rank, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    sample_class = factor(sample_class, levels = sample_classes),
    MarkerSet = factor(MarkerSet, levels = sample_classes, labels = marker_set_labels[sample_classes])
  )

marker_summary_plot <- ggplot2::ggplot(
  marker_summary,
  ggplot2::aes(x = sample_class, y = MeanLog2)
) +
  ggplot2::geom_boxplot(
    outlier.shape = NA,
    width = 0.55,
    linewidth = 0.3,
    fill = "grey92",
    color = "black"
  ) +
  ggplot2::geom_point(
    ggplot2::aes(color = sample_class),
    position = ggplot2::position_jitter(width = 0.12, height = 0),
    size = 1.6,
    alpha = 0.9
  ) +
  ggplot2::facet_wrap(~MarkerSet, scales = "free_y", nrow = 1) +
  ggplot2::scale_color_manual(values = sample_class_colors) +
  ggplot2::labs(x = NULL, y = "Mean log2 intensity") +
  ggplot2::theme_classic(base_size = 8) +
  ggplot2::theme(
    text = ggplot2::element_text(family = "sans", color = "black"),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold", size = 7),
    axis.title = ggplot2::element_text(size = 8, face = "bold"),
    axis.text = ggplot2::element_text(size = 6, color = "black"),
    axis.line = ggplot2::element_line(linewidth = 0.3, color = "black"),
    axis.ticks = ggplot2::element_line(linewidth = 0.3, color = "black"),
    legend.position = "none",
    panel.spacing = grid::unit(1, "lines")
  )

ggplot2::ggsave(
  file.path(saving_dir, "marker_abundance_summary.svg"),
  marker_summary_plot,
  width = 180,
  height = 70,
  units = "mm",
  device = svglite::svglite
)

score_all <- purrr::imap_dfr(marker_sets, function(markers, marker_set) {
  long_by_class %>%
    dplyr::filter(Genes %in% markers) %>%
    dplyr::group_by(sample_class, Sample) %>%
    dplyr::summarise(
      Score = mean(Log2Intensity, na.rm = TRUE),
      DetectedMarkers = sum(!is.na(Log2Intensity)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      MarkerSet = marker_set,
      MarkerN = length(markers)
    )
}) %>%
  dplyr::mutate(
    sample_class = factor(sample_class, levels = sample_classes),
    MarkerSet = factor(MarkerSet, levels = sample_classes, labels = marker_set_labels[sample_classes])
  )

global_stats <- score_all %>%
  dplyr::group_by(MarkerSet) %>%
  dplyr::summarise(
    p_value = stats::kruskal.test(Score ~ sample_class)$p.value,
    .groups = "drop"
  )

pairwise_stats <- score_all %>%
  dplyr::group_by(MarkerSet) %>%
  tidyr::nest() %>%
  dplyr::mutate(
    pairwise = purrr::map(data, function(df) {
      pw <- stats::pairwise.wilcox.test(df$Score, df$sample_class, p.adjust.method = "holm")
      as.data.frame(as.table(pw$p.value), stringsAsFactors = FALSE) %>%
        dplyr::rename(group1 = Var1, group2 = Var2, p_adj = Freq) %>%
        dplyr::filter(!is.na(p_adj))
    })
  ) %>%
  dplyr::select(MarkerSet, pairwise) %>%
  tidyr::unnest(pairwise)

module_score_plot <- ggplot2::ggplot(
  score_all,
  ggplot2::aes(x = sample_class, y = Score)
) +
  ggplot2::geom_boxplot(
    ggplot2::aes(fill = sample_class),
    width = 0.48,
    outlier.shape = NA,
    linewidth = 0.32,
    color = "black",
    alpha = 0.92
  ) +
  ggplot2::geom_point(
    ggplot2::aes(color = sample_class),
    position = ggplot2::position_jitter(width = 0.09, height = 0),
    size = 1.05,
    alpha = 0.55,
    stroke = 0
  ) +
  ggplot2::stat_summary(fun = median, geom = "crossbar", width = 0.42, linewidth = 0.28, color = "black") +
  ggplot2::facet_wrap(~MarkerSet, scales = "free_y", nrow = 1) +
  ggplot2::scale_fill_manual(values = sample_class_colors) +
  ggplot2::scale_color_manual(values = sample_class_colors) +
  ggplot2::labs(x = NULL, y = "Marker module score\n(mean log2 intensity)") +
  ggplot2::theme_classic(base_size = 7) +
  ggplot2::theme(
    text = ggplot2::element_text(family = "Arial", color = "black"),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold", size = 7.2, margin = ggplot2::margin(b = 6)),
    axis.title.y = ggplot2::element_text(size = 7.2, face = "bold", margin = ggplot2::margin(r = 5)),
    axis.text.x = ggplot2::element_text(size = 6.2, color = "black", angle = 35, hjust = 1, vjust = 1),
    axis.text.y = ggplot2::element_text(size = 6.2, color = "black"),
    axis.line = ggplot2::element_line(linewidth = 0.3, color = "black"),
    axis.ticks = ggplot2::element_line(linewidth = 0.3, color = "black"),
    axis.ticks.length = grid::unit(1.5, "mm"),
    panel.spacing = grid::unit(1.1, "lines"),
    legend.position = "none",
    plot.margin = ggplot2::margin(4, 4, 4, 4),
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    panel.background = ggplot2::element_rect(fill = "white", color = NA)
  )

ggplot2::ggsave(
  file.path(saving_dir, "marker_module_scores.svg"),
  module_score_plot,
  width = 180,
  height = 75,
  units = "mm",
  device = svglite::svglite
)

writexl::write_xlsx(
  list(
    marker_abundance_summary = marker_summary,
    marker_module_scores = score_all,
    global_kruskal = global_stats,
    pairwise_wilcox = pairwise_stats,
    rank_marker_rows = rank_data %>%
      dplyr::filter(MarkerSet != "None") %>%
      dplyr::arrange(MarkerSet, Genes, sample_class)
  ),
  path = file.path(saving_dir, "marker_validation_summary.xlsx")
)

message("Rank-abundance and marker-validation QC completed.")
