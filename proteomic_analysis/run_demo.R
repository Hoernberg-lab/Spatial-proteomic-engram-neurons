args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1) sub("^--file=", "", file_arg) else "run_demo.R"
repo_dir <- normalizePath(dirname(script_path), winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(repo_dir, "R", "analysis_labels.R"))) {
  repo_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

source(file.path(repo_dir, "R", "analysis_labels.R"))

input_dir <- file.path(repo_dir, "demo", "input")
expected_dir <- file.path(repo_dir, "demo", "expected_output")
output_dir <- file.path(repo_dir, "demo", "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- read.csv(file.path(input_dir, "demo_sample_metadata.csv"), check.names = FALSE, stringsAsFactors = FALSE)
pg_matrix <- read.csv(file.path(input_dir, "demo_pg_matrix.csv"), check.names = FALSE, stringsAsFactors = FALSE)

parsed_sample_class <- parse_sample_class(metadata$sample_id)
parsed_condition_code <- parse_condition_code(metadata$condition_code)

if (!identical(parsed_sample_class, metadata$sample_class)) {
  stop("Demo sample_class parsing did not match metadata.", call. = FALSE)
}
if (!identical(parsed_condition_code, as.character(metadata$condition_code))) {
  stop("Demo condition_code parsing did not match metadata.", call. = FALSE)
}
if (!identical(normalize_condition(metadata$condition_code), metadata$condition)) {
  stop("Demo condition normalization did not match metadata.", call. = FALSE)
}

sample_cols <- intersect(metadata$sample_id, names(pg_matrix))
if (length(sample_cols) != nrow(metadata)) {
  stop("Demo metadata sample_id values do not match pg matrix sample columns.", call. = FALSE)
}

intensity_mat <- as.matrix(pg_matrix[, sample_cols, drop = FALSE])
storage.mode(intensity_mat) <- "numeric"

result_table <- data.frame(
  sample_id = sample_cols,
  sample_class = metadata$sample_class[match(sample_cols, metadata$sample_id)],
  condition_code = as.character(metadata$condition_code[match(sample_cols, metadata$sample_id)]),
  condition = metadata$condition[match(sample_cols, metadata$sample_id)],
  mean_intensity = as.numeric(colMeans(intensity_mat, na.rm = TRUE)),
  stringsAsFactors = FALSE
)

comparison_cols <- grep("^[A-Za-z.]+\\.", names(pg_matrix), value = TRUE)
comparison_table <- do.call(rbind, lapply(comparison_cols, function(col) {
  split <- regexec("^([A-Za-z.]+)\\.(.+)$", col)
  parts <- regmatches(col, split)[[1]]
  if (length(parts) != 3) {
    stop("Could not split demo comparison column: ", col, call. = FALSE)
  }
  parsed <- parse_comparison_key(parts[3])
  if (is.null(parsed)) {
    stop("Could not parse demo comparison key: ", parts[3], call. = FALSE)
  }
  data.frame(
    metric = parts[2],
    comparison = parts[3],
    canonical_comparison = parsed$name,
    stringsAsFactors = FALSE
  )
}))

write.csv(result_table, file.path(output_dir, "demo_result_table.csv"), row.names = FALSE)
write.csv(comparison_table, file.path(output_dir, "demo_comparison_table.csv"), row.names = FALSE)

expected_result <- read.csv(file.path(expected_dir, "expected_demo_result_table.csv"), check.names = FALSE, stringsAsFactors = FALSE)
expected_comparison <- read.csv(file.path(expected_dir, "expected_demo_comparison_table.csv"), check.names = FALSE, stringsAsFactors = FALSE)
expected_result$condition_code <- as.character(expected_result$condition_code)

if (!isTRUE(all.equal(result_table, expected_result, check.attributes = FALSE))) {
  stop("Demo result table does not match expected output.", call. = FALSE)
}
if (!isTRUE(all.equal(comparison_table, expected_comparison, check.attributes = FALSE))) {
  stop("Demo comparison table does not match expected output.", call. = FALSE)
}

message("Demo completed successfully.")
message("Wrote: ", file.path(output_dir, "demo_result_table.csv"))
message("Wrote: ", file.path(output_dir, "demo_comparison_table.csv"))
