sample_classes <- c("mcherry", "neuropil", "cfos", "neuron")

condition_code_map <- c(
  "1" = "paired_cno",
  "2" = "paired_veh",
  "3" = "unpaired_cno",
  "4" = "unpaired_veh"
)

condition_levels <- unname(condition_code_map)
reference_condition <- "paired_veh"

source_analysis_labels <- function() {
  list(
    sample_classes = sample_classes,
    condition_code_map = condition_code_map,
    condition_levels = condition_levels,
    reference_condition = reference_condition
  )
}

normalize_condition <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("-", "_", x)
  mapped <- unname(condition_code_map[x])
  ifelse(!is.na(mapped), mapped, x)
}

parse_condition_code <- function(x) {
  x <- as.character(x)
  code <- ifelse(x %in% names(condition_code_map), x, NA_character_)
  fallback <- regmatches(x, regexpr("(?<![0-9])[1234](?![0-9])", x, perl = TRUE))
  fallback[fallback == ""] <- NA_character_
  ifelse(!is.na(code), code, fallback)
}

parse_sample_class <- function(sample_id) {
  pattern <- paste0("(?i)(^|[^A-Za-z0-9])(", paste(sample_classes, collapse = "|"), ")(?=$|[^A-Za-z0-9])")
  out <- regmatches(sample_id, regexpr(pattern, sample_id, perl = TRUE))
  out[out == ""] <- NA_character_
  tolower(gsub("^[^A-Za-z0-9]+|[^A-Za-z0-9]+$", "", out))
}
