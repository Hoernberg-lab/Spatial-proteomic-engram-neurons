sample_classes <- c("mcherry", "neuropil", "cfos", "neuron")

sample_class_aliases <- c(
  "bg" = "neuropil",
  "background" = "neuropil",
  "mcherry" = "mcherry",
  "neuropil" = "neuropil",
  "cfos" = "cfos",
  "neuron" = "neuron"
)

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
    sample_class_aliases = sample_class_aliases,
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
  x_chr <- as.character(x)
  out <- rep(NA_character_, length(x_chr))

  direct <- !is.na(x_chr) & x_chr %in% names(condition_code_map)
  out[direct] <- x_chr[direct]

  needs_regex <- !direct & !is.na(x_chr)
  matches <- regexpr("(?<![0-9])[1234](?![0-9])", x_chr[needs_regex], perl = TRUE)
  matched_idx <- which(needs_regex)[matches > 0]
  out[matched_idx] <- regmatches(x_chr[needs_regex], matches)[matches > 0]

  out
}

parse_sample_class <- function(sample_id) {
  sample_id_chr <- as.character(sample_id)
  out <- rep(NA_character_, length(sample_id_chr))
  pattern <- paste0("(?i)(^|[^A-Za-z0-9])(", paste(names(sample_class_aliases), collapse = "|"), ")(?=$|[^A-Za-z0-9])")

  can_match <- !is.na(sample_id_chr)
  matches <- regexpr(pattern, sample_id_chr[can_match], perl = TRUE)
  matched_idx <- which(can_match)[matches > 0]

  if (length(matched_idx) > 0) {
    raw <- regmatches(sample_id_chr[can_match], matches)[matches > 0]
    token <- gsub("^[^A-Za-z0-9]+|[^A-Za-z0-9]+$", "", raw)
    out[matched_idx] <- normalize_sample_class(token)
  }

  out
}

normalize_sample_class <- function(x) {
  x0 <- tolower(as.character(x))
  out <- rep(NA_character_, length(x0))
  known <- !is.na(x0) & x0 %in% names(sample_class_aliases)
  out[known] <- unname(sample_class_aliases[x0[known]])
  out
}

parse_comparison_side <- function(side) {
  side <- as.character(side)
  tokens <- strsplit(side, "_", fixed = TRUE)[[1]]
  tokens <- tokens[nzchar(tokens)]

  if (length(tokens) < 2) {
    return(NULL)
  }

  condition_code <- tokens[length(tokens)]
  if (!condition_code %in% names(condition_code_map)) {
    return(NULL)
  }

  sample_class_token <- tokens[length(tokens) - 1]
  sample_class <- normalize_sample_class(sample_class_token)
  if (is.na(sample_class)) {
    return(NULL)
  }

  analysis_unit <- if (length(tokens) > 2) {
    paste(tokens[seq_len(length(tokens) - 2)], collapse = "_")
  } else {
    NA_character_
  }

  list(
    original = side,
    analysis_unit = analysis_unit,
    sample_class = sample_class,
    condition_code = condition_code,
    condition = unname(condition_code_map[condition_code])
  )
}

format_comparison_side <- function(parsed_side) {
  if (is.null(parsed_side)) {
    return(NA_character_)
  }

  parts <- c(
    parsed_side$analysis_unit,
    parsed_side$sample_class,
    parsed_side$condition
  )
  paste(stats::na.omit(parts), collapse = "_")
}

parse_comparison_key <- function(key) {
  key <- as.character(key)
  sides <- strsplit(key, ".over.", fixed = TRUE)[[1]]
  if (length(sides) != 2) {
    return(NULL)
  }

  case <- parse_comparison_side(sides[1])
  reference <- parse_comparison_side(sides[2])

  if (is.null(case) || is.null(reference)) {
    return(NULL)
  }

  list(
    original = key,
    case = case,
    reference = reference,
    name = paste(format_comparison_side(case), format_comparison_side(reference), sep = "_vs_")
  )
}
