# ================================================================
# Patch tracing_cfos_correlation.r to use normalized raw input files
# ================================================================
# Run this once after normalize_tracing_raw_inputs.R.
# It edits only the input-file lines in tracing_cfos_correlation.r and writes
# a timestamped backup before changing anything.
# ================================================================

script_path <- "c:/Users/topohl/Documents/GitHub/Neha/tracing_cfos_correlation.r"

if (!file.exists(script_path)) {
  stop("Could not find tracing script: ", script_path)
}

x <- readLines(script_path, warn = FALSE)
backup_path <- paste0(script_path, ".backup_", format(Sys.time(), "%Y%m%d_%H%M%S"))
writeLines(x, backup_path)
message("Backup written to: ", backup_path)

replace_line <- function(lines, pattern, replacement) {
  hit <- grep(pattern, lines)
  if (length(hit) != 1) {
    stop("Expected exactly one match for pattern: ", pattern, "; found: ", length(hit))
  }
  lines[hit] <- replacement
  lines
}

x <- replace_line(
  x,
  '^intensity_raw_file\\s*<-',
  'intensity_raw_file <- file.path(input_dir, "raw_input_normalized", "IntensityNoOutliersFin_normalized.xlsx")'
)

x <- replace_line(
  x,
  '^cell_count_raw_file\\s*<-',
  'cell_count_raw_file <- file.path(input_dir, "raw_input_normalized", "CellCountNoOutliersFin_normalized.csv")'
)

# Make the raw reader robust to files where the full region name column is
# called Abbrev/abbrev instead of Abbreviation/abbreviation.
needle <- 'x <- x %>% janitor::clean_names()'
insert <- c(
  '  x <- x %>% janitor::clean_names()',
  '  names(x) <- stringr::str_replace_all(names(x), "^abbrev$", "abbreviation")',
  '  names(x) <- stringr::str_replace_all(names(x), "^abbr$", "abbreviation")',
  '  names(x) <- stringr::str_replace_all(names(x), "^cellcount$", "cell_count")'
)

hit <- grep(fixed = TRUE, needle, x)
if (length(hit) >= 1) {
  # Only patch the first reader occurrence if it has not already been patched.
  already_patched <- any(grepl('str_replace_all\\(names\\(x\\), "\\^abbrev\\$"', x[pmax(1, hit[1] - 2):pmin(length(x), hit[1] + 5)]))
  if (!already_patched) {
    x[hit[1]] <- paste(insert, collapse = "\n")
  }
}

writeLines(x, script_path)
message("Patched tracing script to use normalized raw inputs:")
message("  ", script_path)
message("Now run:")
message('  source("c:/Users/topohl/Documents/GitHub/Neha/tracing_cfos_correlation.r")')
