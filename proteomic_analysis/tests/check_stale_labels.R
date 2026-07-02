args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1) sub("^--file=", "", file_arg) else file.path("tests", "check_stale_labels.R")
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
if (!dir.exists(repo_root)) repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

all_files <- list.files(repo_root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
active_files <- all_files[
  grepl("\\.(r|R|md|Rmd)$", all_files) &
    !grepl("/legacy/", normalizePath(all_files, winslash = "/", mustWork = FALSE)) &
    !grepl("/tests/check_stale_labels\\.R$", normalizePath(all_files, winslash = "/", mustWork = FALSE))
]

read_file <- function(path) paste(readLines(path, warn = FALSE), collapse = "\n")
contents <- setNames(lapply(active_files, read_file), normalizePath(active_files, winslash = "/", mustWork = FALSE))

failures <- character()

add_hits <- function(label, pattern, ignore_case = TRUE) {
  hit <- vapply(contents, function(x) grepl(pattern, x, perl = TRUE, ignore.case = ignore_case), logical(1))
  if (any(hit)) {
    failures <<- c(failures, paste0(label, ": ", paste(names(contents)[hit], collapse = ", ")))
  }
}

add_hits("old con/res/sus condition mapping", '"[123]"\\s*=\\s*"(con|res|sus)"|phenotypes\\s*=\\s*c\\([^)]*"(con|res|sus)"|\\b(con|res|sus)_vs_(con|res|sus)\\b')
add_hits("group-code regex missing code 4", "\\[123\\]")
add_hits("old active sample-class parsing labels", "neuron_soma|neuron_neuropil|sample_class[^\\n]*(microglia|celltype_layer)|celltype_layer[^\\n]*sample_class")
add_hits("drafting-language output names/comments", "\\bNature\\b|publication|manuscript")
add_hits("old plotting helper names", "theme_nature|theme_nature_qc")

if (length(failures) > 0) {
  stop(paste(c("Stale-label audit failed:", failures), collapse = "\n"), call. = FALSE)
}

message("Stale-label audit passed for active code.")
