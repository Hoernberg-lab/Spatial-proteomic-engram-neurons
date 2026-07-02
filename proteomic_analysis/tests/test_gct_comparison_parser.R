source(file.path("R", "analysis_labels.R"))

cases <- c(
  "neuron_3.over.neuron_2" = "neuron_unpaired_cno_vs_neuron_paired_veh",
  "neuron_4.over.neuron_2" = "neuron_unpaired_veh_vs_neuron_paired_veh",
  "neuron_4.over.neuron_3" = "neuron_unpaired_veh_vs_neuron_unpaired_cno",
  "bg_2.over.bg_1" = "neuropil_paired_veh_vs_neuropil_paired_cno",
  "bg_3.over.bg_1" = "neuropil_unpaired_cno_vs_neuropil_paired_cno",
  "bg_4.over.bg_1" = "neuropil_unpaired_veh_vs_neuropil_paired_cno",
  "cfos_1.over.bg_1" = "cfos_paired_cno_vs_neuropil_paired_cno",
  "cfos_2.over.bg_1" = "cfos_paired_veh_vs_neuropil_paired_cno",
  "CA1_so_bg_2.over.CA1_so_bg_1" = "CA1_so_neuropil_paired_veh_vs_CA1_so_neuropil_paired_cno",
  "CA1_so_neuropil_1.over.CA1_so_neuropil_2" = "CA1_so_neuropil_paired_cno_vs_CA1_so_neuropil_paired_veh",
  "mcherry_1.over.mcherry_2" = "mcherry_paired_cno_vs_mcherry_paired_veh"
)

for (key in names(cases)) {
  parsed <- parse_comparison_key(key)
  if (is.null(parsed)) {
    stop("Expected parser success for: ", key, call. = FALSE)
  }
  if (!identical(parsed$name, unname(cases[[key]]))) {
    stop(
      "Unexpected parsed name for ", key,
      "\nExpected: ", unname(cases[[key]]),
      "\nActual:   ", parsed$name,
      call. = FALSE
    )
  }
}

invalid_cases <- c(
  "neuron.over.neuron_2",
  "neuron_5.over.neuron_2",
  "foo_2.over.neuron_2",
  "cfos_1_vs_bg_1"
)

for (key in invalid_cases) {
  if (!is.null(parse_comparison_key(key))) {
    stop("Expected parser failure for: ", key, call. = FALSE)
  }
}

parsed_codes <- parse_condition_code(c("1", "sample_4", NA, "none"))
if (!identical(parsed_codes, c("1", "4", NA_character_, NA_character_))) {
  stop("parse_condition_code must return one value per input.", call. = FALSE)
}

parsed_classes <- parse_sample_class(c("sample_bg_1", "demo_cFos_2", NA, "none"))
if (!identical(parsed_classes, c("neuropil", "cfos", NA_character_, NA_character_))) {
  stop("parse_sample_class must return one value per input.", call. = FALSE)
}

message("GCT comparison parser tests passed.")
