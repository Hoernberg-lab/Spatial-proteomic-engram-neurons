library(readxl)
library(dplyr)
library(writexl)
library(tools)
source(file.path("R", "analysis_labels.R"))

# ---- CONFIGURATION ----
# Set mode: "excel" for sheets in one Excel file, "folder" for all Excel files in a folder
mode <- "folder" # or "folder"

# Define file paths
file_path <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Datasets/gct/data/pg.matrix_filtered_pcaAdjusted_unnormalized.xlsx"
folder_path <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Datasets/gct/data/imputed/grouped"
metadata_path <- ""
output_dir <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Datasets/gct/data/morpheus"

# Read metadata
metadata <- read_excel(metadata_path)
if (!"sample_class" %in% names(metadata)) {
    metadata$sample_class <- parse_sample_class(metadata$sample_id)
}
if (!"condition_code" %in% names(metadata)) {
    source_condition <- if ("condition" %in% names(metadata)) metadata$condition else metadata$sample_id
    metadata$condition_code <- parse_condition_code(source_condition)
}
metadata$condition <- normalize_condition(metadata$condition_code)

# ---- READ DATA FILES ----
sheet_dfs <- list()
if (mode == "excel") {
    # Read all sheets from one Excel file
    sheet_names <- excel_sheets(file_path)
    sheet_dfs <- setNames(
        lapply(sheet_names, function(sheet) read_excel(file_path, sheet = sheet)),
        sheet_names
    )
} else if (mode == "folder") {
    # Read all Excel files in folder (ignore sheets, use file name as key)
    excel_files <- list.files(folder_path, pattern = "\\.xlsx$", full.names = TRUE)
    sheet_dfs <- setNames(
        lapply(excel_files, function(f) read_excel(f)),
        basename(file_path_sans_ext(excel_files))
    )
}

# ---- PROCESS EACH DATA FRAME ----
new_dfs <- lapply(sheet_dfs, function(df) {
    # Rename "T: Protein.Names" to "id"
    # Rename "T: Protein.Names" or "Protein.Names" to "id"
    name_col <- intersect(c("T: Protein.Names", "Protein.Names"), colnames(df))
    if (length(name_col) == 1) {
        colnames(df)[colnames(df) == name_col] <- "id"
        df <- df[, c("id", setdiff(colnames(df), "id"))]
    }
    # Remove columns starting from "T: Protein.Group" or "Protein.Group" and after
    group_col <- intersect(c("T: Protein.Group", "Protein.Group"), colnames(df))
    if (length(group_col) == 1) {
        idx <- which(colnames(df) == group_col)
        if (idx > 1) {
            df <- df[, seq_len(idx - 1)]
        } else {
            # If group_col is the first column, remove all columns (or handle as needed)
            df <- df[, FALSE]
        }
    }
    sample_cols <- intersect(colnames(df), metadata$sample_id)
    if (length(sample_cols) == 0) return(df)
    # Build metadata rows
    meta_rows <- lapply(names(metadata), function(meta_col) {
        row <- rep(NA, ncol(df))
        names(row) <- colnames(df)
        row["id"] <- meta_col
        if (meta_col %in% colnames(metadata) && length(sample_cols) > 0) {
            idxs <- match(sample_cols, metadata$sample_id)
            valid <- !is.na(idxs) & idxs >= 1 & idxs <= nrow(metadata)
            if (any(valid)) {
                row[sample_cols[valid]] <- metadata[[meta_col]][idxs[valid]]
            }
        }
        row
    })
    # ---- Add combined row: sample_class_condition ----
    combined_row <- rep(NA, ncol(df))
    names(combined_row) <- colnames(df)
    combined_row["id"] <- "sample_class_condition"
    if (all(c("sample_class", "condition") %in% colnames(metadata))) {
        idxs <- match(sample_cols, metadata$sample_id)
        valid <- !is.na(idxs) & idxs >= 1 & idxs <= nrow(metadata)
        if (any(valid)) {
            sample_class <- as.character(metadata$sample_class[idxs[valid]])
            condition <- as.character(metadata$condition[idxs[valid]])
            combined <- paste(sample_class, condition, sep = "_")
            combined_row[sample_cols[valid]] <- combined
        }
    }
    meta_df <- as.data.frame(do.call(rbind, meta_rows), stringsAsFactors = FALSE)
    meta_df <- rbind(meta_df, as.data.frame(t(combined_row), stringsAsFactors = FALSE))
    # Combine metadata rows and original data
    final_df <- rbind(meta_df, df)
    # ---- Add phenotypeWithinUnit row after NA replacement ----
    phenotype_row <- rep(NA, ncol(final_df))
    names(phenotype_row) <- colnames(final_df)
    phenotype_row["id"] <- "phenotypeWithinUnit"
    # Use combined_row values for phenotypeWithinUnit
    if (all(c("sample_class", "condition") %in% colnames(metadata))) {
        idxs <- match(sample_cols, metadata$sample_id)
        valid <- !is.na(idxs) & idxs >= 1 & idxs <= nrow(metadata)
        if (any(valid)) {
            sample_class <- as.character(metadata$sample_class[idxs[valid]])
            condition <- as.character(metadata$condition[idxs[valid]])
            combined <- paste(sample_class, condition, sep = "_")
            phenotype_row[sample_cols[valid]] <- combined
        }
    }
    # Insert phenotypeWithinUnit row after sample_class and condition rows
    insert_candidates <- which(final_df$id %in% c("sample_class", "condition", "condition_code"))
    insert_idx <- if (length(insert_candidates) > 0) max(insert_candidates) else nrow(meta_df)
    final_df <- rbind(
        final_df[seq_len(insert_idx), , drop = FALSE],
        as.data.frame(t(phenotype_row), stringsAsFactors = FALSE),
        final_df[(insert_idx + 1):nrow(final_df), , drop = FALSE]
    )
    # ---- Remove columns where "exclude" metadata row is TRUE ----
    exclude_row_idx <- which(final_df$id == "exclude")
    if (length(exclude_row_idx) == 1) {
        exclude_row <- final_df[exclude_row_idx, , drop = FALSE]
        exclude_cols <- which(tolower(as.character(exclude_row)) == "true")
        exclude_cols <- setdiff(exclude_cols, which(colnames(final_df) == "id"))
        if (length(exclude_cols) > 0) {
            final_df <- final_df[, -exclude_cols, drop = FALSE]
        }
    }
    final_df
})

# ---- SAVE EACH DATA FRAME AS EXCEL ----
for (sheet in names(new_dfs)) {
    out_path <- file.path(output_dir, paste0(sheet, "_with_metadata.xlsx"))
    write_xlsx(new_dfs[[sheet]], out_path)
}

# ---- SAVE EACH DATA FRAME AS GCT v1.3 ----
write_gct_v1.3 <- function(df, file, metadata) {
    meta_row_idx <- which(df$id %in% c(names(metadata), "sample_class_condition", "phenotypeWithinUnit"))
    data_rows <- df[-meta_row_idx, , drop = FALSE]
    data_rows <- data_rows[!is.na(data_rows$id) & data_rows$id != "", , drop = FALSE]
    sample_cols <- setdiff(colnames(df), "id")
    sample_cols <- intersect(sample_cols, colnames(data_rows))
    if (length(sample_cols) == 0 || nrow(data_rows) == 0) {
        warning(sprintf("No data to write for file: %s", file))
        return()
    }
    is_numeric_row <- function(row) {
        all(suppressWarnings(!is.na(as.numeric(row[sample_cols]))))
    }
    numeric_rows_idx <- apply(data_rows, 1, is_numeric_row)
    data_rows <- data_rows[numeric_rows_idx, , drop = FALSE]
    for (col in sample_cols) {
        data_rows[[col]] <- as.numeric(data_rows[[col]])
    }
    if (nrow(data_rows) == 0) {
        warning(sprintf("No numeric data to write for file: %s", file))
        return()
    }
    con <- file(file, "wt")
    on.exit(close(con))
    writeLines("#1.3", con)
    writeLines(sprintf("%d\t%d", nrow(data_rows), length(sample_cols)), con)
    writeLines(paste(c("id", sample_cols), collapse = "\t"), con)
    for (i in seq_len(nrow(data_rows))) {
        writeLines(paste(unlist(data_rows[i, c("id", sample_cols)]), collapse = "\t"), con)
    }
}

for (sheet in names(new_dfs)) {
    out_path <- file.path(output_dir, paste0(sheet, "_with_metadata.gct"))
    write_gct_v1.3(new_dfs[[sheet]], out_path, metadata)
}
