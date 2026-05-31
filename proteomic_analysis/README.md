# Proteomic analysis

This directory contains the proteomics analysis code for the `topohl/Neha` project.

The README is intended as an entry point for reproducing and extending the analysis. It documents the expected workflow, data assumptions, output conventions, and quality-control checks. The exact script names should be kept synchronized with the files in this folder as the analysis develops.

## Purpose

The analysis in this folder is designed to process, quality-control, model, and visualize proteomics data. Typical use cases include:

- importing processed protein-intensity matrices and associated sample metadata;
- checking sample-level and protein-level quality metrics;
- normalizing or transforming abundance values where appropriate;
- running differential-abundance or group-comparison analyses;
- summarizing results in tables suitable for downstream interpretation;
- generating exploratory and publication-oriented figures;
- preparing ranked protein lists or gene-symbol tables for enrichment analysis.

## Expected project layout

A recommended structure is:

```text
proteomic_analysis/
├── README.md
├── data/                 # local input data; usually not committed
├── metadata/             # sample annotation files, if separated from data
├── scripts/              # analysis scripts
├── R/                    # reusable R helper functions, if present
├── results/              # generated output tables and figures
└── logs/                 # optional run logs
```

If the repository uses a different layout, prefer the existing layout and update this section accordingly.

## Inputs

The analysis usually requires two input layers:

1. A protein-level abundance matrix, with proteins or protein groups as rows and samples as columns.
2. A sample-metadata table, with one row per sample and variables such as experimental group, batch, sex, tissue, region, layer, or other study-specific annotations.

Recommended minimum metadata columns:

| Column | Meaning |
|---|---|
| `sample_id` | Unique sample identifier matching the abundance matrix column names. |
| `group` | Experimental group or condition. |
| `batch` | Processing or acquisition batch, if applicable. |
| `sex` | Biological sex, if relevant to the design. |
| `region` | Brain region or anatomical unit, if applicable. |
| `layer` | Layer or subregion annotation, if applicable. |

Assumption: raw or large intermediate proteomics files are not tracked in Git unless they are small and non-sensitive. Store large data files externally or add them to `.gitignore`.

## Reproducibility

Before running the analysis, record the R version and package versions:

```r
sessionInfo()
```

Recommended packages for typical R-based proteomics workflows include:

```r
install.packages(c(
  "tidyverse",
  "readxl",
  "openxlsx",
  "limma",
  "ggplot2",
  "pheatmap",
  "janitor",
  "here"
))
```

Use Bioconductor for packages such as `limma` when needed:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("limma")
```

For a fully reproducible project, consider adding an `renv.lock` file:

```r
install.packages("renv")
renv::init()
renv::snapshot()
```

## Suggested workflow

Run scripts from the repository root unless a script explicitly documents another working directory.

A typical workflow is:

1. Import abundance data and metadata.
2. Validate sample-name matching between matrix and metadata.
3. Filter proteins with excessive missingness or low detection.
4. Transform and normalize abundance values if this has not already been done upstream.
5. Perform exploratory QC: missingness, PCA, sample clustering, batch structure, group separation.
6. Fit statistical models for predefined biological contrasts.
7. Export result tables with effect sizes, test statistics, raw p-values, and adjusted p-values.
8. Generate figures for QC, differential-abundance summaries, heatmaps, and enrichment-ready protein lists.

Example command pattern:

```bash
Rscript proteomic_analysis/scripts/01_import_qc.R
Rscript proteomic_analysis/scripts/02_differential_analysis.R
Rscript proteomic_analysis/scripts/03_visualization.R
```

Adjust these commands to the actual script names in this directory.

## Statistical notes

For proteomics group comparisons, avoid relying only on nominal p-values, especially when sample size is small. Report:

- effect size or log2 fold change;
- standard error or confidence interval where available;
- raw p-value;
- multiple-testing-adjusted p-value;
- number of quantified samples per group;
- missingness per protein and per group.

If batch, sex, region, or layer are part of the design, include them explicitly in the model rather than treating all samples as exchangeable. For small-n analyses, interpret enrichment and differential-abundance results as evidence-weighted and exploratory unless independently validated.

## Output conventions

Recommended result folders:

```text
results/
├── qc/
├── tables/
├── figures/
├── enrichment_inputs/
└── logs/
```

Recommended file formats:

- `.xlsx` for human-readable result workbooks;
- `.csv` or `.tsv` for machine-readable tables;
- `.svg` or `.pdf` for vector figures;
- `.txt` or `.rds` for ranked lists or reusable R objects.

## Quality-control checklist

Before interpreting biological results, verify:

- sample IDs match exactly between metadata and abundance matrix;
- group labels and factor levels are correct;
- missingness is not strongly confounded with group or batch;
- PCA or clustering does not reveal sample swaps or obvious outliers;
- batch effects are inspected and, where justified, modeled or corrected;
- multiple-testing correction is applied across the relevant family of tests;
- exported tables contain enough information to reproduce contrasts.

## Version-control recommendations

Commit:

- analysis scripts;
- helper functions;
- small example metadata files if non-sensitive;
- documentation;
- final lightweight result summaries when appropriate.

Do not commit:

- large raw mass-spectrometry files;
- private or identifiable sample metadata;
- temporary cache files;
- machine-specific absolute paths;
- generated figures or tables unless they are intentionally versioned outputs.

## TODO

- Replace placeholder script names in the workflow section with the actual file names.
- Add exact input file names and required metadata columns.
- Add a short description of each analysis script.
- Add expected output file names.
- Add `renv` or another environment lockfile if strict reproducibility is required.
