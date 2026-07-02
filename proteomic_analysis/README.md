# Proteomic analysis

## Overview

This repository contains analysis code for the Neha proteomics workflow.

The goal is to keep preprocessing, quality control, statistical modeling, and enrichment analysis transparent, reproducible, and reusable for related proteomics datasets.

## Canonical Labels

Neha sample classes are:

- `mcherry`
- `neuropil`
- `cfos`
- `neuron`

Experimental condition codes are:

- `1` = `paired_cno`
- `2` = `paired_veh`
- `3` = `unpaired_cno`
- `4` = `unpaired_veh`

Shared R definitions live in [R/analysis_labels.R](R/analysis_labels.R).

## Workflow

The repository is organized by analysis stage:

- `01_preprocessing`: metadata formatting, imputation, and GCT conversion
- `02_id_mapping`: protein identifier mapping
- `03_qc_exploration`: quality-control exploration
- `04_differential_expression_enrichment`: differential expression and pathway enrichment
- `05_celltype_enrichment_EWCE`: EWCE enrichment against external reference cell types

Old Exp9/E9-specific scripts are preserved under `legacy/`.

## Reproducibility

Recommended procedure:

```bash
git clone https://github.com/topohl/Neha.git
cd Neha/proteomic_analysis
```

Record:

```r
sessionInfo()
```

If available, use a locked environment:

```r
renv::restore()
```

## Data

Raw and processed proteomics data may be stored outside this repository. Large intermediate files are intentionally excluded from Git.

## Checks

Run the lightweight terminology audit with:

```bash
Rscript tests/check_stale_labels.R
```
