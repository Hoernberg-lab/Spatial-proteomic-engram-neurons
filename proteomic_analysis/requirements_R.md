# R Requirements

## Tested Environment

- Windows 11 with R 4.5.1
- Ubuntu GitHub Actions runner with R 4.5.1

The lightweight demo in `run_demo.R` uses only base R.

## Full Workflow Packages

Package versions should be recorded with `sessionInfo()` for each full analysis run. The workflow has been checked with the following package families:

| Package | Tested version |
| --- | --- |
| `dplyr` | 1.1.x |
| `tidyr` | 1.3.x |
| `readxl` | 1.4.x |
| `writexl` | 1.5.x |
| `ggplot2` | 3.5.x |
| `patchwork` | 1.3.x |
| `stringr` | 1.5.x |
| `openxlsx` | 4.2.x |
| `tibble` | 3.2.x |
| `svglite` | 2.1.x |
| `pheatmap` | 1.0.x |
| `ggridges` | 0.5.x |
| `ggrepel` | 0.9.x |
| `viridis` | 0.6.x |
| `future` | 1.33.x |
| `future.apply` | 1.11.x |
| `digest` | 0.6.x |

## Bioconductor Packages

The full EWCE workflow requires Bioconductor packages:

- `limma`
- `EWCE`
- `ewceData`
- `org.Mm.eg.db`
- `AnnotationDbi`

These packages are not needed for the lightweight demo or the stale-label/parser checks.
