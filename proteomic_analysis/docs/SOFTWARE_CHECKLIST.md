# Software Checklist

This document maps checklist items to repository files for the Neha proteomics analysis workflow.

| Requirement | Repository location |
| --- | --- |
| Source code | `proteomic_analysis/01_preprocessing/`, `proteomic_analysis/02_id_mapping/`, `proteomic_analysis/03_qc_exploration/`, `proteomic_analysis/04_differential_expression_enrichment/`, `proteomic_analysis/05_celltype_enrichment_EWCE/`, `proteomic_analysis/R/analysis_labels.R` |
| Demo dataset | `proteomic_analysis/demo/input/` |
| README system requirements | `proteomic_analysis/README.md` |
| Dependencies and versions | `proteomic_analysis/requirements_R.md`; exact local versions can be captured with `sessionInfo()` |
| Tested OS and R versions | `proteomic_analysis/README.md`, `proteomic_analysis/requirements_R.md`; CI workflow under `.github/workflows/proteomic-analysis-demo.yml` |
| Non-standard hardware | `proteomic_analysis/README.md` notes that the demo requires none |
| Installation instructions | `proteomic_analysis/README.md` |
| Install time | `proteomic_analysis/README.md` |
| Demo instructions | `proteomic_analysis/README.md`; runnable script is `proteomic_analysis/run_demo.R` |
| Expected output | `proteomic_analysis/demo/expected_output/` |
| Demo runtime | `proteomic_analysis/README.md` |
| Running on own data | `proteomic_analysis/README.md`; configurable script is `proteomic_analysis/01_preprocessing/02_excel_convert.r` |
| Optional full reproduction instructions | `proteomic_analysis/README.md` and stage folders under `proteomic_analysis/` |
| License | `LICENSE` |
| Citation | `CITATION.cff` |
| Code repository link | `https://github.com/topohl/Neha` |

The lightweight demo is designed to run without private shared-drive paths or heavy enrichment analysis.
