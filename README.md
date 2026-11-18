# Variation in Cryptobiosis Potential Among Free-Living Desert Nematode Taxa

Data and R scripts for the salinity dose-response, recovery, long-term cultivation (CE1–CE3), and desiccation (D1–D2) experiments.

## Run the analyses

The analysis was tested with R 4.5.1 on Windows 11. Install the packages listed below. From the repository root, run:

```sh
Rscript R/00_run_all.R
```

The scripts write model results and tables to `Results/` and figures to `Manuscript/`.

## Data sources

- `Data/raw/` contains original experimental source files. Some original files retain historical labels. In the manuscript, ARO.22.05 is *Acrobeloides goughiensis* WWL726, and the cooler ES5 treatment is 16 °C.
- `Data/curation/Cryptobiosis_master_final.xlsx` is the curated workbook used as an analysis input. The pipeline starts from this workbook and does not reconstruct it from the raw files. The cultivation analysis also reads original plate logs. Desiccation analyses read D1 and D2 inputs directly.
- `Data/curation/salinity/recorded_ec50_bootstrap_draws.csv` contains the saved draws used by the default run to reproduce the reported EC50 confidence intervals.

## Tested package versions

| Package | Version |
| --- | --- |
| readxl | 1.4.5 |
| dplyr | 1.1.4 |
| tidyr | 1.3.1 |
| stringr | 1.5.2 |
| purrr | 1.1.0 |
| survival | 3.8.3 |
| ggplot2 | 4.0.1 |
| patchwork | 1.3.2 |
| glmmTMB | 1.1.14 |
| MASS | 7.3.65 |
| readr | 2.1.5 |
| scales | 1.4.0 |
| cowplot | 1.2.0 |
| ggtext | 0.1.2 |
| emmeans | 2.0.4 |
| DHARMa | 0.5.0 |
| png | 0.1.8 |
| jsonlite | 2.0.0 |
| digest | 0.6.37 |

Each run records its environment and execution times in `Results/run_info/`.
