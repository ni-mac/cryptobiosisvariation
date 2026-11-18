# Variation in Cryptobiosis Potential Among Free-Living Desert Nematode Taxa

Data and R scripts for salinity dose-response, recovery, long-term cultivation (CE1–CE3), and desiccation (D1–D2) experiments.

## Run the analyses

The analysis was tested with R 4.5.1 on Windows 11. Install the packages listed below. Download and extract the repository, open a terminal in the extracted folder, and run:

```sh
Rscript R/00_run_all.R
```

Alternatively, set the extracted folder as the working directory in RStudio and run these commands in the R console:

```r
source("R/00_run_all.R")
run_all(".")
```

Running the scripts creates `Results/` for model results and tables, `Manuscript/figures/` for the main figures, and `Manuscript/Supplement/Supplement_figures/` for the supplementary figures. Run information is saved in `Results/run_info/`. Existing output files are overwritten when the analyses are run again.

## Data sources

- `Data/raw/` contains the experimental source files, grouped by experiment.
- `Data/curation/Cryptobiosis_master_final.xlsx` contains the combined observations used for salinity (`Core_master`), recovery (`salinity_recovery_raw_full`), and cultivation (`MAL_daily_long`). Cultivation plate dates are also read from the three original Excel workbooks ending in `20251023.xlsx`.
- The scripts use the prepared master workbook and D2 CSV inputs. Preparing these inputs from the original files is currently a separate step.
- `Data/curation/salinity/recorded_ec50_bootstrap_draws.csv` contains the 2,000 saved bootstrap results used to reproduce the reported EC50 confidence intervals. These are derived results. A fresh bootstrap can be requested through the salinity function and may give slightly different confidence limits.

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
