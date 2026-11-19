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

- `Data/raw/` contains the experimental source files, grouped by experiment. The salinity workbooks contain the supplied adjusted observation sheets used for the analyses.
- `R/00_prepare_data.R` rebuilds the five required analysis tables from these sources. The start script runs this step automatically and saves the tables in `Results/inputs/`.
- `Data/strain_names.csv` maps historical strain identifiers to the manuscript names.
- Reported EC50 confidence intervals are reproduced from the 2,000 saved bootstrap results in `Data/reference/recorded_ec50_bootstrap_draws.csv`, but can also be generated from the code. 

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
