# Variation in Cryptobiosis Potential Among Free-Living Desert Nematode Taxa

Data and R scripts for the salinity dose-response, recovery, long-term cultivation (CE1–CE3), and desiccation (D1–D2) experiments.

## Run the analyses

Install R and the packages listed in `Verification/package_versions.csv`. The recorded R environment is in `Verification/sessionInfo.txt`. From the repository root, run:

```sh
Rscript R/00_run_all.R
```

The scripts write model results and tables to `Results/` and figures to `Manuscript/`.

## Data sources

- `Data/raw/` contains original experimental source files. Some original files retain historical labels. In the manuscript, ARO.22.05 is *Acrobeloides goughiensis* WWL726, and the cooler ES5 treatment is 16 °C.
- `Data/curation/Cryptobiosis_master_final.xlsx` is the curated workbook used as an analysis input. The pipeline starts from this workbook and does not reconstruct it from the raw files. The cultivation analysis also reads original plate logs. Desiccation analyses read D1 and D2 inputs directly.
- `Data/curation/salinity/recorded_ec50_bootstrap_draws.csv` contains the saved draws used by the default run to reproduce the reported EC50 confidence intervals.
