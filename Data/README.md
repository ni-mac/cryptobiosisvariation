# Data guide

`file_inventory.csv` lists every supplied data file and its role in the current pipeline.

## Inputs read by the scripts

- Salinity, recovery and cultivation observations come from `curation/Cryptobiosis_master_final.xlsx`, using `Core_master`, `Mira_raw_full` and `MAL_daily_long`, respectively.
- Cultivation plate dates are read from the three original Excel workbooks ending in `20251023.xlsx` under `raw/cultivation/`.
- D1 counts come from `raw/desiccation/D1/Rena_anhydrobiosis.xlsx`.
- D2 counts come from the combined metadata CSV and the P. shuimeiren technical-count CSV under `raw/desiccation/D2/`.
- The default EC50 confidence intervals use the recorded bootstrap draws under `curation/salinity/`.

## Supporting sources

The remaining experimental files document the observations used to prepare the inputs. The pipeline does not rebuild the complete master workbook or the D2 combined CSV from these source files. Excel and CSV files may contain overlapping observations. Their presence alone does not imply they are exact duplicates.

Original filenames and workbook cells retain historical labels. ARO.22.05 in the salinity experiments is Acrobeloides goughiensis WWL726. The cooler ES5 exposure is 16 °C. Current strain names are listed in `strain_names.csv`.
