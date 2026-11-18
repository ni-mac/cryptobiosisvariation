# Salinity curation and reproducibility

The shared `../Cryptobiosis_master_final.xlsx` is the historical **curated integrated master workbook**, assembled from individual source workbooks, not an untouched raw laboratory workbook. The raw-source archive is retained separately. `Core_master` contains source filename/sheet provenance. The R module imports that master afresh, reconstructs eligible salinity observations and endpoints and refits the selected GLMM. `Mira_raw_full` supplies the independently recomputed recovery statistics.

`recorded_ec50_bootstrap_draws.csv` contains the 2000 historical fitted plate-cluster bootstrap results, copied without changes. These are a **derived reference cache**, not raw input. The old bootstrap used load-balanced parallel processing with random streams assigned to workers, so its exact random sequence depends on scheduling. The default `bootstrap_reference_replay` mode computes all reported interval endpoints directly from these recorded results. A separate `fresh_bootstrap` option refits bootstrap samples with deterministic per-replicate seeds and outputs new results. It must not silently replace the manuscript confidence intervals or be described as an identical replay.


Known explicit curation decisions:

- The cooler ES5 treatment is 16 °C. Historical 15 °C labels are corrected during import. The saved bootstrap file retains its original column header, translated to 16 °C during import.
- Preliminary failed ES5 run is excluded; observation days follow documented per-experiment source mapping.
- Thomas combined affected endpoint includes dead plus dormant counts. Other source-specific mortality/immobility definitions are preserved and recorded.
- The recovery analysis excludes the previously identified ALT16C/17PSU/replicate1 internally inconsistent count (10 recovered versus4 candidates), plus zero-candidate plates for candidate-conditional recovery; the reference picture marks zero-candidate tested doses with crosses.
- Supplementary Figure S2 is rendered from the recomputed recovery counts.


Active analysis and result CSVs use `ES5__16C`. The correction changes the treatment label and preserves the observed counts.
