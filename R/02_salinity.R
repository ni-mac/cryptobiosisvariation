# Salinity analyses for the current manuscript. Sources and curation are documented
# in README.md. No files outside project_root are written.
# bootstrap_reference_replay uses the recorded 2000 fitted cluster-bootstrap draws.
# fresh_bootstrap uses deterministic replicate-specific seeds; it is not expected
# to duplicate the historical load-balanced RNG schedule.
run_salinity <- function(project_root,
                         bootstrap_mode = c("bootstrap_reference_replay", "fresh_bootstrap"),
                         n_boot = 2000L, render_figures = TRUE) {
  bootstrap_mode <- match.arg(bootstrap_mode)
  suppressPackageStartupMessages({
    library(readxl); library(dplyr); library(glmmTMB); library(MASS)
    library(ggplot2); library(patchwork); library(readr); library(scales); library(cowplot)
  })
  root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  cur <- file.path(root, "Data", "reference")
  cache <- file.path(root, "Results", "cache", "salinity")
  stats <- file.path(root, "Results", "statistics", "salinity")
  figures <- file.path(root, "Manuscript", "figures")
  for (folder in c(cache, stats, figures)) dir.create(folder, recursive=TRUE, showWarnings=FALSE)
  input_folder <- file.path(root, "Results", "inputs")
  options(stringsAsFactors=FALSE)
  here <- cache

protocol_map <- tibble::tribble(
  ~Experiment_id, ~rehydration_group, ~metadata_status, ~analysis_role, ~experiment_year, ~protocol_writeup_year,
  "PAP.22.17_salinity_R1", "counted_not_rehydrated", "confirmed", "supporting", NA_integer_, NA_integer_,
  "ES5_salinity_preliminary_R1", "rehydration_not_documented", "failed_preliminary", "exploratory_only", NA_integer_, NA_integer_,
  "PAP.22.17_salinity_R2", "autoclaved_H2O_after_t7", "confirmed", "supporting", NA_integer_, NA_integer_,
  "ARO.22.05_salinity_R1", "autoclaved_H2O_after_t4", "confirmed", "supporting", NA_integer_, NA_integer_,
  "ES5_salinity_R1", "autoclaved_H2O_after_t11", "confirmed_results_not_good", "supporting", 2022L, 2023L,
  "A.max_22C", "M9_after_t14", "confirmed", "primary_salinity_recovery", NA_integer_, NA_integer_,
  "ALT.22.08_16C", "M9_after_t14", "confirmed", "primary_salinity_recovery", NA_integer_, NA_integer_,
  "ALT.22.08_22C", "M9_after_t14", "confirmed", "primary_salinity_recovery", NA_integer_, NA_integer_,
  "ES5_22C", "M9_after_t14", "confirmed", "primary_salinity_recovery", NA_integer_, NA_integer_,
  "PAP.22.17_16C", "M9_after_t14", "confirmed", "primary_salinity_recovery", NA_integer_, NA_integer_,
  "series2023_ALT.22.08_salinity_22C", "counted_not_rehydrated", "confirmed", "supporting", 2023L, NA_integer_,
  "series2023_LAG.22.03_salinity_22C", "counted_not_rehydrated", "confirmed", "supporting", 2023L, NA_integer_,
  "series2023_A.max_salinity_22C", "counted_not_rehydrated", "confirmed", "supporting", 2023L, NA_integer_
)

expected_experiments <- protocol_map$Experiment_id
numeric_cols <- c("Temperature_treatment", "Salinity [PSU]", "Observation_day_raw", "Replicate", "Alive", "Dead", "Total", "Dormant_stage", "Dead_without_dormant")
export_cols <- c(
  "Data_source", "Experiment_id", "Strain", "Species", "reproductive_mode", "Sampling location",
  "Temperature_treatment", "Observation_day_raw", "Salinity [PSU]", "Replicate", "experimental_unit_id",
  "Alive", "Dead", "Dead_without_dormant", "Analysis_dead", "Analysis_not_dead", "Total", "Dormant_stage", "Paralysed_adult",
  "Analysis_endpoint_definition", "rehydration_group", "metadata_status", "experiment_year", "protocol_writeup_year",
  "analysis_role", "Mortality_definition", "Mortality_analysis_eligible", "analysis_included_primary",
  "analysis_exclusion_reason", "mortality_proportion", "survival_proportion", "QC_flag", "Count_consistency_QC",
  "Source_file", "Source_sheet", "Source_date", "notes"
)

as_bool <- function(x) {
  if (is.logical(x)) return(replace(x, is.na(x), FALSE))
  z <- tolower(trimws(as.character(x)))
  z %in% c("true", "t", "1", "yes", "y")
}

collapse_unique <- function(x) paste(unique(x[nzchar(x)]), collapse = "; ")

wilson <- function(k, n, z = qnorm(0.975)) {
  p <- ifelse(n > 0, k / n, NA_real_)
  den <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / den
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  cbind(proportion = p, low = pmax(0, centre - half), high = pmin(1, centre + half))
}


core <- read.csv(file.path(input_folder, "salinity_observations.csv"), check.names = FALSE, na.strings = "")
derived_cols <- c(
  "experimental_unit_id", "rehydration_group", "metadata_status", "analysis_role",
  "experiment_year", "protocol_writeup_year",
  "analysis_included_primary", "analysis_exclusion_reason", "mortality_proportion", "survival_proportion",
  "Analysis_dead", "Analysis_not_dead", "Analysis_endpoint_definition"
)
required <- setdiff(unique(c("Analysis_family", export_cols, numeric_cols)), derived_cols)
missing_cols <- setdiff(required, names(core))
if (length(missing_cols)) stop("Missing source columns: ", paste(missing_cols, collapse = ", "))

sal <- core %>% filter(Analysis_family == "osmobiosis_salinity") %>% left_join(protocol_map, by = "Experiment_id")
unmapped <- setdiff(unique(sal$Experiment_id), expected_experiments)
if (length(unmapped)) stop("Unmapped salinity experiments: ", paste(unmapped, collapse = ", "))
if (any(is.na(sal$rehydration_group))) stop("Protocol mapping produced missing values.")

for (nm in numeric_cols) sal[[nm]] <- suppressWarnings(as.numeric(sal[[nm]]))

# User confirmed 8 Nov 2026: the cool ES5 experiment was conducted at 16 C.
# Correct the historical source label when importing; original workbooks remain preserved.
sal$Temperature_treatment[sal$Experiment_id == "ES5_salinity_R1"] <- 16


# Fail fast if the confirmed metadata are not represented in the final master.
series_2023_ids <- c(
  "series2023_ALT.22.08_salinity_22C",
  "series2023_LAG.22.03_salinity_22C",
  "series2023_A.max_salinity_22C"
)
series_2023 <- sal %>% filter(Experiment_id %in% series_2023_ids)
if (!all(series_2023_ids %in% unique(series_2023$Experiment_id))) stop("One or more confirmed 2023 salinity series are absent.")
if (!all(series_2023$Temperature_treatment == 22)) stop("Confirmed 2023 series must all have treatment temperature 22°C.")
series_2023_design <- series_2023 %>% group_by(Experiment_id) %>%
  summarise(max_psu = max(.data[["Salinity [PSU]"]], na.rm = TRUE), includes_25_psu = any(.data[["Salinity [PSU]"]] == 25), .groups = "drop")
if (!all(series_2023_design$max_psu == 40 & series_2023_design$includes_25_psu)) {
  stop("Confirmed 2023 series must extend to 40 PSU and include 25 PSU.")
}
supporting_es5 <- sal %>% filter(Experiment_id == "ES5_salinity_R1")
if (!nrow(supporting_es5) || !all(supporting_es5$Temperature_treatment == 16)) stop("ES5_salinity_R1 must be recorded at 16°C.")
if (!all(supporting_es5$experiment_year == 2022L & supporting_es5$protocol_writeup_year == 2023L)) {
  stop("ES5_salinity_R1 must be mapped to experiment year 2022 and protocol/write-up year 2023.")
}
lag_rows <- sal %>% filter(Strain == "LAG.22.03")
if (!nrow(lag_rows) || !all(lag_rows$Species == "Plectus pavus")) stop("LAG.22.03 must be recorded as Plectus pavus.")
if (any(grepl("pending", sal$metadata_status, fixed = TRUE))) stop("All salinity metadata must be resolved in the final workflow.")

metadata_validation <- bind_rows(
  series_2023_design %>% transmute(check = paste0(Experiment_id, "_design"), expected = "22°C; max 40 PSU; includes 25 PSU", observed = paste0("22°C; max ", max_psu, " PSU; includes 25 PSU=", includes_25_psu), pass = max_psu == 40 & includes_25_psu),
  tibble(check = "ES5_salinity_R1_metadata", expected = "16°C; experiment 2022; protocol/write-up 2023", observed = "16°C; experiment 2022; protocol/write-up 2023", pass = TRUE),
  tibble(check = "LAG.22.03_species", expected = "Plectus pavus", observed = paste(unique(lag_rows$Species), collapse = "; "), pass = all(lag_rows$Species == "Plectus pavus"))
)
write.csv(metadata_validation, file.path(here, "metadata_validation.csv"), row.names = FALSE, na = "")

# Preserve source counts and use the endpoint prepared for analysis in each
# source workbook. For the A. maximus, ALT.22.08 and PAP.22.17 adjusted sheets this
# is Dead_with_paralysed, stored as Dead during preparation and explicitly named
# in Mortality_definition. This is an operational, not confirmed, mortality endpoint.
sal$Analysis_dead <- sal$Dead
sal$Analysis_not_dead <- sal$Total - sal$Analysis_dead
sal$Analysis_endpoint_definition <- sal$Mortality_definition

# Dose-response days documented in the experiment records and analysis sheets.
salinity_recovery_analysis_days <- list(
  "A.max_22C" = c(4, 7),
  "ALT.22.08_16C" = c(4, 7, 9),
  "ALT.22.08_22C" = c(4, 7, 9),
  "ES5_22C" = c(1, 2, 4),
  "PAP.22.17_16C" = c(4, 7, 9)
)
day_eligible <- vapply(seq_len(nrow(sal)), function(i) {
  allowed <- salinity_recovery_analysis_days[[sal$Experiment_id[i]]]
  is.null(allowed) || sal$Observation_day_raw[i] %in% allowed
}, logical(1))

source_eligible <- as_bool(sal$Mortality_analysis_eligible)
count_ok <- !is.na(sal$Analysis_dead) & !is.na(sal$Total) & sal$Total > 0 & sal$Analysis_dead >= 0 & sal$Analysis_dead <= sal$Total &
  !is.na(sal[["Salinity [PSU]"]]) & !is.na(sal$Observation_day_raw) & !is.na(sal$Replicate)
sal$analysis_included_primary <- source_eligible & count_ok & day_eligible

reason_one <- function(i) {
  z <- character()
  if (!source_eligible[i]) z <- c(z, ifelse(is.na(sal$QC_flag[i]) || sal$QC_flag[i] == "", "source_marked_ineligible", as.character(sal$QC_flag[i])))
  if (!day_eligible[i]) z <- c(z, "outside_documented_dose_response_days")
  if (is.na(sal$Analysis_dead[i]) || is.na(sal$Total[i])) z <- c(z, "missing_analysis_dead_or_total")
  else if (sal$Total[i] <= 0) z <- c(z, "nonpositive_total")
  else if (sal$Analysis_dead[i] < 0 || sal$Analysis_dead[i] > sal$Total[i]) z <- c(z, "analysis_dead_outside_0_total")
  if (is.na(sal[["Salinity [PSU]"]][i]) || is.na(sal$Observation_day_raw[i]) || is.na(sal$Replicate[i])) z <- c(z, "missing_design_variable")
  paste(unique(z), collapse = "; ")
}
sal$analysis_exclusion_reason <- vapply(seq_len(nrow(sal)), reason_one, character(1))
sal$mortality_proportion <- ifelse(sal$analysis_included_primary, sal$Analysis_dead / sal$Total, NA_real_)
sal$survival_proportion <- ifelse(sal$analysis_included_primary, sal$Analysis_not_dead / sal$Total, NA_real_)
sal$experimental_unit_id <- paste0(sal$Experiment_id, "|PSU=", sal[["Salinity [PSU]"]], "|rep=", sal$Replicate)


  write.csv(sal[,export_cols], file.path(cache,"salinity_curated_rows.csv"),row.names=FALSE)
  salinity_rows <- sal

d <- salinity_rows
required <- c(
  "Experiment_id", "Strain", "Species", "reproductive_mode",
  "Sampling location", "Temperature_treatment", "Observation_day_raw",
  "Salinity [PSU]", "Replicate", "experimental_unit_id",
  "Analysis_dead", "Analysis_not_dead", "Total",
  "Dormant_stage", "Paralysed_adult", "Analysis_endpoint_definition",
  "metadata_status", "analysis_role", "Mortality_analysis_eligible",
  "analysis_included_primary", "Source_file", "Source_sheet"
)
missing_columns <- setdiff(required, names(d))
if (length(missing_columns)) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

as_logical_clean <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x %in% c("TRUE", "T", "1")
}

# Keep all documented dose-response observation days, but remove the failed
# preliminary ES5 run and rows excluded by the existing source audit.
keep <- as_logical_clean(d$Mortality_analysis_eligible) &
  as_logical_clean(d$analysis_included_primary) &
  d$analysis_role != "exploratory_only" &
  d$metadata_status != "failed_preliminary" &
  d$Experiment_id != "ES5_salinity_preliminary_R1"

x <- d[keep, , drop = FALSE]

x$endpoint_class <- ifelse(
  grepl("immobile or missing", x$Analysis_endpoint_definition, ignore.case = TRUE),
  "immobile_or_missing_classified_as_dead",
  ifelse(
    grepl("paralysed individuals are counted as dead|Dead_with_paralysed",
          x$Analysis_endpoint_definition, ignore.case = TRUE),
    "dead_or_paralysed",
    ifelse(
      grepl("Original 20-adult cohort", x$Analysis_endpoint_definition,
            ignore.case = TRUE),
      "confirmed_dead_dormant_not_dead",
      "corrected_source_mortality"
    )
  )
)

x$effect_n <- as.numeric(x$Analysis_dead)
x$unaffected_n <- as.numeric(x$Analysis_not_dead)
x$total_n <- as.numeric(x$Total)
x$salinity_PSU <- as.numeric(x[["Salinity [PSU]"]])
x$exposure_day <- as.numeric(x$Observation_day_raw)
x$temperature_C <- as.numeric(x$Temperature_treatment)
x$plate_id <- x$experimental_unit_id
x$experimental_block <- x$Experiment_id

if (any(!is.finite(x$effect_n)) || any(!is.finite(x$unaffected_n)) ||
    any(!is.finite(x$total_n))) stop("Non-finite endpoint counts detected.")
if (any(x$effect_n < 0) || any(x$unaffected_n < 0) || any(x$total_n <= 0)) {
  stop("Endpoint counts outside their permitted ranges.")
}
if (any(x$effect_n + x$unaffected_n != x$total_n)) {
  stop("effect_n + unaffected_n does not equal total_n in every row.")
}

out <- data.frame(
  strain = x$Strain,
  species = x$Species,
  reproductive_mode = x$reproductive_mode,
  sampling_location = x[["Sampling location"]],
  experimental_block = x$experimental_block,
  plate_id = x$plate_id,
  replicate = x$Replicate,
  temperature_C = x$temperature_C,
  exposure_day = x$exposure_day,
  salinity_PSU = x$salinity_PSU,
  effect_n = x$effect_n,
  unaffected_n = x$unaffected_n,
  total_n = x$total_n,
  effect_fraction = x$effect_n / x$total_n,
  dormant_n = x$Dormant_stage,
  paralysed_n = x$Paralysed_adult,
  endpoint_class = x$endpoint_class,
  endpoint_definition = x$Analysis_endpoint_definition,
  source_file = x$Source_file,
  source_sheet = x$Source_sheet,
  stringsAsFactors = FALSE
)

out <- out[order(
  out$strain, out$experimental_block, out$temperature_C,
  out$exposure_day, out$salinity_PSU, out$replicate
), ]


# Recreate exactly the endpoint table used for the model. The original model
# table omitted the two descriptive fields below, which do not enter fitting.
all_time <- out[, setdiff(names(out),c("reproductive_mode","sampling_location"))]
# Preserve the historical row order for optimizer and plotting reproducibility.
all_time <- all_time[order(all_time$strain, all_time$experimental_block,
                          all_time$exposure_day,all_time$salinity_PSU,all_time$replicate),]

is_series2023 <- startsWith(all_time$experimental_block, "series2023_")
all_time <- all_time %>%
  mutate(
    dead_inactive_n = effect_n,
    dead_inactive_n = if_else(is_series2023,
                              effect_n + coalesce(as.numeric(dormant_n), 0),
                              dead_inactive_n),
    active_n = total_n - dead_inactive_n,
    analysis_endpoint = if_else(is_series2023,
                                "confirmed dead + recorded dormant",
                                "source-defined dead/inactive"),
    analysis_endpoint = if_else(endpoint_class == "immobile_or_missing_classified_as_dead",
                                "immobile or missing classified as dead",
                                analysis_endpoint)
  )
stopifnot(all_time$dead_inactive_n >= 0, all_time$active_n >= 0)
write_csv(all_time, file.path(cache, "R_blockwise_dead_cryptobiotic_analysis_ready.csv"))

## Figure 2 input: reconstruct confirmed death and cryptobiosis directly from
## the R-harmonized salinity table.
sal <- salinity_rows
keep <- as_logical_clean(sal$Mortality_analysis_eligible) &
  as_logical_clean(sal$analysis_included_primary) &
  sal$analysis_role != "exploratory_only" &
  sal$metadata_status != "failed_preliminary" &
  sal$Experiment_id != "ES5_salinity_preliminary_R1" &
  sal$Strain %in% c("ALT.22.08", "PAP.22.17") &
  as.numeric(sal$Observation_day_raw) == 4
x <- sal[keep, , drop = FALSE]
num <- function(v) suppressWarnings(as.numeric(v))
x <- x %>% mutate(
  Analysis_dead = num(Analysis_dead), Total = num(Total),
  Dormant_stage = pmax(coalesce(num(Dormant_stage), 0), 0),
  Paralysed_adult = pmax(coalesce(num(Paralysed_adult), 0), 0),
  Dead_without_dormant = num(Dead_without_dormant),
  Temperature_treatment = num(Temperature_treatment),
  Observation_day_raw = num(Observation_day_raw),
  salinity = num(`Salinity [PSU]`),
  definition = coalesce(Analysis_endpoint_definition, ""),
  series2023 = startsWith(Experiment_id, "series2023_"),
  dead_inactive_n = if_else(series2023, Analysis_dead + Dormant_stage, Analysis_dead),
  confirmed_dead_n = if_else(Dead_without_dormant >= 0, Dead_without_dormant, NA_real_)
)
uses_paralysed <- grepl("paralysed individuals are counted as dead|Dead_with_paralysed",
                        x$definition, ignore.case = TRUE)
corrected <- grepl("corrected source summary", x$definition, ignore.case = TRUE)
idx <- uses_paralysed & is.na(x$confirmed_dead_n)
x$confirmed_dead_n[idx] <- x$Analysis_dead[idx] - x$Paralysed_adult[idx]
idx <- corrected & is.na(x$confirmed_dead_n)
x$confirmed_dead_n[idx] <- x$Analysis_dead[idx] - x$Dormant_stage[idx]
idx <- x$series2023 & is.na(x$confirmed_dead_n)
x$confirmed_dead_n[idx] <- x$Analysis_dead[idx]
if (anyNA(x$confirmed_dead_n)) stop("Confirmed deaths could not be reconstructed in R.")
x <- x %>% mutate(
  confirmed_dead_n = pmax(confirmed_dead_n, 0),
  dead_inactive_n = pmax(dead_inactive_n, confirmed_dead_n),
  cryptobiotic_n = dead_inactive_n - confirmed_dead_n,
  affected_n = confirmed_dead_n + cryptobiotic_n
)
detail <- x %>% transmute(
  Strain, Temperature_treatment, Observation_day_raw,
  `Salinity [PSU]` = salinity, Experiment_id, experimental_unit_id, Replicate,
  confirmed_dead_n, inactive_or_unconfirmed_n = cryptobiotic_n, affected_n, Total
)
write_csv(detail, file.path(cache, "R_day4_dead_cryptobiotic_analysis_ready.csv"))
pooled <- detail %>%
  group_by(Strain, Temperature_treatment, `Salinity [PSU]`) %>%
  summarise(
    confirmed_dead_n = sum(confirmed_dead_n),
    inactive_or_unconfirmed_n = sum(inactive_or_unconfirmed_n),
    dead_inactive_n = sum(affected_n), Total = sum(Total),
    observation_days = paste(sort(unique(Observation_day_raw)), collapse = "; "),
    .groups = "drop"
  )
write_csv(pooled, file.path(cache, "R_day4_dead_cryptobiotic_pooled_counts.csv"))


day4_rows <- detail
blockwise <- all_time
write.csv(blockwise,file.path(cache,"blockwise_dead_inactive_analysis_ready.csv"),row.names=FALSE)
# Refit selected global response model from curated master counts, never from
# the saved reference model. The cool ES5 group uses the confirmed 16C label.
d <- blockwise
d$group <- factor(paste0(d$strain,"__",d$temperature_C,"C"))
d$salinity_10 <- d$salinity_PSU/10
day_mean <- mean(d$exposure_day)
d$day_c <- d$exposure_day-day_mean
d$source_block <- factor(d$experimental_block)
d$plate_uid <- factor(paste(d$experimental_block,d$plate_id,sep="__"))
d$active_n <- d$total_n-d$dead_inactive_n
model_formula <- cbind(dead_inactive_n,active_n) ~ group*salinity_10*day_c +
                    (1|source_block)+(1|plate_uid)
message("Salinity: refitting selected GLMM")
model <- glmmTMB(model_formula,data=d,family=binomial("logit"),
                  control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS")))
stopifnot(model$fit$convergence==0L,isTRUE(model$sdr$pdHess))
saveRDS(model,file.path(cache,"global_salinity_model.rds"))
beta <- fixef(model)$cond
V <- vcov(model)$cond
write.csv(data.frame(term=names(beta),estimate=beta,SE=sqrt(diag(V))),
          file.path(stats,"global_model_coefficients.csv"),row.names=FALSE)
write.csv(data.frame(AIC=AIC(model),convergence=model$fit$convergence,pdHess=model$sdr$pdHess),
          file.path(stats,"global_model_validation.csv"),row.names=FALSE)

set.seed(1211)

beta <- fixef(model)$cond
cov_beta <- vcov(model)$cond
draws <- MASS::mvrnorm(2000, mu = beta, Sigma = cov_beta)

fixed_formula <- ~ group * salinity_10 * day_c

marginal_curve <- function(group_name, salinity_values, coefficients) {
  observed_days <- sort(unique(d$exposure_day[d$group == group_name]))
  nd <- expand.grid(
    salinity_PSU = salinity_values,
    exposure_day = observed_days,
    KEEP.OUT.ATTRS = FALSE
  )
  nd$group <- factor(group_name, levels = levels(d$group))
  nd$salinity_10 <- nd$salinity_PSU / 10
  nd$day_c <- nd$exposure_day - day_mean
  mm <- model.matrix(fixed_formula, nd)
  mm <- mm[, names(beta), drop = FALSE]
  probability <- plogis(drop(mm %*% coefficients))
  aggregate(probability, list(salinity_PSU = nd$salinity_PSU), mean)$x
}

marginal_curve_draws <- function(group_name, salinity_values, coefficient_draws) {
  observed_days <- sort(unique(d$exposure_day[d$group == group_name]))
  nd <- expand.grid(
    salinity_PSU = salinity_values,
    exposure_day = observed_days,
    KEEP.OUT.ATTRS = FALSE
  )
  nd$group <- factor(group_name, levels = levels(d$group))
  nd$salinity_10 <- nd$salinity_PSU / 10
  nd$day_c <- nd$exposure_day - day_mean
  mm <- model.matrix(fixed_formula, nd)[, names(beta), drop = FALSE]
  n_grid <- length(salinity_values)
  result <- matrix(0, nrow = n_grid, ncol = nrow(coefficient_draws))
  for (i in seq_along(observed_days)) {
    index <- ((i - 1) * n_grid + 1):(i * n_grid)
    result <- result + plogis(mm[index, , drop = FALSE] %*% t(coefficient_draws)) /
      length(observed_days)
  }
  result
}

find_ec50 <- function(grid, probability) {
  delta <- probability - 0.5
  crossing <- which(delta[-length(delta)] * delta[-1] <= 0)
  if (!length(crossing)) return(NA_real_)
  i <- crossing[1]
  approx(probability[c(i, i + 1)], grid[c(i, i + 1)], xout = 0.5)$y
}

rows <- list()
curve_rows <- list()

for (group_name in levels(d$group)) {
  z <- d[d$group == group_name, ]
  tested_min <- min(z$salinity_PSU)
  tested_max <- max(z$salinity_PSU)
  broad_grid <- seq(-20, max(100, tested_max * 2), length.out = 2401)
  plotted_grid <- seq(tested_min, tested_max, length.out = 301)

  point_probability <- marginal_curve(group_name, broad_grid, beta)
  point_ec50 <- find_ec50(broad_grid, point_probability)

  draw_probabilities <- marginal_curve_draws(group_name, broad_grid, draws)
  draw_ec50 <- apply(
    draw_probabilities, 2,
    function(probability) find_ec50(broad_grid, probability)
  )
  finite_draws <- draw_ec50[is.finite(draw_ec50)]
  interval <- if (length(finite_draws) >= 100) {
    quantile(finite_draws, c(0.025, 0.975), names = FALSE)
  } else {
    c(NA_real_, NA_real_)
  }

  observed <- aggregate(
    cbind(dead_inactive_n, total_n) ~ salinity_PSU,
    data = z,
    FUN = sum
  )
  observed$fraction <- observed$dead_inactive_n / observed$total_n
  crosses_observed <- min(observed$fraction) <= 0.5 && max(observed$fraction) >= 0.5
  in_range <- is.finite(point_ec50) && point_ec50 >= tested_min && point_ec50 <= tested_max

  parts <- strsplit(group_name, "__", fixed = TRUE)[[1]]
  rows[[length(rows) + 1]] <- data.frame(
    strain_temperature = group_name,
    strain = parts[1],
    temperature_C = as.numeric(sub("C$", "", parts[2])),
    observed_days = paste(sort(unique(z$exposure_day)), collapse = "; "),
    source_blocks = paste(sort(unique(z$experimental_block)), collapse = "; "),
    tested_salinity_min_PSU = tested_min,
    tested_salinity_max_PSU = tested_max,
    marginal_EC50_PSU = point_ec50,
    CI_low = interval[1],
    CI_high = interval[2],
    simulation_draws_with_crossing = length(finite_draws),
    observed_crosses_50_percent = crosses_observed,
    EC50_within_tested_range = in_range,
    reportable = crosses_observed && in_range
  )

  pred <- marginal_curve(group_name, plotted_grid, beta)
  plotted_draws <- marginal_curve_draws(group_name, plotted_grid, draws)
  pred_low <- apply(plotted_draws, 1, quantile, probs = 0.025, na.rm = TRUE)
  pred_high <- apply(plotted_draws, 1, quantile, probs = 0.975, na.rm = TRUE)
  curve_rows[[length(curve_rows) + 1]] <- data.frame(
    strain_temperature = group_name,
    salinity_PSU = plotted_grid,
    predicted_dead_inactive = pred,
    CI_low = pred_low,
    CI_high = pred_high
  )
}

results <- do.call(rbind, rows)
curves <- do.call(rbind, curve_rows)
write.csv(
  results,
  file.path(cache, "global_model_marginal_EC50_by_strain_temperature.csv"),
  row.names = FALSE
)
write.csv(
  curves,
  file.path(cache, "global_model_marginal_dose_response_curves.csv"),
  row.names = FALSE
)

print(results, row.names = FALSE)

beta <- fixef(model)$cond
V <- vcov(model)$cond
coef_names <- names(beta)

groups <- c(
  "A. max__22C", "ALT.22.08__16C", "ALT.22.08__22C",
  "ARO.22.05__22C", "ES5__16C", "ES5__22C",
  "LAG.22.03__22C", "PAP.22.17__16C", "PAP.22.17__22C"
)
reference_group <- "A. max__22C"

# Return the four coefficients defining a group's response surface:
# intercept, salinity slope, day slope, and salinity-by-day interaction.
surface_matrix <- function(group) {
  M <- matrix(0, nrow = 4, ncol = length(beta),
              dimnames = list(c("intercept", "salinity", "day", "salinity_day"),
                              coef_names))
  M["intercept", "(Intercept)"] <- 1
  M["salinity", "salinity_10"] <- 1
  M["day", "day_c"] <- 1
  M["salinity_day", "salinity_10:day_c"] <- 1

  if (group != reference_group) {
    M["intercept", paste0("group", group)] <- 1
    M["salinity", paste0("group", group, ":salinity_10")] <- 1
    M["day", paste0("group", group, ":day_c")] <- 1
    M["salinity_day", paste0("group", group, ":salinity_10:day_c")] <- 1
  }
  M
}

wald_pair <- function(group_1, group_2, family, terms = rownames(surface_matrix(group_1))) {
  L <- (surface_matrix(group_1) - surface_matrix(group_2))[terms, , drop = FALSE]
  estimate <- as.vector(L %*% beta)
  LVL <- L %*% V %*% t(L)
  rank_L <- qr(LVL)$rank
  statistic <- as.numeric(t(estimate) %*% MASS::ginv(LVL) %*% estimate)
  data.frame(
    comparison_family = family,
    group_1 = group_1,
    group_2 = group_2,
    chi_square = statistic,
    df = rank_L,
    p_value = pchisq(statistic, df = rank_L, lower.tail = FALSE)
  )
}

temperature_pairs <- list(
  c("ALT.22.08__22C", "ALT.22.08__16C"),
  c("PAP.22.17__22C", "PAP.22.17__16C"),
  c("ES5__22C", "ES5__16C")
)
strain_22 <- c(
  "A. max__22C", "ALT.22.08__22C", "ARO.22.05__22C",
  "ES5__22C", "LAG.22.03__22C", "PAP.22.17__22C"
)
strain_pairs <- combn(strain_22, 2, simplify = FALSE)

salinity_results <- do.call(rbind, c(
  lapply(temperature_pairs, function(z) wald_pair(
    z[1], z[2], "temperature_within_strain", c("salinity", "salinity_day")
  )),
  lapply(strain_pairs, function(z) wald_pair(
    z[1], z[2], "strain_at_22C", c("salinity", "salinity_day")
  ))
))
salinity_results$p_adj_Holm <- ave(
  salinity_results$p_value, salinity_results$comparison_family,
  FUN = function(p) p.adjust(p, method = "holm")
)
salinity_results$significant_Holm_0.05 <- salinity_results$p_adj_Holm < 0.05

write.csv(
  salinity_results,
  file.path(stats, "pairwise_salinity_effect_tests.csv"),
  row.names = FALSE
)

print(salinity_results, row.names = FALSE)

# Retain exact manuscript bootstrap draws as an explicit reference replay.
if (bootstrap_mode == "bootstrap_reference_replay") {
  boot <- read.csv(file.path(cur,"recorded_ec50_bootstrap_draws.csv"),check.names=FALSE)
  names(boot) <- sub("ES5__15C", "ES5__16C", names(boot), fixed=TRUE)
  stopifnot(nrow(boot)==2000L)
} else {
  boot <- local({
    d <- blockwise
    d$group <- factor(paste0(d$strain,"__",d$temperature_C,"C"))
    group_levels <- levels(d$group)
    d$salinity_10 <- d$salinity_PSU/10
    day_mean <- mean(d$exposure_day)
    d$day_c <- d$exposure_day-day_mean
    d$source_block <- factor(d$experimental_block)
    d$plate_uid <- factor(paste(d$experimental_block,d$plate_id,sep="__"))
    d$active_n <- d$total_n-d$dead_inactive_n

grid <- seq(-30, 100, length.out = 1301)
fixed_formula <- ~ group * salinity_10 * day_c
observed_days <- lapply(group_levels, function(g) sort(unique(d$exposure_day[d$group == g])))
names(observed_days) <- group_levels

resample_plates <- function(data) {
  pieces <- list()
  k <- 0L
  for (block in unique(data$experimental_block)) {
    z <- data[data$experimental_block == block, , drop = FALSE]
    plates <- unique(z$plate_uid)
    sampled <- sample(plates, length(plates), replace = TRUE)
    for (j in seq_along(sampled)) {
      zz <- z[z$plate_uid == sampled[j], , drop = FALSE]
      zz$plate_uid <- paste0(block, "__bootstrap_", j)
      k <- k + 1L
      pieces[[k]] <- zz
    }
  }
  out <- do.call(rbind, pieces)
  out$group <- factor(out$group, levels = group_levels)
  out$source_block <- factor(out$source_block, levels = levels(d$source_block))
  out$plate_uid <- factor(out$plate_uid)
  out
}

find_ec50 <- function(probability) {
  delta <- probability - 0.5
  crossing <- which(delta[-length(delta)] * delta[-1] <= 0)
  if (!length(crossing)) return(NA_real_)
  i <- crossing[1]
  approx(probability[c(i, i + 1)], grid[c(i, i + 1)], xout = 0.5)$y
}

ec50_from_fit <- function(model) {
  beta <- fixef(model)$cond
  result <- setNames(rep(NA_real_, length(group_levels)), group_levels)
  for (g in group_levels) {
    days <- observed_days[[g]]
    nd <- expand.grid(salinity_PSU = grid, exposure_day = days,
                      KEEP.OUT.ATTRS = FALSE)
    nd$group <- factor(g, levels = group_levels)
    nd$salinity_10 <- nd$salinity_PSU / 10
    nd$day_c <- nd$exposure_day - day_mean
    mm <- model.matrix(fixed_formula, nd)[, names(beta), drop = FALSE]
    n_grid <- length(grid)
    probability <- rep(0, n_grid)
    for (i in seq_along(days)) {
      index <- ((i - 1) * n_grid + 1):(i * n_grid)
      probability <- probability + plogis(drop(mm[index, , drop = FALSE] %*% beta)) /
        length(days)
    }
    result[g] <- find_ec50(probability)
  }
  result
}

bootstrap_one <- function(i) {
  bd <- resample_plates(d)
  fit <- tryCatch(
    suppressWarnings(glmmTMB(
      model_formula, data = bd, family = binomial(link = "logit"),
      control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
    )),
    error = function(e) NULL
  )
  if (is.null(fit) || !identical(fit$fit$convergence, 0L) || !isTRUE(fit$sdr$pdHess)) {
    return(c(success = 0, setNames(rep(NA_real_, length(group_levels)), group_levels)))
  }
  c(success = 1, ec50_from_fit(fit))
}


    ans <- vector("list", n_boot)
    for (i in seq_len(n_boot)) {
      set.seed(1211L+i) # per-replicate, independent of number/speed of workers
      ans[[i]] <- bootstrap_one(i)
      if(i %% 10L==0L) message("Fresh EC50 bootstrap ",i,"/",n_boot)
    }
    z <- as.data.frame(do.call(rbind,ans),check.names=FALSE)
    z$bootstrap_replicate <- seq_len(n_boot)
    z
  })
}
write.csv(boot,file.path(cache,"ec50_bootstrap_draws.csv"),row.names=FALSE)
point <- read.csv(file.path(cache,"global_model_marginal_EC50_by_strain_temperature.csv"))
ec <- do.call(rbind,lapply(point$strain_temperature,function(g){
  vals <- boot[[g]][boot$success==1 & is.finite(boot[[g]])]
  ci <- quantile(pmax(0,vals),c(.025,.975),names=FALSE,type=7)
  data.frame(strain_temperature=g,
             point_EC50_PSU=point$marginal_EC50_PSU[point$strain_temperature==g],
             nonnegative_bootstrap_CI_low=ci[1], nonnegative_bootstrap_CI_high=ci[2],
             valid_EC50_replicates=length(vals), requested_replicates=nrow(boot),
             bootstrap_mode=bootstrap_mode)
}))
write.csv(ec,file.path(stats,"cluster_bootstrap_marginal_EC50_nonnegative_summary.csv"),row.names=FALSE)

d <- day4_rows

names(d)[names(d) == "Temperature_treatment"] <- "temperature_C"
names(d)[names(d) == "Observation_day_raw"] <- "exposure_day"
names(d)[names(d) == "Salinity [PSU]"] <- "salinity_PSU"
d <- subset(d, exposure_day == 4)
d$temperature <- factor(d$temperature_C, levels = c(16, 22))
d$salinity_10 <- d$salinity_PSU / 10
stopifnot(!anyNA(d), !anyDuplicated(d$experimental_unit_id),
          all(d$confirmed_dead_n >= 0), all(d$inactive_or_unconfirmed_n >= 0),
          all(d$confirmed_dead_n + d$inactive_or_unconfirmed_n <= d$Total))
models <- list(); rows <- list(); counts <- list()
run_test <- function(strain, outcome, column) {
  x <- d[d$Strain == strain, ]
  ranges <- lapply(split(x$salinity_PSU, x$temperature), range)
  low <- max(vapply(ranges, `[`, numeric(1), 1))
  high <- min(vapply(ranges, `[`, numeric(1), 2))
  x <- subset(x, salinity_PSU >= low & salinity_PSU <= high)
  x$success <- x[[column]]; x$failure <- x$Total - x$success
  full <- glm(cbind(success, failure) ~ temperature * salinity_10,
              family = quasibinomial("logit"), data = x)
  additive <- update(full, . ~ temperature + salinity_10)
  salinity_only <- update(full, . ~ salinity_10)
  stopifnot(full$converged, additive$converged, salinity_only$converged,
            full$rank == 4, nobs(full) == nrow(x),
            all(is.finite(coef(summary(full)))))
  a <- anova(additive, full, test = "F")
  j <- anova(salinity_only, full, test = "F")
  # Independent arithmetic cross-check of the nested quasi-F statistic.
  manual <- (deviance(additive) - deviance(full)) / summary(full)$dispersion
  stopifnot(isTRUE(all.equal(unname(a$F[2]), unname(manual), tolerance = 1e-8)))
  key <- paste(strain, outcome, sep = "_")
  models[[key]] <<- list(full = full, additive = additive)
  rows[[key]] <<- x
  counts[[key]] <<- transform(aggregate(cbind(success, failure) ~ temperature_C + salinity_PSU, x, sum),
                              strain = strain, outcome = outcome)
  b <- coef(full); cs <- coef(summary(full))
  data.frame(strain = strain, outcome = outcome, min_PSU = low, max_PSU = high,
    n = nrow(x), F = a$F[2], df1 = a$Df[2], df2 = df.residual(full),
    p_value = a$`Pr(>F)`[2], dispersion = summary(full)$dispersion,
    slope16 = unname(b["salinity_10"]),
    slope22 = unname(b["salinity_10"] + b["temperature22:salinity_10"]),
    interaction_SE = cs["temperature22:salinity_10", "Std. Error"],
    joint_F = j$F[2], joint_p = j$`Pr(>F)`[2])
}
res <- rbind(run_test("ALT.22.08", "Inactive", "inactive_or_unconfirmed_n"),
             run_test("PAP.22.17", "Inactive", "inactive_or_unconfirmed_n"),
             run_test("ALT.22.08", "Dead", "confirmed_dead_n"),
             run_test("PAP.22.17", "Dead", "confirmed_dead_n"))
res$p_adj_Holm <- p.adjust(res$p_value, "holm")
res$joint_p_adj_Holm <- p.adjust(res$joint_p, "holm")
write.csv(res, file.path(stats, "interaction_tests.csv"), row.names = FALSE)
write.csv(do.call(rbind, counts), file.path(stats, "counts_by_dose.csv"), row.names = FALSE)
saveRDS(list(models = models, analysis_rows = rows), file.path(cache, "day4_models_and_rows.rds"))
capture.output(sessionInfo(), file = file.path(stats, "sessionInfo.txt"))
print(res,row.names=FALSE)

d <- day4_rows

names(d)[names(d) == "Temperature_treatment"] <- "temperature_C"
names(d)[names(d) == "Observation_day_raw"] <- "exposure_day"
names(d)[names(d) == "Salinity [PSU]"] <- "salinity_PSU"

# This analysis concerns the composition of affected animals at the only
# observation day shared by all four strain-temperature groups.
d <- d[d$exposure_day == 4, ]
# Restrict the primary comparison to doses represented in every one of the
# four strain-temperature groups (ALT/PAP at 16/22 degrees C).
shared_salinity <- Reduce(
  intersect,
  split(d$salinity_PSU, interaction(d$Strain, d$temperature_C, drop = TRUE))
)
d <- d[d$salinity_PSU %in% shared_salinity, ]
d$affected_n <- d$confirmed_dead_n + d$inactive_or_unconfirmed_n
d <- d[d$affected_n > 0, ]
d$strain <- relevel(factor(d$Strain), ref = "ALT.22.08")
d$temperature <- relevel(factor(d$temperature_C), ref = "16")
d$salinity_10 <- d$salinity_PSU / 10
d$salinity_factor <- factor(d$salinity_PSU)

# Primary adjusted comparison. It asks whether the temperature-related change
# in the dead/inactive composition differs between ALT and PAP after adjustment
# for the four matched salinity levels.
adjusted_full <- glm(
  cbind(confirmed_dead_n, inactive_or_unconfirmed_n) ~
    salinity_factor + strain * temperature,
  data = d,
  family = quasibinomial(link = "logit")
)
adjusted_reduced <- update(adjusted_full, . ~ . - strain:temperature)
adjusted_test <- anova(adjusted_reduced, adjusted_full, test = "F")


composition <- data.frame(F=adjusted_test$F[2],df1=adjusted_test$Df[2],
  df2=df.residual(adjusted_full),dispersion=summary(adjusted_full)$dispersion,
  p_value=adjusted_test$`Pr(>F)`[2])
write.csv(composition,file.path(stats,"day4_dead_vs_inactive_global_interaction.csv"),row.names=FALSE)

m9_ids <- c("A.max_22C", "ALT.22.08_16C", "ALT.22.08_22C", "PAP.22.17_16C")
m9_raw <- read.csv(file.path(input_folder, "recovery_observations.csv"), check.names = FALSE, na.strings = "") %>%
  mutate(
    Day = suppressWarnings(as.numeric(Day)),
    salinity_PSU = suppressWarnings(as.numeric(sub("^([0-9.]+).*", "\\1", as.character(Treatment)))),
    temperature_C = suppressWarnings(as.numeric(Temperature_treatment_C)),
    replicate_source = suppressWarnings(as.numeric(Replicate)),
    recovered_n = suppressWarnings(as.numeric(Rehyadrated)),
    cryptobiotic_candidates_n = rowSums(cbind(
      suppressWarnings(as.numeric(Doramnt_stage)),
      suppressWarnings(as.numeric(Paralaysed_adult)),
      suppressWarnings(as.numeric(`Paralaysed_juvenil`))
    ), na.rm = TRUE),
    initial_n = ifelse(grepl("^ALT\\.22\\.08", Experiment_id), 30, 20)
  ) %>%
  filter(Experiment_id %in% m9_ids, Day == 14)

if (nrow(m9_raw) != 120L || any(is.na(m9_raw$salinity_PSU)) || any(is.na(m9_raw$recovered_n))) {
  stop("M9 day-14 source must contain 120 complete plate records with salinity and recovered counts.")
}
if (!all(table(m9_raw$Experiment_id) == 30L)) stop("Each M9 experiment must contain 30 plates.")
if (any(m9_raw$recovered_n < 0 | m9_raw$recovered_n > m9_raw$initial_n)) stop("M9 recovered counts fall outside 0..initial_n.")

# One original record is internally inconsistent: ten animals were recorded as
# recovered although only four animals had been classified as coiled or
# paralysed before M9 addition. It is retained in the source workbook but
# excluded from candidate-based recovery analyses.
m9_exclusions <- m9_raw %>%
  filter(
    Experiment_id == "ALT.22.08_16C",
    salinity_PSU == 17,
    replicate_source == 1
  ) %>%
  transmute(
    Experiment_id, salinity_PSU, replicate_source,
    cryptobiotic_candidates_n, recovered_n,
    exclusion_reason = "Recovered count exceeds the number of pre-M9 cryptobiotic candidates"
  )
if (nrow(m9_exclusions) != 1L) stop("Expected exactly one predefined inconsistent M9 record.")

m9_raw <- m9_raw %>%
  anti_join(
    m9_exclusions %>% dplyr::select(Experiment_id, salinity_PSU, replicate_source),
    by = c("Experiment_id", "salinity_PSU", "replicate_source")
  )

m9_tested_grid <- m9_raw %>%
  group_by(Experiment_id, Strain = Strain_standardized, temperature_C, salinity_PSU) %>%
  summarise(
    tested_plates = n(),
    cryptobiotic_candidates = sum(cryptobiotic_candidates_n),
    has_cryptobiotic_candidates = cryptobiotic_candidates > 0,
    .groups = "drop"
  )

m9_raw <- m9_raw %>% filter(cryptobiotic_candidates_n > 0)

if (any(m9_raw$recovered_n > m9_raw$cryptobiotic_candidates_n)) {
  stop("M9 recovered counts exceed the number of pre-M9 cryptobiotic candidates.")
}

m9 <- m9_raw %>%
  group_by(Experiment_id, salinity_PSU) %>%
  mutate(replicate_clean = row_number()) %>%
  ungroup() %>%
  transmute(
    Experiment_id,
    Strain = Strain_standardized,
    temperature_C,
    salinity_PSU,
    replicate_clean,
    replicate_source,
    plate_id = paste0(Experiment_id, "|PSU=", salinity_PSU, "|rep=", replicate_clean),
    initial_n,
    cryptobiotic_candidates_n,
    recovered_n,
    recovered_fraction_candidate = recovered_n / cryptobiotic_candidates_n,
    recovery_present = recovered_n > 0,
    exposed = salinity_PSU > 0
  )

m9_dose_summary <- m9 %>%
  group_by(Experiment_id, Strain, temperature_C, salinity_PSU) %>%
  summarise(
    plates = n(),
    plates_with_recovery = sum(recovery_present),
    recovered_animals = sum(recovered_n),
    cryptobiotic_candidates = sum(cryptobiotic_candidates_n),
    recovered_fraction_candidate = recovered_animals / cryptobiotic_candidates,
    .groups = "drop"
  )
m9_ci <- wilson(m9_dose_summary$plates_with_recovery, m9_dose_summary$plates)
m9_dose_summary$recovery_plate_fraction <- m9_ci[, "proportion"]
m9_dose_summary$plate_wilson_95_low <- m9_ci[, "low"]
m9_dose_summary$plate_wilson_95_high <- m9_ci[, "high"]

m9_experiment_summary <- m9 %>%
  group_by(Experiment_id, Strain, temperature_C) %>%
  summarise(
    plates = n(),
    plates_with_recovery = sum(recovery_present),
    recovered_animals = sum(recovered_n),
    exposed_plates = sum(exposed),
    exposed_plates_with_recovery = sum(recovery_present & exposed),
    exposed_recovered_animals = sum(recovered_n[exposed]),
    control_recovered_animals = sum(recovered_n[!exposed]),
    .groups = "drop"
  )

alt_m9 <- m9 %>% filter(grepl("^ALT\\.22\\.08", Experiment_id), exposed)
alt_m9_table <- table(
  factor(alt_m9$temperature_C, levels = c(16, 22)),
  factor(alt_m9$recovery_present, levels = c(FALSE, TRUE))
)
alt_m9_fisher <- fisher.test(alt_m9_table)
alt_m9_by_temperature <- alt_m9 %>%
  group_by(temperature_C) %>%
  summarise(plates = n(), plates_with_recovery = sum(recovery_present), plate_fraction = mean(recovery_present), .groups = "drop")
alt_m9_test <- tibble(
  comparison = "ALT.22.08 exposed plates: M9 recovery occurrence at 16 vs 22 C",
  recovery_plates_16C = alt_m9_by_temperature$plates_with_recovery[alt_m9_by_temperature$temperature_C == 16],
  total_plates_16C = alt_m9_by_temperature$plates[alt_m9_by_temperature$temperature_C == 16],
  recovery_plates_22C = alt_m9_by_temperature$plates_with_recovery[alt_m9_by_temperature$temperature_C == 22],
  total_plates_22C = alt_m9_by_temperature$plates[alt_m9_by_temperature$temperature_C == 22],
  risk_difference_16C_minus_22C = alt_m9_by_temperature$plate_fraction[alt_m9_by_temperature$temperature_C == 16] - alt_m9_by_temperature$plate_fraction[alt_m9_by_temperature$temperature_C == 22],
  fisher_p_value = alt_m9_fisher$p.value
)

write.csv(m9, file.path(stats, "m9_recovery_analysis_ready.csv"), row.names = FALSE, na = "")
write.csv(m9_dose_summary, file.path(stats, "m9_recovery_by_salinity.csv"), row.names = FALSE, na = "")
write.csv(m9_experiment_summary, file.path(stats, "m9_recovery_experiment_summary.csv"), row.names = FALSE, na = "")
write.csv(alt_m9_test, file.path(stats, "m9_ALT_temperature_fisher.csv"), row.names = FALSE, na = "")
write.csv(m9_exclusions, file.path(stats, "m9_recovery_exclusions.csv"), row.names = FALSE, na = "")
write.csv(m9_tested_grid, file.path(stats, "m9_recovery_tested_salinity_grid.csv"), row.names = FALSE, na = "")

if (render_figures) {
  # A clean graphics process avoids session-dependent glyph anti-aliasing after
  # fitting mixed models. The child reads only freshly generated analysis CSVs.
  renderer <- file.path(root, "R", "02_salinity.R")
  rscript <- file.path(R.home("bin"), "Rscript.exe")
  if (!file.exists(rscript)) rscript <- file.path(R.home("bin"), "Rscript")
  status <- system2(rscript,
    args = c("--vanilla", shQuote(renderer), "--render-only", shQuote(root)))
  if (!identical(as.integer(status), 0L)) {
    stop("Clean-process salinity figure rendering failed with status ", status)
  }
}

  writeLines(c("Selected GLMM refit from curated master workbook.",
    paste("EC50 bootstrap mode:",bootstrap_mode),
    "Replay means derived reference draws, not a fresh bootstrap.",
    "ES5 cool exposure and group identifier consistently labelled 16C.",
    "Recovery statistics are recomputed here; artwork is rebuilt by run_supplement_figures."),
    file.path(stats,"REPRODUCTION_MODE.txt"))
  capture.output(sessionInfo(),file=file.path(stats,"sessionInfo.txt"))
  invisible(list(model=model,ec50=ec,interaction_tests=res,composition=composition))
}


# Rendering uses the original plotting-library order in a fresh R process.
# This function neither refits the scientific GLMM nor reads reference images.
render_salinity_figures <- function(project_root) {
  suppressPackageStartupMessages({
    library(ggplot2); library(patchwork); library(dplyr)
    library(readr); library(scales); library(cowplot)
  })
  root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  sys.source(file.path(root,'R/01_figure_style.R'),envir=environment())
  cache <- file.path(root, "Results", "cache", "salinity")
  stats <- file.path(root, "Results", "statistics", "salinity")
  figures <- file.path(root, "Manuscript", "figures")
  dir.create(figures, recursive = TRUE, showWarnings = FALSE)

  here <- figures
  dose <- cache
  red <- "#CC6677"; blue <- "#4477AA"; grey_strip <- "#D9D9D9"

theme_manuscript <- function(base_size = 9) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(colour = "#E6E6E6", linewidth = .3),
      strip.background = element_rect(fill = grey_strip, colour = "#333333", linewidth = .5),
      strip.text = element_text(face = "plain", size = base_size),
      axis.title = element_text(face = "plain"),
      panel.border = element_rect(colour = "#333333", linewidth = .5),
      plot.background = element_rect(fill = "white", colour = NA)
    )
}

save_png <- function(plot, filename, width, height) {
  ggsave(file.path(here, filename), plot, width = width, height = height,
         units = "in", dpi = 300, bg = "white")
}

global_axes <- function(plot, x_label, y_label, left=.065, bottom=.045) {
 ggdraw() +
  draw_plot(plot,x=left,y=bottom,width=1-left,height=1-bottom) +
  draw_label(x_label,x=(1+left)/2,y=.030,fontface='plain',size=14) +
  draw_label(y_label,x=.043,y=(1+bottom)/2,angle=90,fontface='plain',size=14)
}

species_label <- c(
  "ALT.22.08" = "bolditalic(Panagrolaimus)~bold(\"sp.\")~plain(\"(WWL002)\")",
  "ARO.22.05" = "bolditalic(Acrobeloides~goughiensis)~plain(\"(WWL726)\")",
  "PAP.22.17" = "bolditalic(Acrobeloides~tricornis)~plain(\"(WWL709)\")",
  "LAG.22.03" = "bolditalic(Plectus~pavus)~plain(\"(WWL505)\")",
  "ES5" = "bolditalic(Panagrolaimus~einhardi)~plain(\"(ES5)\")",
  "A. max" = "bolditalic(Acrobeloides~maximus)~plain(\"(A. max)\")"
)

# Render Figure 2 before Figure 1: plotmath glyph anti-aliasing depends on
# earlier device/font-cache use in this R environment. This order reproduces
# both manuscript images exactly, without reading reference images.
## Figure 2: dead and cryptobiotic counts on day 4 -------------------------
species_label_fig2 <- c(
  "ALT.22.08" = "bolditalic(Panagrolaimus)~bold(\"sp.\")~plain(\"(WWL002)\")",
  "PAP.22.17" = "bolditalic(Acrobeloides~tricornis)~plain(\"(WWL709)\")"
)

d2 <- read_csv(file.path(cache, "R_day4_dead_cryptobiotic_pooled_counts.csv"), show_col_types = FALSE) %>%
  rename(salinity = `Salinity [PSU]`, temperature = Temperature_treatment,
         dead = confirmed_dead_n, cryptobiotic = inactive_or_unconfirmed_n, total = Total)

d2_replicates <- read_csv(file.path(cache, "R_day4_dead_cryptobiotic_analysis_ready.csv"), show_col_types = FALSE) %>%
  rename(salinity = `Salinity [PSU]`, temperature_C = Temperature_treatment,
         dead = confirmed_dead_n, cryptobiotic = inactive_or_unconfirmed_n, total = Total) %>%
  mutate(temperature = relevel(factor(temperature_C), ref = "16"),
         salinity_10 = salinity / 10)

curve_for <- function(strain, temperature_value, response, label) {
  z <- filter(d2_replicates, Strain == strain) %>%
    mutate(success = .data[[response]], failure = total - success)
  fit <- glm(cbind(success, failure) ~ temperature * salinity_10,
             data = z, family = quasibinomial(link = "logit"))
  grid <- data.frame(
    salinity = seq(min(z$salinity[z$temperature_C == temperature_value]),
                   max(z$salinity[z$temperature_C == temperature_value]), length.out = 250),
    temperature = factor(as.character(temperature_value), levels = levels(z$temperature))
  ) %>% mutate(salinity_10 = salinity / 10)
  pred <- predict(fit, newdata = grid, type = "link", se.fit = TRUE)
  grid$percent <- 100 * plogis(pred$fit)
  grid$lower <- 100 * plogis(pred$fit - qnorm(.975) * pred$se.fit)
  grid$upper <- 100 * plogis(pred$fit + qnorm(.975) * pred$se.fit)
  grid$category <- label
  grid
}

panel_dead_crypto <- function(strain, temperature) {
  z <- filter(d2, Strain == strain, temperature == !!temperature)
  observed <- bind_rows(
    transmute(z, salinity, percent = 100 * cryptobiotic / total, category = "Cryptobiotic"),
    transmute(z, salinity, percent = 100 * dead / total, category = "Dead")
  )
  curves <- bind_rows(curve_for(strain, temperature, "cryptobiotic", "Cryptobiotic"),
                      curve_for(strain, temperature, "dead", "Dead"))
  ggplot() +
    geom_ribbon(data = curves,
                aes(salinity, ymin = lower, ymax = upper, fill = category),
                alpha = .14, colour = NA) +
    geom_line(data = curves, aes(salinity, percent, colour = category), linewidth = .75) +
    geom_point(data = filter(observed, category == "Cryptobiotic"),
               aes(salinity, percent, colour = category), size = 2) +
    geom_point(data = filter(observed, category == "Dead"),
               aes(salinity, percent, colour = category), shape = 21,
               fill = "white", size = 2, stroke = .8) +
    scale_colour_manual(values = c(Cryptobiotic = blue, Dead = red)) +
    scale_fill_manual(values = c(Cryptobiotic = blue, Dead = red), guide = "none") +
    facet_wrap(~factor(species_label_fig2[[strain]]), labeller = label_parsed) +
    coord_cartesian(ylim = c(-3, 103)) +
    labs(x = NULL, y = NULL, colour = NULL,
         title=if(strain=="ALT.22.08") paste0(temperature," °C") else NULL) +
    theme_manuscript(12) +
    theme(
      strip.text = element_text(size = 12, face = "plain"),
      plot.title=element_text(hjust=.5,size=16,face="bold"),
      plot.title.position="panel",
      axis.text = element_text(size = 11),
      legend.text = element_text(size = 11),
      legend.key.width = unit(1.4, "lines")
    )
}

top2 <- panel_dead_crypto("ALT.22.08", 22) + panel_dead_crypto("ALT.22.08", 16)
bottom2 <- panel_dead_crypto("PAP.22.17", 22) + panel_dead_crypto("PAP.22.17", 16)
header22_fig2 <- ggplot() +
  annotate("text", 0, 0, label = "22 °C", fontface = "bold", size = 5.5) + theme_void()
header16_fig2 <- ggplot() +
  annotate("text", 0, 0, label = "16 °C", fontface = "bold", size = 5.5) + theme_void()
fig2_core <- (top2 / bottom2) +
  plot_layout(heights = c(1, 1), guides = "collect") &
  theme(legend.position = "bottom")
fig2_legend <- get_legend(fig2_core)
fig2_core <- fig2_core & theme(legend.position = "none")
fig2 <- ggdraw() +
  draw_plot(fig2_core, x = .08, y = .14, width = .92, height = .86) +
  draw_label("Salinity (PSU)", x = .54, y = .128, size = 14) +
  draw_label("Counts (%)", x = .065, y = .57, angle = 90, size = 14) +
  draw_plot(fig2_legend, x = .30, y = .055, width = .48, height = .055)
save_png(fig2, "Figure_2_mortality_cryptobiosis.png", 9, 7.5)


## Figure 1: global dose-response curves and EC50 intervals -----------------
raw1 <- read_csv(file.path(cache, "R_blockwise_dead_cryptobiotic_analysis_ready.csv"), show_col_types = FALSE) %>%
  mutate(strain_temperature = paste0(strain, "__", format(temperature_C, trim = TRUE), "C"))
points1 <- raw1 %>%
  group_by(strain_temperature, salinity_PSU) %>%
  summarise(dead_inactive_n = sum(dead_inactive_n), total_n = sum(total_n), .groups = "drop") %>%
  mutate(percent = 100 * dead_inactive_n / total_n)
curves1 <- read_csv(file.path(dose, "global_model_marginal_dose_response_curves.csv"), show_col_types = FALSE)
ec1 <- read_csv(file.path(stats, "cluster_bootstrap_marginal_EC50_nonnegative_summary.csv"), show_col_types = FALSE)

panel_ec50 <- function(group, strain) {
  obs <- filter(points1, strain_temperature == group)
  pred <- filter(curves1, strain_temperature == group)
  est <- filter(ec1, strain_temperature == group)
  strip <- species_label[[strain]]
  ggplot() +
    geom_hline(yintercept = 50, linetype = "dashed", colour = "grey55", linewidth = .35) +
    geom_rect(data = est,
              aes(xmin = pmax(min(obs$salinity_PSU), nonnegative_bootstrap_CI_low),
                  xmax = pmin(max(obs$salinity_PSU), nonnegative_bootstrap_CI_high),
                  ymin = -Inf, ymax = Inf), fill = red, alpha = .10) +
    geom_ribbon(data = pred,
                aes(salinity_PSU, ymin = 100 * CI_low, ymax = 100 * CI_high),
                fill = blue, alpha = .15) +
    geom_line(data = pred, aes(salinity_PSU, 100 * predicted_dead_inactive),
              colour = blue, linewidth = .65) +
    geom_vline(data = est, aes(xintercept = point_EC50_PSU),
               colour = red, linetype = "dotted", linewidth = .55) +
    geom_point(data = obs, aes(salinity_PSU, percent), size = 1.6) +
    facet_wrap(~factor(strip), labeller = label_parsed) +
    scale_x_continuous(breaks = seq(0, 50, 10)) +
    coord_cartesian(xlim = c(0, 50), ylim = c(-3, 103)) +
    labs(x = NULL, y = NULL) + theme_manuscript(12) +
    theme(axis.text=element_text(size=11,colour='black'),
          strip.text=element_text(size=12,face='plain'),
          plot.margin=margin(3,4,3,4))
}

blank_panel <- plot_spacer()
layout1 <- list(
  list("ALT.22.08__22C", "ALT.22.08"), list("ALT.22.08__16C", "ALT.22.08"),
  list("ARO.22.05__22C", "ARO.22.05"), NULL,
  list("PAP.22.17__22C", "PAP.22.17"), list("PAP.22.17__16C", "PAP.22.17"),
  list("LAG.22.03__22C", "LAG.22.03"), NULL,
  list("ES5__22C", "ES5"), list("ES5__16C", "ES5"),
  list("A. max__22C", "A. max"), NULL
)
panels1 <- lapply(layout1, function(x) if (is.null(x)) blank_panel else panel_ec50(x[[1]], x[[2]]))
panels1[[1]] <- panels1[[1]] + labs(title="22 °C") + theme(plot.title=element_text(hjust=.5,size=16,face="bold"),plot.title.position="panel")
panels1[[2]] <- panels1[[2]] + labs(title="16 °C") + theme(plot.title=element_text(hjust=.5,size=16,face="bold"),plot.title.position="panel")
header22 <- ggplot() + annotate("text", 0, 0, label = "22 °C", fontface = "bold", size = 5.5) + theme_void()
# ES5 R1 was confirmed to have been run at 16 C (2026-10-16).
# Retain legacy model keys above to reuse the unchanged fitted predictions.
header16 <- ggplot() + annotate("text", 0, 0, label = "16 °C", fontface = "bold", size = 5.5) + theme_void()
fig1 <- wrap_plots(panels1, ncol = 2, guides = "collect") +
  plot_annotation(theme = theme(plot.margin = margin(4, 4, 4, 4))) &
  theme(legend.position = "none")
fig1 <- global_axes(fig1, "Salinity (PSU)", "Dead or cryptobiotic (%)", .065, .045)
save_png(fig1, "Figure_1_dose_response.png", 9, 11)


  invisible(figures)
}

.salinity_cli_args <- commandArgs(trailingOnly = TRUE)
.salinity_script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(.salinity_script_arg) &&
    basename(sub("^--file=", "", .salinity_script_arg[1])) == "02_salinity.R" &&
    length(.salinity_cli_args) && identical(.salinity_cli_args[1], "--render-only")) {
  if (length(.salinity_cli_args) != 2L) stop("Usage: 02_salinity.R --render-only PROJECT_ROOT")
  render_salinity_figures(.salinity_cli_args[2])
}
