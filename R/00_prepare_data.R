# Rebuild the analysis tables from the experimental source files.
# Outputs are the five CSV tables required by the analysis modules.
run_prepare_data <- function(project_root) {
  suppressPackageStartupMessages({
    library(readxl)
    library(dplyr)
    library(tidyr)
  })
  root <- normalizePath(project_root, winslash = '/', mustWork = TRUE)
  output <- file.path(root, 'Results', 'inputs')
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  numeric_value <- function(x) suppressWarnings(as.numeric(as.character(x)))
  column <- function(x, name) {
    if (name %in% names(x)) x[[name]] else rep(NA, nrow(x))
  }
  number_in_label <- function(x) numeric_value(sub('^[^0-9-]*(-?[0-9.]+).*', '\\1', as.character(x)))
  read_source_csv <- function(path) {
    read.csv(path, sep = ';', colClasses = 'character', check.names = FALSE,
             na.strings = c('', 'NA'), fileEncoding = 'UTF-8-BOM')
  }
  save_table <- function(x, name) {
    write.csv(x, file.path(output, name), row.names = FALSE, na = '')
    x
  }
  source_path <- function(path) substring(path, nchar(root) + 2L)

  # Standardize only the salinity columns needed by the analysis and audit.
  salinity_table <- function(x, experiment, strain, species, reproduction,
                             origin, temperature, source, sheet, source_group) {
    data.frame(
      Data_source = source_group, Experiment_id = experiment, Strain = strain,
      Species = species, reproductive_mode = reproduction,
      'Sampling location' = origin, Temperature_treatment = temperature,
      'Salinity [PSU]' = number_in_label(column(x, 'Treatment')),
      Observation_day_raw = number_in_label(column(x, 'Day')),
      Replicate = numeric_value(column(x, 'Replicate')),
      Alive = numeric_value(column(x, 'Alive')),
      Dead = numeric_value(column(x, 'Dead')),
      Total = numeric_value(column(x, 'Alive')) + numeric_value(column(x, 'Dead')),
      Dead_without_dormant = numeric_value(column(x, 'Dead_without_dormant')),
      Dormant_stage = numeric_value(column(x, 'Dormant_stage')),
      Paralysed_adult = numeric_value(column(x, 'Paralaysed_adult')),
      Source_file = source_path(source), Source_sheet = sheet,
      Source_date = as.character(column(x, 'Date')),
      Source_reference = seq_len(nrow(x)) + 1L,
      Mortality_definition = 'Dead uses the corrected source summary column.',
      Mortality_analysis_eligible = TRUE, QC_flag = '', Count_consistency_QC = '',
      notes = '', Analysis_family = 'osmobiosis_salinity', check.names = FALSE
    )
  }

  # Adjusted observation sheets also retain the counts before M9 recovery.
  source_dir <- file.path(root, 'Data/raw/salinity')
  specs <- list(
    list('A_maximus/A_maximus_salinity_22C.xlsx', 'Amax_22°C_adujsted_final',
         'A.max_22C', 'A. max', 'Acrobeloides maximus', 'Senegal', 22, 20),
    list('WWL709/WWL709_salinity_16C.xlsx', 'PAP22.17_16°C_adjusted',
         'PAP.22.17_16C', 'PAP.22.17', 'Acrobeloides tricornis', 'Chile, Paposo', 16, 20),
    list('WWL002/WWL002_salinity_16_22C.xlsx', 'ALT22.08_16°C_adjusted',
         'ALT.22.08_16C', 'ALT.22.08', 'Panagrolaimus sp.', 'Chile, Altiplano', 16, 30),
    list('WWL002/WWL002_salinity_16_22C.xlsx', 'ALT22.08_22°C_adjusted ',
         'ALT.22.08_22C', 'ALT.22.08', 'Panagrolaimus sp.', 'Chile, Altiplano', 22, 30)
  )
  recovery <- list()
  salinity <- list()
  for (spec in specs) {
    path <- file.path(source_dir, spec[[1]])
    x <- as.data.frame(read_excel(path, sheet = spec[[2]]))
    x$Source_row <- seq_len(nrow(x)) + 1L
    x$Experiment_id <- spec[[3]]
    x$Strain_standardized <- spec[[4]]
    x$Temperature_treatment_C <- spec[[7]]
    x$Source_file <- source_path(path)
    x$Source_sheet <- spec[[2]]
    recovery[[spec[[3]]]] <- x
    z <- x[!is.na(numeric_value(x$Day)) & !is.na(number_in_label(x$Treatment)) &
             !is.na(numeric_value(x$Reference)), , drop = FALSE]
    z <- z[order(numeric_value(z$Day), number_in_label(z$Treatment), numeric_value(z$Reference)), ]
    dose_days <- paste(numeric_value(z$Day), number_in_label(z$Treatment))
    if (any(table(dose_days) != 3L)) stop('Expected three source plates per dose/day in ', spec[[3]])
    corrected_replicate <- ave(seq_len(nrow(z)), dose_days, FUN = seq_along)
    y <- salinity_table(z, spec[[3]], spec[[4]], spec[[5]], 'parthenogenetic',
                        spec[[6]], spec[[7]], path, spec[[2]], 'salinity_recovery')
    y$Replicate <- as.numeric(corrected_replicate)
    y$Dead <- numeric_value(z$Dead_with_paralysed)
    y$Total <- spec[[8]]
    y$Alive <- y$Total - y$Dead
    y$Dormant_stage <- numeric_value(column(z, 'Doramnt_stage'))
    y$Source_reference <- numeric_value(z$Reference)
    y$Mortality_analysis_eligible <- y$Observation_day_raw != 14 &
      !is.na(y$Dead) & y$Dead >= 0 & y$Dead <= y$Total
    y$Mortality_definition <- 'Dead uses source column Dead_with_paralysed; paralysed individuals are counted as dead'
    salinity[[spec[[3]]]] <- y
  }

  # ES5 at 22 C uses its corrected summary and three replicates per dose/day.
  path <- file.path(root, 'Data/raw/salinity/ES5/ES5_salinity_22C.csv')
  x <- read_source_csv(path)
  y <- salinity_table(x, 'ES5_22C', 'ES5', 'Panagrolaimus einhardi', 'sexual',
                      'Bornheim, Bonn, Germany', 22, path, 'CSV', 'salinity_recovery')
  y <- y[order(y$Observation_day_raw, y[['Salinity [PSU]']], y$Replicate), ]
  y$Total <- 20
  y$Alive <- 20 - y$Dead
  y$Mortality_analysis_eligible <- y$Observation_day_raw != 14 & !is.na(y$Dead) & y$Dead >= 0 & y$Dead <= 20
  y$Mortality_definition <- 'Dead uses source column Dead; paralysed individuals are counted as dead'
  salinity[['ES5_22C']] <- y

  summaries <- list(
    list('WWL709/WWL709_salinity_22C_R1.csv', 'PAP.22.17_salinity_R1',
         'PAP.22.17', 'Acrobeloides tricornis', 'parthenogenetic', 'Chile, Paposo', 22),
    list('WWL709/WWL709_salinity_22C_R2.csv', 'PAP.22.17_salinity_R2',
         'PAP.22.17', 'Acrobeloides tricornis', 'parthenogenetic', 'Chile, Paposo', 22),
    list('WWL726/WWL726_salinity_22C.csv', 'ARO.22.05_salinity_R1',
         'ARO.22.05', 'Acrobeloides goughiensis', 'asexual', 'Atacama Desert, Chile', 22),
    list('ES5/ES5_salinity_16C.csv', 'ES5_salinity_R1',
         'ES5', 'Panagrolaimus einhardi', 'sexual', 'Bornheim, Bonn, Germany', 16)
  )
  for (spec in summaries) {
    path <- file.path(root, 'Data/raw/salinity', spec[[1]])
    x <- read_source_csv(path)
    y <- salinity_table(x, spec[[2]], spec[[3]], spec[[4]], spec[[5]],
                        spec[[6]], spec[[7]], path, 'CSV', 'supporting')
    if ('Strain' %in% names(x)) {
      preliminary <- x$Strain == 'ES5'
      y$Experiment_id[preliminary] <- 'ES5_salinity_preliminary_R1'
      y$Strain[preliminary] <- 'ES5'
      y$Species[preliminary] <- 'Panagrolaimus einhardi'
      y$reproductive_mode[preliminary] <- 'sexual'
      y[['Sampling location']][preliminary] <- 'Bornheim, Bonn, Germany'
    }
    if ('Dormant_stage?' %in% names(x)) y$Dormant_stage <- numeric_value(x[['Dormant_stage?']])
    y$Mortality_analysis_eligible <- !is.na(y$Total) & y$Total == 20 &
      y$Dead >= 0 & y$Alive >= 0
    if (spec[[2]] == 'ES5_salinity_R1') {
      # Day 11 is a later rehydration check with no original-cohort counts.
      y$Mortality_analysis_eligible[y$Observation_day_raw == 11] <- FALSE
      y$Mortality_definition <- 'Main source summary: immobile or missing animals were counted as dead; the later ad hoc water-rehydration check is exploratory and stored separately.'
    }
    salinity[[spec[[2]]]] <- y
  }

  # The 2023 source counts juveniles separately. Reconstruct the 20-adult
  # cohort only where Adult + Dead + Dormant supports that denominator.
  path <- file.path(root, 'Data/raw/salinity/series2023/WWL002_WWL505_Amax_salinity_22C.xlsx')
  x <- as.data.frame(read_excel(path, sheet = 'ALT_LAG'))
  metadata <- list(
    ALT = c('ALT.22.08', 'Panagrolaimus sp.', 'parthenogenetic', 'Altiplano, Chile'),
    AMAX = c('A. max', 'Acrobeloides maximus', 'parthenogenetic', 'Senegal'),
    LAG = c('LAG.22.03', 'Plectus pavus', 'asexual', 'Atacama Desert, Chile')
  )
  for (label in names(metadata)) {
    z <- x[trimws(x$Species) == label, , drop = FALSE]
    spec <- metadata[[label]]
    id <- paste0('series2023_', gsub(' ', '', spec[[1]], fixed = TRUE), '_salinity_22C')
    y <- salinity_table(z, id, spec[[1]], spec[[2]], spec[[3]], spec[[4]],
                        22, path, 'ALT_LAG', 'series2023')
    alive <- numeric_value(z$Alive)
    dormant <- numeric_value(z$Dormant)
    adult <- numeric_value(z$Adult)
    original_total <- alive + y$Dead + dormant
    reconstruct <- original_total != 20 & adult + y$Dead + dormant == 20
    reconstruct[is.na(reconstruct)] <- FALSE
    y$Alive <- ifelse(reconstruct, adult + dormant, alive + dormant)
    y$Total <- ifelse(reconstruct, 20, original_total)
    # One source plate has irreconcilable cohort counts. Preserve the source
    # workbook and mark its analysis counts missing, as in the reviewed input.
    invalid <- is.na(y$Total) | y$Total != 20
    y$Alive[invalid] <- NA_real_
    y$Dead[invalid] <- NA_real_
    y$Total[invalid] <- NA_real_
    y$Dormant_stage <- dormant
    y$Dead_without_dormant <- numeric_value(z$Dead)
    y$Source_reference <- match(rownames(z), rownames(x)) + 1L
    y$Mortality_analysis_eligible <- !is.na(y$Total) & y$Total == 20 &
      y$Dead >= 0 & y$Alive >= 0 & dormant >= 0
    y$Mortality_definition <- 'Original 20-adult cohort: Dead/20; dormant animals are treated as not dead descriptively but were counted without standardized rehydration. Where offspring increased Alive, Adult + Dead + Dormant reconstructs the original cohort and juveniles are excluded.'
    salinity[[id]] <- y
  }
  core <- bind_rows(salinity)
  core$Count_consistency_QC <- ifelse(core$Alive + core$Dead == core$Total, 'OK', 'mismatch')
  save_table(core, 'salinity_observations.csv')
  # Recovery needs the source counts, labels and original replicate numbers.
  recovery <- bind_rows(lapply(recovery, function(x) {
    data.frame(
      Experiment_id = x$Experiment_id, Strain_standardized = x$Strain_standardized,
      Temperature_treatment_C = x$Temperature_treatment_C,
      Day = numeric_value(x$Day), Treatment = as.character(x$Treatment),
      Replicate = numeric_value(x$Replicate), Rehyadrated = numeric_value(x$Rehyadrated),
      Doramnt_stage = numeric_value(column(x, 'Doramnt_stage')),
      Paralaysed_adult = numeric_value(column(x, 'Paralaysed_adult')),
      Paralaysed_juvenil = numeric_value(column(x, 'Paralaysed_juvenil')),
      Source_file = x$Source_file, Source_sheet = x$Source_sheet,
      Source_row = x$Source_row, check.names = FALSE
    )
  }))
  save_table(recovery, 'recovery_observations.csv')

  # Each DailyLog contains all recorded observation days for its experiment.
  cultivation <- list()
  for (ce in 1:3) {
    directory <- file.path(root, 'Data/raw/cultivation', paste0('CE', ce))
    path <- list.files(directory, pattern = '_daily_log[.]csv$', full.names = TRUE)
    stopifnot(length(path) == 1L)
    wide <- read_source_csv(path)
    long <- pivot_longer(wide, matches('^t[0-9]+$'), names_to = 'Observation_day_raw', values_to = 'Status_raw', cols_vary = 'slowest')
    long <- long[!is.na(long$Status_raw) & nzchar(trimws(long$Status_raw)), ]
    long$Status_raw <- trimws(long$Status_raw)
    long$Experiment_id <- paste0('MAL', ce)
    long$Strain <- if (ce == 3L) 'ALT.22.08' else 'PAP.22.17'
    long$Salinity_PSU <- numeric_value(long[['Treatment(PSU)']])
    long$Replicate <- numeric_value(long$Replicate)
    long$Day <- as.numeric(sub('^t', '', long$Observation_day_raw))
    long$Status_code <- numeric_value(sub('.*?([0-7]).*', '\\1', long$Status_raw))
    long$Source_file <- source_path(path)
    cultivation[[ce]] <- long[, c('Experiment_id', 'Strain', 'Salinity_PSU', 'Replicate',
                                  'Generation', 'Observation_day_raw', 'Status_raw', 'Day',
                                  'Status_code', 'Source_file')]
  }
  save_table(bind_rows(cultivation), 'cultivation_observations.csv')

  # D2: read six complete strains from the original survival worksheet.
  path <- file.path(root, 'Data/raw/desiccation/D2/D2_survival.xlsx')
  cells <- read_excel(path, sheet = 'survival_numbers', col_names = FALSE,
                      col_types = 'text', .name_repair = 'minimal')
  strain_metadata <- list(
    ALT2208 = c('Panagrolaimus sp.', 'asexual', 'Altiplano, Chile'),
    'ES5 (P. einhardi)' = c('Panagrolaimus einhardi', 'sexual', 'Bornheim, Bonn, Germany'),
    JU1365 = c('Panagrolaimus sp.', 'sexual', 'Periyar-Madurai, Tamil Nadu, India'),
    PAP2229 = c('Panagrolaimus nebliphilus', 'sexual', 'Paposo, Chile'),
    PS1159 = c('Panagrolaimus sp.', 'sexual', 'San Francisco, USA'),
    QG3759 = c('Panagrolaimus sp.', 'sexual', 'Ukraine')
  )
  d2 <- list()
  treatment <- NA_character_
  for (i in seq_len(nrow(cells))) {
    if (length(d2) == 36L) break
    label <- as.character(cells[[1]][i])
    if (is.na(label)) next
    if (startsWith(label, 'Acclimated survival')) treatment <- 'Acclimated'
    if (startsWith(label, 'Non-acclimated survival')) treatment <- 'Non_acclimated'
    if (!label %in% names(strain_metadata)) next
    spec <- strain_metadata[[label]]
    for (replicate in 1:3) {
      j <- 2L + (replicate - 1L) * 3L
      alive <- numeric_value(cells[[j]][i])
      total <- numeric_value(cells[[j + 1L]][i])
      if (!is.finite(alive) || !is.finite(total) || total <= 0 || alive < 0 || alive > total) {
        stop('Invalid D2 count at source row ', i, ', replicate ', replicate)
      }
      d2[[length(d2) + 1L]] <- data.frame(
        Strain = label, Species = spec[[1]], reproduction = spec[[2]],
        'Sampling location' = spec[[3]], Treatment = treatment,
        Replicate = paste0('R', replicate), Alive = alive, Total = total,
        Dead = total - alive, Survival_percent = round(100 * alive / total, 2),
        check.names = FALSE
      )
    }
  }
  stopifnot(length(d2) == 36L)
  d2 <- bind_rows(d2)
  stopifnot(!anyDuplicated(d2[c('Strain', 'Treatment', 'Replicate')]),
            all(table(d2$Treatment) == 18L))
  save_table(d2, 'D2_survival.csv')
  path <- file.path(root, 'Data/raw/desiccation/D2/D2_WWL115_technical_counts.csv')
  technical <- read.csv(path, check.names = FALSE)
  # TR1--TR3 are repeated counts of the same 100-animal biological plate.
  stopifnot(nrow(technical) == 6L, all(technical[c('TR1', 'TR2', 'TR3')] >= 0),
            all(technical[c('TR1', 'TR2', 'TR3')] <= 100))
  technical$Alive_sum <- rowSums(technical[c('TR1', 'TR2', 'TR3')])
  technical$Total_counted <- 300
  technical$Dead_sum <- 300 - technical$Alive_sum
  # The original technical-count table records these means to ten decimals.
  technical$Alive_mean_per_100 <- round(technical$Alive_sum / 3, 10)
  technical$Survival_percent <- technical$Alive_mean_per_100
  save_table(technical, 'D2_technical_counts.csv')
  message('Prepared salinity, recovery, cultivation and D2 tables from source files.')
  invisible(output)
}
