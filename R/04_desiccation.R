# Desiccation analysis

run_desiccation <- function(project_root) {
  suppressPackageStartupMessages({
    library(readxl)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
    library(ggtext)
    library(glmmTMB)
    library(emmeans)
  })
  options(stringsAsFactors = FALSE, scipen = 999)
  set.seed(1211)
  raw_folder <- file.path(project_root, "Data", "raw", "desiccation")
  output_folder <- file.path(project_root, "Results", "statistics", "desiccation")
  tables <- file.path(project_root, "Results", "tables", "desiccation")
  cache <- file.path(project_root, "Results", "cache", "desiccation")
  figures <- file.path(project_root, "Manuscript", "figures")
  for (folder in c(output_folder, tables, cache, figures)) dir.create(folder, recursive = TRUE,
    showWarnings = FALSE)

  z95 <- qnorm(0.975)
  wilson <- function(x, n) {
    p <- x/n
    den <- 1 + z95^2/n
    mid <- (p + z95^2/(2 * n))/den
    half <- z95 * sqrt(p * (1 - p)/n + z95^2/(4 * n^2))/den
    data.frame(low = pmax(0, mid - half), high = pmin(1, mid + half))
  }

  summarise_recovery <- function(x) {
    x$experiment <- factor(x$experiment, levels = c("D2", "D1"))
    ans <- x |>
      group_by(experiment, strain, group) |>
      summarise(n_replicates = n(), recovered = sum(recovered), initial = sum(initial), pooled_recovery = recovered/initial,
        replicate_mean = mean(recovery), replicate_sd = sd(recovery), .groups = "drop")
    ans <- bind_cols(ans, wilson(ans$recovered, ans$initial))
    set.seed(1211)
    historical_rng_gap_consumed <- FALSE
    boot <- x |>
      group_by(experiment, strain, group) |>
      group_modify(~{
        if (as.character(.y$experiment) == "D1" && !historical_rng_gap_consumed) {
          for (j in seq_len(6L)) replicate(20000L, sample.int(5L, 5L, replace = TRUE))
          historical_rng_gap_consumed <<- TRUE
        }
        y <- .x$recovery
        b <- replicate(20000, mean(sample(y, length(y), replace = TRUE)))
        tibble(bootstrap_low = unname(quantile(b, 0.025)), bootstrap_high = unname(quantile(b,
          0.975)))
      }) |>
      ungroup()
    left_join(ans, boot, by = c("experiment", "strain", "group"))
  }

  d1_source <- read_excel(file.path(raw_folder, "D1", "D1_survival.xlsx"), sheet = "Nematoden Daten")
  d1_data <- transmute(d1_source, experiment = "D1", strain = gsub(" ", ".", trimws(Strain), fixed = TRUE),
    species = trimws(Species), group = if_else(Treatment == "Desiccation", "Desiccation", "Control"),
    replicate = as.integer(Replicate), recovered = as.integer(Alive), initial = as.integer(Total),
    failure = as.integer(Dead), missing = 0L, recovery = recovered/initial, note = "No preconditioning; 8 h desiccation, 17 h rehydration")
  d1_data$strain <- recode(d1_data$strain, GBG22.06 = "GBG.22.06")

  d2_source <- read.csv(file.path(project_root, "Results/inputs/D2_survival.csv"), check.names = FALSE)
  shuimeiren_technical <- read.csv(file.path(project_root, "Results/inputs/D2_technical_counts.csv"),
    check.names = FALSE) |>
    transmute(Strain, Species, reproduction, `Sampling location` = Sampling_location, Treatment,
      Replicate, Alive = Alive_mean_per_100, Total = 100, Dead = 100 - Alive_mean_per_100, Survival_percent)
  d2_source <- bind_rows(d2_source, shuimeiren_technical)
  d2_data <- transmute(d2_source, experiment = "D2", strain = recode(Strain, ALT2208 = "ALT.22.08",
    PAP2229 = "PAP.22.29", `ES5 (P. einhardi)` = "ES5"), species = Species, group = if_else(Treatment ==
    "Acclimated", "Preconditioned", "Not preconditioned"), replicate = as.integer(sub("R", "", Replicate)),
    recovered = as.numeric(Alive), initial = as.numeric(Total), failure = as.numeric(Dead), missing = 0,
    recovery = recovered/initial, note = "48 h desiccation and 24 h rehydration; groups differ only in 96 h preconditioning")

  all_data <- bind_rows(d1_data, d2_data)
  stopifnot(all(all_data$recovered >= 0), all(all_data$initial > 0), all(all_data$recovered + all_data$failure +
    all_data$missing == all_data$initial))

  all_data <- all_data |>
    mutate(species = case_when(strain == "PAP.22.29" ~ "Panagrolaimus nebliphilus", TRUE ~ species),
      origin = case_when(strain %in% c("ALT.22.08", "ARO.22.05", "LAG.22.03", "PAP.22.17", "PAP.22.29",
        "TDT.23.02", "TDT.23.16") ~ "Atacama", strain %in% c("860.22.04", "GBG.22.06", "PWB.22.01") ~
        "Namib", strain == "ES5" ~ "Germany", strain == "DF5050" ~ "Iceland", strain == "JU1365" ~
        "India", strain == "PS1159" ~ "USA", strain == "QG3759" ~ "Ukraine", TRUE ~ NA_character_),
      atacama_origin = if_else(origin == "Atacama", "Atacama", "Non-Atacama"), desiccated = if_else(experiment ==
        "D1" & group == "Control", 0L, 1L), preconditioned = if_else(experiment == "D2" & group ==
        "Preconditioned", 1L, 0L), protocol_group = case_when(experiment == "D1" & group == "Control" ~
        "D1_control", experiment == "D1" ~ "D1_desiccation_no_preconditioning", experiment ==
        "D2" & group == "Preconditioned" ~ "D2_desiccation_96h_preconditioning", experiment ==
        "D2" ~ "D2_desiccation_no_preconditioning"), biological_unit = "plate/filter", included_joint_model = TRUE,
      survived = recovered, not_survived = initial - recovered, survival_proportion = recovery)

  descriptive_summary <- summarise_recovery(all_data)
  write.csv(all_data, file.path(output_folder, "analysis_data_all_groups.csv"), row.names = FALSE)
  write.csv(descriptive_summary, file.path(tables, "survival_descriptive_bootstrap.csv"), row.names = FALSE)

  model_data_all <- all_data |>
    mutate(survived = round(survived), recovered = round(recovered), failure = initial - recovered,
      recovery = recovered/initial, survival_proportion = recovery, not_survived = initial - survived,
      non_survived = initial - survived, experiment = relevel(factor(experiment), "D2"), strain = factor(strain),
      desiccated = relevel(factor(desiccated), "0"), preconditioned = relevel(factor(preconditioned),
        "0"), atacama_origin = relevel(factor(atacama_origin), "Non-Atacama"))

  model_data <- model_data_all |>
    filter(desiccated == "1") |>
    droplevels()

  stopifnot(nrow(model_data) == 82L, nlevels(model_data$strain) == 12L, all(model_data$survived >=
    0), all(model_data$non_survived >= 0), all(model_data$survived + model_data$non_survived == model_data$initial))

  form_full <- cbind(survived, non_survived) ~ experiment + preconditioned + atacama_origin + (1 |
    strain)

  joint_model <- glmmTMB(form_full, data = model_data, family = betabinomial(link = "logit"))

  model_without_origin <- update(joint_model, . ~ . - atacama_origin)
  model_without_preconditioning <- update(joint_model, . ~ . - preconditioned)
  model_without_protocol <- update(joint_model, . ~ . - experiment)
  binomial_model <- glmmTMB(form_full, data = model_data, family = binomial(link = "logit"))

  lrt_row <- function(reduced, full, term) {
    a <- anova(reduced, full)
    data.frame(term = term, likelihood_ratio_chisq = a$Chisq[2], df = a$`Chi Df`[2], p_value = a$`Pr(>Chisq)`[2])
  }

  term_tests <- bind_rows(lrt_row(model_without_origin, joint_model, "Atacama origin"), lrt_row(model_without_preconditioning,
    joint_model, "Preconditioning"), lrt_row(model_without_protocol, joint_model, "Experimental protocol"))

  model_validation <- data.frame(model = c("Primary beta-binomial GLMM", "Ordinary binomial GLMM"),
    n = c(nobs(joint_model), nobs(binomial_model)), logLik = c(as.numeric(logLik(joint_model)), as.numeric(logLik(binomial_model))),
    df = c(attr(logLik(joint_model), "df"), attr(logLik(binomial_model), "df")), AIC = c(AIC(joint_model),
      AIC(binomial_model)), convergence = c(joint_model$fit$convergence, binomial_model$fit$convergence),
    positive_definite_hessian = c(joint_model$sdr$pdHess, binomial_model$sdr$pdHess)) |>
    mutate(delta_AIC = AIC - min(AIC))

  coefficient_table <- as.data.frame(summary(joint_model)$coefficients$cond)
  coefficient_table$term <- rownames(coefficient_table)
  rownames(coefficient_table) <- NULL
  coefficient_table <- coefficient_table |>
    dplyr::select(term, everything())

  strain_variance <- VarCorr(joint_model)$cond$strain
  variance_table <- data.frame(component = "Strain random intercept", variance = as.numeric(strain_variance[1]),
    standard_deviation = attr(strain_variance, "stddev")[[1]])

  joint_predictions <- emmeans(joint_model, ~preconditioned | experiment, type = "response")
  estimated_survival <- as.data.frame(joint_predictions) |>
    filter(experiment == "D2")

  write.csv(model_data, file.path(output_folder, "analysis_data_desiccated.csv"), row.names = FALSE)
  write.csv(term_tests, file.path(output_folder, "fixed_effect_likelihood_ratio_tests.csv"), row.names = FALSE)
  write.csv(model_validation, file.path(output_folder, "distribution_validation.csv"), row.names = FALSE)
  write.csv(coefficient_table, file.path(output_folder, "fixed_effect_coefficients.csv"), row.names = FALSE)
  write.csv(variance_table, file.path(output_folder, "strain_random_effect.csv"), row.names = FALSE)
  write.csv(estimated_survival, file.path(output_folder, "estimated_survival_by_treatment.csv"), row.names = FALSE)
  saveRDS(joint_model, file.path(cache, "joint_model.rds"))

  for (v in c("experiment", "preconditioned", "strain", "atacama_origin")) model_data[[v]] <- factor(model_data[[v]])
  cat("DESIGN TABLE\n")
  print(with(model_data, table(experiment, preconditioned)))
  cat("NONINTEGER COUNTS\n")
  print(subset(model_data, abs(survived - round(survived)) > 1e-08 | abs(initial - round(initial)) >
    1e-08)[, c("experiment", "strain", "replicate", "survived", "initial")])
  direct_data <- droplevels(subset(model_data, experiment == "D2"))
  cat("DIRECT EXPERIMENT\n")
  print(with(direct_data, table(strain, preconditioned)))
  write.csv(direct_data, file.path(output_folder, "D2_analysis_data.csv"), row.names = FALSE)
  design_matrix <- model.matrix(~experiment + preconditioned + atacama_origin, model_data)
  cat("DESIGN RANK", qr(design_matrix)$rank, "/", ncol(design_matrix), "\n")
  direct_model <- glmmTMB(cbind(survived, non_survived) ~ preconditioned + (1 | strain), data = direct_data, family = betabinomial())
  direct_null <- update(direct_model, . ~ . - preconditioned)
  fixed_strain_model <- glmmTMB(cbind(survived, non_survived) ~ preconditioned + strain, data = direct_data,
    family = betabinomial())
  fixed_strain_null <- update(fixed_strain_model, . ~ . - preconditioned)
  no_strain_model <- update(fixed_strain_model, . ~ . - strain)
  interaction_model <- update(fixed_strain_model, . ~ . + preconditioned:strain)
  d2_models <- list(direct_random = direct_model, direct_random_null = direct_null, direct_fixed = fixed_strain_model,
    direct_fixed_null = fixed_strain_null, no_strain = no_strain_model, strain_interaction = interaction_model)
  validation <- do.call(rbind, lapply(names(d2_models), function(n) {
    x <- d2_models[[n]]
    data.frame(model = n, n = nobs(x), convergence = x$fit$convergence, pdHess = x$sdr$pdHess,
      AIC = AIC(x))
  }))
  write.csv(validation, file.path(output_folder, "model_validation.csv"), row.names = FALSE)
  print(validation)
  lrt <- function(a, b, label) {
    z <- anova(a, b)
    data.frame(test = label, chisq = z$Chisq[2], df = z$`Chi Df`[2], p = z$`Pr(>Chisq)`[2])
  }
  tests <- rbind(lrt(direct_null, direct_model, "D2 preconditioning random strain"), lrt(fixed_strain_null, fixed_strain_model,
    "D2 preconditioning fixed strain"), lrt(no_strain_model, fixed_strain_model, "D2 omnibus fixed strain"),
    lrt(fixed_strain_model, interaction_model, "D2 preconditioning by strain"))
  write.csv(tests, file.path(output_folder, "direct_tests.csv"), row.names = FALSE)
  print(tests)
  b <- fixef(direct_model)$cond["preconditioned1"]
  se <- sqrt(vcov(direct_model)$cond["preconditioned1", "preconditioned1"])
  write.csv(data.frame(OR = exp(b), lower = exp(b - 1.96 * se), upper = exp(b + 1.96 * se)), file.path(output_folder,
    "direct_preconditioning_effect.csv"), row.names = FALSE)
  saveRDS(d2_models, file.path(cache, "D2_models.rds"))
  capture.output(sessionInfo(), file = file.path(output_folder, "sessionInfo.txt"))

  d1_desiccated <- subset(model_data, experiment == "D1")
  d1_desiccated$strain <- factor(d1_desiccated$strain)
  d1_desiccated$proportion <- d1_desiccated$survived/d1_desiccated$initial
  keep <- c("860.22.04", "DF5050", "ES5", "PAP.22.29")
  d1_comparison_data <- droplevels(subset(d1_desiccated, strain %in% keep))
  stopifnot(nrow(d1_desiccated) == 40, nrow(d1_comparison_data) == 20, all(table(d1_desiccated$strain) == 5))
  stat <- function(y, g) {
    n <- length(y)
    k <- nlevels(g)
    mu <- mean(y)
    means <- tapply(y, g, mean)
    nn <- table(g)
    ssb <- sum(nn * (means - mu)^2)
    ssw <- sum((y - means[g])^2)
    (ssb/(k - 1))/(ssw/(n - k))
  }
  set.seed(1211)
  observed <- stat(d1_desiccated$proportion, d1_desiccated$strain)
  permutation_iterations <- 9999L
  permutation_statistics <- replicate(permutation_iterations, stat(sample(d1_desiccated$proportion), d1_desiccated$strain))
  global <- data.frame(test = "All8 strains: plate-level omnibus permutation F", F = observed, permutations = permutation_iterations,
    p = (1 + sum(permutation_statistics >= observed - 1e-12))/(permutation_iterations + 1))
  write.csv(global, file.path(output_folder, "experiment1_all8_global.csv"), row.names = FALSE)
  d1_model <- glmmTMB(cbind(survived, non_survived) ~ strain, data = d1_comparison_data, family = betabinomial())
  d1_null <- update(d1_model, . ~ 1)
  stopifnot(d1_model$fit$convergence == 0, isTRUE(d1_model$sdr$pdHess), d1_null$fit$convergence == 0,
    isTRUE(d1_null$sdr$pdHess))
  marginal_means <- emmeans(d1_model, ~strain, data = d1_comparison_data)
  d1_pairs <- as.data.frame(summary(pairs(marginal_means), infer = c(TRUE, TRUE), adjust = "none"))
  d1_pairs$p_Holm6 <- p.adjust(d1_pairs$p.value, "holm")
  d1_pairs$OR <- exp(d1_pairs$estimate)
  d1_pairs$OR_lower95 <- exp(d1_pairs$asymp.LCL)
  d1_pairs$OR_upper95 <- exp(d1_pairs$asymp.UCL)
  write.csv(d1_pairs, file.path(output_folder, "experiment1_panagrolaimus_model_pairs.csv"), row.names = FALSE)
  write.csv(anova(d1_null, d1_model), file.path(output_folder, "experiment1_panagrolaimus_global.csv"), row.names = FALSE)
  write.csv(data.frame(model = "Experiment1 four Panagrolaimus beta-binomial", convergence = d1_model$fit$convergence,
    pdHess = d1_model$sdr$pdHess, AIC = AIC(d1_model)), file.path(output_folder, "experiment1_model_validation.csv"),
    row.names = FALSE)
  saveRDS(list(model = d1_model, null = d1_null, data = d1_comparison_data), file.path(cache, "D1_models.rds"))

  marginal_means <- emmeans(d2_models$strain_interaction, ~preconditioned | strain)
  within_strain_results <- as.data.frame(summary(pairs(marginal_means, reverse = TRUE), infer = c(TRUE,
    TRUE), adjust = "none"))
  within_strain_results$p_Holm <- p.adjust(within_strain_results$p.value, "holm")
  write.csv(within_strain_results, file.path(tables, "D2_within_strain_preconditioning.csv"), row.names = FALSE)
  marginal_means <- emmeans(d2_models$strain_interaction, ~strain | preconditioned, data = direct_data)
  between_strain_results <- as.data.frame(summary(pairs(marginal_means), infer = c(TRUE, TRUE), adjust = "none"))
  between_strain_results$p_Holm_all42 <- p.adjust(between_strain_results$p.value, "holm")
  between_strain_results$odds_ratio <- exp(between_strain_results$estimate)
  between_strain_results$OR_lower_unadjusted95 <- exp(between_strain_results$asymp.LCL)
  between_strain_results$OR_upper_unadjusted95 <- exp(between_strain_results$asymp.UCL)
  write.csv(between_strain_results, file.path(tables, "D2_between_strain_pairs.csv"), row.names = FALSE)

  split_pair <- function(x) {
    z <- strsplit(as.character(x), " - ", fixed = TRUE)
    data.frame(strain_A = vapply(z, `[`, character(1), 1), strain_B = vapply(z, `[`, character(1),
      2))
  }
  pair_table <- rbind(data.frame(experiment = "D1", preconditioning = "None", split_pair(d1_pairs$contrast),
    OR = d1_pairs$OR, lower95 = d1_pairs$OR_lower95, upper95 = d1_pairs$OR_upper95, p_raw = d1_pairs$p.value, p_Holm = d1_pairs$p_Holm6),
    data.frame(experiment = "D2", preconditioning = ifelse(between_strain_results$preconditioned == 0, "None", "96 h"),
      split_pair(between_strain_results$contrast), OR = between_strain_results$odds_ratio, lower95 = between_strain_results$OR_lower_unadjusted95, upper95 = between_strain_results$OR_upper_unadjusted95,
      p_raw = between_strain_results$p.value, p_Holm = between_strain_results$p_Holm_all42))
  ids <- c(`860.22.04` = "WWL115", ALT.22.08 = "WWL002", PAP.22.29 = "WWL072", PAP.22.17 = "WWL709",
    ARO.22.05 = "WWL726", LAG.22.03 = "WWL505", TDT.23.02 = "WWL701", TDT.23.16 = "WWL702", GBG.22.06 = "WWL703",
    PWB.22.01 = "WWL724")
  map_id <- function(x) ifelse(x %in% names(ids), unname(ids[x]), x)
  pair_table$strain_A <- map_id(pair_table$strain_A)
  pair_table$strain_B <- map_id(pair_table$strain_B)
  write.csv(pair_table, file.path(tables, "S8_strain_pairwise_comparisons.csv"), row.names = FALSE)
  stopifnot(nrow(pair_table) == 48L, sum(pair_table$p_Holm < 0.05 & pair_table$experiment == "D1") ==
    3L, sum(pair_table$p_Holm < 0.05 & pair_table$experiment == "D2") == 9L)

  render_desiccation_figure(project_root, all_data, descriptive_summary)

  if (requireNamespace("DHARMa", quietly = TRUE)) {
    set.seed(1211)
    diagnostics <- do.call(rbind, lapply(list(joint = joint_model, D2_random = d2_models$direct_random,
      D2_interaction = d2_models$strain_interaction), function(model) {
      sim <- DHARMa::simulateResiduals(model, n = 1000, plot = FALSE)
      data.frame(uniformity = DHARMa::testUniformity(sim, plot = FALSE)$p.value, dispersion = DHARMa::testDispersion(sim,
        plot = FALSE)$p.value, zero_inflation = DHARMa::testZeroInflation(sim, plot = FALSE)$p.value)
    }))
    write.csv(diagnostics, file.path(output_folder, "residual_diagnostics.csv"), row.names = TRUE)
  }
  invisible(list(joint = joint_model, D2 = d2_models, pairwise = pair_table, summary = descriptive_summary))
}

render_desiccation_figure <- function(project_root, all_data, descriptive_summary) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(ggtext)
  })
  sys.source(file.path(project_root, "R/01_figure_style.R"), envir = environment())
  figures <- file.path(project_root, "Manuscript/figures")
  set.seed(1211)
  species_labels <- c(`860.22.04` = "<b><i>P. shuimeiren</i></b><br>(WWL115)", DF5050 = "<b><i>P. superbus</i></b><br>(DF5050)",
    ES5 = "<b><i>P. einhardi</i></b><br>(ES5)", GBG.22.06 = "<b><i>Acrobeloides</i> sp.</b><br>(WWL703)",
    PAP.22.29 = "<b><i>P. nebliphilus</i></b><br>(WWL072)", PWB.22.01 = "<b><i>P. tabacum</i></b><br>(WWL724)",
    TDT.23.02 = "<b><i>Acrobeles</i> sp.</b><br>(WWL701)", TDT.23.16 = "<b><i>S. ophioglossa</i></b><br>(WWL702)",
    ALT.22.08 = "<b><i>Panagrolaimus</i> sp.</b><br>(WWL002)", JU1365 = "<b><i>Panagrolaimus</i> sp.</b><br>(JU1365)",
    PS1159 = "<b><i>Panagrolaimus</i> sp.</b><br>(PS1159)", QG3759 = "<b><i>Panagrolaimus</i> sp.</b><br>(QG3759)",
    ARO.22.05 = "<b><i>A. goughiensis</i></b><br>(WWL726)", PAP.22.17 = "<b><i>A. tricornis</i></b><br>(WWL709)",
    LAG.22.03 = "<b><i>P. pavus</i></b><br>(WWL505)")

  label_treatments <- function(data) {
    data |>
      mutate(treatment = case_when(experiment == "D1" & group == "Control" ~ "Control", experiment ==
        "D1" & group == "Desiccation" ~ "Treatment", experiment == "D2" & group == "Not preconditioned" ~
        "Control", experiment == "D2" & group == "Preconditioned" ~ "Treatment", TRUE ~ as.character(group)))
  }
  publication_data <- all_data |>
    filter(!(experiment == "D2" & strain == "860.22.04")) |>
    label_treatments()
  publication_summary <- label_treatments(descriptive_summary)

  strain_orders <- list(D1 = c("ES5", "PAP.22.29", "860.22.04", "DF5050", "TDT.23.02", "TDT.23.16",
    "GBG.22.06", "PWB.22.01"), D2 = c("ES5", "PAP.22.29", "860.22.04", "ALT.22.08", "JU1365", "PS1159",
    "QG3759"))

  make_experiment_panel <- function(experiment_name, colours, shapes) {
    treatment_order <- names(colours)
    ds <- filter(publication_summary, experiment == experiment_name) |>
      mutate(treatment = factor(treatment, levels = treatment_order), strain = factor(strain, levels = strain_orders[[experiment_name]]))
    dp <- filter(publication_data, experiment == experiment_name) |>
      mutate(treatment = factor(treatment, levels = treatment_order), strain = factor(strain, levels = strain_orders[[experiment_name]]))
    dodge <- position_dodge(width = 0.3)
    ggplot(ds, aes(strain, 100 * replicate_mean, colour = treatment, shape = treatment)) + geom_point(data = dp,
      aes(y = 100 * recovery), alpha = 0.38, size = 2.1, position = position_jitterdodge(jitter.width = 0.045,
        jitter.height = 0, dodge.width = 0.3)) + geom_errorbar(aes(ymin = 100 * bootstrap_low,
      ymax = 100 * bootstrap_high), width = 0.12, linewidth = 0.55, position = dodge) + geom_point(size = 4.1,
      position = dodge) + scale_colour_manual(values = colours) + scale_shape_manual(values = shapes) +
      guides(colour = guide_legend(ncol = 2, byrow = TRUE), shape = guide_legend(ncol = 2, byrow = TRUE)) +
      scale_x_discrete(limits = strain_orders[[experiment_name]], labels = species_labels) + scale_y_continuous(breaks = seq(0,
      100, 25), expand = expansion(mult = c(0.03, 0.04))) + coord_cartesian(ylim = c(-3, 103)) +
      labs(x = NULL, y = "Survival (%)", colour = NULL, shape = NULL, tag = if (experiment_name ==
        "D1")
        "A" else "B") + theme_classic(base_size = 12) + theme(plot.background = element_rect(fill = "white",
      colour = NA), panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
      panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.35), panel.grid.minor = element_blank(),
      axis.text.x = ggtext::element_markdown(angle = 32, hjust = 1, size = 10.5, lineheight = 1.12,
        margin = margin(t = 3)), axis.text.y = element_text(size = 11, margin = margin(r = 3)),
      axis.title.y = element_text(size = 14, face = "plain", margin = margin(r = 6)), legend.position = "bottom",
      legend.justification = "center", legend.direction = "horizontal", legend.box = "horizontal",
      legend.text = element_text(size = 11), legend.key.width = grid::unit(1.1, "lines"), plot.tag = element_text(size = 18,
        face = "bold"), plot.margin = margin(4, 5, 0, 5))
  }

  p_D1 <- make_experiment_panel("D1", c(Control = "#666666", Treatment = "#C43C39"), c(Control = 16,
    Treatment = 17))
  p_D2 <- make_experiment_panel("D2", c(Control = "#666666", Treatment = "#C43C39"), c(Control = 16,
    Treatment = 17))
  p_publication_panels <- (p_D1/p_D2) + plot_layout(heights = c(1, 1), guides = "collect") + plot_annotation(theme = theme(plot.margin = margin(5,
    5, 5, 5))) & theme(legend.position = "bottom")

  ggsave(file.path(figures, "Figure_4_desiccation.png"), p_publication_panels, width = 9, height = 6.8,
    dpi = 300)

  invisible(figures)
}
