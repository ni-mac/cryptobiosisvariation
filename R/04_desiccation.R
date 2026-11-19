# Reproducible desiccation analysis (D1 D1, D2 D2); D3 excluded from the current manuscript.
# Raw observations preserved byte-identically; historical codes are internal model IDs.
# Bootstrap and original random seed retained for manuscript figure parity.
run_desiccation <- function(project_root) {
suppressPackageStartupMessages({library(readxl);library(dplyr);library(tidyr);library(ggplot2);library(patchwork);library(ggtext);library(glmmTMB);library(emmeans)})
options(stringsAsFactors=FALSE,scipen=999)
set.seed(1211)
raw <- file.path(project_root,"Data","raw","desiccation")
out <- file.path(project_root,"Results","statistics","desiccation")
tables <- file.path(project_root,"Results","tables","desiccation")
cache <- file.path(project_root,"Results","cache","desiccation")
figures <- file.path(project_root,"Manuscript","figures")
for(d in c(out,tables,cache,figures)) dir.create(d,recursive=TRUE,showWarnings=FALSE)
safe_write_csv <- function(x, path) {
  tryCatch(
    write.csv(x, path, row.names = FALSE),
    error = function(e) warning("Could not update open output file: ", basename(path))
  )
}

z95 <- qnorm(0.975)
wilson <- function(x, n) {
  p <- x / n
  den <- 1 + z95^2 / n
  mid <- (p + z95^2 / (2 * n)) / den
  half <- z95 * sqrt(p * (1 - p) / n + z95^2 / (4 * n^2)) / den
  data.frame(low = pmax(0, mid - half), high = pmin(1, mid + half))
}

summarise_recovery <- function(x) {
  # Preserve the original bootstrap group order with neutral experiment labels.
  x$experiment <- factor(x$experiment, levels = c("D2", "D1"))
  ans <- x |>
    group_by(experiment, strain, group) |>
    summarise(
      n_replicates = n(), recovered = sum(recovered), initial = sum(initial),
      pooled_recovery = recovered / initial,
      replicate_mean = mean(recovery), replicate_sd = sd(recovery),
      .groups = "drop"
    )
  ans <- bind_cols(ans, wilson(ans$recovered, ans$initial))
  set.seed(1211)
  historical_rng_gap_consumed <- FALSE
  boot <- x |>
    group_by(experiment, strain, group) |>
    group_modify(~{
      # The historical group order inserted six D3 groups (five plates each)
      # between D2 and D1. Advance RNG only, without reading D3 data,
      # to retain the accepted D1 bootstrap intervals after D3 removal.
      if (as.character(.y$experiment) == "D1" && !historical_rng_gap_consumed) {
        for (j in seq_len(6L)) replicate(20000L, sample.int(5L, 5L, replace=TRUE))
        historical_rng_gap_consumed <<- TRUE
      }
      y <- .x$recovery
      b <- replicate(20000, mean(sample(y, length(y), replace = TRUE)))
      tibble(bootstrap_low = unname(quantile(b, 0.025)),
             bootstrap_high = unname(quantile(b, 0.975)))
    }) |>
    ungroup()
  left_join(ans, boot, by = c("experiment", "strain", "group"))
}

# D1: five replicate plates per strain and group; no preconditioning.
D10 <- read_excel(file.path(raw, "D1", "D1_survival.xlsx"), sheet = "Nematoden Daten")
D1 <- transmute(
  D10,
  experiment = "D1", strain = gsub(" ", ".", trimws(Strain), fixed = TRUE),
  species = trimws(Species), group = if_else(Treatment == "Desiccation", "Desiccation", "Control"),
  replicate = as.integer(Replicate), recovered = as.integer(Alive),
  initial = as.integer(Total), failure = as.integer(Dead), missing = 0L,
  recovery = recovered / initial,
  note = "No preconditioning; 8 h desiccation, 17 h rehydration"
)
D1$strain <- recode(D1$strain, `TDT.23.02` = "TDT.23.02", `TDT.23.16` = "TDT.23.16",
                      `PAP.22.29` = "PAP.22.29", `PWB.22.01` = "PWB.22.01",
                      `GBG22.06` = "GBG.22.06")

# D2: same drying/rehydration protocol in both groups; only 96 h preconditioning differs.
D20 <- read.csv(file.path(project_root, "Results/inputs/D2_survival.csv"), check.names = FALSE)
shuimeiren_technical <- read.csv(
  file.path(project_root, "Results/inputs/D2_technical_counts.csv"),
  check.names = FALSE
) |>
  transmute(
    Strain, Species, reproduction,
    `Sampling location` = Sampling_location,
    Treatment, Replicate,
    # TR1--TR3 are repeated counts from the same biological plate, not
    # independent animals. Retain one observation per plate and represent the
    # mean of the three 100-animal counts on its original denominator.
    Alive = Alive_mean_per_100,
    Total = 100,
    Dead = 100 - Alive_mean_per_100,
    Survival_percent
  )
D20 <- bind_rows(D20, shuimeiren_technical)
D2 <- transmute(
  D20,
  experiment = "D2", strain = recode(Strain, ALT2208 = "ALT.22.08", PAP2229 = "PAP.22.29",
                                           `ES5 (P. einhardi)` = "ES5"),
  species = Species,
  group = if_else(Treatment == "Acclimated", "Preconditioned", "Not preconditioned"),
  replicate = as.integer(sub("R", "", Replicate)), recovered = as.numeric(Alive),
  initial = as.numeric(Total), failure = as.numeric(Dead), missing = 0,
  recovery = recovered / initial,
  note = "48 h desiccation and 24 h rehydration; groups differ only in 96 h preconditioning"
)

# D3 raw files remain preserved but are not read by this analysis.
all_data <- bind_rows(D1, D2)
stopifnot(all(all_data$recovered >= 0), all(all_data$initial > 0),
          all(all_data$recovered + all_data$failure + all_data$missing == all_data$initial))

# Analysis metadata for the two desiccation protocols in the manuscript.
all_data <- all_data |>
  mutate(
    species = case_when(
      strain == "PAP.22.29" ~ "Panagrolaimus nebliphilus",
      TRUE ~ species
    ),
    origin = case_when(
      strain %in% c("ALT.22.08", "ARO.22.05", "LAG.22.03", "PAP.22.17",
                    "PAP.22.29", "TDT.23.02", "TDT.23.16") ~ "Atacama",
      strain %in% c("860.22.04", "GBG.22.06", "PWB.22.01") ~ "Namib",
      strain == "ES5" ~ "Germany",
      strain == "DF5050" ~ "Iceland",
      strain == "JU1365" ~ "India",
      strain == "PS1159" ~ "USA",
      strain == "QG3759" ~ "Ukraine",
      TRUE ~ NA_character_
    ),
    atacama_origin = if_else(origin == "Atacama", "Atacama", "Non-Atacama"),
    desiccated = if_else(
      (experiment == "D1" & group == "Control") |
        (experiment == "D3" & group == "Control"), 0L, 1L
    ),
    preconditioned = if_else(
      experiment == "D3" |
        (experiment == "D2" & group == "Preconditioned"), 1L, 0L
    ),
    protocol_group = case_when(
      experiment == "D1" & group == "Control" ~ "D1_control",
      experiment == "D1" ~ "D1_desiccation_no_preconditioning",
      experiment == "D2" & group == "Preconditioned" ~ "D2_desiccation_96h_preconditioning",
      experiment == "D2" ~ "D2_desiccation_no_preconditioning",
      experiment == "D3" & group == "Control" ~ "D3_72h_preconditioning_control",
      experiment == "D3" ~ "D3_72h_preconditioning_desiccation"
    ),
    biological_unit = "plate/filter",
    included_joint_model = TRUE,
    survived = recovered,
    not_survived = initial - recovered,
    survival_proportion = recovery
  )

summary <- summarise_recovery(all_data)
write.csv(all_data, file.path(out, "analysis_data_all_groups.csv"), row.names=FALSE)
write.csv(summary, file.path(tables, "survival_descriptive_bootstrap.csv"), row.names=FALSE)

d_all <- all_data |>
  mutate(
    # Historical model input rounded technical-repeat means to integer counts.
    # Preserve that model convention separately from unrounded descriptive data.
    survived = round(survived),
    recovered = round(recovered),
    failure = initial - recovered,
    recovery = recovered / initial,
    survival_proportion = recovery,
    not_survived = initial - survived,
    non_survived = initial - survived,
    experiment = relevel(factor(experiment), "D2"),
    strain = factor(strain),
    desiccated = relevel(factor(desiccated), "0"),
    preconditioned = relevel(factor(preconditioned), "0"),
    atacama_origin = relevel(factor(atacama_origin), "Non-Atacama")
  )

# The inferential endpoint is survival after desiccation and rehydration.
# Non-desiccated controls remain relevant descriptively but are not observations
# of this endpoint and therefore are not included in the joint GLMM.
d <- d_all |>
  filter(desiccated == "1") |>
  droplevels()

stopifnot(
  nrow(d) == 82L,
  nlevels(d$strain) == 12L,
  all(d$survived >= 0), all(d$non_survived >= 0),
  all(d$survived + d$non_survived == d$initial)
)

form_full <- cbind(survived, non_survived) ~
  experiment + preconditioned + atacama_origin + (1 | strain)

# The beta-binomial family is the primary model because plate-level survival
# is overdispersed relative to an ordinary binomial model.
fit <- glmmTMB(form_full, data = d, family = betabinomial(link = "logit"))

# Nested models are inferential checks for prespecified terms, not separate
# analyses of the three experiments.
fit_no_origin <- update(fit, . ~ . - atacama_origin)
fit_no_preconditioning <- update(fit, . ~ . - preconditioned)
fit_no_experiment <- update(fit, . ~ . - experiment)
fit_binomial <- glmmTMB(form_full, data = d, family = binomial(link = "logit"))

lrt_row <- function(reduced, full, term) {
  a <- anova(reduced, full)
  data.frame(
    term = term,
    likelihood_ratio_chisq = a$Chisq[2],
    df = a$`Chi Df`[2],
    p_value = a$`Pr(>Chisq)`[2]
  )
}

term_tests <- bind_rows(
  lrt_row(fit_no_origin, fit, "Atacama origin"),
  lrt_row(fit_no_preconditioning, fit, "Preconditioning"),
  lrt_row(fit_no_experiment, fit, "Experimental protocol")
)

model_validation <- data.frame(
  model = c("Primary beta-binomial GLMM", "Ordinary binomial GLMM"),
  n = c(nobs(fit), nobs(fit_binomial)),
  logLik = c(as.numeric(logLik(fit)), as.numeric(logLik(fit_binomial))),
  df = c(attr(logLik(fit), "df"), attr(logLik(fit_binomial), "df")),
  AIC = c(AIC(fit), AIC(fit_binomial)),
  convergence = c(fit$fit$convergence, fit_binomial$fit$convergence),
  positive_definite_hessian = c(fit$sdr$pdHess, fit_binomial$sdr$pdHess)
) |>
  mutate(delta_AIC = AIC - min(AIC))

coef_table <- as.data.frame(summary(fit)$coefficients$cond)
coef_table$term <- rownames(coef_table)
rownames(coef_table) <- NULL
coef_table <- coef_table |>
  dplyr::select(term, everything())

vc <- VarCorr(fit)$cond$strain
variance_table <- data.frame(
  component = "Strain random intercept",
  variance = as.numeric(vc[1]),
  standard_deviation = attr(vc, "stddev")[[1]]
)

# Model-based survival by preconditioning status within the experiment in
# which both levels were directly observed (D2).
emm <- emmeans(fit, ~ preconditioned | experiment, type = "response")
estimated_survival <- as.data.frame(emm) |>
  filter(experiment == "D2")

write.csv(d, file.path(out, "analysis_data_desiccated.csv"), row.names = FALSE)
write.csv(term_tests, file.path(out, "fixed_effect_likelihood_ratio_tests.csv"), row.names = FALSE)
write.csv(model_validation, file.path(out, "distribution_validation.csv"), row.names = FALSE)
write.csv(coef_table, file.path(out, "fixed_effect_coefficients.csv"), row.names = FALSE)
write.csv(variance_table, file.path(out, "strain_random_effect.csv"), row.names = FALSE)
write.csv(estimated_survival, file.path(out, "estimated_survival_by_treatment.csv"), row.names = FALSE)
saveRDS(fit, file.path(cache, "joint_model.rds"))


for(v in c('experiment','preconditioned','strain','atacama_origin'))d[[v]]<-factor(d[[v]])
cat('DESIGN TABLE\n');print(with(d,table(experiment,preconditioned)))
cat('NONINTEGER COUNTS\n');print(subset(d,abs(survived-round(survived))>1e-8 | abs(initial-round(initial))>1e-8)[,c('experiment','strain','replicate','survived','initial')])
k<-droplevels(subset(d,experiment=='D2'))
cat('DIRECT EXPERIMENT\n');print(with(k,table(strain,preconditioned)))
write.csv(k,file.path(out,'D2_analysis_data.csv'),row.names=FALSE)
mm<-model.matrix(~experiment+preconditioned+atacama_origin,d)
cat('DESIGN RANK',qr(mm)$rank,'/',ncol(mm),'\n')
# Direct within-protocol comparison, preserving random strain effect of old model.
m<-glmmTMB(cbind(survived,non_survived)~preconditioned+(1|strain),data=k,family=betabinomial())
m0<-update(m,.~.-preconditioned)
# Fixed-strain sensitivity and omnibus strain test; no boundary variance test.
f<-glmmTMB(cbind(survived,non_survived)~preconditioned+strain,data=k,family=betabinomial())
f0<-update(f,.~.-preconditioned)
fs0<-update(f,.~.-strain)
fi<-update(f,.~.+preconditioned:strain)
fits<-list(direct_random=m,direct_random_null=m0,direct_fixed=f,direct_fixed_null=f0,no_strain=fs0,strain_interaction=fi)
validation<-do.call(rbind,lapply(names(fits),function(n){x<-fits[[n]];data.frame(model=n,n=nobs(x),convergence=x$fit$convergence,pdHess=x$sdr$pdHess,AIC=AIC(x))}))
write.csv(validation,file.path(out,'model_validation.csv'),row.names=FALSE);print(validation)
lrt<-function(a,b,label){z<-anova(a,b);data.frame(test=label,chisq=z$Chisq[2],df=z$`Chi Df`[2],p=z$`Pr(>Chisq)`[2])}
tests<-rbind(lrt(m0,m,'D2 preconditioning random strain'),lrt(f0,f,'D2 preconditioning fixed strain'),lrt(fs0,f,'D2 omnibus fixed strain'),lrt(f,fi,'D2 preconditioning by strain'))
write.csv(tests,file.path(out,'direct_tests.csv'),row.names=FALSE);print(tests)
b<-fixef(m)$cond['preconditioned1'];se<-sqrt(vcov(m)$cond['preconditioned1','preconditioned1'])
write.csv(data.frame(OR=exp(b),lower=exp(b-1.96*se),upper=exp(b+1.96*se)),file.path(out,'direct_preconditioning_effect.csv'),row.names=FALSE)
saveRDS(fits,file.path(cache,'D2_models.rds'))
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))

r<-subset(d,experiment=='D1');r$strain<-factor(r$strain);r$proportion<-r$survived/r$initial
keep<-c('860.22.04','DF5050','ES5','PAP.22.29')
x<-droplevels(subset(r,strain %in% keep))
stopifnot(nrow(r)==40,nrow(x)==20,all(table(r$strain)==5))
# Omnibus comparison: all eight strains, complete plate-level permutation.
stat<-function(y,g){n<-length(y);k<-nlevels(g);mu<-mean(y);means<-tapply(y,g,mean);nn<-table(g);ssb<-sum(nn*(means-mu)^2);ssw<-sum((y-means[g])^2);(ssb/(k-1))/(ssw/(n-k))}
set.seed(1211);observed<-stat(r$proportion,r$strain)
b<-9999L;sim<-replicate(b,stat(sample(r$proportion),r$strain))
global<-data.frame(test='All8 strains: plate-level omnibus permutation F',F=observed,permutations=b,p=(1+sum(sim>=observed-1e-12))/(b+1))
write.csv(global,file.path(out,'experiment1_all8_global.csv'),row.names=FALSE)
# Direct post-desiccation survival comparisons, not control-normalized effects.
m<-glmmTMB(cbind(survived,non_survived)~strain,data=x,family=betabinomial())
m0<-update(m,.~1)
stopifnot(m$fit$convergence==0,isTRUE(m$sdr$pdHess),m0$fit$convergence==0,isTRUE(m0$sdr$pdHess))
e<-emmeans(m,~strain,data=x)
c<-as.data.frame(summary(pairs(e),infer=c(TRUE,TRUE),adjust='none'))
c$p_Holm6<-p.adjust(c$p.value,'holm');c$OR<-exp(c$estimate);c$OR_lower95<-exp(c$asymp.LCL);c$OR_upper95<-exp(c$asymp.UCL)
write.csv(c,file.path(out,'experiment1_panagrolaimus_model_pairs.csv'),row.names=FALSE)
write.csv(anova(m0,m),file.path(out,'experiment1_panagrolaimus_global.csv'),row.names=FALSE)
write.csv(data.frame(model='Experiment1 four Panagrolaimus beta-binomial',convergence=m$fit$convergence,pdHess=m$sdr$pdHess,AIC=AIC(m)),file.path(out,'experiment1_model_validation.csv'),row.names=FALSE)
saveRDS(list(model=m,null=m0,data=x),file.path(cache,'D1_models.rds'))


e <- emmeans(fits$strain_interaction, ~ preconditioned | strain)
within <- as.data.frame(summary(pairs(e, reverse=TRUE), infer=c(TRUE,TRUE), adjust="none"))
within$p_Holm <- p.adjust(within$p.value, "holm")
write.csv(within, file.path(tables,"D2_within_strain_preconditioning.csv"), row.names=FALSE)
e <- emmeans(fits$strain_interaction, ~ strain | preconditioned, data=k)
b <- as.data.frame(summary(pairs(e), infer=c(TRUE,TRUE), adjust="none"))
b$p_Holm_all42 <- p.adjust(b$p.value,"holm")
b$odds_ratio <- exp(b$estimate)
b$OR_lower_unadjusted95 <- exp(b$asymp.LCL)
b$OR_upper_unadjusted95 <- exp(b$asymp.UCL)
write.csv(b,file.path(tables,"D2_between_strain_pairs.csv"),row.names=FALSE)
# Consolidated S8: current experiment numbers and current strain identifiers.
a <- c
split_pair <- function(x) {
 z <- strsplit(as.character(x)," - ",fixed=TRUE)
 data.frame(strain_A=vapply(z,`[`,character(1),1),strain_B=vapply(z,`[`,character(1),2))
}
pair_table <- rbind(
 data.frame(experiment="D1",preconditioning="None",split_pair(a$contrast),OR=a$OR,lower95=a$OR_lower95,upper95=a$OR_upper95,p_raw=a$p.value,p_Holm=a$p_Holm6),
 data.frame(experiment="D2",preconditioning=ifelse(b$preconditioned==0,"None","96 h"),split_pair(b$contrast),OR=b$odds_ratio,lower95=b$OR_lower_unadjusted95,upper95=b$OR_upper_unadjusted95,p_raw=b$p.value,p_Holm=b$p_Holm_all42))
ids <- c("860.22.04"="WWL115","ALT.22.08"="WWL002","PAP.22.29"="WWL072","PAP.22.17"="WWL709","ARO.22.05"="WWL726","LAG.22.03"="WWL505","TDT.23.02"="WWL701","TDT.23.16"="WWL702","GBG.22.06"="WWL703","PWB.22.01"="WWL724")
map_id <- function(x) ifelse(x %in% names(ids),unname(ids[x]),x)
pair_table$strain_A <- map_id(pair_table$strain_A)
pair_table$strain_B <- map_id(pair_table$strain_B)
write.csv(pair_table,file.path(tables,"S8_strain_pairwise_comparisons.csv"),row.names=FALSE)
stopifnot(nrow(pair_table)==48L,sum(pair_table$p_Holm<.05 & pair_table$experiment=="D1")==3L,sum(pair_table$p_Holm<.05 & pair_table$experiment=="D2")==9L)

render_desiccation_figure(project_root, all_data, summary)

# Residual checks are diagnostics, not manuscript hypothesis tests.
if(requireNamespace("DHARMa",quietly=TRUE)) {
 set.seed(1211)
 diagnostics <- do.call(rbind,lapply(list(joint=fit,D2_random=fits$direct_random,D2_interaction=fits$strain_interaction),function(model){
  sim <- DHARMa::simulateResiduals(model,n=1000,plot=FALSE)
  data.frame(uniformity=DHARMa::testUniformity(sim,plot=FALSE)$p.value,dispersion=DHARMa::testDispersion(sim,plot=FALSE)$p.value,zero_inflation=DHARMa::testZeroInflation(sim,plot=FALSE)$p.value)
 }))
 write.csv(diagnostics,file.path(out,"residual_diagnostics.csv"),row.names=TRUE)
}
invisible(list(joint=fit,D2=fits,pairwise=pair_table,summary=summary))
}

# Render only, without refitting models or changing numerical outputs.
render_desiccation_figure <- function(project_root, all_data, summary) {
 suppressPackageStartupMessages({library(dplyr);library(ggplot2);library(patchwork);library(ggtext)})
 sys.source(file.path(project_root,'R/01_figure_style.R'),envir=environment())
 figures <- file.path(project_root,'Manuscript/figures')
set.seed(1211)
species_labels <- c(
  "860.22.04" = "<b><i>P. shuimeiren</i></b><br>(WWL115)",
  "DF5050" = "<b><i>P. superbus</i></b><br>(DF5050)",
  "ES5" = "<b><i>P. einhardi</i></b><br>(ES5)",
  "GBG.22.06" = "<b><i>Acrobeloides</i> sp.</b><br>(WWL703)",
  "PAP.22.29" = "<b><i>P. nebliphilus</i></b><br>(WWL072)",
  "PWB.22.01" = "<b><i>P. tabacum</i></b><br>(WWL724)",
  "TDT.23.02" = "<b><i>Acrobeles</i> sp.</b><br>(WWL701)",
  "TDT.23.16" = "<b><i>S. ophioglossa</i></b><br>(WWL702)",
  "ALT.22.08" = "<b><i>Panagrolaimus</i> sp.</b><br>(WWL002)",
  "JU1365" = "<b><i>Panagrolaimus</i> sp.</b><br>(JU1365)",
  "PS1159" = "<b><i>Panagrolaimus</i> sp.</b><br>(PS1159)",
  "QG3759" = "<b><i>Panagrolaimus</i> sp.</b><br>(QG3759)",
  "ARO.22.05" = "<b><i>A. goughiensis</i></b><br>(WWL726)",
  "PAP.22.17" = "<b><i>A. tricornis</i></b><br>(WWL709)",
  "LAG.22.03" = "<b><i>P. pavus</i></b><br>(WWL505)"
)

# Preserve manuscript plotting dataset: historical figure omits the six D2
# WWL115 plate-level markers although its means/CI and all models include them.
# Keep explicit for byte/pixel parity; see Documentation/desiccation_provenance.md.
publication_data <- all_data |>
  filter(!(experiment == "D2" & strain == "860.22.04")) |>
  mutate(treatment = case_when(
    experiment == "D1" & group == "Control" ~ "Control",
    experiment == "D1" & group == "Desiccation" ~ "Treatment",
    experiment == "D2" & group == "Not preconditioned" ~ "Control",
    experiment == "D2" & group == "Preconditioned" ~ "Treatment",
    experiment == "D3" & group == "Control" ~ "Control",
    experiment == "D3" & group == "Desiccation after preconditioning" ~
      "Treatment",
    TRUE ~ as.character(group)
  ))

publication_summary <- summary |>
  mutate(treatment = case_when(
    experiment == "D1" & group == "Control" ~ "Control",
    experiment == "D1" & group == "Desiccation" ~ "Treatment",
    experiment == "D2" & group == "Not preconditioned" ~ "Control",
    experiment == "D2" & group == "Preconditioned" ~ "Treatment",
    experiment == "D3" & group == "Control" ~ "Control",
    experiment == "D3" & group == "Desiccation after preconditioning" ~
      "Treatment",
    TRUE ~ as.character(group)
  ))

strain_orders <- list(
  D1 = c("ES5", "PAP.22.29", "860.22.04", "DF5050", "TDT.23.02",
           "TDT.23.16", "GBG.22.06", "PWB.22.01"),
  D2 = c("ES5", "PAP.22.29", "860.22.04", "ALT.22.08", "JU1365", "PS1159", "QG3759"),
  D3 = c("PAP.22.17", "ARO.22.05", "LAG.22.03")
)

make_experiment_panel <- function(experiment_name, colours, shapes) {
  treatment_order <- names(colours)
  ds <- filter(publication_summary, experiment == experiment_name) |>
    mutate(treatment = factor(treatment, levels = treatment_order),
           strain = factor(strain, levels = strain_orders[[experiment_name]]))
  dp <- filter(publication_data, experiment == experiment_name) |>
    mutate(treatment = factor(treatment, levels = treatment_order),
           strain = factor(strain, levels = strain_orders[[experiment_name]]))
  dodge <- position_dodge(width = 0.30)
  ggplot(ds, aes(strain, 100 * replicate_mean, colour = treatment, shape = treatment)) +
    geom_point(data = dp, aes(y = 100 * recovery), alpha = 0.38, size = 2.1,
               position = position_jitterdodge(jitter.width = 0.045,
                                               jitter.height = 0,
                                               dodge.width = 0.30)) +
    geom_errorbar(aes(ymin = 100 * bootstrap_low, ymax = 100 * bootstrap_high),
                  width = 0.12, linewidth = 0.55, position = dodge) +
    geom_point(size = 4.1, position = dodge) +
    scale_colour_manual(values = colours) +
    scale_shape_manual(values = shapes) +
    guides(colour = guide_legend(ncol = 2, byrow = TRUE),
           shape = guide_legend(ncol = 2, byrow = TRUE)) +
    scale_x_discrete(limits = strain_orders[[experiment_name]], labels = species_labels) +
    scale_y_continuous(breaks = seq(0, 100, 25), expand = expansion(mult = c(0.03, 0.04))) +
    coord_cartesian(ylim = c(-3, 103)) +
    labs(x = NULL, y = "Survival (%)", colour = NULL, shape = NULL,
         tag=if(experiment_name=="D1") "A" else "B") +
    theme_classic(base_size = 12) +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
      panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      axis.text.x = ggtext::element_markdown(angle = 32, hjust = 1, size = 10.5,
                                             lineheight = 1.12, margin=margin(t=3)),
      axis.text.y = element_text(size = 11, margin=margin(r=3)),
      axis.title.y = element_text(size = 14, face = "plain", margin=margin(r=6)),
      legend.position = "bottom",
      legend.justification = "center",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.text = element_text(size = 11),
      legend.key.width = grid::unit(1.1, "lines"),
      plot.tag=element_text(size=18,face="bold"),
      plot.margin = margin(4, 5, 0, 5)
    )
}

p_D1 <- make_experiment_panel(
  "D1",
  c("Control" = "#666666", "Treatment" = "#C43C39"),
  c("Control" = 16, "Treatment" = 17)
)
p_D2 <- make_experiment_panel(
  "D2",
  c("Control" = "#666666", "Treatment" = "#C43C39"),
  c("Control" = 16, "Treatment" = 17)
)
p_publication_panels <- (p_D1 / p_D2) +
  plot_layout(heights = c(1, 1), guides = "collect") +
  plot_annotation(theme = theme(plot.margin = margin(5, 5, 5, 5))) &
  theme(legend.position = "bottom")

p_publication <- p_publication_panels

ggsave(file.path(figures, "Figure_4_desiccation.png"), p_publication,
       width = 9, height = 6.8, dpi = 300)




 invisible(figures)
}
