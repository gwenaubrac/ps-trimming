## ---------------------------
##
## Program: PS-Trimming of new antidepressant users
##
## Purpose: The goal is to understand what information PS-trimming gives
## regarding confounding and EMM in the sample. Specifically, we want to:
## - Apply asymmetric trimming to new antidepressant users to estimate IRRs/IRDs/HRs of all-cause mortality
## - Assess the impact of asymmetric trimming (2%, 5%, 10%, 20%, PF) on the ATE, ATT, and ATU
## - Compare these findings with PS-stratified estimates
##
## To achieve this, we will perform the following steps:
## 1. Estimate initial propensity score (PS) and preference score (PF) based on set of covariates
## 2. Stratify sample into deciles of PS and estimate HR/IRR/IRD
## 3. Trim sample using asymmetric trimming at 2%, 5%, 10%, 20%, and 30-70% preference-score trimming
## 4. Re-estimate the PS and compute IPTW, ATU, and ATT weights
## 5. Estimate IRRs, IRDs, and HRs**
## **Bootstraps (R = 2,000) to obtain robust CIs for IRDs and IRRs in a separate program. 
##
## Authors: Gwen Aubrac & Michael Webster-Clark
##
## Date Created: 2025-02-24
##
## ---------------------------
##
## Notes: Incidence rates calculated per 100 person-years.
##
## ---------------------------


#### 0. SET UP ####
library(dplyr)
library(ggplot2)
library(magrittr)
library(survival)
library(writexl)

cohort <- readRDS('Z:/EPI/Protocol 24_004042/Gwen - ipcw_methods/results/antidepressant/all-cause mortality/main/intermediate/cohort_iptw.rds')

cohort %<>%
  select(id,
         trt,
         age_at_entry,
         entry_date,
         itt_exit_date,
         itt_event,
         age_group,
         sex,
         year,
         ethnicity,
         deprivation, 
         anxiety_base, 
         arrhythmia_base, 
         cardiomyopathy_base,
         cerebrovascular_disease_base,
         chronic_renal_disease_base, 
         copd_base, 
         depression_base, 
         epilepsy_base,
         heart_failure_base, 
         hepatic_disease_base, 
         hyperlipidemia_base,
         hypertension_base,
         ischemic_heart_disease_base,
         lvh_base,
         myocardial_infarction_base,
         pacemaker_base,
         pvd_base,
         stroke_base,
         suicidal_ideation_self_harm_base,
         valvular_heart_disease_base) %>% 
  rename(death = itt_event) %>% 
  mutate(follow_up = as.numeric(itt_exit_date - entry_date),
         new_trt = as.factor(if_else(trt == "ssri", 1, 0)))

#### 1. ESTIMATE PS AND PF ####

# ps: propensity score
# ps_marg: marginal propensity score (i.e., prevalence of trt)
# pf: preference score

ps_model <- glm(
  new_trt ~ age_group + sex + year + ethnicity + deprivation + anxiety_base + 
    arrhythmia_base + cardiomyopathy_base +cerebrovascular_disease_base +
    chronic_renal_disease_base + copd_base + depression_base + epilepsy_base +
    heart_failure_base + hepatic_disease_base + hyperlipidemia_base + hypertension_base +
    ischemic_heart_disease_base + lvh_base + myocardial_infarction_base +pacemaker_base +
    pvd_base + stroke_base + suicidal_ideation_self_harm_base + valvular_heart_disease_base,
  family = binomial, 
  data = cohort
)

cohort$ps <- predict(ps_model, type = "response")

ps_model_marg <- glm(
  new_trt ~ 1,
  family = binomial,
  data = cohort
)

cohort$ps_marg <- predict(ps_model_marg, type = "response")

cohort %<>%
  mutate(
    logit_ps = log(ps/(1-ps)),
    logit_prev = log(ps_marg/(1-ps_marg)),
    logit_pf = logit_ps - logit_prev,
    odds_pf = exp(logit_pf),
    pf = odds_pf/(1+odds_pf)
  )

p <- ggplot(data=cohort, aes(x=pf, group=new_trt, fill=new_trt)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_manual(values=c("0"="blue", "1"="red"), 
                    labels=c("0"="SNRI", "1"="SSRI"),
                    name = "Treatment") +
  theme_minimal() + 
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(title="Preference Score Distribution in Untrimmed Population", x = "preference score") 

p

ggsave("C:/Users/gwen.aubrac/Desktop/trimming_res/plots/pf_untrimmed.pdf", 
       plot = p,  
       width = 8, 
       height = 6,
       units = "in")

#### 2. PS-STRATIFY AND ESTIMATE ATE/ATT/ATU ####

## a) Get cutoff points and stratify data
ps_deciles <- quantile(cohort$ps, probs = seq(0.1, 0.9, by = 0.1))
cohort$decile <- cut(cohort$ps, breaks = c(0, ps_deciles, 1), labels = 1:10, include.lowest = TRUE)

test <- cohort %>% group_by(decile) %>% summarize(mean_ps = mean(ps))
ps_deciles

cohort_strat <- split(cohort, cohort$decile)

## b) Estimate HRs 

coxph_res_strat <- data.frame(
  decile = rep(c(1:10), each = 1),
  coef = NA,
  hr = NA,
  hr_lower = NA,
  hr_upper = NA
)

fit_coxph_strat <- function(decile) {
  
  coxph_fit <- coxph(Surv(follow_up, death) ~ new_trt,
                      data = cohort_strat[[decile]])
  
  coxph_res_strat$coef[coxph_res_strat$decile == decile] <- summary(coxph_fit)$coefficients[, "coef"]
  coxph_res_strat$hr[coxph_res_strat$decile == decile] <- exp(summary(coxph_fit)$coefficients[, "coef"])
  coxph_res_strat$hr_lower[coxph_res_strat$decile == decile] <- summary(coxph_fit)$conf.int[, "lower .95"]
  coxph_res_strat$hr_upper[coxph_res_strat$decile == decile] <- summary(coxph_fit)$conf.int[, "upper .95"]
  
  coxph_res_strat <<- coxph_res_strat
}


fit_coxph_strat(decile = 1)
fit_coxph_strat(decile = 2)
fit_coxph_strat(decile = 3)
fit_coxph_strat(decile = 4)
fit_coxph_strat(decile = 5)
fit_coxph_strat(decile = 6)
fit_coxph_strat(decile = 7)
fit_coxph_strat(decile = 8)
fit_coxph_strat(decile = 9)
fit_coxph_strat(decile = 10)

## d) Estimate IRRs and IRDs

ir_res_strat <- data.frame(
  decile = rep(c(1:10), each = 1),
  ir_trt1 = NA,
  ir_trt0 = NA,
  irr = NA,
  ird = NA
)

get_ir_strat <- function(decile) {
  
  # for treated
  ir_all_trt1 <- cohort_strat[[decile]] %>% 
    filter(new_trt == 1) %>% 
    summarize(
      sum_pt = sum(follow_up),
      sum_death = sum(death)
    )
  
  ir_res_strat$ir_trt1[ir_res_strat$decile == decile] <- ir_all_trt1$sum_death/ir_all_trt1$sum_pt*100

  # for untreated
  ir_all_trt0 <- cohort_strat[[decile]] %>% 
    filter(new_trt == 0) %>% 
    summarize(
      sum_pt = sum(follow_up),
      sum_death = sum(death)
    )
  
  ir_res_strat$ir_trt0[ir_res_strat$decile == decile] <- ir_all_trt0$sum_death/ir_all_trt0$sum_pt*100

  ir_res_strat <<- ir_res_strat
}

get_ir_strat(decile = 1)
get_ir_strat(decile = 2)
get_ir_strat(decile = 3)
get_ir_strat(decile = 4)
get_ir_strat(decile = 5)
get_ir_strat(decile = 6)
get_ir_strat(decile = 7)
get_ir_strat(decile = 8)
get_ir_strat(decile = 9)
get_ir_strat(decile = 10)

ir_res_strat %<>%
  mutate(
    irr = ir_trt1/ir_trt0,
    ird = ir_trt1 - ir_trt0
  )

## e) Save results
ir_res_strat %<>%
  mutate_if(is.numeric, round, digits = 4)

coxph_res_strat %<>%
  mutate_if(is.numeric, round, digits = 4)

write_xlsx(coxph_res_strat, path = "C:/Users/gwen.aubrac/Desktop/trimming_res/coxph_res_strat.xlsx")
write_xlsx(ir_res_strat, path = "C:/Users/gwen.aubrac/Desktop/trimming_res/ir_res_strat.xlsx")

#### 3. TRIM SAMPLE ####

## a) Identify percentile cutoffs

perc_upper <- cohort %>%
  filter(new_trt == 0) %>%
  summarise(
    perc_99 = quantile(ps, 0.99, na.rm = TRUE),
    perc_97_5 = quantile(ps, 0.975, na.rm = TRUE),
    perc_95 = quantile(ps, 0.95, na.rm = TRUE),
    perc_90 = quantile(ps, 0.90, na.rm = TRUE)
  )

perc_lower <- cohort %>%
  filter(new_trt == 1) %>%
  summarise(
    perc_1 = quantile(ps, 0.01, na.rm = TRUE),
    perc_2_5 = quantile(ps, 0.025, na.rm = TRUE),
    perc_5 = quantile(ps, 0.05, na.rm = TRUE),
    perc_10 = quantile(ps, 0.1, na.rm = TRUE)
  )

summary(cohort$ps)

cohort %<>%
  mutate(
    perc_99 = round(perc_upper$perc_99, 4),
    perc_97_5 = round(perc_upper$perc_97_5, 4),
    perc_95 = round(perc_upper$perc_95, 4),
    perc_90 = round(perc_upper$perc_90, 4),
    perc_1 = round(perc_lower$perc_1, 4),
    perc_2_5 = round(perc_lower$perc_2_5, 4),
    perc_5 = round(perc_lower$perc_5, 4),
    perc_10 = round(perc_lower$perc_10, 4)
  )

trim_2 <- cohort %>% filter(ps < perc_99 & ps > perc_1)
trim_5 <- cohort %>% filter(ps < perc_97_5 & ps > perc_2_5)
trim_10 <- cohort %>% filter(ps < perc_95 & ps > perc_5)
trim_20 <- cohort %>% filter(ps < perc_90 & ps > perc_10)
trim_pf <- cohort %>% filter(pf > 0.3 & pf < 0.7)

## b) Visualizate PS in trimmed samples (without weighting)

# untrimmed
p <- ggplot(data=cohort, aes(x=ps, group=new_trt, fill=new_trt)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_manual(values=c("0"="blue", "1"="red"), 
                    labels=c("0"="SNRI", "1"="SSRI"),
                    name = "Treatment") +
  theme_minimal() + 
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(title="Propensity Score Distribution in Untrimmed Population", x = "propensity score") 

p

ggsave("C:/Users/gwen.aubrac/Desktop/trimming_res/plots/ps_untrimmed.pdf", 
       plot = p,  
       width = 8, 
       height = 6,
       units = "in")

# 2% trimmed
p <- ggplot(data=trim_2, aes(x=ps, group=new_trt, fill=new_trt)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_manual(values=c("0"="blue", "1"="red"), 
                    labels=c("0"="SNRI", "1"="SSRI"),
                    name = "Treatment") +
  theme_minimal() + 
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(title="Propensity Score Distribution after 2% Asymmetric Trimming", x = "propensity score") 

p

ggsave("C:/Users/gwen.aubrac/Desktop/trimming_res/plots/ps_trim2.pdf", 
       plot = p,  
       width = 8, 
       height = 6,
       units = "in")

# 5% trimmed
p <- ggplot(data=trim_5, aes(x=ps, group=new_trt, fill=new_trt)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_manual(values=c("0"="blue", "1"="red"), 
                    labels=c("0"="SNRI", "1"="SSRI"),
                    name = "Treatment") +
  theme_minimal() + 
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(title="Propensity Score Distribution after 5% Asymmetric Trimming", x = "propensity score") 

p

ggsave("C:/Users/gwen.aubrac/Desktop/trimming_res/plots/ps_trim5.pdf", 
       plot = p,  
       width = 8, 
       height = 6,
       units = "in")


# 10% trimmed
p <- ggplot(data=trim_10, aes(x=ps, group=new_trt, fill=new_trt)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_manual(values=c("0"="blue", "1"="red"), 
                    labels=c("0"="SNRI", "1"="SSRI"),
                    name = "Treatment") +
  theme_minimal() + 
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(title="Propensity Score Distribution after 10% Asymmetric Trimming", x = "propensity score") 

p

ggsave("C:/Users/gwen.aubrac/Desktop/trimming_res/plots/ps_trim10.pdf", 
       plot = p,  
       width = 8, 
       height = 6,
       units = "in")

# 20% trimmed
p <- ggplot(data=trim_20, aes(x=ps, group=new_trt, fill=new_trt)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_manual(values=c("0"="blue", "1"="red"), 
                    labels=c("0"="SNRI", "1"="SSRI"),
                    name = "Treatment") +
  theme_minimal() + 
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(title="Propensity Score Distribution after 20% Asymmetric Trimming", x = "propensity score") 

p

ggsave("C:/Users/gwen.aubrac/Desktop/trimming_res/plots/ps_trim20.pdf", 
       plot = p,  
       width = 8, 
       height = 6,
       units = "in")

# PF trimmed
p <- ggplot(data=trim_pf, aes(x=ps, group=new_trt, fill=new_trt)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_manual(values=c("0"="blue", "1"="red"), 
                    labels=c("0"="SNRI", "1"="SSRI"),
                    name = "Treatment") +
  theme_minimal() + 
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(title="Propensity Score Distribution after Preference Score Trimming", x = "propensity score") 

p

ggsave("C:/Users/gwen.aubrac/Desktop/trimming_res/plots/ps_trim_pf.pdf", 
       plot = p,  
       width = 8, 
       height = 6,
       units = "in")

#### 4. RE-ESTIMATE PS AND COMPUTE ATT/ATU/IPTW WEIGHTS ####

## a) Estimate weights

make_weights <- function(data) {
  new_ps_model <- glm(
    new_trt ~ age_group + sex + year + ethnicity + deprivation + anxiety_base + 
      arrhythmia_base + cardiomyopathy_base +cerebrovascular_disease_base +
      chronic_renal_disease_base + copd_base + depression_base + epilepsy_base +
      heart_failure_base + hepatic_disease_base + hyperlipidemia_base + hypertension_base +
      ischemic_heart_disease_base + lvh_base + myocardial_infarction_base +pacemaker_base +
      pvd_base + stroke_base + suicidal_ideation_self_harm_base + valvular_heart_disease_base,
    family = binomial, 
    data = data
  )
  
  new_data <- data
  new_data$new_ps <- predict(new_ps_model, type = "response")
  new_data$IPTW <- if_else(new_data$new_trt == 1, 1/new_data$new_ps, 1/(1-new_data$new_ps))
  new_data$ATTw <- if_else(new_data$new_trt == 1, 1, new_data$new_ps/(1-new_data$new_ps))
  new_data$ATUw <- if_else(new_data$new_trt == 1, (1-new_data$new_ps)/new_data$new_ps, 1)

  return(new_data)
}

cohort_weights <- make_weights(data = cohort)
trim_2_weights <- make_weights(data = trim_2)
trim_5_weights <- make_weights(data = trim_5)
trim_10_weights <- make_weights(data = trim_10)
trim_20_weights <- make_weights(data = trim_20)
trim_pf_weights <- make_weights(data = trim_pf)

## b) Visualize PS in trimmed samples (with weighting)

# untrimmed ATTw
p <- ggplot(data=cohort_weights, aes(x=ps, group=new_trt, fill=new_trt)) +
  geom_density(adjust=1.5, alpha=.4, (aes(weight = ATTw))) +
  scale_fill_manual(values=c("0"="blue", "1"="red"), 
                    labels=c("0"="SNRI", "1"="SSRI"),
                    name = "Treatment") +
  theme_minimal() + 
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(title="Propensity Score Distribution in Untrimmed Population (ATTw)", x = "propensity score") 

p

ggsave("C:/Users/gwen.aubrac/Desktop/trimming_res/plots/ps_untrimmed_ATTw.pdf", 
       plot = p,  
       width = 8, 
       height = 6,
       units = "in")

# 2% trimming ATTw
p <- ggplot(data=trim_2_weights, aes(x=ps, group=new_trt, fill=new_trt)) +
  geom_density(adjust=1.5, alpha=.4, (aes(weight = ATTw))) +
  scale_fill_manual(values=c("0"="blue", "1"="red"), 
                    labels=c("0"="SNRI", "1"="SSRI"),
                    name = "Treatment") +
  theme_minimal() + 
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(title="Propensity Score Distribution after 2% Asymmetric Trimming (ATTw)", x = "propensity score") 

p

ggsave("C:/Users/gwen.aubrac/Desktop/trimming_res/plots/ps_trim2_ATTw.pdf", 
       plot = p,  
       width = 8, 
       height = 6,
       units = "in")

# 5% trimming ATTw
p <- ggplot(data=trim_5_weights, aes(x=ps, group=new_trt, fill=new_trt)) +
  geom_density(adjust=1.5, alpha=.4, (aes(weight = ATTw))) +
  scale_fill_manual(values=c("0"="blue", "1"="red"), 
                    labels=c("0"="SNRI", "1"="SSRI"),
                    name = "Treatment") +
  theme_minimal() + 
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(title="Propensity Score Distribution after 5% Asymmetric Trimming (ATTw)", x = "propensity score") 

p

ggsave("C:/Users/gwen.aubrac/Desktop/trimming_res/plots/ps_trim5_ATTw.pdf", 
       plot = p,  
       width = 8, 
       height = 6,
       units = "in")

# 10% trimming ATTw
p <- ggplot(data=trim_10_weights, aes(x=ps, group=new_trt, fill=new_trt)) +
  geom_density(adjust=1.5, alpha=.4, (aes(weight = ATTw))) +
  scale_fill_manual(values=c("0"="blue", "1"="red"), 
                    labels=c("0"="SNRI", "1"="SSRI"),
                    name = "Treatment") +
  theme_minimal() + 
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(title="Propensity Score Distribution after 10% Asymmetric Trimming (ATTw)", x = "propensity score") 

p

ggsave("C:/Users/gwen.aubrac/Desktop/trimming_res/plots/ps_trim10_ATTw.pdf", 
       plot = p,  
       width = 8, 
       height = 6,
       units = "in")

# 20% trimming ATTw
p <- ggplot(data=trim_20_weights, aes(x=ps, group=new_trt, fill=new_trt)) +
  geom_density(adjust=1.5, alpha=.4, (aes(weight = ATTw))) +
  scale_fill_manual(values=c("0"="blue", "1"="red"), 
                    labels=c("0"="SNRI", "1"="SSRI"),
                    name = "Treatment") +
  theme_minimal() + 
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(title="Propensity Score Distribution after 20% Asymmetric Trimming (ATTw)", x = "propensity score") 

p

ggsave("C:/Users/gwen.aubrac/Desktop/trimming_res/plots/ps_trim20_ATTw.pdf", 
       plot = p,  
       width = 8, 
       height = 6,
       units = "in")

# pf trimming ATTw
p <- ggplot(data=trim_pf_weights, aes(x=ps, group=new_trt, fill=new_trt)) +
  geom_density(adjust=1.5, alpha=.4, (aes(weight = ATTw))) +
  scale_fill_manual(values=c("0"="blue", "1"="red"), 
                    labels=c("0"="SNRI", "1"="SSRI"),
                    name = "Treatment") +
  theme_minimal() + 
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(title="Propensity Score Distribution after Preference Score Trimming (ATTw)", x = "propensity score") 

p

ggsave("C:/Users/gwen.aubrac/Desktop/trimming_res/plots/ps_trim_pf_ATTw.pdf", 
       plot = p,  
       width = 8, 
       height = 6,
       units = "in")

#### 5. ESTIMATE HR/IRR/IRD IN TRIMMED SAMPLES ####

## a) Hazard ratios

coxph_res_trim <- data.frame(
  sample = rep(c("all", "trim_2", "trim_5", "trim_10", "trim_20", "trim_pf"), each = 3),
  weights = rep(c("IPTW", "ATTw", "ATUw"), times = 3),
  coef = NA,
  hr = NA,
  hr_lower = NA,
  hr_upper = NA
)

fit_trim_coxph <- function(sample_data, sample_name) {
  
  # IPTW
  coxph_IPTW <- coxph(Surv(follow_up, death) ~ new_trt,
                   data = sample_data,
                   weights = sample_data$IPTW)
  
  coxph_res_trim$coef[coxph_res_trim$sample == sample_name & coxph_res_trim$weights == "IPTW"] <- summary(coxph_IPTW)$coefficients[, "coef"]
  coxph_res_trim$hr[coxph_res_trim$sample == sample_name & coxph_res_trim$weights == "IPTW"] <- exp(summary(coxph_IPTW)$coefficients[, "coef"])
  coxph_res_trim$hr_lower[coxph_res_trim$sample == sample_name & coxph_res_trim$weights == "IPTW"] <- summary(coxph_IPTW)$conf.int[, "lower .95"]
  coxph_res_trim$hr_upper[coxph_res_trim$sample == sample_name & coxph_res_trim$weights == "IPTW"] <- summary(coxph_IPTW)$conf.int[, "upper .95"]
  
  # ATTw
  coxph_ATTw <- coxph(Surv(follow_up, death) ~ new_trt,
                      data = sample_data,
                        weights = sample_data$ATTw)
  
  coxph_res_trim$coef[coxph_res_trim$sample == sample_name & coxph_res_trim$weights == "ATTw"] <- summary(coxph_ATTw)$coefficients[, "coef"]
  coxph_res_trim$hr[coxph_res_trim$sample == sample_name & coxph_res_trim$weights == "ATTw"] <- exp(summary(coxph_ATTw)$coefficients[, "coef"])
  coxph_res_trim$hr_lower[coxph_res_trim$sample == sample_name & coxph_res_trim$weights == "ATTw"] <- summary(coxph_ATTw)$conf.int[, "lower .95"]
  coxph_res_trim$hr_upper[coxph_res_trim$sample == sample_name & coxph_res_trim$weights == "ATTw"] <- summary(coxph_ATTw)$conf.int[, "upper .95"]
  
  # ATUw
  coxph_ATUw <- coxph(Surv(follow_up, death) ~ new_trt,
                      data = sample_data,
                      weights = sample_data$ATUw)
  
  coxph_res_trim$coef[coxph_res_trim$sample == sample_name & coxph_res_trim$weights == "ATUw"] <- summary(coxph_ATUw)$coefficients[, "coef"]
  coxph_res_trim$hr[coxph_res_trim$sample == sample_name & coxph_res_trim$weights == "ATUw"] <- exp(summary(coxph_ATUw)$coefficients[, "coef"])
  coxph_res_trim$hr_lower[coxph_res_trim$sample == sample_name & coxph_res_trim$weights == "ATUw"] <- summary(coxph_ATUw)$conf.int[, "lower .95"]
  coxph_res_trim$hr_upper[coxph_res_trim$sample == sample_name & coxph_res_trim$weights == "ATUw"] <- summary(coxph_ATUw)$conf.int[, "upper .95"]
  
  coxph_res_trim <<- coxph_res_trim
  
}

fit_trim_coxph(sample_data = cohort_weights, sample_name = "all")
fit_trim_coxph(sample_data = trim_2_weights, sample_name = "trim_2")
fit_trim_coxph(sample_data = trim_5_weights, sample_name = "trim_5")
fit_trim_coxph(sample_data = trim_10_weights, sample_name = "trim_10")
fit_trim_coxph(sample_data = trim_20_weights, sample_name = "trim_20")
fit_trim_coxph(sample_data = trim_pf_weights, sample_name = "trim_pf")

## b) Incidence rates

# Get weighted person-time and weighted events
# in each trimmed sample
cohort_weights %<>%
  mutate(pt_IPTW = follow_up * IPTW,
         pt_ATTw = follow_up * ATTw, 
         pt_ATUw = follow_up * ATUw,
         death_IPTW = death * IPTW,
         death_ATTw = death * ATTw,
         death_ATUw = death * ATUw)

trim_2_weights %<>%
  mutate(pt_IPTW = follow_up * IPTW,
         pt_ATTw = follow_up * ATTw, 
         pt_ATUw = follow_up * ATUw,
         death_IPTW = death * IPTW,
         death_ATTw = death * ATTw,
         death_ATUw = death * ATUw)

trim_5_weights %<>%
  mutate(pt_IPTW = follow_up * IPTW,
         pt_ATTw = follow_up * ATTw, 
         pt_ATUw = follow_up * ATUw,
         death_IPTW = death * IPTW,
         death_ATTw = death * ATTw,
         death_ATUw = death * ATUw)

trim_10_weights %<>%
  mutate(pt_IPTW = follow_up * IPTW,
         pt_ATTw = follow_up * ATTw, 
         pt_ATUw = follow_up * ATUw,
         death_IPTW = death * IPTW,
         death_ATTw = death * ATTw,
         death_ATUw = death * ATUw)

trim_20_weights %<>%
  mutate(pt_IPTW = follow_up * IPTW,
         pt_ATTw = follow_up * ATTw, 
         pt_ATUw = follow_up * ATUw,
         death_IPTW = death * IPTW,
         death_ATTw = death * ATTw,
         death_ATUw = death * ATUw)

trim_pf_weights %<>%
  mutate(pt_IPTW = follow_up * IPTW,
         pt_ATTw = follow_up * ATTw, 
         pt_ATUw = follow_up * ATUw,
         death_IPTW = death * IPTW,
         death_ATTw = death * ATTw,
         death_ATUw = death * ATUw)

# Separately calculate IR per 100 in treated and untreated

ir_res_trim <- data.frame(
  sample = rep(c("all", "trim_2", "trim_5", "trim_10", "trim_20", "trim_pf"), each = 3),
  weights = rep(c("IPTW", "ATTw", "ATUw"), times = 3),
  ir_trt1 = NA,
  ir_trt0 = NA,
  irr = NA,
  ird = NA
)

get_ir_trim <- function(sample_data, sample_name) {
  
  # for treated
  ir_all_trt1 <- sample_data %>% 
    filter(new_trt == 1) %>% 
    summarize(
      sum_pt_IPTW = sum(pt_IPTW),
      sum_pt_ATTw = sum(pt_ATTw),
      sum_pt_ATUw = sum(pt_ATUw),
      sum_death_IPTW = sum(death_IPTW),
      sum_death_ATTw = sum(death_ATTw),
      sum_death_ATUw = sum(death_ATUw)
    )
  
  ir_res_trim$ir_trt1[ir_res_trim$sample == sample_name & ir_res_trim$weights == "IPTW"] <- ir_all_trt1$sum_death_IPTW/ir_all_trt1$sum_pt_IPTW*100
  ir_res_trim$ir_trt1[ir_res_trim$sample == sample_name & ir_res_trim$weights == "ATTw"] <- ir_all_trt1$sum_death_ATTw/ir_all_trt1$sum_pt_ATTw*100
  ir_res_trim$ir_trt1[ir_res_trim$sample == sample_name & ir_res_trim$weights == "ATUw"] <- ir_all_trt1$sum_death_ATUw/ir_all_trt1$sum_pt_ATUw*100
  
  # for untreated
  ir_all_trt0 <- sample_data %>% 
    filter(new_trt == 0) %>% 
    summarize(
      sum_pt_IPTW = sum(pt_IPTW),
      sum_pt_ATTw = sum(pt_ATTw),
      sum_pt_ATUw = sum(pt_ATUw),
      sum_death_IPTW = sum(death_IPTW),
      sum_death_ATTw = sum(death_ATTw),
      sum_death_ATUw = sum(death_ATUw)
    )
  
  ir_res_trim$ir_trt0[ir_res_trim$sample == sample_name & ir_res_trim$weights == "IPTW"] <- ir_all_trt0$sum_death_IPTW/ir_all_trt0$sum_pt_IPTW*100
  ir_res_trim$ir_trt0[ir_res_trim$sample == sample_name & ir_res_trim$weights == "ATTw"] <- ir_all_trt0$sum_death_ATTw/ir_all_trt0$sum_pt_ATTw*100
  ir_res_trim$ir_trt0[ir_res_trim$sample == sample_name & ir_res_trim$weights == "ATUw"] <- ir_all_trt0$sum_death_ATUw/ir_all_trt0$sum_pt_ATUw*100

  ir_res_trim <<- ir_res_trim
}


get_ir_trim(sample_data = cohort_weights, sample_name = "all")
get_ir_trim(sample_data = trim_2_weights, sample_name = "trim_2")
get_ir_trim(sample_data = trim_5_weights, sample_name = "trim_5")
get_ir_trim(sample_data = trim_10_weights, sample_name = "trim_10")
get_ir_trim(sample_data = trim_20_weights, sample_name = "trim_20")
get_ir_trim(sample_data = trim_pf_weights, sample_name = "trim_pf")

ir_res_trim %<>%
  mutate(
    irr = ir_trt1/ir_trt0,
    ird = ir_trt1 - ir_trt0
  )

## c) Save results
ir_res_trim %<>%
  mutate_if(is.numeric, round, digits = 4)

coxph_res_trim %<>%
  mutate_if(is.numeric, round, digits = 4)

write_xlsx(ir_res_trim, path = "C:/Users/gwen.aubrac/Desktop/trimming_res/ir_res_trim.xlsx")
write_xlsx(coxph_res_trim, path = "C:/Users/gwen.aubrac/Desktop/trimming_res/coxph_res_trim.xlsx")
