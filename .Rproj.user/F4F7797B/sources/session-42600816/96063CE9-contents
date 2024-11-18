#### Analysis script for Flutamide study: REP (REPEAT) calls ###################
###### Bayesian Multilevel models ##################################
############## BWalkenhorst 2024 ###############################################

#### SETUP ####

# clear everything and use garbage collector at start
rm(list = ls(all.names = TRUE)) #includes hidden objects.
gc() 

# load all necessary libraries
library(ggplot2) 
library(tidyverse)
library(tidybayes) 
library(brms)
library(bayestestR) #e.g. diagnostic_posterior
library(bayesplot)
library(ggokabeito) # colour palette
library(emmeans) # emtrends
library(extrafont)# use font_import() on first use

set.seed(23)

REP_data <- readRDS('data/REP_data.rds')

#Generic weakly informative prior: normal(0, 1);
priors <- c(set_prior("normal(0,1)", class = "Intercept"), set_prior("normal(0,1)", class='b'))

# Custom ggplot theme 
theme_clean <- function() {
  theme_minimal(base_family='Calibri') +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = rel(2), hjust = 0.5),
      axis.title = element_text(face = "bold", size = rel(2)),
      axis.text = element_text(face = "bold", size= rel(1.25)),
      strip.text = element_text(face = "bold", size = rel(2), color='white'),
      strip.background = element_rect(fill = "grey80", color = NA),
      legend.title = element_text(face = "bold", size = rel(2)),
      legend.text = element_text(face = 'italic', size = rel(1.5)))
}

#### FUNCTIONS ####
run_summary_on_models <- function(model_list, output_file) {
  sink(output_file)
  for (i in seq_along(model_list)) {
    cat("Summary for model:", i, "\n\n")
    model_summary <- capture.output(summary(model_list[[i]]))
    cat(model_summary, sep = "\n")
    cat("\n\n")
  }
  sink()
}

get_age_vars <- function(){
  age_vars <- c((30- mean(REP_data$REC_AGE_D))/sd(REP_data$REC_AGE_D), 
                (45- mean(REP_data$REC_AGE_D))/sd(REP_data$REC_AGE_D), 
                (60- mean(REP_data$REC_AGE_D))/sd(REP_data$REC_AGE_D), 
                (75- mean(REP_data$REC_AGE_D))/sd(REP_data$REC_AGE_D), 
                (90- mean(REP_data$REC_AGE_D))/sd(REP_data$REC_AGE_D),
                (105- mean(REP_data$REC_AGE_D))/sd(REP_data$REC_AGE_D), 
                (120- mean(REP_data$REC_AGE_D))/sd(REP_data$REC_AGE_D))
  return(age_vars)
}

################################################################################
######################## Call length ###########################################
################################################################################
# 1) Check for non-linearity ####
B_REP_len_A <- brms::brm(formula = BEG_avg_Len ~ AGE_z + (1|LITTER_CODE/ID),
                          data = REP_data, family = lognormal(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                          prior = priors, threads = threading(4), 
                          file="B_REP_len_A")
B_REP_len_A <- add_criterion(B_REP_len_A, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_REP_len_A)
# plot(B_REP_len_A)
# pp_check(B_REP_len_A, ndraws=100)
B_REP_len_W <- brms::brm(formula = BEG_avg_Len ~ WEIGHT_z + (1|LITTER_CODE/ID), 
                            data = REP_data, family = lognormal(link='identity'),
                            chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                            save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                            prior = priors, threads = threading(4),
                            file="B_REP_len_W")
B_REP_len_W <- add_criterion(B_REP_len_W, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_REP_len_W)
# plot(B_REP_len_W)
# pp_check(B_REP_len_W, ndraws=100)
B_REP_len_CN <- brms::brm(formula = BEG_avg_Len ~ COMP_NORM_z + (1|LITTER_CODE/ID) ,
                            data = REP_data, family = lognormal(link='identity'),
                            chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                            save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                            prior = priors, threads = threading(4),
                            file="B_REP_len_CN")
B_REP_len_CN <- add_criterion(B_REP_len_CN, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_REP_len_CN)
# plot(B_REP_len_CN)
# pp_check(B_REP_len_CN, ndraws=100)
B_REP_len_G <- brms::brm(formula = BEG_avg_Len ~ GS_z +(1|LITTER_CODE/ID),
                            data = REP_data, family = lognormal(link='identity'),
                            chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                            save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                            prior = priors, threads = threading(4),
                            file="B_REP_len_G")
B_REP_len_G <- add_criterion(B_REP_len_G, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_REP_len_G)
# plot(B_REP_len_G)
# pp_check(B_REP_len_G, ndraws=100)
B_REP_len_R <- brms::brm(formula = BEG_avg_Len ~ RAIN_z + (1|LITTER_CODE/ID),
                           data = REP_data, family = lognormal(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                           save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                           prior = priors, threads = threading(4),
                           file="B_REP_len_R")
B_REP_len_R <- add_criterion(B_REP_len_R, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_REP_len_R)
# plot(B_REP_len_R)
# pp_check(B_REP_len_R, ndraws=100)

# non linear
B_REP_len_A2 <- brms::brm(formula = BEG_avg_Len ~  AGE_z + I(AGE_z^2) + (1|LITTER_CODE/ID),
                          data = REP_data, family = lognormal(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                          prior = priors, threads = threading(4),
                          file="B_REP_len_A2")
B_REP_len_A2 <- add_criterion(B_REP_len_A2, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_REP_len_A2)
# plot(B_REP_len_A2)
# pp_check(B_REP_len_A2, ndraws=100)
B_REP_len_W2 <- brms::brm(formula = BEG_avg_Len ~ WEIGHT_z + I(WEIGHT_z^2) + (1|LITTER_CODE/ID),
                             data = REP_data, family = lognormal(link='identity'),
                             chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                             save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                             prior = priors, threads = threading(4),
                             file="B_REP_len_W2")
B_REP_len_W2 <- add_criterion(B_REP_len_W2, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_REP_len_W2)
# plot(B_REP_len_W2)
# pp_check(B_REP_len_W2, ndraws=100)
B_REP_len_CN2 <- brms::brm(formula = BEG_avg_Len ~ COMP_NORM_z + I(COMP_NORM_z^2) + (1|LITTER_CODE/ID),
                             data = REP_data, family = lognormal(link='identity'),
                             chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                             save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                             prior = priors, threads = threading(4),
                             file="B_REP_len_CN2")
B_REP_len_CN2 <- add_criterion(B_REP_len_CN2, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_REP_len_CN2)
# plot(B_REP_len_CN2)
# pp_check(B_REP_len_CN2, ndraws=100)
B_REP_len_G2 <- brms::brm(formula = BEG_avg_Len ~ GS_z +I(GS_z^2) + (1|LITTER_CODE/ID),
                             data = REP_data, family = lognormal(link='identity'),
                             chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                             save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                             prior = priors, threads = threading(4),
                             file="B_REP_len_G2")
B_REP_len_G2 <- add_criterion(B_REP_len_G2, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_REP_len_G2)
# plot(B_REP_len_G2)
# pp_check(B_REP_len_G2, ndraws=100)
B_REP_len_R2 <- brms::brm(formula = BEG_avg_Len ~ RAIN_z + I(RAIN_z^2) + (1|LITTER_CODE/ID),
                            data = REP_data, family = lognormal(link='identity'),
                            chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                            save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                            prior = priors, threads = threading(4),
                            file="B_REP_len_R2")
B_REP_len_R2 <- add_criterion(B_REP_len_R2, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_REP_len_R2)
# plot(B_REP_len_R2)
# pp_check(B_REP_len_R2, ndraws=100)

loo(B_REP_len_A, B_REP_len_W, B_REP_len_W2, B_REP_len_CN, 
    B_REP_len_CN2, B_REP_len_A2, B_REP_len_G, B_REP_len_G2, 
    B_REP_len_R, B_REP_len_R2)

# best models
loo(B_REP_len_CN, B_REP_len_A2, B_REP_len_W, 
    B_REP_len_G, B_REP_len_R)

loo_R2(B_REP_len_CN)
loo_R2(B_REP_len_A2)
Loo_R2(B_REP_len_W)
loo_R2(B_REP_len_G)
loo_R2(B_REP_len_R)

# save all summaries of the models:
model_list <-
  list(
  B_REP_len_W, B_REP_len_W2, B_REP_len_CN, 
  B_REP_len_CN2, B_REP_len_A2, B_REP_len_G, B_REP_len_G2, 
  B_REP_len_R, B_REP_len_R2, B_REP_len_A
)
run_summary_on_models(model_list, "REP_len_sum.txt")

# save only relevant models
model_list <- list(
  B_REP_len_CN, B_REP_len_A2, B_REP_len_W, 
  B_REP_len_G, B_REP_len_R)

run_summary_on_models(model_list, "REP_len_best_sum.txt")

rm( B_REP_len_W, B_REP_len_W2, B_REP_len_CN, 
   B_REP_len_CN2, B_REP_len_A, B_REP_len_G, B_REP_len_G2, 
   B_REP_len_R, B_REP_len_R2, model_list, B_REP_len_A2) 

# 2) Define model for REP LENGTH ####
B_REP_len <- brms::brm(formula = BEG_avg_Len ~ TREATMENT + AGE_z + I(AGE_z^2) + SEX + WEIGHT_z + COMP_NORM_z + GS_z + RAIN_z +
                               TREATMENT:AGE_z + TREATMENT:I(AGE_z^2) + TREATMENT:SEX + AGE_z:SEX + I(AGE_z^2):SEX +
                               TREATMENT:AGE_z:SEX + TREATMENT:I(AGE_z^2):SEX + (1|LITTER_CODE/ID),
                             data = REP_data, family = lognormal(link='identity'),
                             chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                             save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                             prior = priors, threads = threading(4),
                             file="B_REP_len")
B_REP_len <- add_criterion(B_REP_len, c("loo", "loo_R2"), moment_match = TRUE)

B_REP_len <- readRDS('models/B_REP_len.rds')

#### RESULTS: REP LENGTH ####
summary(B_REP_len)

plot(B_REP_len)
pp_check(B_REP_len, ndraws = 100)

describe_posterior(
  B_REP_len,
  effects = "all", #fixed vs all (for random effects)
  component = "all",
  rope_range = rope_range(B_REP_len),  
  test = c("p_direction", "p_significance", "rope"),
  centrality = "all",
  dispersion = TRUE
)

loo_R2(B_REP_len) 
bayes_R2(B_REP_len)

performance::variance_decomposition(B_REP_len)


#### EMMs at different ages (given non-linearity) ####
(emm_REP_LEN <- emmeans(B_REP_len, ~ TREATMENT:AGE_z, at=list(AGE_z = get_age_vars())))

pd(emm_REP_LEN) # all 100%
equivalence_test(emm_REP_LEN, range = rope_range(B_REP_len)) 
p_significance(emm_REP_LEN, threshold = rope_range(B_REP_len))

# contrasts BETWEEN treatments
(REP_LEN_contrasts <-pairs(emm_REP_LEN, by = c("AGE_z")))

pd(REP_LEN_contrasts)

equivalence_test(REP_LEN_contrasts, range = rope_range(B_REP_len))

p_significance(REP_LEN_contrasts, threshold = rope_range(B_REP_len))

# Plot emmeans: REP length ####
emm_data <- as.data.frame(emm_REP_LEN)
emm_data$REC_AGE_D <- emm_data$AGE_z * sd(REP_data$REC_AGE_D) + mean(REP_data$REC_AGE_D)
emm_data$REP_len <- exp(emm_data$emmean)
emm_data$HPD_low <- exp(emm_data$lower.HPD)
emm_data$HPD_high <- exp(emm_data$upper.HPD)

emm_data$TREATMENT <- factor(emm_data$TREATMENT, levels = c("CTRL", "SUB", "FLUT"))

# 800* 500: EMMs_REP_LEN
ggplot(emm_data, aes(x = as.factor(REC_AGE_D), y = REP_len, color = TREATMENT)) +
  geom_point(position = position_dodge(width = 0.9), size = 3) +  # Show estimates as points
  geom_errorbar(aes(ymin = HPD_low, ymax = HPD_high), position = position_dodge(width = 0.9), width = 0.75) +
  scale_color_okabe_ito(order= c(2, 1, 3), name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  scale_y_continuous(n.breaks=10)+
  labs(x = "Age (days)", y = "Repeat call length (s)\n") +
  theme_clean()

# 800* 500: EMMs_REP_LEN_ribbon
ggplot(emm_data, aes(x = as.factor(REC_AGE_D), y = REP_len, color = TREATMENT, group = TREATMENT)) +
  geom_line(position = position_dodge(width = 0.9), linewidth = 1) +  # Connect estimates as lines
  geom_ribbon(aes(ymin = HPD_low, ymax = HPD_high, fill = TREATMENT), position = position_dodge(width = 0.9), alpha = 0.3) +  # Ribbon for credible intervals
  scale_color_okabe_ito(order= c(2, 1, 3), name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order= c(2, 1, 3), name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  scale_y_continuous(n.breaks=10) +
  labs(x = "Age (days)", y = "Repeat call length (s)\n") +
  theme_clean()

rm(emm_REP_LEN, emm_data, REP_LEN_contrasts)

### PLOTS: REP LENGTH ####
# coefficients:
posterior_desc <- describe_posterior(
  B_REP_len,
  effects = "fixed",
  component = "all",
  rope_range = rope_range(B_REP_len),
  centrality = "all",
  dispersion = TRUE
)
# drop sigma row
posterior_desc <- posterior_desc[-c(1,23),]

# clean up labels:
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'b_', '')
posterior_desc$TREATMENT <- ifelse(grepl("TREATMENTSUB", posterior_desc$Parameter), "SC",
                             ifelse(grepl("TREATMENTFLUT", posterior_desc$Parameter), "DT", "DC"))
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTFLUT', 'DT')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTSUB', 'SC')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'IAGE_zE2', 'Age^2')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'AGE_z', 'Age')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'SEXM', 'Male')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'WEIGHT_z', 'Weight offset')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'COMP_NORM_z', 'Competition')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'GS_z', 'Group size')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'RAIN_z', 'Monthly rainfall')

custom_order <- c('Monthly rainfall',
                  'Group size',
                  'Competition', 
                  'Weight offset', 
                  'DT:Age^2:Male','SC:Age^2:Male','Age^2:Male',
                  'DT:Age:Male','SC:Age:Male','Age:Male', 
                  'DT:Male','SC:Male','Male', 
                  'DT:Age^2','SC:Age^2',"Age^2", 
                  'DT:Age','SC:Age',"Age", 
                  'DT','SC')#,'Intercept') 

posterior_desc$Parameter <- factor(posterior_desc$Parameter, levels = custom_order)
posterior_desc$TREATMENT <- factor(posterior_desc$TREATMENT, levels = c("DC", "SC", "DT"))
# Coeff_REP_LEN 700*800
ggplot(posterior_desc, aes(y = Parameter, x = Median, xmin = CI_low, xmax = CI_high, color=TREATMENT)) +
  geom_vline(xintercept = 0, color='grey', linetype = 'dotted', linewidth =1)+
  geom_point() +
  geom_errorbarh(height = 0) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT'))+
  labs(x = "Parameter value", y = "Parameters") +
 # ggtitle('REP call length')+
  theme_clean()+
  theme(legend.position="none")

# Ontogeny 30 - 130 plot ####
# get all needed values
(sd_age <- sd(REP_data$REC_AGE_D))#
(mean_age <- mean(REP_data$REC_AGE_D))#

range(REP_data$REC_AGE_D)# 31, 130
rec_age_c <- seq(30, 130, by=1)
age_z_vals <- (rec_age_c - mean(REP_data$REC_AGE_D))/sd(REP_data$REC_AGE_D)
age_z_vals <- seq(min(age_z_vals), max(age_z_vals), length.out = 20)

REP_LEN_pred <- B_REP_len %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(REP_data$TREATMENT),
                                    SEX = levels(REP_data$SEX),
                                    AGE_z = age_z_vals,
                                    WEIGHT_z = mean(REP_data$WEIGHT_z),
                                    COMP_NORM_z = mean(REP_data$COMP_NORM_z),
                                    GS_z = mean(REP_data$GS_z),
                                    RAIN_z = mean(REP_data$RAIN_z)),
               re_formula = NA,  robust=T) 

#unscale AGE_z values:
REP_LEN_pred$REC_AGE_D <- REP_LEN_pred$AGE_z * sd_age + mean_age
# ensure right format
REP_LEN_pred$REP_len <- REP_LEN_pred$.epred
REP_LEN_pred$TREATMENT <- factor(REP_LEN_pred$TREATMENT, levels = c("CTRL", "SUB", "FLUT"))
REP_data$TREATMENT <- factor(REP_data$TREATMENT, levels = c("CTRL", "SUB", "FLUT"))

#800*500: REP_LEN_TA
ggplot(REP_LEN_pred, aes(x = REC_AGE_D, y = REP_len, color = TREATMENT, fill = TREATMENT)) +  
  geom_point(data=REP_data, aes(y=BEG_avg_Len))+
  stat_lineribbon(.width = .95) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order = c(2, 1, 3), alpha = 0.3, name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)", y = "Repeat call length (s)\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Repeat call length (s) \n", n.breaks = 10) +
  theme_clean()

rm(REP_LEN_pred,  B_REP_len)

################################################################################
######################## Call interval length ##################################
################################################################################

#### BAYES MODELS ##############################################################
# 1) Check for non-linearity ####
B_REP_int_A <- brms::brm(formula = BEG_avg_Int ~ AGE_z + (1|LITTER_CODE/ID),
                          data = REP_data, family = lognormal(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                          prior = priors, threads = threading(4),
                          file="B_REP_int_A")
B_REP_int_A <- add_criterion(B_REP_int_A, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_REP_int_A)
# plot(B_REP_int_A)
# pp_check(B_REP_int_A, ndraws=100)
B_REP_int_W <- brms::brm(formula = BEG_avg_Int ~ WEIGHT_z + (1|LITTER_CODE/ID), 
                           data = REP_data, family = lognormal(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                           save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                           prior = priors, threads = threading(4),
                           file="B_REP_int_W")
B_REP_int_W <- add_criterion(B_REP_int_W, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_REP_int_W)
# plot(B_REP_int_W)
# pp_check(B_REP_int_W, ndraws=100)
B_REP_int_CN <- brms::brm(formula = BEG_avg_Int ~ COMP_NORM_z + (1|LITTER_CODE/ID),
                           data = REP_data, family = lognormal(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                           save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                           prior = priors, threads = threading(4),
                           file="B_REP_int_CN")
B_REP_int_CN <- add_criterion(B_REP_int_CN, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_REP_int_CN)
# plot(B_REP_int_CN)
# pp_check(B_REP_int_CN, ndraws=100)
B_REP_int_G <- brms::brm(formula = BEG_avg_Int ~ GS_z + (1|LITTER_CODE/ID),
                            data = REP_data, family = lognormal(link='identity'),
                            chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                            save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                            prior = priors, threads = threading(4),
                            file="B_REP_int_G")
B_REP_int_G <- add_criterion(B_REP_int_G, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_REP_int_G)
# plot(B_REP_int_G)
# pp_check(B_REP_int_G, ndraws=100)
B_REP_int_R <- brms::brm(formula = BEG_avg_Int ~ RAIN_z + (1|LITTER_CODE/ID),
                           data = REP_data, family = lognormal(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                           save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                           prior = priors, threads = threading(4),
                           file="B_REP_int_R")
B_REP_int_R <- add_criterion(B_REP_int_R, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_REP_int_R)
# plot(B_REP_int_R)
# pp_check(B_REP_int_R, ndraws=100)

# non linear
B_REP_int_A2 <- brms::brm(formula = BEG_avg_Int ~  AGE_z + I(AGE_z^2) + (1|LITTER_CODE/ID), 
                           data = REP_data, family = lognormal(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                           save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                           prior = priors, threads = threading(4),
                           file="B_REP_int_A2")
B_REP_int_A2 <- add_criterion(B_REP_int_A2, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_REP_int_A2)
# plot(B_REP_int_A2)
# pp_check(B_REP_int_A2, ndraws=100)
B_REP_int_W2 <- brms::brm(formula = BEG_avg_Int ~ WEIGHT_z + I(WEIGHT_z^2) + (1|LITTER_CODE/ID),
                            data = REP_data, family = lognormal(link='identity'),
                            chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                            save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                            prior = priors, threads = threading(4),
                            file="B_REP_int_W2")
B_REP_int_W2 <- add_criterion(B_REP_int_W2, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_REP_int_W2)
# plot(B_REP_int_W2)
# pp_check(B_REP_int_W2, ndraws=100)
B_REP_int_CN2 <- brms::brm(formula = BEG_avg_Int ~ COMP_NORM_z + I(COMP_NORM_z^2) + (1|LITTER_CODE/ID),
                            data = REP_data, family = lognormal(link='identity'),
                            chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                            save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                            prior = priors, threads = threading(4),
                            file="B_REP_int_CN2")
B_REP_int_CN2 <- add_criterion(B_REP_int_CN2, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_REP_int_CN2)
# plot(B_REP_int_CN2)
# pp_check(B_REP_int_CN2, ndraws=100)
B_REP_int_G2 <- brms::brm(formula = BEG_avg_Int ~ GS_z + I(GS_z^2) + (1|LITTER_CODE/ID),
                             data = REP_data, family = lognormal(link='identity'),
                             chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                             save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                             prior = priors, threads = threading(4),
                             file="B_REP_int_G2")
B_REP_int_G2 <- add_criterion(B_REP_int_G2, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_REP_int_G2)
# plot(B_REP_int_G2)
# pp_check(B_REP_int_G2, ndraws=100)
B_REP_int_R2 <- brms::brm(formula = BEG_avg_Int ~ RAIN_z + I(RAIN_z^2) +(1|LITTER_CODE/ID),
                            data = REP_data, family = lognormal(link='identity'),
                            chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                            save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                            prior = priors, threads = threading(4),
                            file="B_REP_int_R2")
B_REP_int_R2 <- add_criterion(B_REP_int_R2, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_REP_int_R2)
# plot(B_REP_int_R2)
# pp_check(B_REP_int_R2, ndraws=100)

loo(B_REP_int_W, B_REP_int_W2, B_REP_int_CN, B_REP_int_A,
    B_REP_int_CN2, B_REP_int_A2, B_REP_int_G, B_REP_int_G2,
    B_REP_int_R, B_REP_int_R2)

# best models
loo(B_REP_int_CN, B_REP_int_W, B_REP_int_A2,  
    B_REP_int_G2, B_REP_int_R)

loo_R2(B_REP_int_CN)
loo_R2(B_REP_int_A2)
loo_R2(B_REP_int_W)
loo_R2(B_REP_int_G2)
loo_R2(B_REP_int_R)

# save all summaries of the models:
model_list <- list(
  B_REP_int_W, B_REP_int_CN, B_REP_int_W2, 
 B_REP_int_CN2,  B_REP_int_A2,B_REP_int_G, B_REP_int_G2,
 B_REP_int_R, B_REP_int_R2
)
run_summary_on_models(model_list, "REP_int_sum.txt")

# save only relevant models: # use 
model_list <- list(
  B_REP_int_CN, B_REP_int_A2, B_REP_int_W,  
  B_REP_int_G2, B_REP_int_R)

run_summary_on_models(model_list, "REP_int_best_sum.txt")

rm(B_REP_int_W, B_REP_int_CN, B_REP_int_W2, B_REP_int_A,
     B_REP_int_CN2,  B_REP_int_A2,B_REP_int_G, B_REP_int_G2,
     B_REP_int_R, B_REP_int_R2, model_list)

# 2) Define models for REP interval length ####
B_REP_int<- brms::brm(formula = BEG_avg_Int ~ TREATMENT + AGE_z + I(AGE_z^2) + SEX + WEIGHT_z + COMP_NORM_z + GS_z + I(GS_z^2) + RAIN_z + 
                        TREATMENT:AGE_z + TREATMENT:I(AGE_z^2) +  TREATMENT:SEX + AGE_z:SEX + I(AGE_z^2):SEX +
                        TREATMENT:AGE_z:SEX + TREATMENT:I(AGE_z^2):SEX + (1|LITTER_CODE/ID),
                             data = REP_data, family = lognormal(link='identity'),
                             chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                             save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                             prior = priors, threads = threading(4),
                             file="B_REP_int")
B_REP_int <- add_criterion(B_REP_int, c("loo", "loo_R2"), moment_match = TRUE)

B_REP_int <- readRDS('models/B_REP_int.rds')

#### RESULTS: REP interval length ####
summary(B_REP_int)

plot(B_REP_int)
pp_check(B_REP_int, ndraws = 100)

describe_posterior(
  B_REP_int,
  effects = "all", #fixed vs all (for random effects)
  component = "all",
  rope_range = rope_range(B_REP_int),  
  test = c("p_direction", "pnificance", "rope"),
  centrality = "all",
  dispersion = TRUE
)

loo_R2(B_REP_int)
bayes_R2(B_REP_int)

performance::variance_decomposition(B_REP_int)

#### EMMs at different ages (given non-linearity) ####
(emm_REP_INT <- emmeans(B_REP_int, ~ TREATMENT:AGE_z, at=list(AGE_z = get_age_vars())))

pd(emm_REP_INT) 
equivalence_test(emm_REP_INT, range = rope_range(B_REP_int)) 
p_significance(emm_REP_INT, threshold = rope_range(B_REP_int))
# contrasts BETWEEN treatments
(REP_INT_contrasts <-pairs(emm_REP_INT, by = c("AGE_z")))

pd(REP_INT_contrasts)

equivalence_test(REP_INT_contrasts, range = rope_range(B_REP_int))

p_significance(REP_INT_contrasts, threshold = rope_range(B_REP_int))

# Plot emmeans: REP Interval duration ####
emm_data <- as.data.frame(emm_REP_INT)
emm_data$REC_AGE_D <- emm_data$AGE_z * sd(REP_data$REC_AGE_D) + mean(REP_data$REC_AGE_D)
emm_data$REP_INT <- exp(emm_data$emmean)
emm_data$HPD_low <- exp(emm_data$lower.HPD)
emm_data$HPD_high <- exp(emm_data$upper.HPD)

emm_data$TREATMENT <- factor(emm_data$TREATMENT, levels = c("CTRL", "SUB", "FLUT"))

# 700* 500: EMMs_REP_INT
ggplot(emm_data, aes(x = as.factor(REC_AGE_D), y = REP_INT, color = TREATMENT)) +
  geom_point(position = position_dodge(width = 0.9), size = 3) +  # Show estimates as points
  geom_errorbar(aes(ymin = HPD_low, ymax = HPD_high), position = position_dodge(width = 0.9), width = 0.75) +
  scale_color_okabe_ito(order= c(2, 1, 3), name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  scale_y_continuous(n.breaks=10)+
  labs(x = "Age (days)", y = "Repeat call interval duration (s)\n") +
  theme_clean()

# 800* 500: EMMs_REP_INT_ribbon
ggplot(emm_data, aes(x = as.factor(REC_AGE_D), y = REP_INT, color = TREATMENT, group = TREATMENT)) +
  geom_line(position = position_dodge(width = 0.9), linewidth = 1) +  # Connect estimates as lines
  geom_ribbon(aes(ymin = HPD_low, ymax = HPD_high, fill = TREATMENT), position = position_dodge(width = 0.9), alpha = 0.3) +  # Ribbon for credible intervals
  scale_color_okabe_ito(order= c(2, 1, 3), name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order= c(2, 1, 3), name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  scale_y_continuous(n.breaks=10) +
  labs(x = "Age (days)", y = "Repeat call interval duration (s)\n") +
  theme_clean()

rm(emm_REP_INT, emm_data, REP_INT_contrasts)

### EMMs, T*S ####
(emm_REP_INT <- emmeans(B_REP_int, pairwise ~ TREATMENT:SEX))

pd(emm_REP_INT)

equivalence_test(emm_REP_INT, range = rope_range(B_REP_int))

p_significance(emm_REP_INT, threshold = rope_range(B_REP_int)) 

### PLOTS: REP INTERVAL LENGTH ####
# coefficients:
posterior_desc <- describe_posterior(
  B_REP_int,
  effects = "fixed",
  component = "all",
  rope_range = rope_range(B_REP_int),
  test = c("p_direction", "pnificance", "rope"),
  centrality = "all",
  dispersion = TRUE
)
# drop sigma row
posterior_desc <- posterior_desc[-c(24, 1),]

# clean up labels:
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'b_', '')
posterior_desc$TREATMENT <- ifelse(grepl("TREATMENTSUB", posterior_desc$Parameter), "SC",
                                   ifelse(grepl("TREATMENTFLUT", posterior_desc$Parameter), "DT", "DC"))
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTSUB', 'SC')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTFLUT', 'DT')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'IAGE_zE2', 'Age^2')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'AGE_z', 'Age')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'SEXM', 'Male')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'WEIGHT_z', 'Weight offset')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'COMP_NORM_z', 'Competition')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'IGS_zE2', 'Group size^2')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'GS_z', 'Group size')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'RAIN_z', 'Monthly rainfall')

custom_order <- c('DT:Age^2:Male','SC:Age^2:Male','Age^2:Male',
                  'DT:Age:Male','SC:Age:Male','Age:Male', 
                  'DT:Monthly rainfall','SC:Monthly rainfall','Monthly rainfall',
                  'DT:Group size^2','SC:Group size^2','Group size^2',
                  'DT:Group size','SC:Group size','Group size',
                  'DT:Competition','SC:Competition','Competition', 
                  'DT:Weight offset','SC:Weight offset','Weight offset', 
                  'DT:Male','SC:Male','Male', 
                  'DT:Age^2','SC:Age^2',"Age^2", 
                  'DT:Age','SC:Age',"Age", 
                  'DT','SC') 
posterior_desc$Parameter <- factor(posterior_desc$Parameter, levels = custom_order)
posterior_desc$TREATMENT <- factor(posterior_desc$TREATMENT, levels = c("DC", "SC", "DT"))

# Coeff_REP_INT 700*800
ggplot(posterior_desc, aes(y = Parameter, x = Median, xmin = CI_low, xmax = CI_high, color=TREATMENT)) +
  geom_vline(xintercept = 0, color='grey', linetype = 'dotted', linewidth =1)+
  geom_point() +
  geom_errorbarh(height = 0) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT'))+
  labs(x = "Parameter value", y = "Parameters") +
 # ggtitle('REP call interval length')+
  theme_clean()+
  theme(legend.position="none")

# Ontogeny 30 - 130 ####
# get all needed values
sd_age <- sd(REP_data$REC_AGE_D)#
mean_age <- mean(REP_data$REC_AGE_D)#

range(REP_data$REC_AGE_D)# 31, 130
rec_age_c <- seq(30, 130, by=1)
age_z_vals <- (rec_age_c - mean(REP_data$REC_AGE_D))/sd(REP_data$REC_AGE_D)
age_z_vals <- seq(min(age_z_vals), max(age_z_vals), length.out = 20)

REP_INT_pred <- B_REP_int %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(REP_data$TREATMENT),
                                    SEX = levels(REP_data$SEX),
                                    AGE_z = age_z_vals,
                                    WEIGHT_z = mean(REP_data$WEIGHT_z),
                                    COMP_NORM_z = mean(REP_data$COMP_NORM_z),
                                    GS_z = mean(REP_data$GS_z),
                                    RAIN_z = mean(REP_data$RAIN_z)),
               re_formula = NA,  robust=T) 


#unscale AGE_z values:
REP_INT_pred$REC_AGE_D <- REP_INT_pred$AGE_z * sd_age + mean_age
# ensure right format
REP_INT_pred$REP_int <- REP_INT_pred$.epred
REP_INT_pred$TREATMENT <- factor(REP_INT_pred$TREATMENT, levels = c("CTRL", "SUB", "FLUT"))
REP_data$TREATMENT <- factor(REP_data$TREATMENT, levels = c("CTRL", "SUB", "FLUT"))

#800*500: REP_INT_TA
ggplot(REP_INT_pred, aes(x = REC_AGE_D, y = REP_int, color = TREATMENT, fill = TREATMENT)) +  
  geom_point(data=REP_data, aes(y=BEG_avg_Int))+
  stat_lineribbon(.width = .95) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order = c(2, 1, 3), alpha = 0.3, name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)", y = "Repeat call interval duration (s)\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Repeat call interval duration (s) \n", n.breaks = 10) +
  theme_clean()

# #700*500: REP_INT_TS
# Summarise predictions: Calculate median and credible intervals
REP_INT_summary <- REP_INT_pred %>%
  group_by(TREATMENT, SEX) %>%
  summarise(
    median_REP_int = median(REP_int),
    lower_CI = quantile(REP_int, 0.025),
    upper_CI = quantile(REP_int, 0.975)
  )

# Plot: Median estimates with credible intervals
ggplot(REP_INT_summary, aes(x = TREATMENT, y = median_REP_int, group = SEX, shape = SEX, colour=TREATMENT)) +
  geom_point(size = 3, position = position_dodge(width = 0.3)) +  # Shapes for SEX
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.2, 
                position = position_dodge(width = 0.3)) +  # Credible intervals
  scale_shape_manual(values = c(16, 17)) +  # Use shapes for SEX (e.g., circles and triangles)
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  scale_x_discrete(labels = c('DC', 'SC', 'DT')) +
  scale_y_continuous(n.breaks = 10)+
  labs(
    x = "Maternal treatment",
    y = "Repeat call interval duration (s)\n",
    shape = "Sex"
  ) +
  theme_clean() 

rm(REP_INT_pred, B_REP_int)

################################################################################
######################## Call rate #############################################
################################################################################

#### BAYES MODELS ##############################################################
# 1) Check for non-linearity ####
B_REP_rat_A <- brms::brm(formula = BEG_rate ~ AGE_z + (1|LITTER_CODE/ID),
                          data = REP_data, family = lognormal(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                          prior = priors, threads = threading(4),
                          file="B_REP_rat_A")
B_REP_rat_A <- add_criterion(B_REP_rat_A, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_REP_rat_TA)
# plot(B_REP_rat_TA)
# pp_check(B_REP_rat_TA, ndraws=100) # slight left shift
B_REP_rat_W <- brms::brm(formula = BEG_rate ~ WEIGHT_z + (1|LITTER_CODE/ID), 
                           data = REP_data, family = lognormal(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                           save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                           prior = priors, threads = threading(4),
                           file="B_REP_rat_W")
B_REP_rat_W <- add_criterion(B_REP_rat_W, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_REP_rat_W)
# plot(B_REP_rat_W)
# pp_check(B_REP_rat_W, ndraws=100)
B_REP_rat_CN <- brms::brm(formula = BEG_rate ~ COMP_NORM_z + (1|LITTER_CODE/ID),
                           data = REP_data, family = lognormal(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                           save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                           prior = priors, threads = threading(4),
                           file="B_REP_rat_CN")
B_REP_rat_CN <- add_criterion(B_REP_rat_CN, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_REP_rat_CN)
# plot(B_REP_rat_CN)
# pp_check(B_REP_rat_CN, ndraws=100)
B_REP_rat_G <- brms::brm(formula = BEG_rate ~ GS_z + (1|LITTER_CODE/ID),
                           data = REP_data, family = lognormal(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                           save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                           prior = priors, threads = threading(4),
                           file="B_REP_rat_G")
B_REP_rat_G <- add_criterion(B_REP_rat_G, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_REP_rat_G)
# plot(B_REP_rat_G)
# pp_check(B_REP_rat_G, ndraws=100)
B_REP_rat_R <- brms::brm(formula = BEG_rate ~ RAIN_z + (1|LITTER_CODE/ID),
                           data = REP_data, family = lognormal(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                           save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                           prior = priors, threads = threading(4),
                           file="B_REP_rat_R")
B_REP_rat_R <- add_criterion(B_REP_rat_R, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_REP_rat_R)
# plot(B_REP_rat_R)
# pp_check(B_REP_rat_R, ndraws=100)

# non linear
B_REP_rat_A2 <- brms::brm(formula = BEG_rate ~ AGE_z + I(AGE_z^2) + (1|LITTER_CODE/ID), 
                          data = REP_data, family = lognormal(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                          prior = priors, threads = threading(4),
                          file="B_REP_rat_A2")
B_REP_rat_A2 <- add_criterion(B_REP_rat_A2, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_REP_rat_A2)
# plot(B_REP_rat_A2)
# pp_check(B_REP_rat_A2, ndraws=100)
B_REP_rat_W2 <- brms::brm(formula = BEG_rate ~ WEIGHT_z + I(WEIGHT_z^2) + (1|LITTER_CODE/ID),
                            data = REP_data, family = lognormal(link='identity'),
                            chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                            save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                            prior = priors, threads = threading(4),
                            file="B_REP_rat_W2")
B_REP_rat_W2 <- add_criterion(B_REP_rat_W2, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_REP_rat_W2)
# plot(B_REP_rat_W2)
# pp_check(B_REP_rat_W2, ndraws=100)
B_REP_rat_CN2 <- brms::brm(formula = BEG_rate ~ COMP_NORM_z + I(COMP_NORM_z^2) + (1|LITTER_CODE/ID),
                            data = REP_data, family = lognormal(link='identity'),
                            chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                            save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                            prior = priors, threads = threading(4),
                            file="B_REP_rat_CN2")
B_REP_rat_CN2 <- add_criterion(B_REP_rat_CN2, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_REP_rat_CN2)
# plot(B_REP_rat_CN2)
# pp_check(B_REP_rat_CN2, ndraws=100)
B_REP_rat_G2 <- brms::brm(formula = BEG_rate ~ GS_z + I(GS_z^2) + (1|LITTER_CODE/ID),
                             data = REP_data, family = lognormal(link='identity'),
                             chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                             save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                             prior = priors, threads = threading(4),
                             file="B_REP_rat_G2")
B_REP_rat_G2 <- add_criterion(B_REP_rat_G2, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_REP_rat_G2)
# plot(B_REP_rat_G2)
# pp_check(B_REP_rat_G2, ndraws=100)
B_REP_rat_R2 <- brms::brm(formula = BEG_rate ~ RAIN_z + I(RAIN_z^2) + (1|LITTER_CODE/ID),
                            data = REP_data, family = lognormal(link='identity'),
                            chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                            save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                            prior = priors, threads = threading(4),
                            file="B_REP_rat_R2")
B_REP_rat_R2 <- add_criterion(B_REP_rat_R2, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_REP_rat_R2)
# plot(B_REP_rat_R2)
# pp_check(B_REP_rat_R2, ndraws=100)

loo(B_REP_rat_A, B_REP_rat_W, B_REP_rat_W2, B_REP_rat_CN, B_REP_rat_CN2,
    B_REP_rat_A2, B_REP_rat_G, B_REP_rat_G2, B_REP_rat_R, B_REP_rat_R2)

# best models
loo(B_REP_rat_CN, B_REP_rat_A2, B_REP_rat_W, 
    B_REP_rat_G, B_REP_rat_R)

loo_R2(B_REP_rat_CN)
loo_R2(B_REP_rat_A2)
loo_R2(B_REP_rat_W)
loo_R2(B_REP_rat_G)
loo_R2(B_REP_rat_R)

# save all summaries of the models:
model_list <- list(
  B_REP_rat_W, B_REP_rat_W2, B_REP_rat_CN, B_REP_rat_CN2, B_REP_rat_A,
  B_REP_rat_A2, B_REP_rat_G, B_REP_rat_G2, B_REP_rat_R, B_REP_rat_R2
)
run_summary_on_models(model_list, "REP_rat_sum.txt")

# save only relevant models
model_list <- list(
  B_REP_rat_CN, B_REP_rat_A2, B_REP_rat_W, 
  B_REP_rat_G, B_REP_rat_R)

run_summary_on_models(model_list, "REP_rat_best_sum.txt")

rm(B_REP_rat_W, B_REP_rat_W2, B_REP_rat_CN, B_REP_rat_CN2,
   B_REP_rat_A, B_REP_rat_G, B_REP_rat_G2, B_REP_rat_R, B_REP_rat_R2, model_list) 

# 2) Define models for REP rate ####
# use TAG2, TACN, TA2, TAW, TAR, TAS
B_REP_rat <- brms::brm(formula = BEG_rate ~ TREATMENT + AGE_z + I(AGE_z^2) + SEX + WEIGHT_z + COMP_NORM_z + GS_z + RAIN_z + 
                         TREATMENT:AGE_z + TREATMENT:I(AGE_z^2) +  TREATMENT:SEX + AGE_z:SEX + I(AGE_z^2):SEX +
                         TREATMENT:AGE_z:SEX + TREATMENT:I(AGE_z^2):SEX + (1|LITTER_CODE/ID),
                             data = REP_data, family = lognormal(link='identity'),
                             chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                             save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                             prior = priors, threads = threading(4),
                             file="B_REP_rat")
B_REP_rat <- add_criterion(B_REP_rat, c("loo", "loo_R2"), moment_match = TRUE)

B_REP_rat <- readRDS('models/B_REP_rat.rds')

#### RESULTS: REP rate ####
summary(B_REP_rat)

plot(B_REP_rat)
pp_check(B_REP_rat, ndraws = 100)

describe_posterior(
  B_REP_rat,
  effects = "all", #fixed vs all (for random effects)
  component = "all",
  rope_range = rope_range(B_REP_rat),  
  test = c("p_direction", "p_significance", "rope"),
  centrality = "all",
  dispersion = TRUE
)

loo_R2(B_REP_rat) 
bayes_R2(B_REP_rat)

performance::variance_decomposition(B_REP_rat)

#### EMMs at different ages (given non-linearity) ####
(emm_REP_RAT <- emmeans(B_REP_rat, ~ TREATMENT:AGE_z, at=list(AGE_z = get_age_vars())))

pd(emm_REP_RAT)
equivalence_test(emm_REP_RAT, range = rope_range(B_REP_rat)) 
p_significance(emm_REP_RAT, threshold = rope_range(B_REP_rat))

# contrasts BETWEEN treatments
(REP_RAT_contrasts <-pairs(emm_REP_RAT, by = c("AGE_z")))

pd(REP_RAT_contrasts)

equivalence_test(REP_RAT_contrasts, range = rope_range(B_REP_rat))

p_significance(REP_RAT_contrasts, threshold = rope_range(B_REP_rat))

# Plot emmeans: REP RATE ####
emm_data <- as.data.frame(emm_REP_RAT)
emm_data$REC_AGE_D <- emm_data$AGE_z * sd(REP_data$REC_AGE_D) + mean(REP_data$REC_AGE_D)
emm_data$REP_RAT <- exp(emm_data$emmean)
emm_data$HPD_low <- exp(emm_data$lower.HPD)
emm_data$HPD_high <- exp(emm_data$upper.HPD)

emm_data$TREATMENT <- factor(emm_data$TREATMENT, levels = c("CTRL", "SUB", "FLUT"))

# 700* 500: EMMs_REP_RAT
ggplot(emm_data, aes(x = as.factor(REC_AGE_D), y = REP_RAT, color = TREATMENT)) +
  geom_point(position = position_dodge(width = 0.9), size = 3) +  # Show estimates as poRATs
  geom_errorbar(aes(ymin = HPD_low, ymax = HPD_high), position = position_dodge(width = 0.9), width = 0.75) +
  scale_color_okabe_ito(order= c(2, 1, 3), name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  scale_y_continuous(n.breaks=10)+
  labs(x = "Age (days)", y = "Repeat call rate (calls/second)\n") +
  theme_clean()

# 800* 500: EMMs_REP_RAT_ribbon
ggplot(emm_data, aes(x = as.factor(REC_AGE_D), y = REP_RAT, color = TREATMENT, group = TREATMENT)) +
  geom_line(position = position_dodge(width = 0.9), linewidth = 1) +  # Connect estimates as lines
  geom_ribbon(aes(ymin = HPD_low, ymax = HPD_high, fill = TREATMENT), position = position_dodge(width = 0.9), alpha = 0.3) +  # Ribbon for credible RATervals
  scale_color_okabe_ito(order= c(2, 1, 3), name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order= c(2, 1, 3), name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  scale_y_continuous(n.breaks=10) +
  labs(x = "Age (days)", y = "Repeat call rate (calls/second)\n") +
  theme_clean()

rm(emm_REP_RAT, emm_data, REP_RAT_contrasts)

### PLOTS: REP RATE ####
# coefficients:
posterior_desc <- describe_posterior(
  B_REP_rat,
  effects = "fixed",
  component = "all",
  rope_range = rope_range(B_REP_rat),
  centrality = "all",
  dispersion = TRUE
)
# drop sigma row
posterior_desc <- posterior_desc[-c(23, 1),]

# clean up labels:
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'b_', '')
posterior_desc$TREATMENT <- ifelse(grepl("TREATMENTSUB", posterior_desc$Parameter), "SC",
                                   ifelse(grepl("TREATMENTFLUT", posterior_desc$Parameter), "DT", "DC"))
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTSUB', 'SC')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTFLUT', 'DT')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'IAGE_zE2', 'Age^2')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'AGE_z', 'Age')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'SEXM', 'Male')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'WEIGHT_z', 'Weight offset')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'COMP_NORM_z', 'Competition')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'GS_z', 'Group size')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'RAIN_z', 'Monthly rainfall')

custom_order <- c('Monthly rainfall',
                  'Group size',
                  'Competition', 
                  'Weight offset', 
                  'DT:Age^2:Male','SC:Age^2:Male','Age^2:Male',
                  'DT:Age:Male','SC:Age:Male','Age:Male', 
                  'DT:Male','SC:Male','Male', 
                  'DT:Age^2','SC:Age^2',"Age^2", 
                  'DT:Age','SC:Age',"Age", 
                  'DT','SC') 

posterior_desc$Parameter <- factor(posterior_desc$Parameter, levels = custom_order)
posterior_desc$TREATMENT <- factor(posterior_desc$TREATMENT, levels = c("DC", "SC", "DT"))

# Coeff_REP_RAT 700*800
ggplot(posterior_desc, aes(y = Parameter, x = Median, xmin = CI_low, xmax = CI_high, color=TREATMENT)) +
  geom_vline(xintercept = 0, color='grey', linetype = 'dotted', linewidth =1)+
  geom_point() +
  geom_errorbarh(height = 0) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT'))+
  labs(x = "Parameter value", y = "Parameters") +
 # ggtitle('REP call rate')+
  theme_clean()+
  theme(legend.position="none")

rm(posterior_desc)

# Ontogeny 30 - 130 
# get all needed values
(sd_age <- sd(REP_data$REC_AGE_D))#
(mean_age <- mean(REP_data$REC_AGE_D))#

range(REP_data$REC_AGE_D)# 31, 130
rec_age_c <- seq(30, 130, by=1)
age_z_vals <- (rec_age_c - mean(REP_data$REC_AGE_D))/sd(REP_data$REC_AGE_D)
age_z_vals <- seq(min(age_z_vals), max(age_z_vals), length.out = 20)

REP_RAT_pred <- B_REP_rat %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(REP_data$TREATMENT),
                                    SEX = levels(REP_data$SEX),
                                    AGE_z = age_z_vals,
                                    WEIGHT_z = mean(REP_data$WEIGHT_z),
                                    COMP_NORM_z = mean(REP_data$COMP_NORM_z),
                                    GS_z = mean(REP_data$GS_z),
                                    RAIN_z = mean(REP_data$RAIN_z)),
               re_formula = NA,  robust=T) 

#unscale AGE_z values:
REP_RAT_pred$REC_AGE_D <- REP_RAT_pred$AGE_z * sd_age + mean_age
# ensure right format
REP_RAT_pred$REP_rate <- REP_RAT_pred$.epred
REP_RAT_pred$TREATMENT <- factor(REP_RAT_pred$TREATMENT, levels = c("CTRL", "SUB", "FLUT"))
REP_data$TREATMENT <- factor(REP_data$TREATMENT, levels = c("CTRL", "SUB", "FLUT"))

#800*500: REP_RAT_TA
ggplot(REP_RAT_pred, aes(x = REC_AGE_D, y = REP_rate, color = TREATMENT, fill = TREATMENT)) +  
  geom_point(data=REP_data, aes(y=BEG_rate))+
  stat_lineribbon(.width = .95) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order = c(2, 1, 3), alpha = 0.3, name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)", y = "Repeat call rate (calls/s)\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Repeat call rate (calls/s) \n", n.breaks = 10) +
  theme_clean()

rm(REP_RAT_pred)

# Cleanup ####
rm(priors, REP_data)
