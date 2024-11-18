#### Analysis script for Flutamide study: DIG (DIGGING) calls ###################
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

#Generic weakly informative prior: normal(0, 1);
priors <- c(set_prior("normal(0,1)", class = "Intercept"), set_prior("normal(0,1)", class='b'))

DIG_data <- readRDS('data/DIG_data.rds')

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


################################################################################
######################## Call length ###########################################
################################################################################
#### BAYES MODELS ##############################################################
# 1) Check for non-linearity ####
B_DIG_len_A <- brms::brm(formula = DIG_avg_Len ~ AGE_z + (1|LITTER_CODE/ID),
                           data = DIG_data, family = lognormal(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                           save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                           prior = priors, threads = threading(4),
                           file="B_DIG_len_A")
B_DIG_len_A <- add_criterion(B_DIG_len_A, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_DIG_len_A)
# plot(B_DIG_len_A)
# pp_check(B_DIG_len_A, ndraws=100)
B_DIG_len_W <- brms::brm(formula = DIG_avg_Len ~  WEIGHT_z + (1|LITTER_CODE/ID), 
                            data = DIG_data, family = lognormal(link='identity'),
                            chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                            save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                            prior = priors, threads = threading(4),
                            file="B_DIG_len_W")
B_DIG_len_W <- add_criterion(B_DIG_len_W, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_DIG_len_W)
# plot(B_DIG_len_W)
# pp_check(B_DIG_len_W, ndraws=100)
B_DIG_len_CN <- brms::brm(formula = DIG_avg_Len ~ COMP_NORM_z + (1|LITTER_CODE/ID) ,
                            data = DIG_data, family = lognormal(link='identity'),
                            chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                            save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                            prior = priors, threads = threading(4),
                            file="B_DIG_len_CN")
B_DIG_len_CN <- add_criterion(B_DIG_len_CN, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_DIG_len_CN)
# plot(B_DIG_len_CN)
# pp_check(B_DIG_len_CN, ndraws=100)
B_DIG_len_G <- brms::brm(formula = DIG_avg_Len ~ GS_z +(1|LITTER_CODE/ID),
                            data = DIG_data, family = lognormal(link='identity'),
                            chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                            save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                            prior = priors, threads = threading(4),
                            file="B_DIG_len_G")
B_DIG_len_G <- add_criterion(B_DIG_len_G, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_DIG_len_G)
# plot(B_DIG_len_G)
# pp_check(B_DIG_len_G, ndraws=100)
B_DIG_len_R <- brms::brm(formula = DIG_avg_Len ~ RAIN_z + (1|LITTER_CODE/ID),
                           data = DIG_data, family = lognormal(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                           save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                           prior = priors, threads = threading(4),
                           file="B_DIG_len_R")
B_DIG_len_R <- add_criterion(B_DIG_len_R, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_DIG_len_R)
# plot(B_DIG_len_R)
# pp_check(B_DIG_len_R, ndraws=100)

# non linear
B_DIG_len_A2 <- brms::brm(formula = DIG_avg_Len ~  AGE_z + I(AGE_z^2) + (1|LITTER_CODE/ID),
                          data = DIG_data, family = lognormal(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                          prior = priors, threads = threading(4),
                          file="B_DIG_len_A2")
B_DIG_len_A2 <- add_criterion(B_DIG_len_A2, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_DIG_len_A2)
# plot(B_DIG_len_A2)
# pp_check(B_DIG_len_A2, ndraws=100)
B_DIG_len_W2 <- brms::brm(formula = DIG_avg_Len ~ WEIGHT_z + I(WEIGHT_z^2) + (1|LITTER_CODE/ID),
                             data = DIG_data, family = lognormal(link='identity'),
                             chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                             save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                             prior = priors, threads = threading(4),
                             file="B_DIG_len_W2")
B_DIG_len_W2 <- add_criterion(B_DIG_len_W2, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_DIG_len_W2)
# plot(B_DIG_len_W2)
# pp_check(B_DIG_len_W2, ndraws=100)
B_DIG_len_CN2 <- brms::brm(formula = DIG_avg_Len ~  COMP_NORM_z + I(COMP_NORM_z^2) + (1|LITTER_CODE/ID),
                             data = DIG_data, family = lognormal(link='identity'),
                             chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                             save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                             prior = priors, threads = threading(4),
                             file="B_DIG_len_CN2")
B_DIG_len_CN2 <- add_criterion(B_DIG_len_CN2, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_DIG_len_CN2)
# plot(B_DIG_len_CN2)
# pp_check(B_DIG_len_CN2, ndraws=100)
B_DIG_len_G2 <- brms::brm(formula = DIG_avg_Len ~ GS_z +I(GS_z^2) + (1|LITTER_CODE/ID),
                             data = DIG_data, family = lognormal(link='identity'),
                             chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                             save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                             prior = priors, threads = threading(4),
                             file="B_DIG_len_G2")
B_DIG_len_G2 <- add_criterion(B_DIG_len_G2, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_DIG_len_G2)
# plot(B_DIG_len_G2)
# pp_check(B_DIG_len_G2, ndraws=100)
B_DIG_len_R2 <- brms::brm(formula = DIG_avg_Len ~ RAIN_z +I(RAIN_z^2) + (1|LITTER_CODE/ID),
                            data = DIG_data, family = lognormal(link='identity'),
                            chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                            save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                            prior = priors, threads = threading(4),
                            file="B_DIG_len_R2")
B_DIG_len_R2 <- add_criterion(B_DIG_len_R2, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_DIG_len_R2)
# plot(B_DIG_len_R2)
# pp_check(B_DIG_len_R2, ndraws=100)

loo(B_DIG_len_W, B_DIG_len_W2, B_DIG_len_CN, 
    B_DIG_len_CN2, B_DIG_len_A, B_DIG_len_A2, 
    B_DIG_len_G, B_DIG_len_G2, 
    B_DIG_len_R, B_DIG_len_R2)

# best models
loo(B_DIG_len_CN, B_DIG_len_A, B_DIG_len_W,  
    B_DIG_len_G, B_DIG_len_R)

loo_R2(B_DIG_len_CN)

loo_R2(B_DIG_len_A)
loo_R2(B_DIG_len_W)
loo_R2(B_DIG_len_G)
loo_R2(B_DIG_len_R)

# save all summaries of the models:
model_list <-
  list(
  B_DIG_len_W, B_DIG_len_W2, B_DIG_len_CN, B_DIG_len_A2,
  B_DIG_len_CN2, B_DIG_len_A2, B_DIG_len_G, B_DIG_len_G2, 
  B_DIG_len_R, B_DIG_len_R2
)
run_summary_on_models(model_list, "DIG_LEN_sum.txt")

# save only relevant models
model_list <- list(
  B_DIG_len_CN, B_DIG_len_A, B_DIG_len_W, 
  B_DIG_len_G, B_DIG_len_R)

run_summary_on_models(model_list, "DIG_LEN_best_sum.txt")

rm(B_DIG_len_W, B_DIG_len_W2, B_DIG_len_CN, 
   B_DIG_len_CN2, B_DIG_len_A, B_DIG_len_A2, B_DIG_len_G, B_DIG_len_G2, 
   B_DIG_len_R, B_DIG_len_R2, model_list) 

# 2) Define models for DIG LENGTH ####
B_DIG_len <- brms::brm(formula = DIG_avg_Len ~ TREATMENT + AGE_z + SEX + WEIGHT_z + COMP_NORM_z + GS_z + RAIN_z +
                               TREATMENT:AGE_z + TREATMENT:SEX + AGE_z:SEX + TREATMENT:AGE_z:SEX + (1|LITTER_CODE/ID),
                             data = DIG_data, family = lognormal(link='identity'),
                             chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                             save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                             prior = priors, threads = threading(4),
                             file="B_DIG_len")
B_DIG_len <- add_criterion(B_DIG_len, c("loo", "loo_R2"), moment_match = TRUE)

B_DIG_len <- readRDS('models/B_DIG_len.rds')

#### RESULTS: DIG LENGTH ####
summary(B_DIG_len)

plot(B_DIG_len)
pp_check(B_DIG_len, ndraws = 100)

describe_posterior(
  B_DIG_len,
  effects = "all", #fixed vs all (for random effects)
  component = "all",
  rope_range = rope_range(B_DIG_len),  
  test = c("p_direction", "p_significance", "rope"),
  centrality = "all",
  dispersion = TRUE
)

loo_R2(B_DIG_len) 

bayes_R2(B_DIG_len)

performance::variance_decomposition(B_DIG_len)

#### EMMs: TA ####
(emm_DIG_LEN <- emtrends(B_DIG_len, pairwise ~ TREATMENT, var='AGE_z'))

pd(emm_DIG_LEN) 

equivalence_test(emm_DIG_LEN, range = rope_range(B_DIG_len)) 

p_significance(emm_DIG_LEN, threshold = rope_range(B_DIG_len))


rm(emm_DIG_LEN)

### TS ####
(emm_DIG_LEN <- emmeans(B_DIG_len, pairwise ~ TREATMENT:SEX))

pd(emm_DIG_LEN) 

equivalence_test(emm_DIG_LEN, range = rope_range(B_DIG_len)) 

p_significance(emm_DIG_LEN, threshold = rope_range(B_DIG_len))

rm(emm_DIG_LEN)

### PLOTS: DIG LENGTH ####
# coefficients:
posterior_desc <- describe_posterior(
  B_DIG_len,
  effects = "fixed",
  component = "all",
  rope_range = rope_range(B_DIG_len),
  test = c("p_direction", "p_significance", "rope"),
  centrality = "all",
  dispersion = TRUE
)
# drop sigma row
posterior_desc <- posterior_desc[-c(17, 1),]

# clean up labels:
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'b_', '')
posterior_desc$TREATMENT <- ifelse(grepl("TREATMENTSUB", posterior_desc$Parameter), "SC",
                             ifelse(grepl("TREATMENTFLUT", posterior_desc$Parameter), "DT", "DC"))
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTFLUT', 'DT')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTSUB', 'SC')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'AGE_z', 'Age')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'SEXM', 'Male')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'WEIGHT_z', 'Weight offset')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'COMP_NORM_z', 'Competition')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'GS_z', 'Group size')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'RAIN_z', 'Monthly rainfall')

custom_order <- c(
                  'Monthly rainfall',
                  'Group size',
                  'Competition', 
                  'Weight offset', 
                  'DT:Age:Male','SC:Age:Male','Age:Male', 
                  'DT:Male','SC:Male','Male', 
                  'DT:Age','SC:Age',"Age", 
                  'DT','SC') 
posterior_desc$Parameter <- factor(posterior_desc$Parameter, levels = custom_order)
posterior_desc$TREATMENT <- factor(posterior_desc$TREATMENT, levels = c("DC", "SC", "DT"))
# Coeff_DIG_LEN 700*800
ggplot(posterior_desc, aes(y = Parameter, x = Median, xmin = CI_low, xmax = CI_high, color=TREATMENT)) +
  geom_vline(xintercept = 0, color='grey', linetype = 'dotted', linewidth =1)+
  geom_point() +
  geom_errorbarh(height = 0) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT'))+
  labs(x = "Parameter value", y = "Parameters") +
 # ggtitle('DIG call length')+
  theme_clean()+
  theme(legend.position="none")

# Ontogeny 30 - 130 
# get all needed values
(sd_age <- sd(DIG_data$REC_AGE_D))#25.35114
(mean_age <- mean(DIG_data$REC_AGE_D))#86.1989

range(DIG_data$REC_AGE_D)# 31, 130
rec_age_c <- seq(30, 130, by=1)
age_z_vals <- (rec_age_c - mean(DIG_data$REC_AGE_D))/sd(DIG_data$REC_AGE_D)
age_z_vals <- seq(min(age_z_vals), max(age_z_vals), length.out = 20)

DIG_LEN_pred <- B_DIG_len %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(DIG_data$TREATMENT),
                                    SEX = levels(DIG_data$SEX),
                                    AGE_z = age_z_vals,
                                    WEIGHT_z = mean(DIG_data$WEIGHT_z),
                                    COMP_NORM_z = mean(DIG_data$COMP_NORM_z),
                                    GS_z = mean(DIG_data$GS_z),
                                    RAIN_z = mean(DIG_data$RAIN_z)),
            re_formula = NA,  robust=T)

#unscale AGE_z values:
DIG_LEN_pred$REC_AGE_D <- DIG_LEN_pred$AGE_z * sd_age + mean_age
# ensure right format
DIG_LEN_pred$DIG_len <- DIG_LEN_pred$.epred
DIG_LEN_pred$TREATMENT <- factor(DIG_LEN_pred$TREATMENT, levels = c("CTRL", "SUB", "FLUT"))
DIG_data$TREATMENT <- factor(DIG_data$TREATMENT, levels = c("CTRL", "SUB", "FLUT"))

#800*500: DIG_LEN_TA
ggplot(DIG_LEN_pred, aes(x = REC_AGE_D, y = DIG_len, color = TREATMENT, fill = TREATMENT)) +  
  geom_point(data=DIG_data, aes(y=DIG_avg_Len))+
  stat_lineribbon(.width = .95) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order = c(2, 1, 3), alpha = 0.3, name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)", y = "Digging call length (s)\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Digging call length (s) \n", n.breaks = 10) +
  theme_clean()

# #700*500: DIG_LEN_TS
# Summarise predictions: Calculate median and credible intervals
DIG_LEN_summary <- DIG_LEN_pred %>%
  group_by(TREATMENT, SEX) %>%
  summarise(
    median_DIG_len = median(DIG_len),
    lower_CI = quantile(DIG_len, 0.025),
    upper_CI = quantile(DIG_len, 0.975)
  )

# Plot: Median estimates with credible intervals
ggplot(DIG_LEN_summary, aes(x = TREATMENT, y = median_DIG_len, group = SEX, shape = SEX, colour=TREATMENT)) +
  geom_point(size = 3, position = position_dodge(width = 0.3)) +  # Shapes for SEX
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.2, 
                position = position_dodge(width = 0.3)) +  # Credible intervals
  scale_shape_manual(values = c(16, 17)) +  # Use shapes for SEX (e.g., circles and triangles)
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  scale_x_discrete(labels = c('DC', 'SC', 'DT')) +
  labs(
    x = "Maternal treatment",
    y = "Digging call length (s)\n",
    shape = "Sex"
  ) +
  theme_clean() 

rm(DIG_LEN_pred, B_DIG_len)

################################################################################
######################## Call interval length ##################################
################################################################################

#### BAYES MODELS ##############################################################
# 1) Check for non-linearity ####
B_DIG_int_A <- brms::brm(formula = DIG_avg_Int ~ AGE_z + (1|LITTER_CODE/ID),
                          data = DIG_data, family = lognormal(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                          prior = priors, threads = threading(4),
                          file="B_DIG_int_A")
B_DIG_int_A <- add_criterion(B_DIG_int_A, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_DIG_int_A)
# plot(B_DIG_int_A)
# pp_check(B_DIG_int_A, ndraws=100)
B_DIG_int_W <- brms::brm(formula = DIG_avg_Int ~ WEIGHT_z + (1|LITTER_CODE/ID), 
                           data = DIG_data, family = lognormal(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                           save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                           prior = priors, threads = threading(4),
                           file="B_DIG_int_W")
B_DIG_int_W <- add_criterion(B_DIG_int_W, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_DIG_int_W)
# plot(B_DIG_int_W)
# pp_check(B_DIG_int_W, ndraws=100)
B_DIG_int_CN <- brms::brm(formula = DIG_avg_Int ~ COMP_NORM_z + (1|LITTER_CODE/ID),
                           data = DIG_data, family = lognormal(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                           save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                           prior = priors, threads = threading(4),
                           file="B_DIG_int_CN")
B_DIG_int_CN <- add_criterion(B_DIG_int_CN, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_DIG_int_CN)
# plot(B_DIG_int_CN)
# pp_check(B_DIG_int_CN, ndraws=100)
B_DIG_int_G <- brms::brm(formula = DIG_avg_Int ~ GS_z + (1|LITTER_CODE/ID),
                            data = DIG_data, family = lognormal(link='identity'),
                            chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                            save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                            prior = priors, threads = threading(4),
                            file="B_DIG_int_G")
B_DIG_int_G <- add_criterion(B_DIG_int_G, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_DIG_int_G)
# plot(B_DIG_int_G)
# pp_check(B_DIG_int_G, ndraws=100)
B_DIG_int_R <- brms::brm(formula = DIG_avg_Int ~ RAIN_z + (1|LITTER_CODE/ID),
                           data = DIG_data, family = lognormal(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                           save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                           prior = priors, threads = threading(4),
                           file="B_DIG_int_R")
B_DIG_int_R <- add_criterion(B_DIG_int_R, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_DIG_int_R)
# plot(B_DIG_int_R)
# pp_check(B_DIG_int_R, ndraws=100)

# non linear
B_DIG_int_A2 <- brms::brm(formula = DIG_avg_Int ~ AGE_z + I(AGE_z^2) + (1|LITTER_CODE/ID), 
                          data = DIG_data, family = lognormal(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                          prior = priors, threads = threading(4),
                          file="B_DIG_int_A2")
B_DIG_int_A2 <- add_criterion(B_DIG_int_A2, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_DIG_int_A2)
# plot(B_DIG_int_A2)
# pp_check(B_DIG_int_A2, ndraws=100)
B_DIG_int_W2 <- brms::brm(formula = DIG_avg_Int ~ WEIGHT_z + I(WEIGHT_z^2) + (1|LITTER_CODE/ID),
                            data = DIG_data, family = lognormal(link='identity'),
                            chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                            save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                            prior = priors, threads = threading(4),
                            file="B_DIG_int_W2")
B_DIG_int_W2 <- add_criterion(B_DIG_int_W2, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_DIG_int_W2)
# plot(B_DIG_int_W2)
# pp_check(B_DIG_int_W2, ndraws=100)
B_DIG_int_CN2 <- brms::brm(formula = DIG_avg_Int ~ COMP_NORM_z + I(COMP_NORM_z^2) + (1|LITTER_CODE/ID),
                            data = DIG_data, family = lognormal(link='identity'),
                            chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                            save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                            prior = priors, threads = threading(4),
                            file="B_DIG_int_CN2")
B_DIG_int_CN2 <- add_criterion(B_DIG_int_CN2, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_DIG_int_CN2)
# plot(B_DIG_int_CN2)
# pp_check(B_DIG_int_CN2, ndraws=100)
B_DIG_int_G2 <- brms::brm(formula = DIG_avg_Int ~ GS_z + I(GS_z^2) + (1|LITTER_CODE/ID),
                             data = DIG_data, family = lognormal(link='identity'),
                             chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                             save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                             prior = priors, threads = threading(4),
                             file="B_DIG_int_G2")
B_DIG_int_G2 <- add_criterion(B_DIG_int_G2, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_DIG_int_G2)
# plot(B_DIG_int_G2)
# pp_check(B_DIG_int_G2, ndraws=100)
B_DIG_int_R2 <- brms::brm(formula = DIG_avg_Int ~ RAIN_z +I(RAIN_z^2) + (1|LITTER_CODE/ID),
                            data = DIG_data, family = lognormal(link='identity'),
                            chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                            save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                            prior = priors, threads = threading(4),
                            file="B_DIG_int_R2")
B_DIG_int_R2 <- add_criterion(B_DIG_int_R2, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_DIG_int_R2)
# plot(B_DIG_int_R2)
# pp_check(B_DIG_int_R2, ndraws=100)

loo(B_DIG_int_W, B_DIG_int_W2, B_DIG_int_CN, B_DIG_int_A2,
    B_DIG_int_CN2, B_DIG_int_A, B_DIG_int_G, B_DIG_int_G2,
    B_DIG_int_R, B_DIG_int_R2)

# best models
loo(B_DIG_int_CN2, B_DIG_int_W, B_DIG_int_A,  
    B_DIG_int_G, B_DIG_int_R)

loo_R2(B_DIG_int_CN2)
loo_R2(B_DIG_int_A)
loo_R2(B_DIG_int_W)
loo_R2(B_DIG_int_G)
loo_R2(B_DIG_int_R)

# save all summaries of the models:
model_list <- list(
B_DIG_int_W, B_DIG_int_CN, B_DIG_int_W2, 
 B_DIG_int_CN2,  B_DIG_int_A, B_DIG_int_G, B_DIG_int_G2,
 B_DIG_int_R, B_DIG_int_R2, B_DIG_int_A2
)
run_summary_on_models(model_list, "DIG_int_sum.txt")

# save only relevant models: # use 
model_list <- list(
  B_DIG_int_CN2, B_DIG_int_A, B_DIG_int_W,  
  B_DIG_int_G, B_DIG_int_R)

run_summary_on_models(model_list, "DIG_int_best_sum.txt")

rm(B_DIG_int_W, B_DIG_int_CN, B_DIG_int_W2, 
     B_DIG_int_CN2,  B_DIG_int_A2, B_DIG_int_G, B_DIG_int_G2,
     B_DIG_int_R, B_DIG_int_R2, B_DIG_int_A, model_list)

# 2) Define models for DIG interval length ####
B_DIG_int <- brms::brm(formula = DIG_avg_Int ~ TREATMENT + AGE_z + SEX + WEIGHT_z + COMP_NORM_z + I(COMP_NORM_z^2) + GS_z + RAIN_z +
                               TREATMENT:AGE_z + TREATMENT:SEX + AGE_z:SEX + TREATMENT:AGE_z:SEX + (1|LITTER_CODE/ID),
                             data = DIG_data, family = lognormal(link='identity'),
                             chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                             save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                             prior = priors, threads = threading(4),
                             file="B_DIG_int")
B_DIG_int <- add_criterion(B_DIG_int, c("loo", "loo_R2"), moment_match = TRUE)

B_DIG_int <- readRDS('models/B_DIG_int.rds')

#### RESULTS: DIG interval length ####
summary(B_DIG_int)

plot(B_DIG_int)
pp_check(B_DIG_int, ndraws = 100)

describe_posterior(
  B_DIG_int,
  effects = "all", #fixed vs all (for random effects)
  component = "all",
  rope_range = rope_range(B_DIG_int),  
  test = c("p_direction", "p_significance", "rope"),
  centrality = "all",
  dispersion = TRUE
)

loo_R2(B_DIG_int) 
bayes_R2(B_DIG_int)

performance::variance_decomposition(B_DIG_int)


### EMMS: DIG interval length: TAS ####
(emm_DIG_INT <- emtrends(B_DIG_int, pairwise ~ TREATMENT:SEX, var = 'AGE_z'))

pd(emm_DIG_INT) 

equivalence_test(emm_DIG_INT, range = rope_range(B_DIG_int)) 

p_significance(emm_DIG_INT, threshold = rope_range(B_DIG_int)) 

rm(emm_DIG_INT)

### PLOTS: DIG INTERVAL DURATION ####
# coefficients:
posterior_desc <- describe_posterior(
  B_DIG_int,
  effects = "fixed",
  component = "all",
  rope_range = rope_range(B_DIG_int),
  test = c("p_direction", "p_significance", "rope"),
  centrality = "all",
  dispersion = TRUE
)
# drop sigma row
posterior_desc <- posterior_desc[-c(18,1),]

# clean up labels:
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'b_', '')
posterior_desc$TREATMENT <- ifelse(grepl("TREATMENTSUB", posterior_desc$Parameter), "SC",
                                   ifelse(grepl("TREATMENTFLUT", posterior_desc$Parameter), "DT", "DC"))
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTFLUT', 'DT')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTSUB', 'SC')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'AGE_z', 'Age')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'SEXM', 'Male')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'WEIGHT_z', 'Weight offset')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'COMP_NORM_z', 'Competition')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'ICompetitionE2', 'Competition2')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'GS_z', 'Group size')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'RAIN_z', 'Monthly rainfall')

custom_order <- c('Monthly rainfall',
                  'Group size',
                  'Competition2',
                  'Competition', 
                  'Weight offset', 
                  'DT:Age:Male','SC:Age:Male','Age:Male', 
                  'DT:Male','SC:Male','Male', 
                  'DT:Age','SC:Age',"Age", 
                  'DT','SC') 
posterior_desc$Parameter <- factor(posterior_desc$Parameter, levels = custom_order)
posterior_desc$TREATMENT <- factor(posterior_desc$TREATMENT, levels = c("DC", "SC", "DT"))

# Coeff_DIG_INT 700*800
ggplot(posterior_desc, aes(y = Parameter, x = Median, xmin = CI_low, xmax = CI_high, color=TREATMENT)) +
  geom_vline(xintercept = 0, color='grey', linetype = 'dotted', linewidth =1)+
  geom_point() +
  geom_errorbarh(height = 0) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT'))+
  labs(x = "Parameter value", y = "Parameters") +
 # ggtitle('DIG call interval duration')+
  theme_clean()+
  theme(legend.position="none")

# Ontogeny 30 - 130 
# get all needed values
(sd_age <- sd(DIG_data$REC_AGE_D))#
(mean_age <- mean(DIG_data$REC_AGE_D))#

range(DIG_data$REC_AGE_D)# 31, 130
rec_age_c <- seq(30, 130, by=1)
age_z_vals <- (rec_age_c - mean(DIG_data$REC_AGE_D))/sd(DIG_data$REC_AGE_D)
age_z_vals <- seq(min(age_z_vals), max(age_z_vals), length.out = 20)

DIG_INT_pred <- B_DIG_int %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(DIG_data$TREATMENT),
                                    SEX = levels(DIG_data$SEX),
                                    AGE_z = age_z_vals,
                                    WEIGHT_z = mean(DIG_data$WEIGHT_z),
                                    COMP_NORM_z = mean(DIG_data$COMP_NORM_z),
                                    GS_z = mean(DIG_data$GS_z),
                                    RAIN_z = mean(DIG_data$RAIN_z)),
              re_formula = NA,  robust=T)

#unscale AGE_z values:
DIG_INT_pred$REC_AGE_D <- DIG_INT_pred$AGE_z * sd_age + mean_age
# ensure right format
DIG_INT_pred$DIG_int <- DIG_INT_pred$.epred
DIG_INT_pred$TREATMENT <- factor(DIG_INT_pred$TREATMENT, levels = c("CTRL", "SUB", "FLUT"))
DIG_data$TREATMENT <- factor(DIG_data$TREATMENT, levels = c("CTRL", "SUB", "FLUT"))

#800*500: DIG_INT_TAS
ggplot(DIG_INT_pred, aes(x = REC_AGE_D, y = DIG_int, color = TREATMENT, fill = TREATMENT)) +  
  geom_point(data=DIG_data, aes(y=DIG_avg_Int))+
  stat_lineribbon(.width = .95) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order = c(2, 1, 3), alpha = 0.3, name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)", y = "Digging call interval duration (s)\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Digging call interval duration (s) \n", n.breaks = 10) +
  theme_clean()+
  facet_wrap(~SEX)

rm(DIG_INT_pred, B_DIG_int)

################################################################################
######################## Call rate #############################################
################################################################################

#### BAYES MODELS ##############################################################
# 1) Check for non-linearity
B_DIG_rat_A <- brms::brm(formula = DIG_rate ~ AGE_z + (1|LITTER_CODE/ID),
                          data = DIG_data, family = lognormal(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                          prior = priors, threads = threading(4),
                          file="B_DIG_rat_A")
B_DIG_rat_A <- add_criterion(B_DIG_rat_A, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_DIG_rat_A)
# plot(B_DIG_rat_A)
# pp_check(B_DIG_rat_A, ndraws=100) # slight left shift
B_DIG_rat_W <- brms::brm(formula = DIG_rate ~ WEIGHT_z + (1|LITTER_CODE/ID), 
                           data = DIG_data, family = lognormal(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                           save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                           prior = priors, threads = threading(4),
                           file="B_DIG_rat_W")
B_DIG_rat_W <- add_criterion(B_DIG_rat_W, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_DIG_rat_W)
# plot(B_DIG_rat_W)
# pp_check(B_DIG_rat_W, ndraws=100)
B_DIG_rat_CN <- brms::brm(formula = DIG_rate ~ COMP_NORM_z + (1|LITTER_CODE/ID),
                           data = DIG_data, family = lognormal(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                           save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                           prior = priors, threads = threading(4),
                           file="B_DIG_rat_CN")
B_DIG_rat_CN <- add_criterion(B_DIG_rat_CN, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_DIG_rat_CN)
# plot(B_DIG_rat_CN)
# pp_check(B_DIG_rat_CN, ndraws=100)
B_DIG_rat_G <- brms::brm(formula = DIG_rate ~ GS_z + (1|LITTER_CODE/ID),
                           data = DIG_data, family = lognormal(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                           save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                           prior = priors, threads = threading(4),
                           file="B_DIG_rat_G")
B_DIG_rat_G <- add_criterion(B_DIG_rat_G, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_DIG_rat_G)
# plot(B_DIG_rat_G)
# pp_check(B_DIG_rat_G, ndraws=100)
B_DIG_rat_R <- brms::brm(formula = DIG_rate ~ RAIN_z + (1|LITTER_CODE/ID),
                           data = DIG_data, family = lognormal(link='identity'),
                           chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                           save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                           prior = priors, threads = threading(4),
                           file="B_DIG_rat_R")
B_DIG_rat_R <- add_criterion(B_DIG_rat_R, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_DIG_rat_R)
# plot(B_DIG_rat_R)
# pp_check(B_DIG_rat_R, ndraws=100)

# non linear
B_DIG_rat_A2 <- brms::brm(formula = DIG_rate ~ TREATMENT + AGE_z + I(AGE_z^2) + TREATMENT:AGE_z + 
                            TREATMENT:I(AGE_z^2) + (1|LITTER_CODE/ID), 
                          data = DIG_data, family = lognormal(link='identity'),
                          chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                          save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                          prior = priors, threads = threading(4),
                          file="B_DIG_rat_A2")
B_DIG_rat_A2 <- add_criterion(B_DIG_rat_A2, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_DIG_rat_A2)
# plot(B_DIG_rat_A2)
# pp_check(B_DIG_rat_A2, ndraws=100)
B_DIG_rat_W2 <- brms::brm(formula = DIG_rate ~ WEIGHT_z + I(WEIGHT_z^2) + (1|LITTER_CODE/ID),
                            data = DIG_data, family = lognormal(link='identity'),
                            chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                            save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                            prior = priors, threads = threading(4),
                            file="B_DIG_rat_W2")
B_DIG_rat_W2 <- add_criterion(B_DIG_rat_W2, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_DIG_rat_W2)
# plot(B_DIG_rat_W2)
# pp_check(B_DIG_rat_W2, ndraws=100)
B_DIG_rat_CN2 <- brms::brm(formula = DIG_rate ~ COMP_NORM_z + I(COMP_NORM_z^2) + (1|LITTER_CODE/ID),
                            data = DIG_data, family = lognormal(link='identity'),
                            chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                            save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                            prior = priors, threads = threading(4),
                            file="B_DIG_rat_CN2")
B_DIG_rat_CN2 <- add_criterion(B_DIG_rat_CN2, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_DIG_rat_CN2)
# plot(B_DIG_rat_CN2)
# pp_check(B_DIG_rat_CN2, ndraws=100)
B_DIG_rat_G2 <- brms::brm(formula = DIG_rate ~ GS_z +I(GS_z^2) + (1|LITTER_CODE/ID),
                             data = DIG_data, family = lognormal(link='identity'),
                             chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                             save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                             prior = priors, threads = threading(4),
                             file="B_DIG_rat_G2")
B_DIG_rat_G2 <- add_criterion(B_DIG_rat_G2, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_DIG_rat_G2)
# plot(B_DIG_rat_G2)
# pp_check(B_DIG_rat_G2, ndraws=100)
B_DIG_rat_R2 <- brms::brm(formula = DIG_rate ~ RAIN_z +I(RAIN_z^2) + (1|LITTER_CODE/ID),
                            data = DIG_data, family = lognormal(link='identity'),
                            chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                            save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                            prior = priors, threads = threading(4),
                            file="B_DIG_rat_R2")
B_DIG_rat_R2 <- add_criterion(B_DIG_rat_R2, c("loo", "loo_R2"), moment_match = TRUE)
# summary(B_DIG_rat_R2)
# plot(B_DIG_rat_R2)
# pp_check(B_DIG_rat_R2, ndraws=100)

loo(B_DIG_rat_W, B_DIG_rat_W2, B_DIG_rat_CN, B_DIG_rat_CN2, B_DIG_rat_A2,
    B_DIG_rat_A, B_DIG_rat_G, B_DIG_rat_G2, B_DIG_rat_R, B_DIG_rat_R2)

# best models
loo(B_DIG_rat_CN2, B_DIG_rat_A, B_DIG_rat_W,  
    B_DIG_rat_G, B_DIG_rat_R)

loo_R2(B_DIG_rat_CN2)
loo_R2(B_DIG_rat_A)
loo_R2(B_DIG_rat_W)
loo_R2(B_DIG_rat_G)
loo_R2(B_DIG_rat_R)

# save all summaries of the models:
model_list <- list(
  B_DIG_rat_W, B_DIG_rat_W2, B_DIG_rat_CN, B_DIG_rat_CN2, B_DIG_rat_A2,
  B_DIG_rat_A, B_DIG_rat_G, B_DIG_rat_G2, B_DIG_rat_R, B_DIG_rat_R2
)
run_summary_on_models(model_list, "DIG_rat_sum.txt")

# save only relevant models
model_list <- list(
  B_DIG_rat_CN2, B_DIG_rat_A, B_DIG_rat_W,  
  B_DIG_rat_G, B_DIG_rat_R)

run_summary_on_models(model_list, "DIG_rat_best_sum.txt")

rm(B_DIG_rat_W, B_DIG_rat_W2, B_DIG_rat_CN, B_DIG_rat_CN2, B_DIG_rat_A,
   B_DIG_rat_A2, B_DIG_rat_G, B_DIG_rat_G2, B_DIG_rat_R, B_DIG_rat_R2, model_list) 

# 2) Define models for DIG rate ####
B_DIG_rat <- brms::brm(formula = DIG_rate ~ TREATMENT + AGE_z + SEX + WEIGHT_z + COMP_NORM_z + I(COMP_NORM_z^2) + GS_z + RAIN_z +
                               TREATMENT:AGE_z + TREATMENT:SEX + AGE_z:SEX + TREATMENT:AGE_z:SEX + (1|LITTER_CODE/ID),
                             data = DIG_data, family = lognormal(link='identity'),
                             chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 25, adapt_delta=0.99),
                             save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                             prior = priors, threads = threading(4),
                             file="B_DIG_rat")
B_DIG_rat <- add_criterion(B_DIG_rat, c("loo", "loo_R2"), moment_match = TRUE)

B_DIG_rat <- readRDS('models/B_DIG_rat.rds')

#### RESULTS: DIG rate ####
summary(B_DIG_rat)

plot(B_DIG_rat)
pp_check(B_DIG_rat, ndraws = 100) # shifted to left!

describe_posterior(
  B_DIG_rat,
  effects = "all", #fixed vs all (for random effects)
  component = "all",
  rope_range = rope_range(B_DIG_rat),  
  test = c("p_direction", "p_significance", "rope"),
  centrality = "all",
  dispersion = TRUE
)

loo_R2(B_DIG_rat) 
bayes_R2(B_DIG_rat)

performance::variance_decomposition(B_DIG_rat)


#### EMMs TA ####
(emm_DIG_RAT <- emtrends(B_DIG_rat, pairwise ~ TREATMENT, var='AGE_z'))

pd(emm_DIG_RAT) 


equivalence_test(emm_DIG_RAT, range = rope_range(B_DIG_rat)) 

p_significance(emm_DIG_RAT, threshold = rope_range(B_DIG_rat))

rm(emm_DIG_RAT)

### TS ####
(emm_DIG_RAT <- emmeans(B_DIG_rat, pairwise ~ TREATMENT:SEX))

pd(emm_DIG_RAT) 

equivalence_test(emm_DIG_RAT, range = rope_range(B_DIG_rat)) 

p_significance(emm_DIG_RAT, threshold = rope_range(B_DIG_rat))

rm(emm_DIG_RAT)

### PLOTS: DIG RATE ####
# coefficients:
posterior_desc <- describe_posterior(
  B_DIG_rat,
  effects = "fixed",
  component = "all",
  rope_range = rope_range(B_DIG_rat),
  centrality = "all",
  dispersion = TRUE
)
# drop sigma row
posterior_desc <- posterior_desc[-c(18, 1),]

# clean up labels:
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'b_', '')
posterior_desc$TREATMENT <- ifelse(grepl("TREATMENTSUB", posterior_desc$Parameter), "SC",
                                   ifelse(grepl("TREATMENTFLUT", posterior_desc$Parameter), "DT", "DC"))
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTFLUT', 'DT')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTSUB', 'SC')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'AGE_z', 'Age')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'SEXM', 'Male')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'WEIGHT_z', 'Weight offset')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'COMP_NORM_z', 'Competition')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'ICompetitionE2', 'Competition2')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'GS_z', 'Group size')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'RAIN_z', 'Monthly rainfall')

custom_order <- c('Monthly rainfall',
                  'Group size',
                  'Competition2', 
                  'Competition', 
                  'Weight offset',
                  'DT:Age:Male','SC:Age:Male','Age:Male', 
                  'DT:Male','SC:Male','Male', 
                  'DT:Age','SC:Age',"Age", 
                  'DT','SC') 
posterior_desc$Parameter <- factor(posterior_desc$Parameter, levels = custom_order)
posterior_desc$TREATMENT <- factor(posterior_desc$TREATMENT, levels = c("DC", "SC", "DT"))

# Coeff_DIG_RAT 700*800
ggplot(posterior_desc, aes(y = Parameter, x = Median, xmin = CI_low, xmax = CI_high, color=TREATMENT)) +
  geom_vline(xintercept = 0, color='grey', linetype = 'dotted', linewidth =1)+
  geom_point() +
  geom_errorbarh(height = 0) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT'))+
  labs(x = "Parameter value", y = "Parameters") +
 # ggtitle('DIG call rate')+
  theme_clean()+
  theme(legend.position="none")

# Ontogeny 30 - 130 ####
# get all needed values
sd_age <- sd(DIG_data$REC_AGE_D)#
mean_age <- mean(DIG_data$REC_AGE_D)#

range(DIG_data$REC_AGE_D)# 31, 130
rec_age_c <- seq(30, 130, by=1)
age_z_vals <- (rec_age_c - mean(DIG_data$REC_AGE_D))/sd(DIG_data$REC_AGE_D)
age_z_vals <- seq(min(age_z_vals), max(age_z_vals), length.out = 20)

DIG_RAT_pred <- B_DIG_rat %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(DIG_data$TREATMENT),
                                    SEX = levels(DIG_data$SEX),
                                    AGE_z = age_z_vals,
                                    WEIGHT_z = mean(DIG_data$WEIGHT_z),
                                    COMP_NORM_z = mean(DIG_data$COMP_NORM_z),
                                    GS_z = mean(DIG_data$GS_z),
                                    RAIN_z = mean(DIG_data$RAIN_z)),
             re_formula = NA,  robust=T)

#unscale AGE_z values:
DIG_RAT_pred$REC_AGE_D <- DIG_RAT_pred$AGE_z * sd_age + mean_age
# ensure right format
DIG_RAT_pred$DIG_rate <- DIG_RAT_pred$.epred
DIG_RAT_pred$TREATMENT <- factor(DIG_RAT_pred$TREATMENT, levels = c("CTRL", "SUB", "FLUT"))
DIG_data$TREATMENT <- factor(DIG_data$TREATMENT, levels = c("CTRL", "SUB", "FLUT"))

#800*500: DIG_RAT_TA
ggplot(DIG_RAT_pred, aes(x = REC_AGE_D, y = DIG_rate, color = TREATMENT, fill = TREATMENT)) +  
  geom_point(data=DIG_data, aes(y=DIG_rate))+
  stat_lineribbon(.width = .95) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order = c(2, 1, 3), alpha = 0.3, name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)", y = "Digging call rate (calls/s)\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Digging call rate (calls/s) \n", n.breaks = 10) +
  theme_clean()

# #700*500: DIG_RAT_TS
# Summarise predictions: Calculate median and credible intervals
DIG_RAT_summary <- DIG_RAT_pred %>%
  group_by(TREATMENT, SEX) %>%
  summarise(
    median_DIG_rat = median(DIG_rate),
    lower_CI = quantile(DIG_rate, 0.025),
    upper_CI = quantile(DIG_rate, 0.975)
  )

# Plot: Median estimates with credible intervals
ggplot(DIG_RAT_summary, aes(x = TREATMENT, y = median_DIG_rat, group = SEX, shape = SEX, colour=TREATMENT)) +
  geom_point(size = 3, position = position_dodge(width = 0.3)) +  # Shapes for SEX
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.2, 
                position = position_dodge(width = 0.3)) +  # Credible intervals
  scale_shape_manual(values = c(16, 17)) +  # Use shapes for SEX (e.g., circles and triangles)
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\ntreatment", labels = c('DC', 'SC', 'DT')) +
  scale_x_discrete(labels = c('DC', 'SC', 'DT')) +
  scale_y_continuous(n.breaks = 10)+
  labs(
    x = "Maternal treatment",
    y = "Digging call rate (calls/second)\n",
    shape = "Sex"
  ) +
  theme_clean() 

rm(DIG_RAT_pred, B_DIG_rat)

# Cleanup ####
rm(DIG_data, priors)
