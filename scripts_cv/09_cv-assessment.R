# CREATED: 2025-05-01
# Prediction Model:  GAD


# https://www.tmwr.org/pre-proc-table
# decorrelate
# use PCA




# LIBRARIES ----
# Core
library(tidymodels)
tidymodels_prefer()

# Post-processing
library(yardstick)

# Helper packages
library(readr)       # for importing data
library(vip)         # for variable importance plots




# IMPORT DATA ----
# Results from original training data
preg_train <- read_rds("data_objects_cv/preg_dat_train_cv7_qval05.rds")
# length(unique(preg_dat$study_id))


# DATA SPLITTING AND RESAMPLING ----
set.seed(123)
# splits <- initial_split(preg_dat, prop = 0.80)
# preg_tr <- training(splits)
# preg_te <- testing(splits)

#resampling for parameter estimates
preg_rs <- vfold_cv(preg_train, v = 10)  #preg_tr



# MODEL SPECIFICATION ----

# * ELASTIC NET ----
elastic_mod <- linear_reg(
    penalty = tune(),
    mixture = tune()) %>% 
  set_engine("glmnet") %>% 
  set_mode("regression")

param_set <- extract_parameter_set_dials(elastic_mod)
# dials::penalty()
# dials::mixture()




# RECIPE ----
lr_rec <-
  recipe(ga_ad ~ ., data = preg_train) %>% 
  # step_rm(study_id, timepoint_day, visit) %>%
  step_rm(study_id, timepoint_day, visit, preg_smoke, mat_age_years, dm_16,
          CD8T, NK, Bcell, Mono, Gran) %>%
  step_dummy(all_nominal_predictors(), -all_outcomes()) %>%
  step_zv(all_predictors()) %>% 
  step_normalize(all_predictors())
  
# summary(lr_rec)
# 
# preg_prep <- lr_rec %>%
#   prep() %>% 
#   juice()
# 
# preg_prep %>% glimpse()



# WORKFLOW ----
elastic_wrk <-
  workflow() %>% 
  add_model(elastic_mod) %>% 
  add_recipe(lr_rec)



# TUNE ----
# 1. Automatic
elastic_grid_filling <- grid_space_filling(
  penalty(),
  mixture(),
  size = 500)

elastic_grid_random <- grid_random(
  penalty(),
  mixture(),
  size = 500)


#2. Manual specification of parameter range
param_set_1 <- param_set %>% 
  update(penalty= penalty(range = c(-2, 0)))  #log10 scale

elastic_grid_random <- grid_random(param_set_1, size = 1000)
summary(elastic_grid_random)




# Train and tune the model
set.seed(234)
elastic_res <- tune_grid(
  elastic_wrk,
  resamples = preg_rs,
  grid      = elastic_grid_random
  )

elastic_res %>% 
  collect_metrics()



# Plot multiple metrics
elastic_res %>%
  collect_metrics() %>%
  ggplot(aes(penalty, mean, color = .metric)) +
  # ggplot(aes(mixture, mean, color = .metric)) +
  # geom_errorbar(aes(
  #   ymin = mean - std_err,
  #   ymax = mean + std_err
  # ),
  # alpha = 0.5
  # ) +
  geom_line(linewidth = 1.5) +
  facet_wrap(~.metric, scales = "free", nrow = 2) +
  scale_x_log10() +
  theme(legend.position = "right")

top_models <-
  elastic_res %>% 
  show_best(metric = "rsq", n = 40) %>% 
  arrange(penalty) 
top_models %>% print(n = 40)


#Choose model parameter set
elastic_best_rsq <- elastic_res %>% 
  select_best(metric = "rsq")

elastic_err_rsq <- elastic_res %>% 
  select_by_one_std_err(metric = "rsq", desc(mixture))

elastic_best_rmse <- elastic_res %>% 
  select_best(metric = "rmse")






# lr_best <- 
#   lr_res %>% 
#   collect_metrics() %>% 
#   arrange(penalty) %>% 
#   slice(12)
# lr_best



elastic_workflow_final <- finalize_workflow(
  elastic_wrk,
  elastic_best_rsq
)


#Fit on training with selected parameters on training
fit_train <- 
  elastic_workflow_final %>% 
  fit(data = preg_train)



# Test on CV hold out sample
preg_test <- read_rds("data_objects_cv/preg_dat_test_cv7_qval05.rds")


test_pred <- predict(fit_train, new_data = preg_test) %>% 
  bind_cols(preg_test)

metrics(test_pred, truth = ga_ad, estimate = .pred)

cor(test_pred$ga_ad, test_pred$.pred)
plot(test_pred$ga_ad, test_pred$.pred)






# Final Elastic Net fit on all data
# n.b. This is not done yet. Must be tuned on an rset created from preg_dat
#estimate coefficients of model

met <- metric_set(rmse, mae, rsq, ccc)
ctrl <- control_resamples(verbose = TRUE)


fit_train <- 
  elastic_workflow_final %>% 
  fit(data = preg_train)

fit_train_coef <-
  fit_train %>%
  extract_fit_parsnip() %>%
  tidy()

fit_train_coef %>% 
  print(n = Inf)

elastic_fit <-
  elastic_workflow_final %>% 
  fit(data = rbind(preg_train, preg_test))


elastic_fit %>% 
  # pluck(".workflow", 1) %>%   
  extract_fit_parsnip() %>% 
  vip(num_features = 50, 
      aesthetics = list(color = 'blue', fill = "grey50", size = 0.3))








# last_fit(
#   lr_workflow_final,
#   splits
# ) %>%
#   collect_metrics()

#cg22501942; H19
#cg01402255; SLC6A3 dopamine transporter
#cg05859845; PRDM16 associated with maternal weight in first trimester

#cg22409100; AHRR
#cg15742777; PRSS23; linked to early pregnancy
#cg06837308; CYP1A1; https://www.sciencedirect.com/science/article/abs/pii/S1382668917303137
#cg08817867; GPR15; https://www.sciencedirect.com/science/article/abs/pii/S0022399920308886; smokers
#cg11263420; MYO1G; maternal smoking
#cg25960191; SLC7A11; placental expression
#cg07227049; GNG12; maternal smoking
#cg23054181; FAM169B
#cg02035751; GPR15; inflammation; smoking

#cg26233331; ZNF20



