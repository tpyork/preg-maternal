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
# preg_dat <- read_rds("data_objects/preg_dat_25cov_qval10_.rds")
preg_dat <- read_rds("data_objects/preg_dat_25cov_qval01.rds")
# length(unique(preg_dat$study_id))


# DATA SPLITTING AND RESAMPLING ----
set.seed(123)
# splits <- initial_split(preg_dat, prop = 0.80)
# 
# preg_tr <- training(splits)
# preg_te <- testing(splits)

#resampling
preg_rs <- vfold_cv(preg_dat, v = 10)  #preg_tr



# MODEL SPECIFICATION ----

# * Penalized Linear Regression - LASSO ----
lasso_mod <- linear_reg(penalty = tune(), mixture = 1) %>%
  set_engine("glmnet")


# * XGBoost ----
xgb_mod <- boost_tree(
    trees = 300,    #1000
    tree_depth = tune(), min_n = tune(),
    loss_reduction = tune(),                     ## first three: model complexity
    sample_size = tune(), mtry = tune(),         ## randomness
    learn_rate = tune()                          ## step size
  ) %>%
  set_engine("xgboost") %>%
  set_mode("regression")


# * MARS ----
mars_mod <- mars(
    mode         = "regression",
    engine       = "earth",
    num_terms    = NULL,
    prod_degree  = NULL,
    prune_method = NULL
  )


# * ELASTIC NET ----
elastic_mod <- linear_reg(
    penalty = tune(),
    mixture = tune()) %>% 
  set_engine("glmnet")



# RECIPE ----
lr_rec <-
  recipe(ga_ad ~ ., data = preg_dat) %>% 
  # step_rm(study_id, timepoint_day, visit) %>%
  step_rm(study_id, timepoint_day, visit, preg_smoke, mat_age_years, dm_16,
          CD8T, NK, Bcell, Mono, Gran) %>%
  step_dummy(all_nominal_predictors(), -all_outcomes()) %>%
  step_zv(all_predictors()) %>% 
  step_normalize(all_predictors())
  
# summary(lr_recipe)
# 
# preg_prep <- lr_recipe %>%
#   prep() %>% 
#   juice()
# 
# preg_prep %>% glimpse()



# WORKFLOW ----
lasso_wrk <- 
  workflow() %>% 
  add_model(lasso_mod) %>% 
  add_recipe(lr_rec)

xgb_wrk <- 
  workflow() %>% 
  add_model(xgb_mod) %>% 
  add_recipe(lr_rec)

mars_wrk <- 
  workflow() %>% 
  add_model(mars_mod) %>% 
  add_recipe(lr_rec)

elastic_wrk <-
  workflow() %>% 
  add_model(elastic_mod) %>% 
  add_recipe(lr_rec)



# TUNE ----
lasso_grid <- tibble(penalty = seq(0.001, 10, length.out = 200))
lasso_grid %>% print(n = Inf)


xgb_grid <- grid_space_filling(
  tree_depth(),
  min_n(),
  loss_reduction(),
  sample_size = sample_prop(),
  finalize(mtry(), preg_dat),
  learn_rate(),
  size = 50
)


elastic_grid <- grid_space_filling(
  penalty(),
  mixture(),
  size = 500)



# Train and tune the model
lasso_res <- 
  lasso_wrk %>% 
  tune_grid(resamples = preg_rs, #rset object
            grid      = lasso_grid,
            control   = control_grid(save_pred = TRUE), #needed to get data for ROC curve
            metrics   = metric_set(rmse, rsq))


set.seed(234)
xgb_res <- tune_grid(
  xgb_wrk,
  resamples = preg_rs,
  grid      = xgb_grid,
  control   = control_grid(save_pred = TRUE)
)


elastic_res <- tune_grid(
  elastic_wrk,
  resamples = preg_rs,
  grid      = elastic_grid
  )


lasso_res %>% 
  collect_metrics() %>% 
  print(n = Inf)

xgb_res %>% 
  collect_metrics %>% 
  filter(.metric == 'rsq') %>% 
  filter(mean > 0)

elastic_res %>% 
  collect_metrics()



# Plot multiple metrics
lasso_res %>%
  collect_metrics() %>%
  ggplot(aes(penalty, mean, color = .metric)) +
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

# Plot single metric
lr_plot <- 
  lr_res %>% 
  collect_metrics() %>% 
  filter(.metric == "rmse") %>% 
  ggplot(aes(x = penalty, y = mean)) + 
  geom_point() + 
  geom_line() + 
  ylab("RMSE") +
  scale_x_log10(labels = scales::label_number())
lr_plot 

# XGB Plot
xgb_res %>%
  collect_metrics() %>%
  filter(.metric == "rsq") %>%
  select(mean, mtry:sample_size) %>%
  pivot_longer(mtry:sample_size,
               values_to = "value",
               names_to = "parameter"
  ) %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_point(alpha = 0.8, show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "rsq")


elastic_res %>%
  collect_metrics() %>%
  ggplot(aes(penalty, mean, color = .metric)) +
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
  lasso_res %>% 
  show_best(metric = "rsq", n = 40) %>% 
  arrange(penalty) 
top_models %>% print(n = 40)


top_models <-
  elastic_res %>% 
  show_best(metric = "rsq", n = 40) %>% 
  arrange(penalty) 
top_models %>% print(n = 40)


elastic_best <- elastic_res %>% 
  select_best(metric = "rmse")


lowest_rmse <- lasso_res %>%
  select_best(metric = "rmse")

highest_rsq <- lasso_res %>% 
  select_best(metric = "rsq")

# lowest_rmse <- elastic_res %>%
#   select_best(metric = "rmse")



# Manual selection
lasso_best <- 
  lasso_res %>% 
  collect_metrics() %>% 
  arrange(penalty) %>% 
  slice(79)
lasso_best

lasso_workflow_final <- finalize_workflow(
  lasso_wrk,
  lasso_best
  # lowest_rmse
)


elastic_workflow_final <- finalize_workflow(
  elastic_wrk,
  elastic_best
  # lowest_rmse
)



# FIT ----
met <- metric_set(rmse, mae, rsq, ccc)
ctrl <- control_resamples(verbose = TRUE)


preg_res <- 
  elastic_workflow_final %>% 
  fit_resamples(resamples = preg_rs, 
                metrics   = met,
                control   = ctrl)

collect_metrics(preg_res) 


preg_fit <-
  elastic_workflow_final %>% 
    fit(data = preg_dat)


preg_fit_coef <-
  preg_fit %>%
    extract_fit_parsnip() %>%
    tidy()

preg_fit_coef %>% 
    print(n = Inf)

# readr::write_rds(preg_fit_coef, file = "data_objects/preg_fit_coef.rds")



skimr::skim(preg_fit_coef)
summary(preg_fit_coef$estimate[-1])

which(preg_fit_coef$estimate[-1] == max(preg_fit_coef$estimate[-1]))

preg_fit_coef[381,]




# variable importance plot
preg_fit %>% 
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


# EXTRACT DATA FROM RSET ----

split_1 <- preg_rs$splits[[1]]

analysis(split_1)
assessment(split_1)


all_splits <- purrr::map(preg_rs$splits, ~ list(
  train = analysis(.x),
  test = assessment(.x)
))

all_splits[[10]]$test %>% class()

